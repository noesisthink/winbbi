package gobigwig

import (
	"bytes"
	"compress/zlib"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
)

// bwSummaryOnDisk 对应 zoom data 的磁盘格式
type bwSummaryOnDisk struct {
	ChromId    uint32
	Start      uint32
	End        uint32
	ValidCount uint32
	MinVal     float32
	MaxVal     float32
	SumData    float32
	SumSquares float32
}

// bwSummary 内存中的 summary 结构
type bwSummary struct {
	ChromId    uint32
	Start      uint32
	End        uint32
	ValidCount uint32
	MinVal     float32
	MaxVal     float32
	SumData    float32
	SumSquares float32
}

// bwGetBestZoom 选择最合适的zoom level
// desiredReduction: 期望的reduction level（每个summary代表多少个碱基）
// 返回zoom level的索引，如果没有合适的返回 -1
func bwSelectBestZoomLevel(zhdr *bwZoomHdr_t, desiredReduction uint32) int {
	if zhdr == nil || desiredReduction <= 1 || len(zhdr.Level) == 0 {
		return -1
	}
	bestIdx := -1
	var closestDiff uint32 = ^uint32(0) // max uint32
	for i, level := range zhdr.Level {
		if desiredReduction >= level {
			diff := desiredReduction - level
			if diff < closestDiff {
				closestDiff = diff
				bestIdx = i
			}
		}
	}
	return bestIdx
}



// bwGetBestZoomClosest 返回最接近目标缩放因子的层级（允许超过）
func bwGetBestZoomClosest(zhdr *bwZoomHdr_t, desiredReduction uint32) int {
	if zhdr == nil || len(zhdr.Level) == 0 {
		return -1
	}

	bestIdx := -1
	closestDiff := uint32(math.MaxUint32)

	for i, level := range zhdr.Level {
		diff := uint32(math.Abs(float64(int64(desiredReduction) - int64(level))))
		if diff < closestDiff {
			closestDiff = diff
			bestIdx = i
		}
	}

	return bestIdx
}

// bwReadZoomIndex 读取指定zoom level的索引
func bwReadZoomIndex(fp *bigWigFile_t, indexOffset uint64) (*bwRTree_t, error) {
	if indexOffset == 0 {
		return nil, errors.New("invalid index offset")
	}

	// 使用已有的 readRTreeIdx 函数
	idx, err := readRTreeIdx(fp, indexOffset)
	if err != nil {
		return nil, err
	}

	// 读取根节点
	root, err := bwGetRTreeNode(fp, idx.RootOffset)
	if err != nil {
		return nil, err
	}
	idx.Root = root

	return idx, nil
}

// bwGetSummariesInRegion 从指定zoom level获取区间内的summaries
func bwGetSummariesInRegion(fp *bigWigFile_t, zoomIdx int, chrom string, start, end uint32) ([]*bwSummary, error) {
	if fp.Hdr == nil || len(fp.Hdr.ZoomHdrs) == 0 {
		return nil, errors.New("no zoom headers available")
	}

	zhdr := fp.Hdr.ZoomHdrs[0]
	if zoomIdx < 0 || zoomIdx >= len(zhdr.Level) {
		return nil, fmt.Errorf("invalid zoom index: %d", zoomIdx)
	}

	tid := bwGetTid(fp, chrom)
	if tid == ^uint32(0) {
		return nil, fmt.Errorf("chromosome not found: %s", chrom)
	}

	// 读取或使用缓存的索引
	var zoomTree *bwRTree_t
	var err error
	
	if zhdr.Idx[zoomIdx] == nil {
		zoomTree, err = bwReadZoomIndex(fp, zhdr.IndexOffset[zoomIdx])
		if err != nil {
			return nil, err
		}
		zhdr.Idx[zoomIdx] = zoomTree
	} else {
		zoomTree = zhdr.Idx[zoomIdx]
	}

	// 查找重叠的数据块
	blocks := walkRTreeNodes(fp, zoomTree.Root, tid, start, end)
	if blocks == nil || blocks.N == 0 {
		return nil, nil
	}

	// 读取并解析summaries
	summaries := []*bwSummary{}
	compressed := fp.Hdr.bufsize > 0

	for i := uint64(0); i < blocks.N; i++ {
		// 定位到数据块
		if bwSetPos(fp, blocks.Offset[i]) != 0 {
			return nil, errors.New("failed to seek to data block")
		}

		// 读取数据
		compBuf := make([]byte, blocks.Size[i])
		n, err := fp.URL.Read(compBuf)
		if err != nil || n != int(blocks.Size[i]) {
			return nil, fmt.Errorf("failed to read data block: %v", err)
		}

		var data []byte
		if compressed {
			data, err = decompressZlibSimple(compBuf)
			if err != nil {
				return nil, fmt.Errorf("failed to decompress: %v", err)
			}
		} else {
			data = compBuf
		}

		// 解析summaries
		// 每个summary的大小是32字节
		summarySize := 32
		numSummaries := len(data) / summarySize

		for j := 0; j < numSummaries; j++ {
			offset := j * summarySize
			if offset+summarySize > len(data) {
				break
			}

			sum := &bwSummary{
				ChromId:    binary.LittleEndian.Uint32(data[offset : offset+4]),
				Start:      binary.LittleEndian.Uint32(data[offset+4 : offset+8]),
				End:        binary.LittleEndian.Uint32(data[offset+8 : offset+12]),
				ValidCount: binary.LittleEndian.Uint32(data[offset+12 : offset+16]),
				MinVal:     math.Float32frombits(binary.LittleEndian.Uint32(data[offset+16 : offset+20])),
				MaxVal:     math.Float32frombits(binary.LittleEndian.Uint32(data[offset+20 : offset+24])),
				SumData:    math.Float32frombits(binary.LittleEndian.Uint32(data[offset+24 : offset+28])),
				SumSquares: math.Float32frombits(binary.LittleEndian.Uint32(data[offset+28 : offset+32])),
			}

			// 过滤出在查询范围内且染色体匹配的summaries
			if sum.ChromId == tid && sum.Start < end && sum.End > start {
				summaries = append(summaries, sum)
			}
		}
	}

	return summaries, nil
}

// bwGetValuesFromZoom 使用指定的zoom level获取区间的值（带详细调试输出）
// summaryType: "mean", "max", "min", "coverage", "sum"
func bwGetValuesFromZoom(fp *bigWigFile_t, zoomIdx int, chrom string, start, end uint32, numBins int, summaryType string) ([]float32, error) {
	summaries, err := bwGetSummariesInRegion(fp, zoomIdx, chrom, start, end)
	if err != nil {
		return nil, err
	}
	values := make([]float32, numBins)
	for i := range values {
		values[i] = float32(math.NaN())
	}
	if len(summaries) == 0 {
		return values, nil
	}
	binSize := float64(end-start) / float64(numBins)
	for i := 0; i < numBins; i++ {
		binStart := start + uint32(float64(i)*binSize)
		binEnd := start + uint32(float64(i+1)*binSize)
		var sumData float64
		var validCount uint32
		var minVal float32 = float32(math.Inf(1))
		var maxVal float32 = float32(math.Inf(-1))
		overlapCount := 0
		// 找到与当前bin重叠的summaries
		for _, sum := range summaries {
			if sum.End <= binStart || sum.Start >= binEnd {
				continue
			}
			overlapStart := max32(sum.Start, binStart)
			overlapEnd := min32(sum.End, binEnd)
			overlap := overlapEnd - overlapStart
			if overlap == 0 {
				continue
			}
			sumWidth := sum.End - sum.Start
			overlapFactor := float64(overlap) / float64(sumWidth)

			validCount += uint32(float64(sum.ValidCount) * overlapFactor)
			sumData += float64(sum.SumData) * overlapFactor
			if sum.MaxVal > maxVal {
				maxVal = sum.MaxVal
			}
			if sum.MinVal < minVal {
				minVal = sum.MinVal
			}
			overlapCount++
		}

		// 根据summaryType计算最终值
		if validCount > 0 {
			switch summaryType {
			case "mean", "average":
				values[i] = float32(sumData / float64(validCount))
			case "max", "maximum":
				values[i] = maxVal
			case "min", "minimum":
				values[i] = minVal
			case "coverage":
				covFactor := float64(numBins) / float64(end-start)
				values[i] = float32(covFactor * float64(validCount))
			case "sum":
				values[i] = float32(sumData)
			default:
				values[i] = float32(sumData / float64(validCount))
			}
		}
	}

	return values, nil
}

// bwGetValuesAutoZoom 自动选择合适的zoom level并获取值
func bwGetValuesAutoZoom(fp *bigWigFile_t, chrom string, start, end uint32, numBins int, summaryType string) ([]float32, error) {
	if fp.Hdr == nil || len(fp.Hdr.ZoomHdrs) == 0 || fp.Hdr.ZoomHdrs[0] == nil {
		// 没有zoom数据，使用原始数据
		return bwGetValuesFromRaw(fp, chrom, start, end, numBins, summaryType)
	}

	// 计算期望的reduction level
	baseSize := end - start
	desiredReduction := baseSize / uint32(numBins)
	if desiredReduction < 2 {
		desiredReduction = 2
	}

	// 选择最佳zoom level
	zhdr := fp.Hdr.ZoomHdrs[0]
	bestIdx := bwSelectBestZoomLevel(zhdr, desiredReduction)

	if bestIdx >= 0 {
		// 使用zoom level
		return bwGetValuesFromZoom(fp, bestIdx, chrom, start, end, numBins, summaryType)
	}

	// 如果没有合适的zoom level，使用原始数据
	return bwGetValuesFromRaw(fp, chrom, start, end, numBins, summaryType)
}

// bwGetValuesFromRaw 从原始数据获取值（无zoom）
func bwGetValuesFromRaw(fp *bigWigFile_t, chrom string, start, end uint32, numBins int, summaryType string) ([]float32, error) {
	intervals := bwGetOverlappingIntervals(fp, chrom, start, end)
	if intervals == nil || intervals.L == 0 {
		values := make([]float32, numBins)
		for i := range values {
			values[i] = float32(math.NaN())
		}
		return values, nil
	}

	values := make([]float32, numBins)
	for i := range values {
		values[i] = float32(math.NaN())
	}

	binSize := float64(end-start) / float64(numBins)

	for i := 0; i < numBins; i++ {
		binStart := start + uint32(float64(i)*binSize)
		binEnd := start + uint32(float64(i+1)*binSize)

		var sumData float64
		var count uint32
		var minVal float32 = float32(math.Inf(1))
		var maxVal float32 = float32(math.Inf(-1))

		for j := uint32(0); j < intervals.L; j++ {
			if intervals.End[j] <= binStart || intervals.Start[j] >= binEnd {
				continue
			}

			overlapStart := max32(intervals.Start[j], binStart)
			overlapEnd := min32(intervals.End[j], binEnd)
			overlap := overlapEnd - overlapStart

			if overlap > 0 {
				count += overlap
				sumData += float64(intervals.Value[j]) * float64(overlap)

				if intervals.Value[j] > maxVal {
					maxVal = intervals.Value[j]
				}
				if intervals.Value[j] < minVal {
					minVal = intervals.Value[j]
				}
			}
		}

		if count > 0 {
			switch summaryType {
			case "mean", "average":
				values[i] = float32(sumData / float64(count))
			case "max", "maximum":
				values[i] = maxVal
			case "min", "minimum":
				values[i] = minVal
			case "coverage":
				covFactor := float64(numBins) / float64(end-start)
				values[i] = float32(covFactor * float64(count))
			case "sum":
				values[i] = float32(sumData)
			default:
				values[i] = float32(sumData / float64(count))
			}
		}
	}

	return values, nil
}

// 辅助函数
func decompressZlibSimple(compBuf []byte) ([]byte, error) {
	r, err := zlib.NewReader(bytes.NewReader(compBuf))
	if err != nil {
		return nil, err
	}
	defer r.Close()

	var buf bytes.Buffer
	_, err = io.Copy(&buf, r)
	if err != nil && err != io.EOF {
		return nil, err
	}
	return buf.Bytes(), nil
}

func max32(a, b uint32) uint32 {
	if a > b {
		return a
	}
	return b
}

func min32(a, b uint32) uint32 {
	if a < b {
		return a
	}
	return b
}


