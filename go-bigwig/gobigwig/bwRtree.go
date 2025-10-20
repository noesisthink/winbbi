package gobigwig

import (
	"bytes"
	"compress/zlib"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
)

func decompressZlibDebug(compBuf []byte) ([]byte, error) {
	r, err := zlib.NewReader(bytes.NewReader(compBuf))
	if err != nil {
		return nil, err
	}
	defer r.Close()

	var buf bytes.Buffer
	_, err = io.Copy(&buf, r) // io.Copy 会尽量读完整
	if err != nil {
		return nil, err
	}

	// fmt.Printf("[DEBUG] 解压：压缩大小=%d, 解压后大小=%d\n", len(compBuf), n)
	return buf.Bytes(), nil
}


// func decompressZlib(compBuf []byte) ([]byte, error) {
// 	r, err := zlib.NewReader(bytes.NewReader(compBuf))
// 	if err != nil {
// 		return nil, err
// 	}
// 	defer r.Close()

// 	// 使用 bytes.Buffer 读取所有解压数据，而不是固定大小
// 	var buf bytes.Buffer
// 	_, err = buf.ReadFrom(r)
// 	if err != nil && err != io.EOF {
// 		return nil, err
// 	}
	
// 	return buf.Bytes(), nil
// }



func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func readRTreeIdx(fp *bigWigFile_t, offset uint64) (*bwRTree_t, error) {
	var magic uint32
	// 定位到 indexOffset 或 offset
	if offset == 0 {
		if bwSetPos(fp, fp.Hdr.indexoffset) != 0 {
			return nil, errors.New("failed to seek to index offset")
		}
	} else {
		if bwSetPos(fp, offset) != 0 {
			return nil, fmt.Errorf("failed to seek to offset %d", offset)
		}
	}

	// 读取并校验 magic number
	if n, err := bwRead(&magic, 4, 1, fp); err != nil || n != 1 {
		return nil, fmt.Errorf("failed to read magic number: %v", err)
	}
	if int(magic) != IDX_MAGIC {
		fmt.Fprintf(os.Stderr, "[readRTreeIdx] Mismatch in the magic number!\n")
		return nil, errors.New("invalid magic number")
	}

	node := &bwRTree_t{}

	// 顺序读取字段
	if _, err := bwRead(&node.BlockSize, 4, 1, fp); err != nil {
		return nil, err
	}
	if _, err := bwRead(&node.NItems, 8, 1, fp); err != nil {
		return nil, err
	}
	if _, err := bwRead(&node.ChrIdxStart, 4, 1, fp); err != nil {
		return nil, err
	}
	if _, err := bwRead(&node.BaseStart, 4, 1, fp); err != nil {
		return nil, err
	}
	if _, err := bwRead(&node.ChrIdxEnd, 4, 1, fp); err != nil {
		return nil, err
	}
	if _, err := bwRead(&node.BaseEnd, 4, 1, fp); err != nil {
		return nil, err
	}
	if _, err := bwRead(&node.IdxSize, 8, 1, fp); err != nil {
		return nil, err
	}
	if _, err := bwRead(&node.NItemsPerSlot, 4, 1, fp); err != nil {
		return nil, err
	}
	// padding
	var padding uint32
	if _, err := bwRead(&padding, 4, 1, fp); err != nil {
		return nil, err
	}

	node.RootOffset = bwTell(fp)
	return node, nil
}

// bwGetRTreeNode 读取一个 R 树节点
// 如果 offset 为 0，则读取根节点
func bwGetRTreeNode(fp *bigWigFile_t, offset uint64) (*bwRTreeNode_t, error) {
	var err error
	var padding uint8

	// 定位到节点偏移
	if offset != 0 {
		if bwSetPos(fp, offset) != 0 {
			return nil, fmt.Errorf("failed to seek to offset %d", offset)
		}
	} else {
		if bwSetPos(fp, fp.Idx.RootOffset) != 0 {
			return nil, fmt.Errorf("failed to seek to root offset")
		}
	}

	node := &bwRTreeNode_t{}

	// 读取节点类型和子节点数量
	if _, err = bwRead(&node.IsLeaf, 1, 1, fp); err != nil {
		return nil, err
	}
	if _, err = bwRead(&padding, 1, 1, fp); err != nil {
		return nil, err
	}
	if _, err = bwRead(&node.NChildren, 2, 1, fp); err != nil {
		return nil, err
	}

	n := int(node.NChildren)

	// 分配 slice
	node.ChrIdxStart = make([]uint32, n)
	node.BaseStart = make([]uint32, n)
	node.ChrIdxEnd = make([]uint32, n)
	node.BaseEnd = make([]uint32, n)
	node.DataOffset = make([]uint64, n)

	if node.IsLeaf != 0 {
		node.Size = make([]uint64, n)
	} else {
		node.Child = make([]*bwRTreeNode_t, n)
	}

	// 读取每个子节点的信息
	for i := 0; i < n; i++ {
		if _, err = bwRead(&node.ChrIdxStart[i], 4, 1, fp); err != nil {
			return nil, err
		}
		if _, err = bwRead(&node.BaseStart[i], 4, 1, fp); err != nil {
			return nil, err
		}
		if _, err = bwRead(&node.ChrIdxEnd[i], 4, 1, fp); err != nil {
			return nil, err
		}
		if _, err = bwRead(&node.BaseEnd[i], 4, 1, fp); err != nil {
			return nil, err
		}
		if _, err = bwRead(&node.DataOffset[i], 8, 1, fp); err != nil {
			return nil, err
		}

		if node.IsLeaf != 0 {
			if _, err = bwRead(&node.Size[i], 8, 1, fp); err != nil {
				return nil, err
			}
		}
	}

	return node, nil
}

// overlapsLeaf 查找叶子节点中与指定区间重叠的数据块
func overlapsLeaf(node *bwRTreeNode_t, tid, start, end uint32) *bwOverlapBlock_t {
	if node == nil || node.IsLeaf == 0 {
		return nil
	}

	o := &bwOverlapBlock_t{}

	// 1. 统计重叠数量
	for i := 0; i < int(node.NChildren); i++ {
		if tid < node.ChrIdxStart[i] || tid > node.ChrIdxEnd[i] {
			continue
		}

		if node.ChrIdxStart[i] != node.ChrIdxEnd[i] {
			if tid == node.ChrIdxStart[i] {
				if node.BaseStart[i] >= end {
					break
				}
			} else if tid == node.ChrIdxEnd[i] {
				if node.BaseEnd[i] <= start {
					continue
				}
			}
		} else {
			if node.BaseStart[i] >= end || node.BaseEnd[i] <= start {
				continue
			}
		}
		o.N++
	}

	if o.N == 0 {
		return o // 没有重叠，返回空结构
	}

	// 2. 分配切片
	o.Offset = make([]uint64, o.N)
	o.Size = make([]uint64, o.N)

	// 3. 填充 Offset 和 Size
	idx := 0
	for i := 0; i < int(node.NChildren); i++ {
		if tid < node.ChrIdxStart[i] || tid > node.ChrIdxEnd[i] {
			continue
		}

		if node.ChrIdxStart[i] != node.ChrIdxEnd[i] {
			if tid == node.ChrIdxStart[i] {
				if node.BaseStart[i] >= end {
					continue
				}
			} else if tid == node.ChrIdxEnd[i] {
				if node.BaseEnd[i] <= start {
					continue
				}
			}
		} else {
			if node.BaseStart[i] >= end || node.BaseEnd[i] <= start {
				continue
			}
		}

		o.Offset[idx] = node.DataOffset[i]
		o.Size[idx] = node.Size[i]
		idx++
		if idx >= int(o.N) {
			break
		}
	}

	if idx != int(o.N) {
		fmt.Println("[overlapsLeaf] Mismatch between number of overlaps calculated and found!")
		return nil
	}

	return o
}

// mergeOverlapBlocks 合并两个 BwOverlapBlock，返回合并后的块
// 如果 b1 或 b2 为 nil，会返回另一个
func mergeOverlapBlocks(b1, b2 *bwOverlapBlock_t) *bwOverlapBlock_t {
	if b2 == nil || b2.N == 0 {
		return b1
	}
	if b1 == nil || b1.N == 0 {
		return b2
	}

	// 合并 offset 和 size
	b1.Offset = append(b1.Offset, b2.Offset...)
	b1.Size = append(b1.Size, b2.Size...)
	b1.N += b2.N

	// Go 不需要手动释放 b2，GC 会处理
	return b1
}

// overlapsNonLeaf 在非叶子节点上查找与 [start, end) 区间重叠的数据块
func overlapsNonLeaf(fp *bigWigFile_t, node *bwRTreeNode_t, tid, start, end uint32) *bwOverlapBlock_t {
	output := &bwOverlapBlock_t{}

	for i := uint16(0); i < node.NChildren; i++ {
		// 如果染色体索引在子节点范围之外，跳过
		if tid < node.ChrIdxStart[i] || tid > node.ChrIdxEnd[i] {
			continue
		}

		// 处理跨染色体的子节点
		if node.ChrIdxStart[i] != node.ChrIdxEnd[i] {
			if tid == node.ChrIdxStart[i] && node.BaseStart[i] >= end {
				continue
			} else if tid == node.ChrIdxEnd[i] && node.BaseEnd[i] <= start {
				continue
			}
		} else {
			if end <= node.BaseStart[i] || start >= node.BaseEnd[i] {
				continue
			}
		}

		if node.Child[i] == nil {
			child, err := bwGetRTreeNode(fp, node.DataOffset[i])
			if err != nil {
				return nil
			}
			node.Child[i] = child
		}

		var nodeBlocks *bwOverlapBlock_t
		if node.Child[i].IsLeaf != 0 {
			nodeBlocks = overlapsLeaf(node.Child[i], tid, start, end)
		} else {
			nodeBlocks = overlapsNonLeaf(fp, node.Child[i], tid, start, end)
		}

		if nodeBlocks == nil {
			return nil
		}

		output = mergeOverlapBlocks(output, nodeBlocks)
		if output == nil {
			return nil
		}
	}

	return output
}

// walkRTreeNodes 遍历 R 树节点，返回重叠的数据块
// 如果发生错误返回 nil
func walkRTreeNodes(bw *bigWigFile_t, root *bwRTreeNode_t, tid, start, end uint32) *bwOverlapBlock_t {
	if root.IsLeaf != 0 {
		return overlapsLeaf(root, tid, start, end)
	}
	return overlapsNonLeaf(bw, root, tid, start, end)
}

// bwGetTid 返回染色体名对应的索引 (tid)
// 如果未找到返回 ^uint32(0)，即 0xFFFFFFFF
func bwGetTid(fp *bigWigFile_t, chrom string) uint32 {
	if chrom == "" {
		return ^uint32(0) // -1 的无符号表示
	}
	for i := uint32(0); i < uint32(fp.Cl.NKeys); i++ {
		if fp.Cl.Chrom[i] == chrom {
			return i
		}
	}
	return ^uint32(0)
}

func bwGetOverlappingBlocks(fp *bigWigFile_t, chrom string, start, end uint32) *bwOverlapBlock_t {
	tid := bwGetTid(fp, chrom)
	if tid == ^uint32(0) { // 未找到染色体
		fmt.Fprintf(os.Stderr, "[bwGetOverlappingBlocks] Non-existent contig: %s\n", chrom)
		return nil
	}
	// 如果索引尚未加载，则读取 R 树索引
	if fp.Idx == nil {
		var err error
		fp.Idx, err = readRTreeIdx(fp, fp.Hdr.indexoffset)
		if err != nil {
			return nil
		}
	}
	// 如果根节点为空，则读取根节点
	if fp.Idx.Root == nil {
		var err error
		fp.Idx.Root, err = bwGetRTreeNode(fp, 0)
		if err != nil {
			return nil
		}
	}
	// 遍历 R 树查找重叠的数据块
	return walkRTreeNodes(fp, fp.Idx.Root, tid, start, end)
}

// bwFillDataHdr 从字节切片 b 填充数据块头信息到 hdr
func bwFillDataHdr(hdr *bwDataHeader_t, b []byte) error {
	if len(b) < 24 { // 最少需要 24 字节才能包含所有字段
		return fmt.Errorf("byte slice too short: %d", len(b))
	}

	hdr.Tid = binary.LittleEndian.Uint32(b[0:4])
	hdr.Start = binary.LittleEndian.Uint32(b[4:8])
	hdr.End = binary.LittleEndian.Uint32(b[8:12])
	hdr.Step = binary.LittleEndian.Uint32(b[12:16])
	hdr.Span = binary.LittleEndian.Uint32(b[16:20])
	hdr.Type = b[20]                                  // uint8
	hdr.NItems = binary.LittleEndian.Uint16(b[22:24]) // uint16，注意字节偏移

	return nil
}

// pushIntervals 对应 C 里的 pushIntervals
func pushIntervals(o *bwOverlappingIntervals_t, start, end uint32, value float32) *bwOverlappingIntervals_t {
	if o.L+1 >= o.M {
		newM := roundup(o.L + 1)

		newStart := make([]uint32, newM)
		copy(newStart, o.Start)
		o.Start = newStart

		newEnd := make([]uint32, newM)
		copy(newEnd, o.End)
		o.End = newEnd

		newValue := make([]float32, newM)
		copy(newValue, o.Value)
		o.Value = newValue

		o.M = newM
	}

	o.Start[o.L] = start
	o.End[o.L] = end
	o.Value[o.L] = value
	o.L++

	return o
}

// pushBBIntervals 对应 C 版本
func pushBBIntervals(o *bbOverlappingEntries_t, start, end uint32, str string, withString bool) *bbOverlappingEntries_t {
	if o.L+1 >= o.M {
		newM := roundup(o.L + 1)

		newStart := make([]uint32, newM)
		copy(newStart, o.Start)
		o.Start = newStart

		newEnd := make([]uint32, newM)
		copy(newEnd, o.End)
		o.End = newEnd

		if withString {
			newStr := make([]string, newM)
			copy(newStr, o.Str)
			o.Str = newStr
		}

		o.M = newM
	}

	o.Start[o.L] = start
	o.End[o.L] = end
	if withString {
		o.Str[o.L] = bwStrdup(str)
	}
	o.L++

	return o
}


func bwStrdup(s string) string {
	return string([]byte(s)) // Go 中直接复制
}

func bwGetOverlappingIntervalsCore(fp *bigWigFile_t, o *bwOverlapBlock_t, tid, ostart, oend uint32) *bwOverlappingIntervals_t {
	if o == nil || o.N == 0 {
		// fmt.Println("[DEBUG] 没有重叠块")
		return &bwOverlappingIntervals_t{}
	}

	// fmt.Printf("[DEBUG] 处理 %d 个重叠块\n", o.N)
	output := &bwOverlappingIntervals_t{}
	compressed := fp.Hdr.bufsize > 0

	for i := uint64(0); i < o.N; i++ {
		// fmt.Printf("\n[DEBUG] === 块 %d/%d ===\n", i+1, o.N)
		// fmt.Printf("[DEBUG] 偏移: %d, 大小: %d\n", o.Offset[i], o.Size[i])

		// 定位到数据块
		if out := bwSetPos(fp, o.Offset[i]); out != 0 {
			// fmt.Fprintf(os.Stderr, "[ERROR] 定位失败\n")
			return nil
		}

		// 读取数据
		compBuf := make([]byte, o.Size[i])
		n, err := fp.URL.Read(compBuf)
		if err != nil || n != int(o.Size[i]) {
			// fmt.Fprintf(os.Stderr, "[ERROR] 读取失败: %v\n", err)
			return nil
		}

		var uncompressed []byte
		if compressed {
			uncompressed, err = decompressZlibDebug(compBuf)
			if err != nil {
				fmt.Fprintf(os.Stderr, "[ERROR] 解压失败: %v\n", err)
				return nil
			}
		} else {
			uncompressed = compBuf
		}

		if len(uncompressed) < 24 {
			// fmt.Fprintf(os.Stderr, "[ERROR] 数据太短\n")
			return nil
		}

		hdr := bwDataHeader_t{}
		if err := bwFillDataHdr(&hdr, uncompressed); err != nil {
			// fmt.Fprintf(os.Stderr, "[ERROR] 解析头失败: %v\n", err)
			return nil
		}

		// fmt.Printf("[DEBUG] 数据头:\n")
		// fmt.Printf("  Tid=%d (查询tid=%d)\n", hdr.Tid, tid)
		// fmt.Printf("  Start=%d, End=%d\n", hdr.Start, hdr.End)
		// fmt.Printf("  Type=%d (1=bedGraph, 2=variableStep, 3=fixedStep)\n", hdr.Type)
		// fmt.Printf("  Step=%d, Span=%d\n", hdr.Step, hdr.Span)
		// fmt.Printf("  NItems=%d\n", hdr.NItems)
// 		fmt.Printf("raw bytes: %v\n", uncompressed[:24])
// fmt.Printf("tid=%d start=%d end=%d step=%d span=%d type=%d nItems=%d\n",
//     binary.LittleEndian.Uint32(uncompressed[0:4]),
//     binary.LittleEndian.Uint32(uncompressed[4:8]),
//     binary.LittleEndian.Uint32(uncompressed[8:12]),
//     binary.LittleEndian.Uint32(uncompressed[12:16]),
//     binary.LittleEndian.Uint32(uncompressed[16:20]),
//     uncompressed[20],
//     binary.LittleEndian.Uint16(uncompressed[22:24]),
// )
// fmt.Printf("raw type bytes: %v\n", uncompressed[20:24])
		if hdr.Tid != tid {
			// fmt.Printf("[DEBUG] 染色体不匹配，跳过\n")
			continue
		}

		p := uncompressed[24:]
		// fmt.Printf("[DEBUG] 数据部分大小: %d 字节\n", len(p))
		// fmt.Printf("[DEBUG] 前32字节数据: % x\n", p[:min(32, len(p))])

		start := hdr.Start
		itemsAdded := 0

		for j := uint16(0); j < hdr.NItems; j++ {
			var end uint32
			var value float32

			switch hdr.Type {
			case 1: // bedGraph
				if len(p) < 12 {
					// fmt.Printf("[DEBUG] bedGraph 数据不足, 结束循环\n")
					break
				}
				start = binary.LittleEndian.Uint32(p[0:4])
				end = binary.LittleEndian.Uint32(p[4:8])
				value = math.Float32frombits(binary.LittleEndian.Uint32(p[8:12]))
				p = p[12:]

			case 2: // variableStep
				if len(p) < 8 {
					// fmt.Printf("[DEBUG] variableStep 数据不足, 结束循环\n")
					break
				}
				start = binary.LittleEndian.Uint32(p[0:4])
				end = start + hdr.Span
				value = math.Float32frombits(binary.LittleEndian.Uint32(p[4:8]))
				p = p[8:]

			case 3: // fixedStep
				if len(p) < 4 {
					// fmt.Printf("[DEBUG] fixedStep 数据不足, 结束循环\n")
					break
				}
				start += hdr.Step
				end = start + hdr.Span
				value = math.Float32frombits(binary.LittleEndian.Uint32(p[0:4]))
				p = p[4:]

			default:
				// fmt.Printf("[DEBUG] 未知类型: %d\n", hdr.Type)
				return nil
			}

			// 跳过不在查询范围的区间
			if end <= ostart || start >= oend {
				continue
			}

			output = pushIntervals(output, start, end, value)
			itemsAdded++
			// if j < 3 {
			// 	// fmt.Printf("[DEBUG]   Item %d: [%d, %d) = %.4f\n", j, start, end, value)
			// }
		}

		// fmt.Printf("[DEBUG] 本块添加了 %d 个区间到结果\n", itemsAdded)
	}

	// fmt.Printf("[DEBUG] 总共返回 %d 个区间\n", output.L)
	return output
}


func bytesToUint32Slice(b []byte) []uint32 {
    n := len(b) / 4
    u32s := make([]uint32, n)
    for i := 0; i < n; i++ {
        u32s[i] = binary.LittleEndian.Uint32(b[i*4 : (i+1)*4])
    }
    return u32s
}

func bwGetOverlappingIntervals(fp *bigWigFile_t, chrom string, start, end uint32) *bwOverlappingIntervals_t {
	tid := bwGetTid(fp, chrom)
	if tid == ^uint32(0) { // tid == -1 的情况
		return nil
	}

	blocks := bwGetOverlappingBlocks(fp, chrom, start, end)
	if blocks == nil {
		return nil
	}
	output := bwGetOverlappingIntervalsCore(fp, blocks, tid, start, end)
	return output
}

func bwOverlappingIntervalsIterator(fp *bigWigFile_t, chrom string, start, end, blocksPerIteration uint32) *bwOverlapIterator_t {
	var output *bwOverlapIterator_t
	tid := bwGetTid(fp, chrom)
	if tid == ^uint32(0) { // tid == -1
		return output
	}

	output = &bwOverlapIterator_t{
		Bw:                 fp,
		Tid:                tid,
		Start:              start,
		End:                end,
		BlocksPerIteration: blocksPerIteration,
	}

	blocks := bwGetOverlappingBlocks(fp, chrom, start, end)
	output.Blocks = blocks

	if blocks != nil {
		n := blocks.N
		if n > uint64(blocksPerIteration) {
			blocks.N = uint64(blocksPerIteration)
		}
		output.Intervals = bwGetOverlappingIntervalsCore(fp, blocks, tid, start, end)
		blocks.N = n
		output.Offset = uint64(blocksPerIteration)
	}

	output.Data = output.Intervals
	return output
}

func bwIteratorNext(iter *bwOverlapIterator_t) *bwOverlapIterator_t {
	if iter == nil || iter.Blocks == nil {
		return nil
	}
	// 类型断言
	blocks, ok := iter.Blocks.(*bwOverlapBlock_t)
	if !ok || blocks == nil {
		return nil
	}
	// 释放上一次迭代的结果
	if iter.Intervals != nil {
		iter.Intervals = nil
	}
	if iter.Entries != nil {
		iter.Entries = nil
	}
	iter.Data = nil
	if iter.Offset < blocks.N {
		// 保存原始值
		n := blocks.N
		// 计算本次迭代的块数量
		var currentN uint64
		if iter.Offset+uint64(iter.BlocksPerIteration) > n {
			currentN = n - iter.Offset
		} else {
			currentN = uint64(iter.BlocksPerIteration)
		}
		// 截取本次迭代的块
		currentBlocks := &bwOverlapBlock_t{
			N:      currentN,
			Offset: blocks.Offset[iter.Offset : iter.Offset+currentN],
			Size:   blocks.Size[iter.Offset : iter.Offset+currentN],
		}
		// 获取区间或条目
		if iter.Bw.Type == 0 {
			iter.Intervals = bwGetOverlappingIntervalsCore(iter.Bw, currentBlocks, iter.Tid, iter.Start, iter.End)
			iter.Data = iter.Intervals
		} 
		iter.Offset += uint64(iter.BlocksPerIteration)
		// 检查是否出错
		if iter.Intervals == nil && iter.Entries == nil {
			return nil
		}
	}
	return iter
}

func bwGetValues(fp *bigWigFile_t, chrom string, start, end uint32, includeNA bool) *bwOverlappingIntervals_t {
	intermediate := bwGetOverlappingIntervals(fp, chrom, start, end)
	if intermediate == nil {
		return nil
	}
	output := &bwOverlappingIntervals_t{}
	if output == nil {
		return nil
	}
	if includeNA {
		// 每个位点都返回一个值
		length := end - start
		output.Value = make([]float32, length)
		for i := range output.Value {
			output.Value[i] = float32(math.NaN())
		}
		for i := uint32(0); i < intermediate.L; i++ {
			for j := intermediate.Start[i]; j < intermediate.End[i]; j++ {
				if j < start || j >= end {
					continue
				}
				output.Value[j-start] = intermediate.Value[i]
			}
		}
		output.L = length
	} else {
		// 只返回实际有值的位置
		var n uint32 = 0
		for i := uint32(0); i < intermediate.L; i++ {
			if intermediate.Start[i] < start {
				intermediate.Start[i] = start
			}
			if intermediate.End[i] > end {
				intermediate.End[i] = end
			}
			n += intermediate.End[i] - intermediate.Start[i]
		}
		output.L = n
		output.Start = make([]uint32, n)
		output.Value = make([]float32, n)

		idx := uint32(0)
		for i := uint32(0); i < intermediate.L; i++ {
			for j := intermediate.Start[i]; j < intermediate.End[i]; j++ {
				if j < start || j >= end {
					continue
				}
				output.Start[idx] = j
				output.Value[idx] = intermediate.Value[i]
				idx++
			}
		}
	}
	return output
}

// bwReadIndex 读取指定 offset 的 RTree 索引，如果 offset 为 0，则读取值的索引
// 返回 nil 表示出错
func bwReadIndex(fp *bigWigFile_t, offset uint64) *bwRTree_t {
	idx, err := readRTreeIdx(fp, offset)
	if err != nil || idx == nil {
		return nil
	}

	// 读取根节点
	root, err := bwGetRTreeNode(fp, idx.RootOffset)
	if err != nil || root == nil {
		return nil
	}
	idx.Root = root
	return idx
}