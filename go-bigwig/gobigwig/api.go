package gobigwig

/*
#cgo CFLAGS: -I.
#cgo LDFLAGS: -lm
#include <stdint.h>
#include <stdlib.h>
#include <windows.h>  // 引入Windows头文件，使用Windows原生类型

// Windows兼容：用Windows原生类型替代C99类型（避免uint64_t未定义）
struct CBWFileInfo {
    WORD    Version;           // 对应uint16_t（2字节）
    WORD    NLevels;           // 对应uint16_t
    WORD    FieldCount;        // 对应uint16_t
    WORD    DefinedFieldCount; // 对应uint16_t
    DWORD   Bufsize;           // 对应uint32_t（4字节）
    ULONGLONG Extensionoffset; // 对应uint64_t（8字节，Windows原生）
    ULONGLONG NBasesCovered;   // 对应uint64_t
    double  MinVal;            // 保持不变（8字节）
    double  MaxVal;
    double  SumData;
    double  SumSquared;
};
*/
import "C"

import (
	"fmt"
	"math"
	"runtime"
	"sync"
	"unsafe"
)

// -------------------------- 你原有结构体（保持不变） --------------------------
type Bigwig_file_out struct {
	bf_fp *bigWigFile_t
	Info  FileInfo_bw_out
}

type FileInfo_bw_out struct {
	Version           uint16
	NLevels           uint16
	FieldCount        uint16
	DefinedFieldCount uint16
	Bufsize           uint32
	Extensionoffset   uint64
	NBasesCovered     uint64
	MinVal            float64
	MaxVal            float64
	SumData           float64
	SumSquared        float64
}

// -------------------------- 你原有核心方法（仅修正1行错误） --------------------------
func OpenBigWig(fname string) (*Bigwig_file_out, error) {
	// 1. 检查是否是 BigWig 文件
	isBw, err := bwisBigWig(fname)
	if err != nil {
		return nil, fmt.Errorf("检查文件格式失败: %w", err)
	}
	if !isBw {
		return nil, fmt.Errorf("不是有效的 BigWig 文件")
	}
	// 2. 打开文件
	url, err := Open(fname)
	if err != nil {
		return nil, fmt.Errorf("打开文件失败: %w", err)
	}
	fp := &bigWigFile_t{
		URL:     url,
		IsWrite: false,
		Type:    0, // 0 = BigWig
	}
	// 3. 读取文件头
	if err := bwHdrRead(fp); err != nil {
		url.Close()
		return nil, fmt.Errorf("读取文件头失败: %w", err)
	}
	// 4. 读取染色体列表
	cl, err := bwReadchromList(fp)
	if err != nil {
		url.Close()
		return nil, fmt.Errorf("读取染色体列表失败: %w", err)
	}
	fp.Cl = cl
	// 5. 读取索引
	idx := bwReadIndex(fp, 0)
	if idx == nil {
		url.Close()
		return nil, fmt.Errorf("读取索引失败")
	}
	fp.Idx = idx

	fbo := FileInfo_bw_out{
		Version:           fp.Hdr.version,
		NLevels:           fp.Hdr.nLevels,
		FieldCount:        fp.Hdr.fieldCount,
		DefinedFieldCount: fp.Hdr.definedFieldCount,
		Bufsize:           fp.Hdr.bufsize,
		Extensionoffset:   fp.Hdr.extensionoffset,
		NBasesCovered:     fp.Hdr.NBasesCovered,
		MinVal:            fp.Hdr.MinVal,
		MaxVal:            fp.Hdr.MaxVal,
		SumData:           fp.Hdr.SumData,
		SumSquared:        fp.Hdr.SumSquared,
	}

	return &Bigwig_file_out{
		bf_fp: fp,
		Info:  fbo,
	}, nil
}

func CloseBigWig(fp *Bigwig_file_out) {
	if fp.bf_fp != nil && fp.bf_fp.URL != nil {
		fp.bf_fp.URL.Close()
	}
}

func (fp *Bigwig_file_out) ReadBigWigSignal(chrom string, start int, end int) []float32 {
	start_uint32 := uint32(start)
	end_uint32 := uint32(end)
	blocksPerIteration := uint32(10) // 每次处理10个块
	iter := bwOverlappingIntervalsIterator(fp.bf_fp, chrom, start_uint32, end_uint32, blocksPerIteration)
	if iter == nil {
		fmt.Println("创建迭代器失败")
		return nil
	}
	output_float32 := []float32{}
	// 迭代所有数据块
	for iter.Data != nil {
		intervals := iter.Intervals
		if intervals != nil {
			output_float32 = append(output_float32, intervals.Value...)
		}
		iter = bwIteratorNext(iter)
	}
	return output_float32
}

func (fp *Bigwig_file_out) Getmeta_hdr() {
	fmt.Println("\n--- 文件头信息 ---")
	getmeta_hdr(fp.bf_fp)
	fmt.Println("\n--- 染色体信息 ---")
}

func (fp *Bigwig_file_out) GetVersion() uint16          { return fp.Info.Version }
func (fp *Bigwig_file_out) GetNLevels() uint16          { return fp.Info.NLevels }
func (fp *Bigwig_file_out) GetFieldCount() uint16       { return fp.Info.FieldCount }
func (fp *Bigwig_file_out) GetDefinedFieldCount() uint16 { return fp.Info.DefinedFieldCount }
func (fp *Bigwig_file_out) GetBufsize() uint32          { return fp.Info.Bufsize }
func (fp *Bigwig_file_out) GetExtensionOffset() uint64  { return fp.Info.Extensionoffset }
func (fp *Bigwig_file_out) GetNBasesCovered() uint64    { return fp.Info.NBasesCovered }
func (fp *Bigwig_file_out) GetMinVal() float64          { return fp.Info.MinVal }
func (fp *Bigwig_file_out) GetMaxVal() float64          { return fp.Info.MaxVal }
func (fp *Bigwig_file_out) GetSumData() float64         { return fp.Info.SumData }
func (fp *Bigwig_file_out) GetSumSquared() float64      { return fp.Info.SumSquared }

func (fp *Bigwig_file_out) PrintZoomInfo() {
	if fp.bf_fp.Hdr == nil || len(fp.bf_fp.Hdr.ZoomHdrs) == 0 {
		fmt.Println("No zoom levels available")
		return
	}

	zhdr := fp.bf_fp.Hdr.ZoomHdrs[0]
	fmt.Println("=== Zoom Levels ===")
	for i, level := range zhdr.Level {
		fmt.Printf("Level %d: reduction=%d, indexOffset=%d, dataOffset=%d\n",
			i, level, zhdr.IndexOffset[i], zhdr.DataOffset[i])
	}
}

// ZoomSelector 定义了选择 zoom level 的函数类型
type ZoomSelector func(zhdr *bwZoomHdr_t, desiredReduction uint32) int

// BWOptions_Zoom 表示 zoom 层级选择和取值的参数
type BWOptions_Zoom struct {
	NumBins        int          // 输出分辨率（输出多少个bin）
	SummaryType    string       // 求值方式，如 "mean" / "max"
	IndexZoomModel ZoomSelector // zoom选择策略函数
}

func (fp *Bigwig_file_out) GetZoomValues(
	chrom string,
	start int,
	end int,
	numBins int,
	useClosest bool,
	desiredReduction int,
) []float32 {

	opts := BWOptions_Zoom{
		NumBins:     numBins,
		SummaryType: "mean",
	}

	if useClosest {
		opts.IndexZoomModel = bwGetBestZoomClosest
	} else {
		opts.IndexZoomModel = bwSelectBestZoomLevel
	}

	if len(fp.bf_fp.Hdr.ZoomHdrs) == 0 {
		fmt.Println("no zoom headers available")
		return nil
	}

	zhdr := fp.bf_fp.Hdr.ZoomHdrs[0]
	// 核心修正：删除 &zhdr 中的 &，直接传入 zhdr（单层指针）
	zoomIdx := opts.IndexZoomModel(zhdr, uint32(desiredReduction))
	if zoomIdx < 0 {
		fmt.Printf("no suitable zoom level found for desiredReduction=%d\n", desiredReduction)
		return nil
	}

	values, err := bwGetValuesFromZoom(
		fp.bf_fp, zoomIdx, chrom,
		uint32(start), uint32(end),
		opts.NumBins, opts.SummaryType,
	)
	if err != nil {
		fmt.Printf("failed to read zoom data: %v\n", err)
		return nil
	}

	// 并行替换 NaN 为 0
	n := len(values)
	if n == 0 {
		return values
	}

	numCPU := runtime.NumCPU()
	chunkSize := (n + numCPU - 1) / numCPU
	var wg sync.WaitGroup

	for i := 0; i < numCPU; i++ {
		startIdx := i * chunkSize
		endIdx := startIdx + chunkSize
		if endIdx > n {
			endIdx = n
		}

		if startIdx >= n {
			break
		}

		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()
			for j := s; j < e; j++ {
				if math.IsNaN(float64(values[j])) {
					values[j] = 0
				}
			}
		}(startIdx, endIdx)
	}

	wg.Wait()
	return values
}

// -------------------------- C绑定接口（Windows兼容，无变化） --------------------------

// 1. 打开文件（返回句柄，失败返回0）
//export BigWigOpen
func BigWigOpen(fname *C.char) C.uintptr_t {
	if fname == nil {
		fmt.Println("BigWigOpen: 文件名不能为空")
		return 0
	}
	goFname := C.GoString(fname)
	fp, err := OpenBigWig(goFname)
	if err != nil {
		fmt.Printf("BigWigOpen: 打开失败: %v\n", err)
		return 0
	}
	// Go指针转句柄返回给Python
	return C.uintptr_t(uintptr(unsafe.Pointer(fp)))
}

// 2. 关闭文件（释放资源）
//export BigWigClose
func BigWigClose(handle C.uintptr_t) {
	if handle == 0 {
		return
	}
	fp := (*Bigwig_file_out)(unsafe.Pointer(uintptr(handle)))
	CloseBigWig(fp)
}

// 3. 读取原始信号（返回float数组指针+长度，Python侧需释放内存）
//export BigWigReadSignal
func BigWigReadSignal(
	handle C.uintptr_t,
	chrom *C.char,
	start C.int,
	end C.int,
	outLen *C.int,
) *C.float {
	// 参数校验
	if handle == 0 || chrom == nil || outLen == nil || start < 0 || end <= start {
		if outLen != nil {
			*outLen = 0
		}
		return nil
	}

	// 类型转换
	fp := (*Bigwig_file_out)(unsafe.Pointer(uintptr(handle)))
	goChrom := C.GoString(chrom)
	goVals := fp.ReadBigWigSignal(goChrom, int(start), int(end))

	// 处理返回结果
	if goVals == nil || len(goVals) == 0 {
		*outLen = 0
		return nil
	}
	*outLen = C.int(len(goVals))

	// 分配C内存并拷贝数据（Python侧需调用BigWigFree释放）
	cVals := (*C.float)(C.malloc(C.size_t(len(goVals)) * C.sizeof_float))
	if cVals == nil {
		fmt.Println("BigWigReadSignal: 内存分配失败")
		*outLen = 0
		return nil
	}
	goBuf := (*[1 << 30]C.float)(unsafe.Pointer(cVals))[:len(goVals):len(goVals)]
	for i, v := range goVals {
		goBuf[i] = C.float(v)
	}

	return cVals
}

// 4. 获取Zoom缩放数据（并行NaN转0）
//export BigWigGetZoomValues
func BigWigGetZoomValues(
	handle C.uintptr_t,
	chrom *C.char,
	start C.int,
	end C.int,
	numBins C.int,
	useClosest C.int,
	desiredReduction C.int,
	outLen *C.int,
) *C.float {
	// 参数校验
	if handle == 0 || chrom == nil || outLen == nil || start < 0 || end <= start || numBins <= 0 || desiredReduction < 1 {
		if outLen != nil {
			*outLen = 0
		}
		return nil
	}

	// 类型转换
	fp := (*Bigwig_file_out)(unsafe.Pointer(uintptr(handle)))
	goChrom := C.GoString(chrom)
	goVals := fp.GetZoomValues(
		goChrom, int(start), int(end),
		int(numBins), useClosest != 0, int(desiredReduction),
	)

	// 处理返回结果
	if goVals == nil || len(goVals) == 0 {
		*outLen = 0
		return nil
	}
	*outLen = C.int(len(goVals))

	// 分配C内存并拷贝数据
	cVals := (*C.float)(C.malloc(C.size_t(len(goVals)) * C.sizeof_float))
	if cVals == nil {
		fmt.Println("BigWigGetZoomValues: 内存分配失败")
		*outLen = 0
		return nil
	}
	goBuf := (*[1 << 30]C.float)(unsafe.Pointer(cVals))[:len(goVals):len(goVals)]
	for i, v := range goVals {
		goBuf[i] = C.float(v)
	}

	return cVals
}

// 5. 获取文件元信息（Windows兼容：使用Windows原生结构体类型）
//export BigWigGetInfo
func BigWigGetInfo(handle C.uintptr_t, info *C.struct_CBWFileInfo) C.int {
	if handle == 0 || info == nil {
		return -1 // 失败返回-1
	}
	// 句柄转回Go结构体
	fp := (*Bigwig_file_out)(unsafe.Pointer(uintptr(handle)))
	goInfo := fp.Info

	// 赋值到Windows兼容的C结构体（类型对应：Go → Windows C类型）
	info.Version = C.WORD(goInfo.Version)
	info.NLevels = C.WORD(goInfo.NLevels)
	info.FieldCount = C.WORD(goInfo.FieldCount)
	info.DefinedFieldCount = C.WORD(goInfo.DefinedFieldCount)
	info.Bufsize = C.DWORD(goInfo.Bufsize)
	info.Extensionoffset = C.ULONGLONG(goInfo.Extensionoffset)
	info.NBasesCovered = C.ULONGLONG(goInfo.NBasesCovered)
	info.MinVal = C.double(goInfo.MinVal)
	info.MaxVal = C.double(goInfo.MaxVal)
	info.SumData = C.double(goInfo.SumData)
	info.SumSquared = C.double(goInfo.SumSquared)

	return 0 // 成功返回0
}

// 6. 释放C内存（Python侧必须调用，避免内存泄漏）
//export BigWigFree
func BigWigFree(ptr unsafe.Pointer) {
	if ptr != nil {
		C.free(ptr)
	}
}