package gobigwig

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"unsafe"
)

const (
	GOBIGWIG_VERSION  = 0.1
	BIGWIG_MAGIC      = 0x888FFC26
	CIRTREE_MAGIC     = 0x78CA8C91
	IDX_MAGIC         = 0x2468ace0
	DEFAULT_nCHILDREN = 64
	DEFAULT_BLOCKSIZE = 32768
)


type bwStatsType struct {
	doesNotExist int
	mean         int
	average      int
	stdev        int
	dev          int
	max          int
	min          int
	cov          int
	coverage     int
	sum          int
}

func newIotabwStatsType() bwStatsType {
	return bwStatsType{
		doesNotExist: -1, /*!< This does nothing 无操作 */
		mean:         0,  /*!< The mean value  平均值*/
		average:      0,  /*!< The mean value 平均值 */
		stdev:        1,  /*!< The standard deviation of the values 标准差  */
		dev:          1,  /*!< The standard deviation of the values 标准差 */
		max:          2,  /*!< The maximum value 最大值 */
		min:          3,  /*!< The minimum value 最小值 */
		cov:          4,  /*!< The number of bases covered 覆盖碱基数 */
		coverage:     4,  /*!<The number of bases covered 覆盖碱基数 */
		sum:          5,  /*!< The sum of per-base values 每碱基值的总和 */
	}
}

//Should hide this from end users
/*!
 * @brief BigWig 文件有多个“缩放（zoom）层级”，每个层级都有自己独立的头部信息（header）。
 *        这个结构体用于保存这些缩放层级的头部。
 *
 * 注意（N.B.）：在文件的实际二进制表示中，`level` 和 `dataOffset` 之间有 4 个字节的填充（padding）。
 */
type bwZoomHdr_t struct {
	// 每个缩放层级的编号（0, 1, 2, ...）
	Level []uint32

	// 每个层级对应的数据在文件中的偏移地址（字节为单位）
	// 目前 libBigWig 实际上并没有用到这个字段
	DataOffset []uint64

	// 每个层级的索引在文件中的偏移地址（用于定位索引）
	IndexOffset []uint64

	// 每个缩放层级对应的一棵 R 树索引
	Idx []*bwRTree_t
}

/*!
 * @brief The header section of a bigWig file.
 *
 * Some of the values aren't currently used for anything. Others may optionally not exist.
 */
type bigWigHdr_t struct {
	version           uint16         /**<The version information of the file.*/
	nLevels           uint16         /**<The number of "zoom" levels.*/
	ctoffset          uint64         /**<The offset to the on-disk chromosome tree list.*/
	dataOffset        uint64         /**<The on-disk offset to the first block of data.*/
	indexoffset       uint64         /**<The on-disk offset to the data index.*/
	fieldCount        uint16         /**<Total number of fields.*/
	definedFieldCount uint16         /**<Number of fixed-format BED fields.*/
	sqloffset         uint64         /**<The on-disk offset to an SQL string. This is unused.*/
	summaryoffset     uint64         /**<If there's a summary, this is the offset to it on the disk.*/
	bufsize           uint32         /**<The compression buffer size (if the data is compressed).*/
	extensionoffset   uint64         /**<Unused*/
	ZoomHdrs          []*bwZoomHdr_t /**<Pointers to the header for each zoom level.*/
	//total Summary
	NBasesCovered uint64  /**<The total bases covered in the file.*/
	MinVal        float64 /**<The minimum value in the file.*/
	MaxVal        float64 /**<The maximum value in the file.*/
	SumData       float64 /**<The sum of all values in the file.*/
	SumSquared    float64 /**<The sum of the squared values in the file.*/
}

//Should probably replace this with a hash
/*!
 * @brief Holds the chromosomes and their lengths
 */
// chromList holds the chromosomes and their lengths
type chromList struct {
	NKeys int64    // The number of chromosomes
	Chrom []string // A list of chromosome names
	Len   []uint32 // The lengths of each chromosome
}

// bwLL is a linked list of R-tree nodes
type bwLL struct {
	Node *bwRTreeNode_t // pointer to R-tree node
	Next *bwLL          // pointer to next linked list element
}

// BwZoomBuffer represents a zoom buffer linked list
type bwZoomBuffer_t struct {
	P    unsafe.Pointer  // pointer to any type (void* equivalent)
	L    uint32          // length or some value
	M    uint32          // another value
	Next *bwZoomBuffer_t // pointer to next buffer
}


type bigWigFile_t struct {
	URL         *URL             // 一个指针，可以处理本地或远程文件（包含缓冲区）
	Hdr         *bigWigHdr_t     // 文件头信息
	Cl          *chromList       // 染色体名称列表（顺序即 ID）
	Idx         *bwRTree_t       // 整个数据集的索引
	WriteBuffer *bwWriteBuffer_t // 写入时使用的缓冲区
	IsWrite     bool             // false: 以读取模式打开，true: 以写入模式打开
	Type        int              // 0: bigWig 文件，1: bigBed 文件
}

type bwWriteBuffer_t struct {
	NBlocks          uint64           // 已写入的块数
	BlockSize        uint32           // 最大子节点数
	NEntries         uint64           // 已处理的条目数，用于第一个 contig 和计算缩放级别
	RunningWidthSum  uint64           // 第一个 contig 条目宽度的累计和，用于计算缩放级别
	Tid              uint32           // 当前处理的 TID
	Start            uint32           // 块的起始位置
	End              uint32           // 块的结束位置
	Span             uint32           // 条目的跨度（如适用）
	Step             uint32           // 步长（如适用）
	LType            uint8            // 最近添加条目的类型
	L                uint32           // 当前 p 的大小，结合类型决定包含条目数
	P                unsafe.Pointer   // 大小为 hdr.bufSize 的缓冲区
	FirstIndexNode   *bwLL            // 链表中的第一个索引节点
	CurrentIndexNode *bwLL            // 链表中的最后一个索引节点
	FirstZoomBuffer  **bwZoomBuffer_t // 叶子节点链表的第一个节点
	LastZoomBuffer   **bwZoomBuffer_t // 叶子节点链表的最后一个节点
	NNodes           *uint64          // 每个缩放级别的叶子节点数，用于确定重复级别
	CompressPsz      uint32           // 压缩缓冲区大小（uLongf 一般映射为 uint32）
	CompressP        unsafe.Pointer   // 压缩缓冲区，大小为 CompressPsz
}

type bwOverlappingIntervals_t struct {
	L     uint32    // 当前持有的区间数量
	M     uint32    // 最大可容纳的区间/值数量
	Start []uint32  // 起始位置（0-based, 半开区间）
	End   []uint32  // 结束位置（0-based, 半开区间）
	Value []float32 // 每个区间对应的值
}

// bbOverlappingEntries_t 保存区间与字符串的关联
type bbOverlappingEntries_t struct {
	L     uint32   // 当前持有的区间数量
	M     uint32   // 结构体可容纳的最大值/区间数量
	Start []uint32 // 起始位置（基于 0 的半开区间）
	End   []uint32 // 结束位置（基于 0 的半开区间）
	Str   []string // 与每个条目关联的字符串
}

// bwOverlapIterator_t 用于迭代 bigWig 或 bigBed 文件的记录
type bwOverlapIterator_t struct {
	Bw                 *bigWigFile_t             // 指向 bigWig/bigBed 文件
	Tid                uint32                    // 染色体/contig ID
	Start              uint32                    // 查询区间的起始位置
	End                uint32                    // 查询区间的结束位置
	Offset             uint64                    // 块的偏移量
	BlocksPerIteration uint32                    // 每次迭代使用的块数
	WithString         int                       // 对于 bigBed entries，是否返回字符串
	Blocks             interface{}               // 重叠的块
	Intervals          *bwOverlappingIntervals_t // 重叠的区间（或 nil）
	Entries            *bbOverlappingEntries_t   // 重叠的条目（或 nil）
	Data               interface{}               // 指向 Intervals 或 Entries，用于判断是否继续迭代
}

type bwRTreeNode_t struct {
	IsLeaf      uint8            // 是否为叶子节点
	_           [1]byte          // 1 字节填充(padding)
	NChildren   uint16           // 节点的子节点数量
	ChrIdxStart []uint32         // 每个子节点起始染色体索引列表
	BaseStart   []uint32         // 每个子节点起始位置列表
	ChrIdxEnd   []uint32         // 每个子节点结束染色体索引列表
	BaseEnd     []uint32         // 每个子节点结束位置列表
	DataOffset  []uint64         // 叶子节点: 数据块偏移; 枝节点: 子节点偏移
	Size        []uint64         // 仅叶子节点: 数据块大小
	Child       []*bwRTreeNode_t // 仅枝节点: 子节点列表
}

type bwRTree_t struct {
	BlockSize     uint32         // 每个节点最多的子节点数量
	NItems        uint64         // R 树指向的数据块总数（冗余信息）
	ChrIdxStart   uint32         // 第一个描述的染色体索引
	BaseStart     uint32         // ChrIdxStart 上的第一个有值的位置
	ChrIdxEnd     uint32         // 最后一个有条目的染色体索引
	BaseEnd       uint32         // ChrIdxEnd 上的最后一个有值的位置
	IdxSize       uint64         // 索引在文件中的偏移（名字可能误导，其实是偏移）
	NItemsPerSlot uint32         // 每个槽的条目数，总是 1
	_padding      [4]byte        // 文件中存在 4 字节填充
	RootOffset    uint64         // R 树根节点在文件中的偏移（冗余）
	Root          *bwRTreeNode_t // 指向根节点的指针
}

// bwOverlapBlock 保存与某个区间重叠的数据块信息
type bwOverlapBlock_t struct {
	N      uint64   // 重叠的数据块数量，可能为 0
	Offset []uint64 // 每个数据块在文件中的偏移
	Size   []uint64 // 每个数据块在文件中的大小（字节）
}

// BwDataHeader 表示某个数据块的头部信息
type bwDataHeader_t struct {
	Tid    uint32 // 染色体 ID
	Start  uint32 // 数据块起始位置
	End    uint32 // 数据块结束位置
	Step   uint32 // 数据值的步长
	Span   uint32 // 每个数据值的跨度
	Type   uint8  // 数据块类型：1=bedGraph，2=variable step，3=fixed step
	NItems uint16 // 数据块中值的数量
}

func roundup(v uint32) uint32 {
	v--
	v |= v >> 1
	v |= v >> 2
	v |= v >> 4
	v |= v >> 8
	v |= v >> 16
	v++
	return v
}

// Seek to a given position, always from the beginning of the file
// Return 0 on success and -1 on error
// To do, use the return code of urlSeek() in a more useful way.
func bwSetPos(fp *bigWigFile_t, pos uint64) int {
	_, err := fp.URL.Seek(int64(pos), io.SeekStart)
	if err != nil {
		return -1
	}
	return 0
}

func bwTell(fp *bigWigFile_t) uint64 {
	pos, _ := fp.URL.Seek(0, io.SeekCurrent)
	return uint64(pos)
}

func bwRead(data any, sz, nmemb int, fp *bigWigFile_t) (int, error) {
	buf := make([]byte, sz*nmemb)
	n, err := fp.URL.Read(buf)
	if err != nil {
		return 0, err
	}
	if n != len(buf) {
		return 0, io.ErrUnexpectedEOF
	}
	// 按小端解析
	return nmemb, binary.Read(bytes.NewReader(buf), binary.LittleEndian, data)
}



// gobigwig/bigwig.go
func (hdr *bigWigHdr_t) getIndexOffset() uint64{
    return hdr.indexoffset
}

// gobigwig/bigwig.go
func (hdr *bigWigHdr_t) setIndexOffset(offset uint64) {
    hdr.indexoffset = offset
}

// IsBigWig 检查文件是否为 BigWig 文件
func bwisBigWig(fname string) (bool, error) {
	url, err := Open(fname)
	if err != nil {
		return false, err
	}
	defer url.Close()
	var magic uint32
	buf := make([]byte, 4)
	n, err := url.Read(buf)
	if err != nil && err != io.EOF {
		return false, err
	}
	if n != 4 {
		return false, nil
	}
	// 按小端解析
	if err := binary.Read(bytes.NewReader(buf), binary.LittleEndian, &magic); err != nil {
		return false, err
	}
	return magic == BIGWIG_MAGIC, nil
}


func bwReadZoomHdrs(r io.Reader, nLevels uint16) (*bwZoomHdr_t, error) {
	if nLevels == 0 {
		return nil, nil
	}
	zhdr := &bwZoomHdr_t{
		Level:       make([]uint32, nLevels),
		DataOffset:  make([]uint64, nLevels),
		IndexOffset: make([]uint64, nLevels),
		Idx:         make([]*bwRTree_t, nLevels),
	}
	var padding uint32
	for i := uint16(0); i < nLevels; i++ {
		// 读取 level
		if err := binary.Read(r, binary.LittleEndian, &zhdr.Level[i]); err != nil {
			return nil, err
		}
		// 读取 padding
		if err := binary.Read(r, binary.LittleEndian, &padding); err != nil {
			return nil, err
		}
		// 读取 dataOffset
		if err := binary.Read(r, binary.LittleEndian, &zhdr.DataOffset[i]); err != nil {
			return nil, err
		}
		// 读取 indexOffset
		if err := binary.Read(r, binary.LittleEndian, &zhdr.IndexOffset[i]); err != nil {
			return nil, err
		}
	}
	return zhdr, nil
}


func bwHdrRead(bw *bigWigFile_t) error {
	if bw.IsWrite {
		return nil
	}

	bw.Hdr = &bigWigHdr_t{}

	// 读取 magic
	var magic uint32
	if err := binary.Read(bw.URL.rs, binary.LittleEndian, &magic); err != nil {
		return fmt.Errorf("[bwHdrRead] failed to read magic: %w", err)
	}
	if magic != BIGWIG_MAGIC {
		return errors.New("[bwHdrRead] invalid magic number")
	}

	// 顺序读取文件头字段
	fields := []interface{}{
		&bw.Hdr.version,
		&bw.Hdr.nLevels,
		&bw.Hdr.ctoffset,
		&bw.Hdr.dataOffset,
		&bw.Hdr.indexoffset,
		&bw.Hdr.fieldCount,
		&bw.Hdr.definedFieldCount,
		&bw.Hdr.sqloffset,
		&bw.Hdr.summaryoffset,
		&bw.Hdr.bufsize,
		&bw.Hdr.extensionoffset,
	}

	for _, f := range fields {
		if err := binary.Read(bw.URL.rs, binary.LittleEndian, f); err != nil {
			bw.Hdr = nil
			return fmt.Errorf("[bwHdrRead] failed to read header field: %w", err)
		}
	}

	// 读取 zoom headers
	if bw.Hdr.nLevels > 0 {
		zoomHdrs, err := bwReadZoomHdrs(bw.URL.rs, bw.Hdr.nLevels)
		if err != nil {
			bw.Hdr = nil
			return fmt.Errorf("[bwHdrRead] failed to read zoom headers: %w", err)
		}
		bw.Hdr.ZoomHdrs = []*bwZoomHdr_t{zoomHdrs}
	}

	// 读取 summary 信息
	if bw.Hdr.summaryoffset > 0 {
		if _, err := bw.URL.rs.Seek(int64(bw.Hdr.summaryoffset), io.SeekStart); err != nil {
			bw.Hdr = nil
			return fmt.Errorf("[bwHdrRead] failed to seek summary: %w", err)
		}

		summaryFields := []interface{}{
			&bw.Hdr.NBasesCovered,
			&bw.Hdr.MinVal,
			&bw.Hdr.MaxVal,
			&bw.Hdr.SumData,
			&bw.Hdr.SumSquared,
		}
		for _, f := range summaryFields {
			if err := binary.Read(bw.URL.rs, binary.LittleEndian, f); err != nil {
				bw.Hdr = nil
				return fmt.Errorf("[bwHdrRead] failed to read summary: %w", err)
			}
		}
	}

	// 设置压缩标志
	bw.URL.IsCompressed = bw.Hdr.bufsize > 0
	return nil
}

func getmeta_hdr(bw *bigWigFile_t) int {
	if bw == nil || bw.Hdr == nil {
		fmt.Println("BigWig header is nil")
		return -1
	}

	hdr := bw.Hdr

	fmt.Println("=== BigWig Header ===")
	fmt.Printf("version: %d\n", hdr.version)
	fmt.Printf("nLevels: %d\n", hdr.nLevels)
	// fmt.Printf("ctoffset: %d\n", hdr.ctoffset)
	// fmt.Printf("dataOffset: %d\n", hdr.dataOffset)
	// fmt.Printf("indexoffset: %d\n", hdr.indexoffset)
	fmt.Printf("fieldCount: %d\n", hdr.fieldCount)
	fmt.Printf("definedFieldCount: %d\n", hdr.definedFieldCount)
	// fmt.Printf("sqloffset: %d\n", hdr.sqloffset)
	// fmt.Printf("summaryoffset: %d\n", hdr.summaryoffset)
	fmt.Printf("bufsize: %d\n", hdr.bufsize)
	fmt.Printf("extensionoffset: %d\n", hdr.extensionoffset)

	// Zoom headers
	
	if hdr.ZoomHdrs != nil && len(hdr.ZoomHdrs) > 0 {
		zh := hdr.ZoomHdrs[0] // 取第一个 zoom header
		for i := 0; i < len(zh.Level); i++ {
			fmt.Printf("Zoom Level %d: level=%d\n",
				i, zh.Level[i])
		}
	}
	// Summary 信息
	fmt.Println("--- Summary ---")
	fmt.Printf("NBasesCovered: %d\n", hdr.NBasesCovered)
	fmt.Printf("MinVal: %f\n", hdr.MinVal)
	fmt.Printf("MaxVal: %f\n", hdr.MaxVal)
	fmt.Printf("SumData: %f\n", hdr.SumData)
	fmt.Printf("SumSquared: %f\n", hdr.SumSquared)
	return 0
}

func readChromLeaf(bw *bigWigFile_t, cl *chromList, valueSize uint32) (uint64, error) {
	var nVals uint16
	if _, err := bwRead(&nVals, 2, 1, bw); err != nil {
		return 0, fmt.Errorf("failed to read nVals: %w", err)
	}
	chrom := make([]byte, valueSize)
	for i := 0; i < int(nVals); i++ {
		// 读取染色体名
		if n, err := bw.URL.Read(chrom); err != nil || n != int(valueSize) {
			return 0, fmt.Errorf("failed to read chrom name: %w", err)
		}
		// 读取索引
		var idx uint32
		if _, err := bwRead(&idx, 4, 1, bw); err != nil {
			return 0, fmt.Errorf("failed to read chrom index: %w", err)
		}
		// 读取长度
		var length uint32
		if _, err := bwRead(&length, 4, 1, bw); err != nil {
			return 0, fmt.Errorf("failed to read chrom length: %w", err)
		}
		// 存入 chromList
		name := string(bytes.Trim(chrom, "\x00"))
		if int(idx) >= len(cl.Chrom) {
			// 如果 chromList 还不够大，扩容
			newSize := int(idx) + 1
			cl.Chrom = append(cl.Chrom, make([]string, newSize-len(cl.Chrom))...)
			cl.Len = append(cl.Len, make([]uint32, newSize-len(cl.Len))...)
		}
		cl.Chrom[idx] = name
		cl.Len[idx] = length
	}
	return uint64(nVals), nil
}

func readChromNonLeaf(bw *bigWigFile_t, cl *chromList, keySize uint32) (uint64, error) {
	var nVals uint16
	if _, err := bwRead(&nVals, 2, 1, bw); err != nil {
		return 0, fmt.Errorf("failed to read nVals: %w", err)
	}

	rv := uint64(0)
	currentPos := bwTell(bw)
	previous := currentPos + uint64(keySize)

	for i := 0; i < int(nVals); i++ {
		if bwSetPos(bw, previous) != 0 {
			return 0, fmt.Errorf("failed to seek to previous position %d", previous)
		}

		var offset uint64
		if _, err := bwRead(&offset, 8, 1, bw); err != nil {
			return 0, fmt.Errorf("failed to read offset: %w", err)
		}

		if bwSetPos(bw, offset) != 0 {
			return 0, fmt.Errorf("failed to seek to child offset %d", offset)
		}

		// 递归读取下级节点
		nRead, err := readChromBlock(bw, cl, keySize)
		if err != nil {
			return 0, fmt.Errorf("readChromBlock failed: %w", err)
		}
		rv += nRead

		previous += uint64(8 + keySize)
	}

	return rv, nil
}

func readChromBlock(bw *bigWigFile_t, cl *chromList, keySize uint32) (uint64, error) {
	var isLeaf, padding uint8

	// 读取 isLeaf 字节
	if _, err := bwRead(&isLeaf, 1, 1, bw); err != nil {
		return 0, fmt.Errorf("failed to read isLeaf: %w", err)
	}

	// 读取 padding 字节
	if _, err := bwRead(&padding, 1, 1, bw); err != nil {
		return 0, fmt.Errorf("failed to read padding: %w", err)
	}

	// 判断叶子节点
	if isLeaf != 0 {
		return readChromLeaf(bw, cl, keySize)
	}
	return readChromNonLeaf(bw, cl, keySize)
}

func ShowChromosomes(bw *bigWigFile_t) error {
	if bw == nil || bw.Cl == nil {
		return fmt.Errorf("invalid BigWig file or chromosome list is nil")
	}

	fmt.Println("=== Chromosome Information ===")
	for i := 0; i < len(bw.Cl.Chrom); i++ {
		chrom := bw.Cl.Chrom[i]
		length := uint32(0)
		if i < len(bw.Cl.Len) {
			length = bw.Cl.Len[i]
		}
		fmt.Printf("Chrom %d: %-10s Length: %d\n", i, chrom, length)
	}
	fmt.Printf("Total Chromosomes: %d\n", len(bw.Cl.Chrom))
	return nil
}

func bwReadchromList(bw *bigWigFile_t) (*chromList, error) {
	if bw == nil || bw.Hdr == nil {
		return nil, errors.New("invalid BigWig file")
	}
	if bw.IsWrite {
		return nil, errors.New("file opened in write mode")
	}

	// 定位到 chrom tree 偏移位置
	if out := bwSetPos(bw,uint64(bw.Hdr.ctoffset)); out != 0 {
		return nil, errors.New("chromList_setpos error")
	}

	cl := &chromList{}

	var magic, keySize, valueSize, itemsPerBlock uint32
	var itemCount uint64

	// 读取魔数
	if err := binary.Read(bw.URL.rs, binary.LittleEndian, &magic); err != nil {
		return nil, err
	}
	if magic != CIRTREE_MAGIC {
		return nil, errors.New("invalid CIRTREE_MAGIC")
	}

	// 顺序读取字段
	if err := binary.Read(bw.URL.rs, binary.LittleEndian, &itemsPerBlock); err != nil {
		return nil, err
	}
	if err := binary.Read(bw.URL.rs, binary.LittleEndian, &keySize); err != nil {
		return nil, err
	}
	if err := binary.Read(bw.URL.rs, binary.LittleEndian, &valueSize); err != nil {
		return nil, err
	}
	if err := binary.Read(bw.URL.rs, binary.LittleEndian, &itemCount); err != nil {
		return nil, err
	}

	cl.NKeys = int64(itemCount)
	cl.Chrom = make([]string, itemCount)
	cl.Len = make([]uint32, itemCount)

	// 跳过两个 magic（占位）
	if err := binary.Read(bw.URL.rs, binary.LittleEndian, &magic); err != nil {
		return nil, err
	}
	if err := binary.Read(bw.URL.rs, binary.LittleEndian, &magic); err != nil {
		return nil, err
	}

	// 读取染色体树块
	rv, err := readChromBlock(bw, cl, keySize)
	if err != nil {
		return nil, err
	}
	if rv != itemCount {
		return nil, errors.New("chromosome count mismatch")
	}

	return cl, nil
}

