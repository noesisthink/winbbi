import ctypes
import os
import platform

class BigWigReader:
    def __init__(self):
        """初始化：自动识别 Windows 系统，加载 winbbi.dll"""
        self.lib = self._load_shared_library()
        self._bind_types()
        self.file_handle = 0

    def _load_shared_library(self):
        """加载 Windows 平台的 winbbi.dll（仅支持 Windows 64位）"""
        # 仅支持 Windows 系统（库名 winbbi 聚焦 Windows 平台）
        if platform.system() != "Windows":
            raise OSError("winbbi 仅支持 Windows 64位系统")
        if platform.machine() != "AMD64":
            raise OSError("winbbi 仅支持 64位 Windows 系统")

        # 构建 DLL 路径（包内内置 DLL）
        lib_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "lib",
            "windows",
            "amd64"
        )
        dll_path = os.path.join(lib_dir, "winbbi.dll")

        # 检查 DLL 是否存在
        if not os.path.exists(dll_path):
            raise FileNotFoundError(
                f"未找到 winbbi.dll: {dll_path}\n"
                "请确认库文件已正确打包，或重新安装 winbbi"
            )

        # 加载 DLL
        try:
            return ctypes.CDLL(dll_path)
        except Exception as e:
            raise RuntimeError(f"加载 winbbi.dll 失败: {e}\n可能是缺少 MinGW 依赖（libgcc_s_seh-1.dll）")

    def _bind_types(self):
        """绑定 C 接口类型，避免内存错误"""
        # 定义元信息结构体（与 Go 侧 C 结构体一致）
        class CBWFileInfo(ctypes.Structure):
            _fields_ = [
                ("Version", ctypes.c_ushort),
                ("NLevels", ctypes.c_ushort),
                ("FieldCount", ctypes.c_ushort),
                ("DefinedFieldCount", ctypes.c_ushort),
                ("Bufsize", ctypes.c_uint32),
                ("Extensionoffset", ctypes.c_uint64),
                ("NBasesCovered", ctypes.c_uint64),
                ("MinVal", ctypes.c_double),
                ("MaxVal", ctypes.c_double),
                ("SumData", ctypes.c_double),
                ("SumSquared", ctypes.c_double),
            ]
        self.CBWFileInfo = CBWFileInfo

        # 1. BigWigOpen: 打开文件
        self.lib.BigWigOpen.argtypes = (ctypes.c_char_p,)
        self.lib.BigWigOpen.restype = ctypes.c_uint64

        # 2. BigWigClose: 关闭文件
        self.lib.BigWigClose.argtypes = (ctypes.c_uint64,)
        self.lib.BigWigClose.restype = None

        # 3. BigWigReadSignal: 读取原始信号
        self.lib.BigWigReadSignal.argtypes = (
            ctypes.c_uint64, ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_int)
        )
        self.lib.BigWigReadSignal.restype = ctypes.POINTER(ctypes.c_float)

        # 4. BigWigGetZoomValues: 读取缩放信号
        self.lib.BigWigGetZoomValues.argtypes = (
            ctypes.c_uint64, ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)
        )
        self.lib.BigWigGetZoomValues.restype = ctypes.POINTER(ctypes.c_float)

        # 5. BigWigGetInfo: 获取元信息
        self.lib.BigWigGetInfo.argtypes = (ctypes.c_uint64, ctypes.POINTER(CBWFileInfo))
        self.lib.BigWigGetInfo.restype = ctypes.c_int

        # 6. BigWigFree: 释放 C 内存
        self.lib.BigWigFree.argtypes = (ctypes.c_void_p,)
        self.lib.BigWigFree.restype = None

    def open(self, bw_file_path):
        """打开 BigWig 文件
        Args:
            bw_file_path: str - BigWig 文件路径（绝对路径/相对路径，支持中文）
        Raises:
            RuntimeError: 已打开文件或打开失败
            FileNotFoundError: 文件不存在
        """
        if self.file_handle != 0:
            raise RuntimeError("已打开一个文件，请先调用 close() 关闭")
        if not os.path.exists(bw_file_path):
            raise FileNotFoundError(f"文件不存在: {bw_file_path}")

        # 处理路径编码（支持中文路径）
        abs_path = os.path.abspath(bw_file_path)
        c_path = ctypes.c_char_p(abs_path.encode("utf-8"))
        self.file_handle = self.lib.BigWigOpen(c_path)

        if self.file_handle == 0:
            raise RuntimeError(f"打开 BigWig 文件失败: {abs_path}\n可能是文件损坏或非标准 BigWig 格式")

    def close(self):
        """关闭文件，释放资源（必须调用）"""
        if self.file_handle != 0:
            self.lib.BigWigClose(self.file_handle)
            self.file_handle = 0

    def get_metadata(self):
        """获取文件元信息
        Returns:
            dict - 包含版本、缩放层级、信号统计等关键信息
        """
        if self.file_handle == 0:
            raise RuntimeError("文件未打开，请先调用 open()")
        
        info = self.CBWFileInfo()
        ret = self.lib.BigWigGetInfo(self.file_handle, ctypes.byref(info))
        if ret != 0:
            raise RuntimeError("获取文件元信息失败")

        return {
            "版本": info.Version,
            "缩放层级数": info.NLevels,
            "覆盖碱基数": info.NBasesCovered,
            "最小信号值": info.MinVal,
            "最大信号值": info.MaxVal,
            "信号总和": info.SumData,
            "信号平方和": info.SumSquared,
            "缓冲区大小": info.Bufsize,
        }

    def read_raw_signal(self, chrom, start, end):
        """读取指定区间的原始信号（稀疏存储，无信号返回空列表）
        Args:
            chrom: str - 染色体名称（如 "chr1"，大小写敏感）
            start: int - 起始位置（基因组数据通常从 1 开始）
            end: int - 结束位置（必须大于 start）
        Returns:
            list[float] - 原始信号值列表
        """
        if self.file_handle == 0:
            raise RuntimeError("文件未打开，请先调用 open()")
        if start < 0 or end <= start:
            raise ValueError(f"无效区间: start={start}（≥0）, end={end}（>start）")

        c_chrom = ctypes.c_char_p(chrom.encode("utf-8"))
        out_len = ctypes.c_int(0)
        data_ptr = self.lib.BigWigReadSignal(
            self.file_handle, c_chrom, start, end, ctypes.byref(out_len)
        )

        try:
            if out_len.value <= 0 or data_ptr is None:
                return []  # 无信号区间返回空列表
            return [data_ptr[i] for i in range(out_len.value)]
        finally:
            # 必须释放 C 内存，避免泄漏
            if data_ptr is not None:
                self.lib.BigWigFree(data_ptr)

    def read_zoom_signal(self, chrom, start, end, num_bins=100, use_closest=True, desired_reduction=100):
        """读取缩放后的信号（并行 NaN 转 0，适合快速可视化）
        Args:
            chrom: str - 染色体名称
            start: int - 起始位置
            end: int - 结束位置
            num_bins: int - 输出分辨率（bin 数量，默认 100）
            use_closest: bool - 是否选择最接近的 zoom 层级（默认 True）
            desired_reduction: int - 目标缩放比（默认 100）
        Returns:
            list[float] - 缩放后的信号值列表（长度 = num_bins）
        """
        if self.file_handle == 0:
            raise RuntimeError("文件未打开，请先调用 open()")
        if start < 0 or end <= start or num_bins <= 0 or desired_reduction < 1:
            raise ValueError(
                f"无效参数: start≥0, end>start, num_bins≥1, desired_reduction≥1\n"
                f"当前参数: start={start}, end={end}, num_bins={num_bins}, desired_reduction={desired_reduction}"
            )

        c_chrom = ctypes.c_char_p(chrom.encode("utf-8"))
        out_len = ctypes.c_int(0)
        data_ptr = self.lib.BigWigGetZoomValues(
            self.file_handle, c_chrom, start, end,
            num_bins, 1 if use_closest else 0, desired_reduction,
            ctypes.byref(out_len)
        )

        try:
            if out_len.value <= 0 or data_ptr is None:
                return [0.0] * num_bins  # 无数据时返回全 0，避免索引错误
            return [data_ptr[i] for i in range(out_len.value)]
        finally:
            if data_ptr is not None:
                self.lib.BigWigFree(data_ptr)

    def __del__(self):
        """析构函数：自动关闭文件（防止用户忘记调用 close()）"""
        self.close()

    def __enter__(self):
        """支持 with 语句（推荐用法）"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """with 语句结束时自动关闭文件"""
        self.close()