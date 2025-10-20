# 暴露核心类，简化用户导入
from ._binding import BigWigReader

# 库元信息（后续发布时可更新）
__version__ = "0.1.0"
__author__ = "lihua"
__description__ = "Windows 平台高性能 BigWig 文件读取库（基于 Go 实现）"