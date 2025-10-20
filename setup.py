from setuptools import setup, find_packages
import os

# 读取 README.md 作为长描述
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# 定义打包配置
setup(
    name="winbbi",  # PyPI 库名（必须唯一，已改为 winbbi）
    version="0.1.0",
    author="lihua",
    author_email="3209414882@qq.com",
    description="Windows 平台高性能 BigWig 文件读取库（基于 Go 实现，支持原始信号和缩放信号）",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/noesisthink/winbbi",  # 可选，推荐填写
    packages=find_packages(),  # 自动识别 winbbi 包
    package_data={
        # 打包时包含内置的 DLL 文件（关键！）
        "winbbi": [
            "lib/windows/amd64/winbbi.dll"
        ]
    },
    classifiers=[
        # 项目分类（帮助用户在 PyPI 搜索到你的库）
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Operating System :: Microsoft :: Windows :: Windows 11",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    python_requires=">=3.7",  # 支持 Python 3.7+
    install_requires=[],  # 无第三方依赖（仅用标准库）
    keywords="bigwig, windows, bioinformatics, genomics, go,pythonbigwig,winbbi,winbigwig",  # 搜索关键词
)