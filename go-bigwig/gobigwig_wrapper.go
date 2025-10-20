package main

// 导入你的 gobigwig 包（触发 C 绑定接口的导出）
import _ "go-bigwig/gobigwig"

// main 函数是编译入口（空实现，因为我们只需要共享库，不需要可执行文件）
func main() {}