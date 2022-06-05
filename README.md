一些注意事项

1. 本地安装：
把 getLCAI_1.0.0.zip 放到当前工作目录下，运行下面的命令：
install.packages("./getLCAI_1.0.0.zip", repos = NULL, type = "win.binary")

2. 输入数据：
    如果是文件，必须是 制表符 分割的文本文件
    如果是R数据，必须是data.frame

3. 输出：
输出的是一个list，运行getlcai()时，注意将其赋值给一个变量，如：

out = getlcai(......)