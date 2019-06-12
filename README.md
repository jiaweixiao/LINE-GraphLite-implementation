# LINE-GraphLite-implementation
A [GraphLite][1] implementation of LINE in Tang's paper [LINE: Large-scale information network embedding][2]

GraphLite is a lightweight graph computation platform in C/C++.  
LINE is a large-scale information network embedding method.


运行方式：
    1. 测试数据为youtube_links.txt，请下载并参考[LINE论文源码][3]，将源码放入experiment文件夹  
    2. 运行Input/mapto.py脚本，将数据顶点名重新映射为id∈{0,...,n-1}，n为顶点数量，前两行分别写入顶点数量和边数量  
	3. 运行bin/hash-partitioner.pl进行数据划分  
	4. 编译并运行example/Line.cc程序，相关参数需提前修改源码中的宏定义  
	5. 运行Input/mapback.py脚本，顶点id返回顶点名  
    6. 结果测评程序见experiment文件下下的LINE论文源码     

[1]: https://github.com/schencoding/GraphLite
[2]: https://arxiv.org/abs/1503.03578
[3]: https://github.com/tangjianpku/LINE

