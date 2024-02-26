### 记录

![alt text](image.png)

距离约束和体积约束，这里需要把模型处理成四面体的格式（参考后面教程）

![alt](image-1.png)

xpbd约束求解方法，相比与pbd，在增加了红色alpha项，使得调参更方便。alpha越大表示越硬。其中采用的是导数形式，所以alpha为0为最大。

![alt text](image-2.png)

体积约束求解方法

> 代码中注意substep步骤设置带一些，不然计算紊乱出问题。




