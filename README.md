# MolCraft

## Description

TODO list

首先读取小分子文件（ADDMOL），计算几何中心，实现可以根据三维旋转矩阵进行旋转。

读取POSCAR，根据用户的输入范围x, y, z在分数坐标下限定填充范围，并且获取填充个数。

随机产生插入位置，以及旋转角度。然后进行尝试填充，利用范德华半径初步进行判断(VdW.ini)，需要考虑周期性。

如果尝试失败的次数过多，ERROR输出。

输出最终结构