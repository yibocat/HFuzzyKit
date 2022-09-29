# HFuzzyKit

### 介绍

`HFuzzyKit` 是一个在开发中的模糊逻辑工具包，包含了基本的模糊处理功能。

`HFuzzyKit` 包含了基本模糊数和模糊元素逻辑，目前实现了对偶犹豫模糊的大部分模糊逻辑功能。模糊和区间值模糊仅包含模糊数的定义和 size 运算。

基本逻辑运算基于阿基米德 t-norms ，包括基于 Algebraic t-norms 和 Einstein t-norms 的模糊逻辑基本运算

具体包括模糊功能见底部。

### 文档

正在制作中......

### 依赖的包

Numpy https://numpy.org/

Matplotlib https://matplotlib.org/

### 导入(仅 DHFES)

导入对偶犹豫模糊元素包

```python
from DHFES import *
```

导入对偶犹豫模糊元素生成器

```python
from DHFES.DHFEGenerator import *
```

### 示例

以对偶犹豫费马模糊元素为例，创建一个 DHFFE

```python
from DHFES import *
x1 = DHFFE([0.8,0.6],[0.5,0.2])
x1
>>
    The cardinal number of  MD is 2
    The cardinal number of NMD is 2

    DHFFE:{ MD: [0.8 0.6],
            NMD:[0.5 0.2] }

```

计算两个对偶犹豫费马模糊元素的代数范数加

```python
from DHFES import *
x1 = DHFFE([0.8,0.6],[0.5,0.2])
x2 = DHFFE([0.6,0.5],[0.3])
DHFFE_Algebraic_Plus(x1,x2)
>>
    The cardinal number of  MD is 4
    The cardinal number of NMD is 2
    
    DHFFE:{ MD: [0.8515 0.8306 0.7277 0.6797],
            NMD:[0.15 0.06] }
```

创建一个对偶犹豫模糊元素生成器，隶属函数和非隶属函数选择高斯函数

```python
from DHFES.DHFEGenerator import *
import numpy as np
from matplotlib import pyplot as plt

x = DHFEGenerator(3)  ## 创建一个对偶犹豫费马模糊元素生成器

x.set_MF('gaussmf')		## 隶属函数为高斯函数
x.set_NMF('gaussmf')	## 非隶属函数为高斯函数
x.set_MF_Num(3)			## 隶属函数个数为 3 个
x.set_NMF_Num(4)		## 非隶属函数个数为 4 个
x.set_MF_parameters([(4,1.36),(2.23,2.07),(4,3)])							## 设置隶属函数的参数
x.set_NMF_parameters([(7.31,5.47),(7.98,3.28),(5.44,1.69),(8.26,6.22)])		## 设置非隶属函数参数

x.set_Variable(0,10,100)		## 设置自变量范围: 0-1 间隔 100
x.generator_function()			## 开始生成

x.generator_DHFE(6,4)			## 生成隶属自变量为 6，非隶属自变量为 4 的 DHFFE
>>
    The cardinal number of  MD is 3
    The cardinal number of NMD is 4

    DHFFE:{ MD: [0.3391 0.1896 0.6654],
            NMD:[0.4233 0.4271 0.69   0.3769] }
```

随机生成一个隶属度或非隶属度数量在 30 个以内的 DHFFE

```python
from DHFES import *
randomDHFE(3,30)
>>
    The cardinal number of  MD is 9
    The cardinal number of NMD is 2

    DHFFE:{ MD: [0.6125 0.583  0.8336 0.3329 0.8221 0.4339 0.4671 0.7599 0.5337],
            NMD:[0.1878 0.559 ] }
```

### 联系

邮箱: yibocat@yeah.net

### 更新中......

下一步，制作一个文档手册

### 附(`HFuzzypyKit`包含详细功能)

1. 模糊数的定义和基本四则运算
   1. 直觉模糊数
   2. 毕达哥拉斯模糊数
   3. 费马模糊数
   4. 区间值直觉模糊数
   5. 区间值毕达哥拉斯模糊数
   6. 区间值费马模糊数
2. 交互模糊数定义和基本四则运算
   1. 交互直觉模糊
   2. 交互毕达哥拉斯模糊
   3. 交互费马模糊
   4. 交互区间值直觉模糊
   5. 交互区间值毕达哥拉斯模糊
   6. 交互区间值费马模糊
3. 对偶犹豫模糊元素的定义和基本四则运算以及交并运算
   1. 对偶犹豫模糊元素
   2. 对偶犹豫毕达哥拉斯模糊模糊元素
   3. 对偶犹豫费马模糊元素
4. 对偶犹豫模糊的相关系数度量，距离公式，模糊熵测度以及随机生成对偶犹豫模糊元素
   1. 对偶犹豫模糊元素相关系数 C1 
   2. 对偶犹豫模糊元素相关系数 C2 (与 C1 相比，分母采用最大值形式)
   3. 对偶犹豫模糊标准距离公式
   4. 对偶犹豫模糊广义 Hausdorff 公式
   5. 对偶犹豫模糊元素模糊熵
   6. 随机生成对偶犹豫模糊元素
5. 隶属函数与非隶属函数
   1. 基本 Sigmoid 函数
   2. 三角函数
   3. Z 函数
   4. 梯形函数
   5. S 函数
   6. 高斯函数
   7. 双高斯结合函数
   8. 广义贝尔函数
6. 对偶犹豫隶属函数生成器
7. 基于隶属函数的对偶犹豫模糊元素生成器

test