import numpy as np
from matplotlib import pyplot as plt

class CustomMemshipFC(object):
    '''
        自建隶属函数生成器
        ========================================================================================
        说明：
            MemshipFC 隶属函数生成器包含了 8 种基本的隶属函数，而该生成器可以创建自定义隶属函数和非隶属函数
            注意：自建函数时，使用如下方法
            ---------------------------------
            |   def func_test(x,*p):        |
            |       return p[0]*x + p[1]    |
            ---------------------------------
            参数以 list 类型传入，表示为参数列表，然后在函数中使用列表元素的形式表示参数
        ========================================================================================
        属性:
            ArbFunc: 表示自定义隶属函数方法，类型为‘function’
            parameter: 表示函数的参数列表
                parameter=[[2,4,5],[1,4,7],[2,6,3]] 表示为：有三个隶属函数，且参数列表分别为 [2,4,5],[1,4,7],[2,6,3] 的隶属函数
            numFunc: 表示自定义隶属函数的个数，例如对于一个犹豫模糊集，定义了 3 个隶属函数，那么该值为 3

        方法:
            setvariable(self,start,end,linspace):  设置自变量范围和间隔
            _min_generateMF(self):                 计算 self 设置好的自变量范围和参数下的范围内函数最小值，目的是被隶属度函数减去以标准化隶属度函数
            _generateMF(self,x):                   生成隶属度集合
            MF_Plot(self):                         画出设置的隶属度函数曲线，方便观察
            calculate_MD(self,x):                  计算 x 的隶属度，返回一个长度为 numFunc 的 array 类型的数组，即可能隶属度或非隶属度
    '''


    parameter = []
    numFunc = 0
    ArbFunc = None                       ## 自建函数

    _variable_start = 0
    _variable_end = 1
    _linspace = 100

    def __init__(self,ArbFunc,parameter,numFunc):
        assert len(parameter) == numFunc,\
        'ERROR!The number of parameters is not equal to numFunc! numFunc=%d means you should create %d sets of parameters'%(numFunc,numFunc)
        assert hasattr(ArbFunc,'__call__') or ArbFunc == None, 'ERROR: The membership function has to be a function!'
        for para in parameter:
            assert type(para) == list, 'ERROR:The self-bulit function\'s parameter type is error! Parameter type should be list!'

        self.ArbFunc = ArbFunc
        self.parameter = np.asarray(parameter)
        self.numFunc = numFunc

    def __repr__(self):
        """
            Print information of membership function.
        """
        s = ''
        for i in range(self.numFunc):
            s+='Function %d: %s, parameter:%s \n'%(i+1,self.ArbFunc.__name__,self.parameter[i])
        return s
        
    def setvariable(self,start,end,linspace):
        """
            设置自变量变化范围，返回一个 start-end ，间隔为 linspace 的 array
            该函数用来生成隶属函数实例的空间范围
        """
        self._variable_start = start
        self._variable_end = end
        self._linspace = linspace

    def _generateMF(self,x):
        '''
            生成函数的隶属度，组成一个列表集合
        '''
        y = []
        for i in range(self.numFunc):
            y.append(self.ArbFunc(x, *self.parameter[i])-self._min_generateMF()[i])
            # y.append(self.ArbFunc(x, *self.parameter[i]))
        return y

    def _min_generateMF(self):
        '''
            计算隶属函数在当前参数和当前自变量范围下的最小值
            该函数的作用是将隶属函数沿 y 轴方向向下平移最小值个单位，保证隶属函数的值 <= 1
            可以理解为该函数的系数
        '''
        min_mf = []
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        for i in range(self.numFunc):
            min_mf.append(min(self.ArbFunc(x, *self.parameter[i])))
        return min_mf

    def _max_generateMF(self):
        '''
            计算隶属函数在自变量范围内的最大值
        '''
        max_mf = []
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        y = self._generateMF(x)

        for i in range(self.numFunc):
            max_mf.append(max(y[i]))
        return max_mf

    def get_min_generateMF(self):
        """
            获取自变量范围内各个隶属函数的最小值
            注意：这里的最小值是原隶属函数的最小值，即还没有向下平移过的隶属度最小值
                平移后的隶属度最小值均为 0!
        """
        return self._min_generateMF()  

    def get_max_generateMF(self):
        """
            获取自变量范围内各个隶属函数的最大值
        """
        return self._max_generateMF()

    def MF_Plot(self,st):
        """
            画出隶属度曲线图，方便设置可能隶属度
        """
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        y = self._generateMF(x)
        plt.figure(figsize=(8,5))
        p = []
        for j in range(self.numFunc):
            plt.plot(x,y[j],label=st+': '+self.ArbFunc.__name__+'_%d'%(j+1))
        plt.grid(linestyle='-.')
        plt.legend()
        plt.show()

    def calculate_MD(self,x):
        """
            通过隶属度方程，计算隶属度集合
        """
        y = self._generateMF(x)
        return np.array(y)