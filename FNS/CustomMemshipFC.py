import numpy as np
from matplotlib import pyplot as plt

class CustomMemshipFC(object):
    '''
        自定义隶属函数生成器
    '''

    parameter = []
    ArbFunc = None

    _variable_start = 0
    _variable_end = 1
    _linspace = 100

    def __init__(self,ArbFunc,parameter):
        assert type(parameter) == list, 'ERROR:The self-bulit function\'s parameter type is error! Parameter type should be list!'

        self.ArbFunc = ArbFunc
        self.parameter = np.asarray(parameter)

    def __repr__(self):
        '''
            打印隶属函数信息
        '''
        return 'Function : %s, parameters: %s'%(self.ArbFunc.__name__,self.parameter)
    
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
        return self.ArbFunc(x, *self.parameter)-self._min_generateMF()
        # return self.ArbFunc(x, *self.parameter)

    def _min_generateMF(self):
        '''
            计算隶属函数在当前参数和当前自变量范围下的最小值
            该函数的作用是将隶属函数沿 y 轴方向向下平移最小值个单位，保证隶属函数的值 <= 1
            可以理解为该函数的系数
        '''
        min_mf = []
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        min_mf=min(self.ArbFunc(x, *self.parameter))
        return min_mf
    
    def _max_generateMF(self):
        '''
            计算隶属函数在自变量范围内的最大值
        '''
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        return max(self._generateMF(x))

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
        plt.plot(x,y,label=st+': '+self.ArbFunc.__name__)
        plt.grid(linestyle='-.')
        # plt.legend()
        plt.show()

    def calculate_MD(self,x):
        """
            通过隶属度方程，计算隶属度集合
        """
        y = self._generateMF(x)
        return np.array(y)