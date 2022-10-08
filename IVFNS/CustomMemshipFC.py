import numpy as np
from matplotlib import pyplot as plt

class CustomMemshipFC(object):

    ArbFunc = None
    upp_para = []
    low_para = []

    _variable_start = 0
    _variable_end = 1
    _linspace = 100

    def __init__(self,ArbFunc,low_para,upp_para):
        assert hasattr(ArbFunc,'__call__') or ArbFunc == None, 'ERROR: The membership function has to be a function!'
        assert len(upp_para) == len(low_para), 'ERROR: Arguments to the same function must be equal.'
        self.ArbFunc = ArbFunc
        self.low_para = np.asarray(low_para)
        self.upp_para = np.asarray(upp_para)

    def __repr__(self):
        '''
            打印隶属函数信息
        '''
        return 'Function: %s\nLower limit membership parameter: '%self.ArbFunc.__name__ + str(self.low_para) +'\nUpper limit membership parameter: ' + str(self.upp_para)

    def setvariable(self,start,end,linspace):
        """
            设置自变量变化范围，返回一个 start-end ，间隔为 linspace 的 array
            该函数用来生成隶属函数实例的空间范围
        """
        self._variable_start = start
        self._variable_end = end
        self._linspace = linspace

    def _generateMF(self,x):
        y_low = self.ArbFunc(x,*self.low_para)-self._min_generateMF()[0]
        y_upp = self.ArbFunc(x,*self.upp_para)-self._min_generateMF()[1]
        return [y_low,y_upp]

    def _min_generateMF(self):
        '''
            计算隶属函数在当前参数和当前自变量范围下的最小值
            该函数的作用是将隶属函数沿 y 轴方向向下平移最小值个单位，保证隶属函数的值 <= 1
            可以理解为该函数的系数
        '''
        min_mf_low,min_mf_upp = [],[]
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        min_mf_low = min(self.ArbFunc(x, *self.low_para))
        min_mf_upp = min(self.ArbFunc(x, *self.upp_para))
        return [min_mf_low,min_mf_upp]

    def _max_generateMF(self):
        """
            计算上下限隶属函数在自变量范围内的最大值
        """
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        max_mf_low=max(self._generateMF(x)[0])
        max_mf_upp=max(self._generateMF(x)[1])
        return [max_mf_low,max_mf_upp]

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
        # for j in range(1):
        plt.plot(x,y[0],label=st+':Lower limit of membership')
        plt.plot(x,y[1],label=st+':Upper limit of membership')
        plt.grid(linestyle='-.')
        plt.legend()
        plt.show()
    
    def calculate_MD(self,x):
        '''
            通过上限隶属度函数和下限隶属度函数计算隶属度
        '''
        y = self._generateMF(x)
        if y[0]>y[1]:
            y[0],y[1] = y[1],y[0]
        return np.array(y)