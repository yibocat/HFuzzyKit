from membershipFun import *
import numpy as np
from matplotlib import pyplot as plt

class MemshipFC(object):
    '''
        区间值模糊数的隶属函数生成器
    '''
    memFunc = ''
    upp_para = []
    low_para = []

    _variable_start = 0
    _variable_end = 1
    _linspace = 100

    def __init__(self,memFunc,low_para,upp_para):
        """
            Set membership conditions：
            ===============================================
            Parameter:
                sigmf:   (a,b)     ->[(a,b)]     ->tuple /
                trimf:   [a,b,c]   ->[[a,b,c]]   ->list  /
                zmf:     (a,b)     ->[(a,b)]     ->tuple /
                trapmf:  [a,b,c,d] ->[[a,b,c,d]] ->list
                smf:     (a,b)     ->[(a,b)]     ->tuple /
                gaussmf: (a,b)     ->[(a,b)]     ->tuple /
                gauss2mf:(a,b,c,d) ->[(a,b,c,d)] ->tuple /
                gbellmf: (a,b,c)   ->[(a,b,c)]   ->tuple /
        """
        assert len(upp_para) == len(low_para), 'ERROR: Arguments to the same function must be equal.'
        if memFunc == 'sigmf' or memFunc == 'zmf' or memFunc == 'smf' or memFunc == 'gaussmf':
            assert type(upp_para) == tuple and type(low_para) == tuple, 'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
            assert len(upp_para) == 2 and len(low_para) == 2,'ERROR:The number of the %s function is error!The number of parameters should be 2'%memFunc
        elif memFunc == 'gauss2mf':
            assert type(upp_para) == tuple and type(low_para) == tuple,'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
            assert len(upp_para) == 4 and len(low_para) == 4,'ERROR:The number of the %s function is error!The number of parameters should be 4'%memFunc
        elif memFunc == 'gbellmf':
            assert type(upp_para) == tuple and type(low_para) == tuple,'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
            assert len(upp_para) == 3 and len(low_para) == 3,'ERROR:The number of the %s function is error!The number of parameters should be 3'%memFunc
        elif memFunc == 'trimf':
            assert type(upp_para) == list and type(low_para) == list,'ERROR:The %s function\'s parameter type is error!Parameter type should be list!'%memFunc
            assert len(upp_para) == 3 and len(low_para) == 3,'ERROR:The number of the %s function is error!The number of parameters should be 3'%memFunc
        elif memFunc == 'trapmf':
            assert type(upp_para) == list and type(low_para) == list,'ERROR:The %s function\'s parameter type is error!Parameter type should be list!'%memFunc
            assert len(upp_para) == 4 and len(low_para) == 4,'ERROR:The number of the %s function is error!The number of parameters should be 4'%memFunc
        else:
            assert memFunc == '','ERROR!Wrong function name!'

        self.memFunc = memFunc
        self.upp_para = upp_para
        self.low_para = low_para

    def __repr__(self):
        '''
            打印隶属函数信息
        '''
        return 'Function: %s\nLower limit membership parameter: '%self.memFunc + str(self.low_para) +'\nUpper limit membership parameter: ' + str(self.upp_para)

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
            生成函数隶属度,首先生成隶属度下限，再生成隶属度上限
        '''

        y_low = 0
        if self.memFunc == 'sigmf':
            y_low=sigmf(x,self.low_para[0],self.low_para[1])-self._min_generateMF()[0]
        if self.memFunc == 'zmf':
            y_low=zmf(x,self.low_para[0],self.low_para[1])-self._min_generateMF()[0]
        if self.memFunc == 'smf':
            y_low=smf(x,self.low_para[0],self.low_para[1])-self._min_generateMF()[0]
        if self.memFunc == 'gaussmf':
            y_low=gaussmf(x,self.low_para[0],self.low_para[1])-self._min_generateMF()[0]
        if self.memFunc == 'gauss2mf':
            y_low=gauss2mf(x,self.low_para[0],self.low_para[1],self.low_para[2],self.low_para[3])-self._min_generateMF()[0]
        if self.memFunc == 'gbellmf':
            y_low=gbellmf(x,self.low_para[0],self.low_para[1],self.low_para[2])-self._min_generateMF()[0]
        if self.memFunc == 'trimf':
            y_low=trimf(x,self.low_para)-self._min_generateMF()[0]
        if self.memFunc == 'trapmf':
            y_low=trapmf(x,self.low_para)-self._min_generateMF()[0]

        y_upp = 0
        if self.memFunc == 'sigmf':
            y_upp=sigmf(x,self.upp_para[0],self.upp_para[1])-self._min_generateMF()[1]
        if self.memFunc == 'zmf':
            y_upp=zmf(x,self.upp_para[0],self.upp_para[1])-self._min_generateMF()[1]
        if self.memFunc == 'smf':
            y_upp=smf(x,self.upp_para[0],self.upp_para[1])-self._min_generateMF()[1]
        if self.memFunc == 'gaussmf':
            y_upp=gaussmf(x,self.upp_para[0],self.upp_para[1])-self._min_generateMF()[1]
        if self.memFunc == 'gauss2mf':
            y_upp=gauss2mf(x,self.upp_para[0],self.upp_para[1],self.upp_para[2],self.upp_para[3])-self._min_generateMF()[1]
        if self.memFunc == 'gbellmf':
            y_upp=gbellmf(x,self.upp_para[0],self.upp_para[1],self.upp_para[2])-self._min_generateMF()[1]
        if self.memFunc == 'trimf':
            y_upp=trimf(x,self.upp_para)-self._min_generateMF()
        if self.memFunc == 'trapmf':
            y_upp=trapmf(x,self.upp_para)-self._min_generateMF()

        # assert y_upp.all>y_low.all, 'ERROR: The upper membership limit must be greater than the lower membership limit. Select a new set of parameters.'
        return [y_low,y_upp]

    def _min_generateMF(self):
        '''
            计算隶属函数在当前自变量范围下的最小值
            该函数的作用是将隶属函数沿 与轴向下平移最小值个单位，保证隶属度函数的值 <= 1
            可以理解为函数的系数
            先计算下限最小值，再计算上限最小值，最后返回 list
        '''
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)

        min_mf_low=0
        if self.memFunc == 'sigmf':
            min_mf_low=min(sigmf(x,self.low_para[0],self.low_para[1]))
        if self.memFunc == 'zmf':
            min_mf_low=min(zmf(x,self.low_para[0],self.low_para[1]))
        if self.memFunc == 'smf':
            min_mf_low=min(smf(x,self.low_para[0],self.low_para[1]))
        if self.memFunc == 'gaussmf':
            min_mf_low=min(gaussmf(x,self.low_para[0],self.low_para[1]))
        if self.memFunc == 'gauss2mf':
            min_mf_low=min(gauss2mf(x,self.low_para[0],self.low_para[1],self.low_para[2],self.low_para[3]))
        if self.memFunc == 'gbellmf':
            min_mf_low=min(gbellmf(x,self.low_para[0],self.low_para[1],self.low_para[2]))
        if self.memFunc == 'trimf':
            min_mf_low=min(trimf(x,self.low_para))
        if self.memFunc == 'trapmf':
            min_mf_low=min(trapmf(x,self.low_para))

        min_mf_upp=0
        if self.memFunc == 'sigmf':
            min_mf_upp=min(sigmf(x,self.upp_para[0],self.upp_para[1]))
        if self.memFunc == 'zmf':
            min_mf_upp=min(zmf(x,self.upp_para[0],self.upp_para[1]))
        if self.memFunc == 'smf':
            min_mf_upp=min(smf(x,self.upp_para[0],self.upp_para[1]))
        if self.memFunc == 'gaussmf':
            min_mf_upp=min(gaussmf(x,self.upp_para[0],self.upp_para[1]))
        if self.memFunc == 'gauss2mf':
            min_mf_upp=min(gauss2mf(x,self.upp_para[0],self.upp_para[1],self.upp_para[2],self.upp_para[3]))
        if self.memFunc == 'gbellmf':
            min_mf_upp=min(gbellmf(x,self.upp_para[0],self.upp_para[1],self.upp_para[2]))
        if self.memFunc == 'trimf':
            min_mf_upp=min(trimf(x,self.upp_para))
        if self.memFunc == 'trapmf':
            min_mf_upp=min(trapmf(x,self.upp_para))

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