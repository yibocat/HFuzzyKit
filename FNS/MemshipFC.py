from membershipFun import *
import numpy as np
from matplotlib import pyplot as plt

class MemshipFC(object):
    memFunc = ''
    parameter = []

    _variable_start = 0
    _variable_end = 1
    _linspace = 100

    def __init__(self,memFunc,parameter):
        """
            Set membership conditions：
            ===============================================
            Parameter:
                sigmf:   (a,b)      ->tuple /
                trimf:   [a,b,c]    ->list  /
                zmf:     (a,b)      ->tuple /
                trapmf:  [a,b,c,d]  ->list  /
                smf:     (a,b)      ->tuple /
                gaussmf: (a,b)      ->tuple /
                gauss2mf:(a,b,c,d)  ->tuple /
                gbellmf: (a,b,c)    ->tuple /
        """
        if memFunc == 'sigmf' or memFunc == 'zmf' or memFunc == 'smf' or memFunc == 'gaussmf':
            assert type(parameter) == tuple, 'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
            assert len(parameter) == 2,'ERROR:The number of the %s function is error!The number of parameters should be 2'%memFunc
        elif memFunc == 'gauss2mf':
            assert type(parameter) == tuple,'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
            assert len(parameter) == 4,'ERROR:The number of the %s function is error!The number of parameters should be 4'%memFunc
        elif memFunc == 'gbellmf':
            assert type(parameter) == tuple,'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
            assert len(parameter) == 3,'ERROR:The number of the %s function is error!The number of parameters should be 3'%memFunc
        elif memFunc == 'trimf':
            assert type(parameter) == list,'ERROR:The %s function\'s parameter type is error!Parameter type should be list!'%memFunc
            assert len(parameter) == 3,'ERROR:The number of the %s function is error!The number of parameters should be 3'%memFunc
        elif memFunc == 'trapmf':
            assert type(parameter) == list,'ERROR:The %s function\'s parameter type is error!Parameter type should be list!'%memFunc
            assert len(parameter) == 4,'ERROR:The number of the %s function is error!The number of parameters should be 4'%memFunc
        else:
            assert memFunc == '','ERROR!Wrong function name!'

        self.memFunc = memFunc
        self.parameter = parameter

    def __repr__(self):
        '''
            打印隶属函数信息
        '''
        return 'Function: %s, parameter: %s'%(self.memFunc,self.parameter)

    def setvariable(self,start,end,linspace):
        """
            设置自变量变化范围，返回一个 start-end ，间隔为 linspace 的 array
            该函数用来生成隶属函数实例的空间范围
        """
        self._variable_start = start
        self._variable_end = end
        self._linspace = linspace

    # def _generateMF(self,x):
    #     '''
    #         生成函数隶属度
    #     '''
    #     y = 0
    #     if self.memFunc == 'sigmf':
    #         y=sigmf(x,self.parameter[0],self.parameter[1])-self._min_generateMF()
    #     if self.memFunc == 'zmf':
    #         y=zmf(x,self.parameter[0],self.parameter[1])-self._min_generateMF()
    #     if self.memFunc == 'smf':
    #         y=smf(x,self.parameter[0],self.parameter[1])-self._min_generateMF()
    #     if self.memFunc == 'gaussmf':
    #         y=gaussmf(x,self.parameter[0],self.parameter[1])-self._min_generateMF()
    #     if self.memFunc == 'gauss2mf':
    #         y=gauss2mf(x,self.parameter[0],self.parameter[1],self.parameter[2],self.parameter[3])-self._min_generateMF()
    #     if self.memFunc == 'gbellmf':
    #         y=gbellmf(x,self.parameter[0],self.parameter[1],self.parameter[2])-self._min_generateMF()
    #     if self.memFunc == 'trimf':
    #         y=trimf(x,self.parameter)-self._min_generateMF()
    #     if self.memFunc == 'trapmf':
    #         y=trapmf(x,self.parameter)-self._min_generateMF()
    #     return y

    def _generateMF(self,x):
        '''
            生成函数隶属度
        '''
        y = 0
        if self.memFunc == 'sigmf':
            y=sigmf(x,self.parameter[0],self.parameter[1])
        if self.memFunc == 'zmf':
            y=zmf(x,self.parameter[0],self.parameter[1])
        if self.memFunc == 'smf':
            y=smf(x,self.parameter[0],self.parameter[1])
        if self.memFunc == 'gaussmf':
            y=gaussmf(x,self.parameter[0],self.parameter[1])
        if self.memFunc == 'gauss2mf':
            y=gauss2mf(x,self.parameter[0],self.parameter[1],self.parameter[2],self.parameter[3])
        if self.memFunc == 'gbellmf':
            y=gbellmf(x,self.parameter[0],self.parameter[1],self.parameter[2])
        if self.memFunc == 'trimf':
            y=trimf(x,self.parameter)
        if self.memFunc == 'trapmf':
            y=trapmf(x,self.parameter)
        return y


    def _min_generateMF(self):
        """
            计算隶属度函数在当前参数和当前自变量范围下的最小值
            该函数的作用是将隶属函数沿 y 轴方向向下平移最小值个单位，保证隶属度函数的值<=1
            可以理解为函数的系数
        """
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        min_mf=0
        if self.memFunc == 'sigmf':
            min_mf=min(sigmf(x,self.parameter[0],self.parameter[1]))
        if self.memFunc == 'zmf':
            min_mf=min(zmf(x,self.parameter[0],self.parameter[1]))
        if self.memFunc == 'smf':
            min_mf=min(smf(x,self.parameter[0],self.parameter[1]))
        if self.memFunc == 'gaussmf':
            min_mf=min(gaussmf(x,self.parameter[0],self.parameter[1]))
        if self.memFunc == 'gauss2mf':
            min_mf=min(gauss2mf(x,self.parameter[0],self.parameter[1],self.parameter[2],self.parameter[3]))
        if self.memFunc == 'gbellmf':
            min_mf=min(gbellmf(x,self.parameter[0],self.parameter[1],self.parameter[2]))
        if self.memFunc == 'trimf':
            min_mf=min(trimf(x,self.parameter))
        if self.memFunc == 'trapmf':
            min_mf=min(trapmf(x,self.parameter))
        return min_mf

    def _max_generateMF(self):
        """
            计算隶属函数在自变量范围内的最大值
        """
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        max_mf=max(self._generateMF(x))
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
        plt.plot(x,y,label=st+': '+self.memFunc)
        plt.grid(linestyle='-.')
        # plt.legend()
        plt.show()
        
    def calculate_MD(self,x):
        """
            通过隶属度方程，计算隶属度集合
        """
        y = self._generateMF(x)
        return np.array(y)