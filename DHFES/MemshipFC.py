from DHFES.membershipFun import *
import numpy as np
from matplotlib import pyplot as plt

class MemshipFC(object):
    """
        隶属度或非隶属度函数类，将各种隶属度函数整理成一个集合类，方便计算对偶犹豫模糊可能隶属度和可能非隶属度
        =====================================================================================
        属性：
            memFunc: 表示隶属度函数名称，共有 8 种隶属度函数，分别是：
                sigmf(x,a,b)          基本 Sigmoid 激活函数      a 表示偏差或偏移量，b表示激活区域(越大越平坦)
                trimf(x,[a,b,c])      三角函数                  [a,b,c]表示三角的三点组成的数组，满足 a<=b<=c
                zmf(x,a,b)            Z-函数                    形状如 Z，a 表示函数做变化，b 表示函数右变化，满足 a<=b
                trapmf(x,[a,b,c,d])   梯形函数                  [a,b,c,d]表示梯形的四个点，满足 a<=b<=c<=d
                smf(x,a,b)            S-函数                    a 表示从 0 开始爬升的点，b 表示趋于 1 的稳定的点
                gaussmf(x,a,b)        高斯函数                  a 表示均值或中心值，b 表示标准差的高斯参数
                gauss2mf(x,a,b,c,d)   双结合高斯函数             a 表示第一个高斯函数的均值，b 表示第一个高斯函数标准差的高斯参数，c 和 d 分别表示第二个高斯函数的参数
                gbellmf(x,a,b,c)      广义贝尔函数               a 表示宽度，b 表示斜率，c 表示偏差或中心点
                
                注意：memFunc 只能是一种隶属函数，不支持多种函数，即 memFunc 不为数组
            parameter: 表示函数的参数数组
                parameter[(a,b),(a,b)] 为 tuple 类型的有2个的参数数组，适用 sigmf,zmf,smf,gaussmf
                parameter[[a,b,c],[a,b,c]] 为 list 类型的有3个参数的参数数组，适用 trimf
            numFunc: 表示函数个个数，也表示可能隶属度的个数。例如一个对偶犹豫模糊元素的隶属度有 3 个可能的值，则意味着 numFunc = 3
            
            _variable_start,_variable_end,_linspace: 表示函数的自变量范围，start 为起始自变量值，end 为终止自变量，linspace 为数间隔。
                默认为(0,1,100)，表示自变量范围 0-1，总共有 100 个数。
                属性为私有属性，外部不可调用，但可以通过 setvariable() 方法修改
        方法:
            setvariable(self,start,end,linspace):  设置自变量范围和间隔
            _min_generateMF(self):                 计算 self 设置好的自变量范围和参数下的范围内函数最小值，目的是被隶属度函数减去以标准化隶属度函数
            _generateMF(self,x):                   生成隶属度集合
            MF_Plot(self):                         画出设置的隶属度函数曲线，方便观察
            calculate_MD(self,x):                  计算 x 的隶属度，返回一个长度为 numFunc 的 array 类型的数组，即可能隶属度或非隶属度
    """
    memFunc = ''
    parameter = []
    numFunc = 0
    
    _variable_start = 0
    _variable_end = 1
    _linspace = 100
    
    def __init__(self,memFunc,parameter,numFunc):
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
        assert len(parameter) == numFunc,\
        'ERROR!The number of parameters is not equal to numFunc! numFunc=%d means you should create %d sets of parameters'%(numFunc,numFunc)
        if memFunc == 'sigmf' or memFunc == 'zmf' or memFunc == 'smf' or memFunc == 'gaussmf':
            for para in parameter:
                assert type(para) == tuple, 'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
                assert len(para) == 2,'ERROR:The number of the %s function is error!The number of parameters should be 2'%memFunc
        elif memFunc == 'gauss2mf':
            for para in parameter:
                assert type(para) == tuple,'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
                assert len(para) == 4,'ERROR:The number of the %s function is error!The number of parameters should be 4'%memFunc
        elif memFunc == 'gbellmf':
            for para in parameter:
                assert type(para) == tuple,'ERROR:The %s function\'s parameter type is error!Parameter type should be tuple!'%memFunc
                assert len(para) == 3,'ERROR:The number of the %s function is error!The number of parameters should be 3'%memFunc
        elif memFunc == 'trimf':
            for para in parameter:
                assert type(para) == list,'ERROR:The %s function\'s parameter type is error!Parameter type should be list!'%memFunc
                assert len(para) == 3,'ERROR:The number of the %s function is error!The number of parameters should be 3'%memFunc
        elif memFunc == 'trapmf':
            for para in parameter:
                assert type(para) == list,'ERROR:The %s function\'s parameter type is error!Parameter type should be list!'%memFunc
                assert len(para) == 4,'ERROR:The number of the %s function is error!The number of parameters should be 4'%memFunc
        else:
            assert memFunc == '','ERROR!Wrong function name!'
        
        self.memFunc = memFunc
        self.parameter = parameter
        self.numFunc = numFunc
        
    def __repr__(self):
        """
            Print information of membership function.
        """
        s = ''
        for i in range(self.numFunc):
            s+='Function %d: %s, parameter:%s \n'%(i+1,self.memFunc,self.parameter[i])
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
        """
            生成各函数的隶属度，组成一个列表集合
        """
        y = []
        for i in range(self.numFunc):
            if self.memFunc == 'sigmf':
                y.append(sigmf(x,self.parameter[i][0],self.parameter[i][1])-self._min_generateMF()[i])
            if self.memFunc == 'zmf':
                y.append(zmf(x,self.parameter[i][0],self.parameter[i][1])-self._min_generateMF()[i])
            if self.memFunc == 'smf':
                y.append(smf(x,self.parameter[i][0],self.parameter[i][1])-self._min_generateMF()[i])
            if self.memFunc == 'gaussmf':
                y.append(gaussmf(x,self.parameter[i][0],self.parameter[i][1])-self._min_generateMF()[i])
            if self.memFunc == 'gauss2mf':
                y.append(gauss2mf(x,self.parameter[i][0],self.parameter[i][1],self.parameter[i][2],self.parameter[i][3])-self._min_generateMF()[i])
            if self.memFunc == 'gbellmf':
                y.append(gbellmf(x,self.parameter[i][0],self.parameter[i][1],self.parameter[i][2])-self._min_generateMF()[i])
            if self.memFunc == 'trimf':
                y.append(trimf(x,self.parameter[i])-self._min_generateMF()[i])
            if self.memFunc == 'trapmf':
                y.append(trapmf(x,self.parameter[i])-self._min_generateMF()[i])
        return y
        
    def _min_generateMF(self):
        """
            计算隶属度函数在当前参数和当前自变量范围下的最小值
            该函数的作用是将隶属函数沿 y 轴方向向下平移最小值个单位，保证隶属度函数的值<=1
            可以理解为函数的系数
        """
        min_mf = []
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        for i in range(self.numFunc):
            if self.memFunc == 'sigmf':
                min_mf.append(min(sigmf(x,self.parameter[i][0],self.parameter[i][1])))
            if self.memFunc == 'zmf':
                min_mf.append(min(zmf(x,self.parameter[i][0],self.parameter[i][1])))
            if self.memFunc == 'smf':
                min_mf.append(min(smf(x,self.parameter[i][0],self.parameter[i][1])))
            if self.memFunc == 'gaussmf':
                min_mf.append(min(gaussmf(x,self.parameter[i][0],self.parameter[i][1])))
            if self.memFunc == 'gauss2mf':
                min_mf.append(min(gauss2mf(x,self.parameter[i][0],self.parameter[i][1],self.parameter[i][2],self.parameter[i][3])))
            if self.memFunc == 'gbellmf':
                min_mf.append(min(gbellmf(x,self.parameter[i][0],self.parameter[i][1],self.parameter[i][2])))
            if self.memFunc == 'trimf':
                min_mf.append(min(trimf(x,self.parameter[i])))
            if self.memFunc == 'trapmf':
                min_mf.append(min(trapmf(x,self.parameter[i])))
        return min_mf    
        
    def _max_generateMF(self):
        """
            计算隶属函数在自变量范围内的最大值
        """
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        y = self._generateMF(x)
        
        max_mf = []
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
            plt.plot(x,y[j],label=st+': '+self.memFunc+'_%d'%(j+1))
        plt.grid(linestyle='-.')
        plt.legend()
        
    def calculate_MD(self,x):
        """
            通过隶属度方程，计算隶属度集合
        """
        y = self._generateMF(x)
        return np.array(y)
