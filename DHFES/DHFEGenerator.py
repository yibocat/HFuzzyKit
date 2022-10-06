# import numpy as np
from DHFES import DHIFE,DHPFE,DHFFE
from DHFES.MemshipFC import *

class DHFEGenerator(object):
    """
        对偶犹豫模糊元素生成器
        =====================================================================================
        属性：
            qrung: 表示 q 级正交模糊，只可以取 1,2,3 三个值，分别表示直觉模糊，毕达哥拉斯模糊和费马模糊
            _variable_start,_variable_end,_linspace 为三个私有属性，表示自变量范围和自变量间隔
            MFunc: 隶属度函数，str 类型，即 8 种隶属度方法名
            NMFunc: 非隶属度函数，str
            MF_parameter: 隶属度参数，list 类型，表示使用的隶属度方法的参数。
            NMF_parameter: 非隶属度参数，list 类型，表示使用的非隶属度方法的参数。
            numMFC: 隶属度函数的个数
            numNMFC: 非隶属度函数的个数

            mf: 生成的隶属度函数集合，MemshipFC 类型
            nmf: 生成的非隶属度函数集合，MemshipFC 类型
        =====================================================================================
        方法:
            set_MF: 设置隶属函数的函数名
            set_NMF: 设置非隶属函数函数名
            set_MF_Num: 设置隶属函数的个数
            set_NMF_Num: 设置非隶属函数的个数
            set_MF_parameters: 设置隶属函数的参数
            set_NMF_parameters: 设置非隶属函数的参数

            set_Variable: 设置自变量范围和间隔
            generator_function: 生成隶属度和非隶属度函数
            MF_NMF_Plot: 画出隶属度和非隶属度的曲线图
            generator_DHFE: 生成 DHFE
        =====================================================================================
        步骤:
            1. 先进行设置，即隶属度函数选择、个数设置和参数设置
            2. 设置自变量范围和间隔
            3. 生成隶属度和非隶属度函数
            4. 生成 DHFE
            注意: 步骤不可乱
    """
    
    qrung = 0
    _variable_start = 0
    _variable_end = 1
    _linspace = 100
    
    def __init__(self,q):
        assert q==1 or q==2 or q==3,'The q-rung setting of DHFE is wrong, please reset it.\n'+'q=1: DHIFE; q=2,DHPFE; q=3:DHFFE'
        self.qrung = q
        
        self.MFunc = ''
        self.NMFunc = ''

        self.MF_parameter = []
        self.NMF_parameter = []

        self.numMFC = 0
        self.numNMFC = 0
        
        self.mf  = MemshipFC(self.MFunc,self.MF_parameter,self.numMFC)
        self.nmf = MemshipFC(self.NMFunc,self.NMF_parameter,self.numNMFC)
    
    def __repr__(self):
        return 'Membership function:\n'+str(self.mf) +'\n'+'Non-Membership function:\n'+str(self.nmf)
    
    def set_MF(self,MFunc):
        """
            设置隶属度函数
            MFunc: Str 类型，表示隶属函数种类
        """
        assert MFunc=='sigmf' or MFunc=='trimf' or MFunc=='zmf' or MFunc=='smf' or MFunc=='gaussmf' or MFunc=='gauss2mf' or MFunc=='gbellmf' or MFunc=='trapmf'\
        ,'ERROR! Wrong membership function!'
        self.MFunc = MFunc
        
    def set_NMF(self,NMFunc):
        """
            设置非隶属度函数
            NMFunc: Str 类型，表示非隶属函数种类
        """
        assert NMFunc=='sigmf' or NMFunc=='trimf' or NMFunc=='zmf' or NMFunc=='smf' or NMFunc=='gaussmf' or NMFunc=='gauss2mf' or NMFunc=='gbellmf' or NMFunc=='trapmf'\
        ,'ERROR! Wrong membership function!'
        self.NMFunc = NMFunc

    def set_MF_Num(self,n):
        """
            设置隶属度函数的个数
        """
        self.numMFC = n
    
    def set_NMF_Num(self,n):
        """
            设置非隶属度函数个数
        """
        self.numNMFC = n
    
    def set_MF_parameters(self,parameter):
        """
            设置隶属度函数的参数
            parameter: list 类型，表示参数列表
        """
        self.MF_parameter = parameter
        assert len(self.MF_parameter) == self.numMFC,'The number of MFCs has not been set or does not match the number of parameters.'
    
    def set_NMF_parameters(self,parameter):
        """
            设置非隶属度函数的参数
            parameter: list 类型，表示参数列表
        """
        self.NMF_parameter = parameter
        assert len(self.NMF_parameter) == self.numNMFC,'The number of NMFCs has not been set or does not match the number of parameters.'
    
    def set_Variable(self,start,end,linspace):
        self._variable_start = start
        self._variable_end = end
        self._linspace = linspace
    
    def generator_function(self):
        """
            隶属度函数与非隶属度函数生成器
            首先设置参数
            其次设置环境
        """
        assert self.MFunc!='' and self.NMFunc!='' and self.MF_parameter!=[] and self.NMF_parameter!=[] and self.numMFC!=0 and self.numNMFC!=0\
        ,'Membership function or parameter or number of function has not been set! Please set the membership and non-membership function first.'
        self.mf  = MemshipFC(self.MFunc,self.MF_parameter,self.numMFC)
        self.nmf = MemshipFC(self.NMFunc,self.NMF_parameter,self.numNMFC)
        
        self.mf.setvariable(self._variable_start,self._variable_end,self._linspace)
        self.nmf.setvariable(self._variable_start,self._variable_end,self._linspace)
        
    def MF_NMF_Plot(self):
        """
            画图，分别画出隶属度和非隶属度的曲线图
        """
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        mf = self.mf
        nmf = self.nmf
        
        mf.MF_Plot('Membership func')
        nmf.MF_Plot('Non-Membership func')
    
    def generator_DHFE(self,x,y):
        """
            生成对偶犹豫模糊元素
            x 表示隶属函数的自变量取值
            y 表示非隶属函数的自变量取值
        """
        assert self._variable_start<=x<=self._variable_end, 'The independent variable x is not in the range of %d and %d'%(self._variable_start,self._variable_end)
        assert self._variable_start<=y<=self._variable_end, 'The independent variable y is not in the range of %d and %d'%(self._variable_start,self._variable_end)

        if self.qrung == 1:
            newDHFE = DHIFE([],[])
        elif self.qrung == 2:
            newDHFE = DHPFE([],[])
        elif self.qrung == 3:
            newDHFE = DHFFE([],[])
        
        MD  = self.mf.calculate_MD(x)
        NMD = self.nmf.calculate_MD(y)
        assert np.max(MD)**self.qrung+np.max(NMD)**self.qrung<=1,'The MD^'+str(self.qrung)+'+NMD^'+str(self.qrung)+'<=1 and >=0. Please reset the parameters'
        assert np.min(MD)**self.qrung+np.min(NMD)**self.qrung>=0,'The MD^'+str(self.qrung)+'+NMD^'+str(self.qrung)+'<=1 and >=0. Please reset the parameters'
        newDHFE.MD  = self.mf.calculate_MD(x)
        newDHFE.NMD = self.nmf.calculate_MD(y)
        
        return newDHFE