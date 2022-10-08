# import numpy as np
from DHFES import DHIFE,DHPFE,DHFFE
from DHFES.MemshipFC import *
from DHFES.CustomMemshipFC import *

class DHFEGenerator(object):
    '''
        对偶犹豫模糊元素生成器
        说明：
            这是一个对偶犹豫模糊元素生成器，利用对偶犹豫模糊隶属度函数生成器生成，有两种方式：
                1. 使用内建的 8 种隶属函数
                2. 自定义隶属函数
            属性 customFunc 来区分使用哪种隶属函数
        ===========================================================================
        属性：
            qrung: 表示创建哪种对偶犹豫模糊元素。若 qrung=1 表示对偶犹豫模糊元素；qrung=2 表示对偶犹豫毕达哥拉斯模糊元素；qrung=3 表示对偶犹豫费马模糊元素
            customFunc 表示自定义函数开关，用来选择哪种隶属函数。若为 False 则使用内建的 8 种隶属函数；若为 True 则使用自定义隶属函数
            _variable_start,_variable_end,_linspace 为三个私有属性，表示自变量范围和自变量间隔
            MF_parameter 表示隶属度参数列表
            NMF_parameter 表示非隶属度参数列表
            numMFC 表示隶属函数个数
            numNMFC 表示非隶属度函数个数
            MFunc, NMFunc:  当 customFunc 为 False 时，MFunc 和 NMFunc 表示内建隶属函数的函数名，‘str’类型
                            当 customFunc 为 True 时，MFunc 和 NMFunc 表示自定义函数本身，‘Function’ 类型

            mf 表示生成的隶属函数集合，MemshipFC 类型或 CustomMemshipFC 类型
            nmf 表示生成的非隶属函数集合，MemshipFC 类型或 CustomMemshipFC 类型
        方法：
            MFgeneratorSet(self,MFunc,MFnum,MFparas) 表示隶属函数设置
            NMFgeneratorSet(self,NMFunc,NMFnum,NMFparas) 表示非隶属函数设置
            set_Variable: 设置自变量范围和间隔
            generator_function: 生成隶属度和非隶属度函数
            MF_NMF_Plot: 画出隶属度和非隶属度的曲线图
            generator_DHFE: 生成 DHFE
        步骤：
            1. 先初始化，创建一个对偶犹豫模糊元素生成器
            2. 设置隶属度函数和非隶属度函数，使用 MFgeneratorSet 和 NMFgeneratorSet 方法
            3. 设置自变量范围和间隔
            4. 生成 DHFE
            注意：步骤不可打乱
    '''
    
    qrung = 0
    customFunc = False       ## 自定义函数开关，False 表示使用内建函数，True 表示使用自定义函数，默认为 False
    
    _variable_start = 0
    _variable_end = 1
    _linspace = 100
    
    def __init__(self,q,customFunc=False):
        assert q==1 or q==2 or q==3,'The q-rung setting of DHFE is wrong, please reset it.\n'+'q=1: DHIFE; q=2,DHPFE; q=3:DHFFE'
        self.qrung = q

        self.MF_parameter = []
        self.NMF_parameter = []

        self.numMFC = 0
        self.numNMFC = 0

        self.customFunc = customFunc
        if self.customFunc == False:
            self.MFunc = ''
            self.NMFunc = ''
            self.mf  = MemshipFC(self.MFunc,self.MF_parameter,self.numMFC)
            self.nmf = MemshipFC(self.NMFunc,self.NMF_parameter,self.numNMFC)
        else:
            self.MFunc = None
            self.NMFunc = None
            self.mf  = CustomMemshipFC(self.MFunc,self.MF_parameter,self.numMFC)
            self.nmf = CustomMemshipFC(self.NMFunc,self.NMF_parameter,self.numNMFC)
    
    def __repr__(self):
        return 'Membership function:\n'+str(self.mf) +'\n'+'Non-Membership function:\n'+str(self.nmf)
    
    def MFgeneratorSet(self,MFunc,MFnum,MFparas):
        '''
            隶属函数设置:
            先判断是否为自定义函数，查看 customFunc 属性
            MFunc:  当不是自定义函数时，str 类型，表示函数的名称
                    当是自定义函数时，function 类型，表示函数本身
            MFnum:  表示隶属函数的个数
            MFparas:表示参数列表
        '''
        if self.customFunc == False:
            assert MFunc=='sigmf' or MFunc=='trimf' or MFunc=='zmf' or MFunc=='smf' or MFunc=='gaussmf' or MFunc=='gauss2mf' or MFunc=='gbellmf' or MFunc=='trapmf',\
            'ERROR! Wrong membership function!'
            self.MFunc = MFunc
            self.numMFC = MFnum
            self.MF_parameter = MFparas
            assert len(self.MF_parameter) == self.numMFC,'The number of MFCs has not been set or does not match the number of parameters.'
        else:
            self.MFunc = MFunc
            self.numMFC = MFnum
            assert hasattr(self.MFunc,'__call__'),'ERROR:The MFunc is not a function!'
            self.MF_parameter = MFparas
            assert len(self.MF_parameter) == self.numMFC,'The number of MFCs has not been set or does not match the number of parameters.'

    def NMFgeneratorSet(self,NMFunc,NMFnum,NMFparas):
        '''
            非隶属函数设置:
            先判断是否为自定义函数，查看 customFunc 属性
            NMFunc:  当不是自定义函数时，str 类型，表示函数的名称
                     当是自定义函数时，function 类型，表示函数本身
            NMFnum:  表示非隶属函数的个数
            NMFparas:表示参数列表
        '''
        if self.customFunc == False:
            assert  NMFunc=='sigmf' or NMFunc=='trimf' or NMFunc=='zmf' or NMFunc=='smf' or NMFunc=='gaussmf' or NMFunc=='gauss2mf' or NMFunc=='gbellmf' or NMFunc=='trapmf',\
            'ERROR! Wrong membership function!'
            self.NMFunc = NMFunc
            self.numNMFC = NMFnum
            self.NMF_parameter = NMFparas
            assert len(self.NMF_parameter) == self.numNMFC,'The number of MFCs has not been set or does not match the number of parameters.'
        else:
            self.NMFunc = NMFunc
            self.numNMFC = NMFnum
            assert hasattr(self.NMFunc,'__call__'),'ERROR:The NMFunc is not a function!'
            self.NMF_parameter = NMFparas
            assert len(self.NMF_parameter) == self.numNMFC,'The number of MFCs has not been set or does not match the number of parameters.'

    def set_Variable(self,start,end,linspace):
        '''
            设置自变量表示范围
            start 表示开始
            end 表示结束
            linspace 表示间隔

            该函数可以用来调整图像的显示范围
        '''
        self._variable_start = start
        self._variable_end = end
        self._linspace = linspace
    
    def generator_function(self):
        """
            隶属度函数与非隶属度函数生成器
            首先设置参数
            其次设置环境
        """

        if self.customFunc == False:
            assert self.MFunc!='' and self.NMFunc!='' and self.MF_parameter!=[] and self.NMF_parameter!=[] and self.numMFC!=0 and self.numNMFC!=0\
            ,'Membership function or parameter or number of function has not been set! Please set the membership and non-membership function first.'
            self.mf  = MemshipFC(self.MFunc,self.MF_parameter,self.numMFC)
            self.nmf = MemshipFC(self.NMFunc,self.NMF_parameter,self.numNMFC)
        else:
            assert hasattr(self.MFunc,'__call__') and hasattr(self.NMFunc,'__call__'),'ERROR:The MFunc and NMFunc are not function type!'
            self.mf  = CustomMemshipFC(self.MFunc,self.MF_parameter,self.numMFC)
            self.nmf = CustomMemshipFC(self.NMFunc,self.NMF_parameter,self.numNMFC)
        
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