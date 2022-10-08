from FNS import IFN,PFN,FFN
from FNS.MemshipFC import *
from FNS.CustomMemshipFC import *

class FNGenerator(object):
    qrung = 0
    customFunc = False

    _variable_start = 0
    _variable_end = 1
    _linspace = 100

    def __init__(self,q,customFunc=False):
        assert q==1 or q==2 or q==3,'The q-rung setting of DHFE is wrong, please reset it.\n'+'q=1: DHIFE; q=2,DHPFE; q=3:DHFFE'
        self.qrung = q

        self.MF_parameter = []
        self.NMF_parameter = []

        self.customFunc = customFunc
        if self.customFunc == False:
            self.MFunc = ''
            self.NMFunc = ''
            self.mf  = MemshipFC(self.MFunc,self.MF_parameter)
            self.nmf = MemshipFC(self.NMFunc,self.NMF_parameter)
        else:
            self.MFunc = None
            self.NMFunc = None
            self.mf = CustomMemshipFC(self.MFunc,self.MF_parameter)
            self.nmf = CustomMemshipFC(self.NMFunc,self.NMF_parameter)

    def __repr__(self):
        return 'Membership function:\n'+str(self.mf) +'\n'+'Non-Membership function:\n'+str(self.nmf)

    def MFgeneratorSet(self,MFunc,MFparas):
        '''
            隶属函数设置:
            先判断是否为自定义函数，查看 customFunc 属性
            MFunc:  当不是自定义函数时，str 类型，表示函数的名称
                    当是自定义函数时，function 类型，表示函数本身
            MFparas:表示参数
        '''
        if self.customFunc == False:
            assert MFunc=='sigmf' or MFunc=='trimf' or MFunc=='zmf' or MFunc=='smf' or MFunc=='gaussmf' or MFunc=='gauss2mf' or MFunc=='gbellmf' or MFunc=='trapmf',\
            'ERROR! Wrong membership function!'
            self.MFunc = MFunc
            self.MF_parameter = MFparas
        else:
            self.MFunc = MFunc
            self.MF_parameter = MFparas
            assert hasattr(self.MFunc,'__call__'),'ERROR:The MFunc is not a function!'
    
    def NMFgeneratorSet(self,NMFunc,NMFparas):
        '''
            非隶属函数设置:
            先判断是否为自定义函数，查看 customFunc 属性
            NMFunc:  当不是自定义函数时，str 类型，表示函数的名称
                     当是自定义函数时，function 类型，表示函数本身
            NMFparas:表示参数
        '''
        if self.customFunc == False:
            assert NMFunc=='sigmf' or NMFunc=='trimf' or NMFunc=='zmf' or NMFunc=='smf' or NMFunc=='gaussmf' or NMFunc=='gauss2mf' or NMFunc=='gbellmf' or NMFunc=='trapmf',\
            'ERROR! Wrong membership function!'
            self.NMFunc = NMFunc
            self.NMF_parameter = NMFparas
        else:
            self.NMFunc = NMFunc
            self.NMF_parameter = NMFparas
            assert hasattr(self.NMFunc,'__call__'),'ERROR:The NMFunc is not a function!'
    
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
            assert self.MFunc!='' and self.NMFunc!='' and self.MF_parameter!=[] and self.NMF_parameter!=[],\
            'Membership function or parameter or number of function has not been set! Please set the membership and non-membership function first.'
            self.mf  = MemshipFC(self.MFunc,self.MF_parameter)
            self.nmf = MemshipFC(self.NMFunc,self.NMF_parameter)
        else:
            assert hasattr(self.MFunc,'__call__') and hasattr(self.NMFunc,'__call__'),'ERROR:The MFunc and NMFunc are not function type!'
            self.mf  = CustomMemshipFC(self.MFunc,self.MF_parameter)
            self.nmf = CustomMemshipFC(self.NMFunc,self.NMF_parameter)
        
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

    def generator_FN(self,x,y):
        """
            生成模糊元素
            x 表示隶属函数的自变量取值
            y 表示非隶属函数的自变量取值
        """
        assert self._variable_start<=x<=self._variable_end, 'The independent variable x is not in the range of %d and %d'%(self._variable_start,self._variable_end)
        assert self._variable_start<=y<=self._variable_end, 'The independent variable y is not in the range of %d and %d'%(self._variable_start,self._variable_end)

        if self.qrung == 1:
            newFN = IFN(0,0)
        elif self.qrung == 2:
            newFN = PFN(0,0)
        elif self.qrung == 3:
            newFN = FFN(0,0)
        
        MD  = self.mf.calculate_MD(x)
        NMD = self.nmf.calculate_MD(y)
        assert np.max(MD)**self.qrung+np.max(NMD)**self.qrung<=1,'The MD^'+str(self.qrung)+'+NMD^'+str(self.qrung)+'<=1 and >=0. Please reset the parameters'
        assert np.min(MD)**self.qrung+np.min(NMD)**self.qrung>=0,'The MD^'+str(self.qrung)+'+NMD^'+str(self.qrung)+'<=1 and >=0. Please reset the parameters'
        newFN.MD  = self.mf.calculate_MD(x)
        newFN.NMD = self.nmf.calculate_MD(y)
        
        return newFN