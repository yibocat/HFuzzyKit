from IVFNS import IVIFN,IVPFN,IVFFN
from IVFNS.MemshipFC import *
from IVFNS.CustomMemshipFC import *

class IVFNGenerator(object):
    qrung = 0
    customFunc = False

    _variable_start = 0
    _variable_end = 1
    _linspace = 100

    def __init__(self,q,customFunc=False):
        assert q==1 or q==2 or q==3,'The q-rung setting of DHFE is wrong, please reset it.\n'+'q=1: DHIFE; q=2,DHPFE; q=3:DHFFE'
        self.qrung = q

        self.upp_MF_para = []
        self.low_MF_para = []
        self.upp_NMF_para = []
        self.low_NMF_para = []

        self.customFunc = customFunc

        if self.customFunc == False:
            self.MFunc = ''
            self.NMFunc = ''
            self.mf  = MemshipFC(self.MFunc,self.low_MF_para,self.upp_MF_para)
            self.nmf = MemshipFC(self.NMFunc,self.low_NMF_para,self.upp_NMF_para)
        else:
            self.MFunc = None
            self.NMFunc = None
            self.mf = CustomMemshipFC(self.MFunc,self.low_MF_para,self.upp_MF_para)
            self.nmf = CustomMemshipFC(self.NMFunc,self.low_NMF_para,self.upp_NMF_para)

    def __repr__(self):
        if self.mf.ArbFunc == None or self.nmf.ArbFunc == None:
            return 'Please set the membership function and parameters.'
        else:
            return 'Membership function:\n'+str(self.mf) +'\n'+'Non-Membership function:\n'+str(self.nmf)

    def MFgeneratorSet(self,MFunc,low_MF_para,upp_MF_para):
        if self.customFunc == False:
            assert MFunc=='sigmf' or MFunc=='trimf' or MFunc=='zmf' or MFunc=='smf' or MFunc=='gaussmf' or MFunc=='gauss2mf' or MFunc=='gbellmf' or MFunc=='trapmf',\
            'ERROR! Wrong membership function!'
            self.MFunc = Mfunc
            self.low_MF_para = low_MF_para
            self.upp_MF_para = upp_MF_para
        else:
            self.MFunc = MFunc
            self.low_MF_para = low_MF_para
            self.upp_MF_para = upp_MF_para
            assert hasattr(self.MFunc,'__call__'),'ERROR:The MFunc is not a function!'

    def NMFgeneratorSet(self,NMFunc,low_NMF_para,upp_NMF_para):
        if self.customFunc == False:
            assert NMFunc=='sigmf' or NMFunc=='trimf' or NMFunc=='zmf' or NMFunc=='smf' or NMFunc=='gaussmf' or NMFunc=='gauss2mf' or NMFunc=='gbellmf' or NMFunc=='trapmf',\
            'ERROR! Wrong membership function!'
            self.NMFunc = NMFunc
            self.low_NMF_para = low_NMF_para
            self.upp_NMF_para = upp_NMF_para
        else:
            self.NMFunc = NMFunc
            self.low_NMF_para = low_NMF_para
            self.upp_NMF_para = upp_NMF_para
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
            assert self.MFunc!='' and self.NMFunc!='' and self.low_MF_para!=[] and self.upp_MF_para!=[] and self.low_NMF_para!=[] and self.upp_NMF_para!=[],\
            'Membership function or parameter or number of function has not been set! Please set the membership and non-membership function first.'
            self.mf = MemshipFC(self.MFunc,self.low_MF_para,self.upp_MF_para)
            self.nmf = MemshipFC(self.NMFunc,self.low_NMF_para,self.upp_NMF_para)
        else:
            assert hasattr(self.MFunc,'__call__') and hasattr(self.NMFunc,'__call__'),'ERROR:The MFunc and NMFunc are not function type!'
            self.mf = CustomMemshipFC(self.MFunc,self.low_MF_para,self.upp_MF_para)
            self.nmf = CustomMemshipFC(self.NMFunc,self.low_NMF_para,self.upp_NMF_para)

        self.mf.setvariable(self._variable_start,self._variable_end,self._linspace)
        self.nmf.setvariable(self._variable_start,self._variable_end,self._linspace)

    def MF_NMF_Plot(self):
        x = np.linspace(self._variable_start,self._variable_end,self._linspace)
        mf = self.mf
        nmf = self.nmf

        mf.MF_Plot('Membership func')
        nmf.MF_Plot('Non-Membership func')

    def generator_IVFN(self,x,y):
        """
            生成模糊元素
            x 表示隶属函数的自变量取值
            y 表示非隶属函数的自变量取值
        """
        assert self._variable_start<=x<=self._variable_end, 'The independent variable x is not in the range of %d and %d'%(self._variable_start,self._variable_end)
        assert self._variable_start<=y<=self._variable_end, 'The independent variable y is not in the range of %d and %d'%(self._variable_start,self._variable_end)

        if self.qrung == 1:
            newIVFN = IVIFN(0,0,0,0)
        elif self.qrung == 2:
            newIVFN = IVPFN(0,0,0,0)
        elif self.qrung == 3:
            newIVFN = IVFFN(0,0,0,0)

        MD  = self.mf.calculate_MD(x)
        NMD = self.nmf.calculate_MD(y)

        assert np.max(MD)**self.qrung+np.max(NMD)**self.qrung<=1,'The MD^'+str(self.qrung)+'+NMD^'+str(self.qrung)+'<=1 and >=0. Please reset the parameters'
        assert np.min(MD)**self.qrung+np.min(NMD)**self.qrung>=0,'The MD^'+str(self.qrung)+'+NMD^'+str(self.qrung)+'<=1 and >=0. Please reset the parameters'

        newIVFN.MDL = self.mf.calculate_MD(x)[0]
        newIVFN.MDU = self.mf.calculate_MD(x)[1]
        newIVFN.NMDL = self.nmf.calculate_MD(y)[0]
        newIVFN.NMDU = self.nmf.calculate_MD(y)[1]

        return newIVFN

