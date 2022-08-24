# from InteractionFNs import *
# from FNs import *
import math
from math import pow
import Archimedean
from Archimedean import *

## Dual Hesitant Intuitionistic Fuzzy Elements
class DHIFE:
    def __init__(self,MD,NMD):
        try:
            assert (not MD or not NMD) or (max(MD)<=1 and max(NMD)<=1 and min(MD)>=0 and min(NMD)>=0) and (0<=max(MD)+max(NMD)<=1 and 0<=min(MD)+min(NMD)<=1)
            self.MD = MD
            self.NMD = NMD
        except AssertionError as er:
            print("ERROE:Construction failed!\n"+
                  "DHIFE max(MD)+max(NMD) and min(MD)+min(NMD) must be in interval[0,1]!\n"+
                  "It will return a variable of type 'NoneType'!",er)
            self = None
    
    def show(self):
        mds,nmds = [],[]
        for md in self.MD:
            mds.append(float(format(md,'.4f')))
        for nmd in self.NMD:
            nmds.append(float(format(nmd,'.4f')))
        print('The cardinal number of  MD is %d'%len(mds))
        print('The cardinal number of NMD is %d'%len(nmds))
        print('\nDHIFE:{\n' + 'MD: '+ str(mds) + ',\n' + 'NMD:' + str(nmds) +'}')

    def getQRung(self):
        return 1
        
    def isEmpty(self):
        ## 判断隶属度和非隶属度是否全为空
        if not self.MD and not self.NMD:
            return True
        else:
            return False
    def isEmpty_half(self):
        ## 判断隶属度和非隶属度是否有空
        if not self.MD or not self.NMD:
            return True
        else:
            return False
        
    def comp(self):
        ## 对偶犹豫模糊元素的补
        newDHFE = DHIFE([],[])
        if not self.MD and self.NMD:
            newDHFE.MD = []
            for nmd in self.NMD:
                newDHFE.NMD.append(1-nmd)
        elif self.MD and not self.NMD:
            newDHFE.NMD = []
            for md in self.MD:
                newDHFE.MD.append(1-md)
        else:
            newDHFE.MD = self.NMD
            newDHFE.NMD = self.MD
        return newDHFE
    
    def DHFEs_Qsort(self):
        ## 对偶犹豫模糊元素隶属度与非隶属度排序
        ## q 表示 q-rung
        ## 从小到大排序
        q = self.getQRung()

        if q==1:
            dhfe = DHIFE([],[])
        elif q==2:
            dhfe = DHPFE([],[])
        elif q==3:
            dhfe = DHFFE([],[])

        for i in sorted(self.MD):
            dhfe.MD.append(i)
        for j in sorted(self.NMD):
            dhfe.NMD.append(j)
        return dhfe
    
    def DHFEs_Qsort_reve(self):
        ## 对偶犹豫模糊元素隶属度与非隶属度排序
        ## q 表示 q-rung
        ## 从小到大排序
        q = self.getQRung()

        if q==1:
            dhfe = DHIFE([],[])
        elif q==2:
            dhfe = DHPFE([],[])
        elif q==3:
            dhfe = DHFFE([],[])

        for i in sorted(self.MD,reverse=True):
            dhfe.MD.append(i)
        for j in sorted(self.NMD,reverse=True):
            dhfe.NMD.append(j)
        return dhfe

    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(in_algebraic_tau(l*algebraic_tau(md)))
        for nmd in self.NMD:
            newDHIFE.NMD.append(in_algebraic_s(l*algebraic_s(nmd)))
        return newDHIFE
    def Algebraic_Times(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(in_algebraic_s(l*algebraic_s(md)))
        for nmd in self.NMD:
            newDHIFE.NMD.append(in_algebraic_tau(l*algebraic_tau(nmd)))
        return newDHIFE
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(in_einstein_tau(l*einstein_tau(md)))
        for nmd in self.NMD:
            newDHIFE.NMD.append(in_einstein_s(l*einstein_s(nmd)))
        return newDHIFE
    def Einstein_Times(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(in_einstein_s(l*einstein_s(md)))
        for nmd in self.NMD:
            newDHIFE.NMD.append(in_einstein_tau(l*einstein_tau(nmd)))
        return newDHIFE
    
    ## Fast Algebraic Basic Operations
    def Fast_Algebraic_Power(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(pow(md,l))
        for nmd in self.NMD:
            newDHIFE.NMD.append(1-pow(1-nmd,l))
        return newDHIFE

    def Fast_Algebraic_Times(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(1-pow(1-md,l))
        for nmd in self.NMD:
            newDHIFE.NMD.append(pow(nmd,l))
        return newDHIFE
    
    ## Fast Einstein Basic Operations
    def Fast_Einstein_Power(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append((2*pow(md,l))/(pow(2-md,l)+pow(md,l)))
        for nmd in self.NMD:
            newDHIFE.NMD.append((pow(1+nmd,l)-pow(1-nmd,l))/(pow(1+nmd,l)+pow(1-nmd,l)))
        return newDHIFE
    def Fast_Einstein_Times(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append((pow(1+md,l)-pow(1-md,l))/(pow(1+md,l)+pow(1-md,l)))
        for nmd in self.NMD:
            newDHIFE.NMD.append((2*pow(nmd,l))/(pow(2-nmd,l)+pow(nmd,l)))
        return newDHIFE