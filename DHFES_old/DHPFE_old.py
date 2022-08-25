# from InteractionFNs import *
# from FNs import *
import math
from math import pow
import Archimedean
from Archimedean import *

## Dual Hesitant Pythagorean Fuzzy Elements
class DHPFE:
    qrung = 2
    def __init__(self,MD,NMD):
        try:
            assert (not MD or not NMD) or (max(MD)<=1 and max(NMD)<=1 and min(MD)>=0 and min(NMD)>=0) and (0<=max(MD)**2+max(NMD)**2<=1 and 0<=min(MD)**2+min(NMD)**2<=1)
            self.MD = MD
            self.NMD = NMD
        except AssertionError as er:
            print("ERROE:Construction failed!\n"+
                  "DHPFE max(MD)^2+max(NMD)^2 and min(MD)^2+min(NMD)^2 must be in interval[0,1]!\n"+
                  "It will return a variable of type 'NoneType'!",er)
            self = None
    
    def __repr__(self):
        mds,nmds = [],[]
        for md in self.MD:
            mds.append(float(format(md,'.4f')))
        for nmd in self.NMD:
            nmds.append(float(format(nmd,'.4f')))
        print('The cardinal number of  MD is %d'%len(mds))
        print('The cardinal number of NMD is %d'%len(nmds))
        if len(mds)>50 or len(nmds)>50:
            return '\nDHPFE:{' + ' MD: '+ str(mds[:50]) + ',\n' + '        NMD:' + str(nmds[:50]) +' }'
        else:
            return '\nDHPFE:{' + ' MD: '+ str(mds) + ',\n' + '        NMD:' + str(nmds) +' }'
    
    def isEmpty(self):
        ## 判断隶属度和非隶属度是否为空
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
        ## 对偶犹豫毕达哥拉斯模糊元素的补
        newDHFE = DHPFE([],[])
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
        ## 快速排序
        dhfe = DHPFE([],[])

        for i in sorted(self.MD):
            dhfe.MD.append(i)
        for j in sorted(self.NMD):
            dhfe.NMD.append(j)
        return dhfe
    
    def DHFEs_Qsort_reve(self):
        ## 对偶犹豫模糊元素隶属度与非隶属度排序
        ## q 表示 q-rung
        ## 从大到小排序
        dhfe = DHPFE([],[])

        for i in sorted(self.MD,reverse=True):
            dhfe.MD.append(i)
        for j in sorted(self.NMD,reverse=True):
            dhfe.NMD.append(j)
        return dhfe
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(in_algebraic_tau(l*algebraic_tau(pow(md,2))),1/2))
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(in_algebraic_s(l*algebraic_s(pow(nmd,2))),1/2))
        return newDHPFE
    def Algebraic_Times(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(in_algebraic_s(l*algebraic_s(pow(md,2))),1/2))
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(in_algebraic_tau(l*algebraic_tau(pow(nmd,2))),1/2))
        return newDHPFE
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(in_einstein_tau(l*einstein_tau(pow(md,2))),1/2))
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(in_einstein_s(l*einstein_s(pow(nmd,2))),1/2))
        return newDHPFE
    def Einstein_Times(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(in_einstein_s(l*einstein_s(pow(md,2))),1/2))
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(in_einstein_tau(l*einstein_tau(pow(nmd,2))),1/2))
        return newDHPFE
    
    ## Fast Algebraic Basic Operations
    def Fast_Algebraic_Power(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(md,l))
        for nmd in self.NMD:
            newDHPFE.NMD.append(
                pow(1-pow(1-nmd**2,l),1/2)
            )
        return newDHPFE
    def Fast_Algebraic_Times(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(
                pow(1-pow(1-md**2,l),1/2)
            )
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(nmd,l))
        return newDHPFE
    
    ## Fast Einstein Basic Operations
    def Fast_Einstein_Power(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(
                pow((2*pow(md**2,l))/(pow(2-md**2,l)+pow(md**2,l)),1/2)
            )
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(
                (pow(1+nmd**2,l)-pow(1-nmd**2,l))/(pow(1+nmd**2,l)+pow(1-nmd**2,l)),1/2)
            )
        return newDHPFE
    def Fast_Einstein_Times(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(
                (pow(1+md**2,l)-pow(1-md**2,l))/(pow(1+md**2,l)+pow(1-md**2,l)),1/2)
            )
        for nmd in self.NMD:
            newDHPFE.NMD.append(
                pow((2*pow(nmd**2,l))/(pow(2-nmd**2,l)+pow(nmd**2,l)),1/2)
            )
        return newDHPFE