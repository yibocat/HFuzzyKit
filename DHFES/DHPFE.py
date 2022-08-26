import DHFES.Archimedean
from DHFES.Archimedean import *

class DHPFE:
    qrung = 2
    def __init__(self,MD,NMD):
        MD = np.asarray(MD)
        NMD = np.asarray(NMD)
        assert (MD.size==0 or NMD.size==0) or (max(MD)<=1 and max(NMD)<=1 and min(MD)>=0 and min(NMD)>=0) and (0<=max(MD)**2+max(NMD)**2<=1 and 0<=min(MD)**2+min(NMD)**2<=1),"ERROE:Construction failed! DHPFE max(MD)^2+max(NMD)^2 and min(MD)^2+min(NMD)^2 must be in interval[0,1]!"
        self.MD = MD
        self.NMD = NMD

    def __repr__(self):
        print('The cardinal number of  MD is %d'%len(self.MD))
        print('The cardinal number of NMD is %d'%len(self.NMD))
        if len(self.MD)>50 or len(self.NMD)>50:
            return '\nDHPFE:{' + ' MD: '+ str(self.MD) + ',\n' + '        NMD:' + str(self.NMD) +' }\n'
        else:
            return '\nDHPFE:{' + ' MD: '+ str(self.MD) + ',\n' + '        NMD:' + str(self.NMD) +' }\n'

    def isEmpty(self):
        ## 判断隶属度和非隶属度是否为空
        if self.MD.size==0 and self.NMD.size==0:
            return True
        else:
            return False
        
    def isEmpty_half(self):
        ## 判断隶属度和非隶属度是否有空
        if self.MD.size==0 or self.NMD.size==0:
            return True
        else:
            return False

    def comp(self):
        ## 对偶犹豫费马模糊元素的补
        newDHFE = DHPFE([],[])
        if self.MD.size==0 and self.NMD.size!=0:
            newDHFE.MD = np.array([])
            newDHFE.NMD = 1-self.NMD
        elif self.MD.size!=0 and self.NMD.size==0:
            newDHFE.NMD = np.array([])
            newDHFE.MD = 1-self.MD
        else:
            newDHFE.MD = self.NMD
            newDHFE.NMD = self.MD
        return newDHFE

    def DHFEs_Qsort(self):
        ## 对偶犹豫模糊元素升序排序
        self.MD.sort()
        self.NMD.sort()
        return self

    def DHFEs_Qsort_reve(self):
        ## 对偶犹豫模糊元素降序排序
        self.MD = np.asarray(sorted(x.MD,reverse=True))
        self.NMD = np.asarray(sorted(x.NMD,reverse=True))
        return self

    def Algebraic_Power(self,l):
        '''
            代数范数幂次运算
            Algebraic norms power operations
            l 表示 l 次幂
        '''
        newDHFE = DHPFE([],[])
        newDHFE.MD = in_algebraic_tau(l*algebraic_tau(self.MD**2))**(1/2)
        newDHFE.NMD = in_algebraic_s(l*algebraic_s(self.NMD**2))**(1/2)
        return newDHFE

    def Algebraic_Times(self,l):
        '''
            代数范数倍运算
            Algebraic norms multiple operations
            l 表示倍数
        '''
        newDHFE = DHPFE([],[])
        newDHFE.MD = in_algebraic_s(l*algebraic_s(self.MD**2))**(1/2)
        newDHFE.NMD = in_algebraic_tau(l*algebraic_tau(self.NMD**2))**(1/2)
        return newDHFE

    def Einstein_Power(self,l):
        '''
            爱因斯坦范数幂运算
            l 表示幂次方
        '''
        newDHFE = DHPFE([],[])
        newDHFE.MD = in_einstein_tau(l*einstein_tau(self.MD**2))**(1/2)
        newDHFE.NMD = in_einstein_s(l*einstein_s(self.NMD**2))**(1/2)
        return newDHFE

    def Einstein_Times(self,l):
        '''
            爱因斯坦倍运算
            l 表示倍数
        '''
        newDHFE = DHPFE([],[])
        newDHFE.MD = in_einstein_s(l*einstein_s(self.MD**2))**(1/2)
        newDHFE.NMD = in_einstein_tau(l*einstein_tau(self.NMD**2))**(1/2)
        return newDHFE

    def Fast_Algebraic_Power(self,l):
        '''
            代数范数快速幂次运算
        '''
        newDHFE = DHPFE([],[])
        newDHFE.MD = self.MD**l
        newDHFE.NMD = (1-(1-self.NMD**2)**l)**(1/2)
        return newDHFE

    def Fast_Algebraic_Times(self,l):
        '''
            代数范数快速倍运算
        '''
        newDHFE = DHPFE([],[])
        newDHFE.MD = (1-(1-self.MD**2)**l)**(1/2)
        newDHFE.NMD = self.NMD**l
        return newDHFE

    def Fast_Einstein_Power(self,l):
        '''
            爱因斯坦快速幂运算
        '''
        newDHFE = DHPFE([],[])
        newDHFE.MD = ((2*(self.MD**2)**l)/((2-self.MD**2)**l+(self.MD**2)**l))**(1/2)
        newDHFE.NMD = (((1+self.NMD**2)**l-(1-self.NMD**2)**l)/((1+self.NMD**2)**l+(1-self.NMD**2)**l))**(1/2)
        return newDHFE

    def Fast_Einstein_Times(self,l):
        '''
            爱因斯坦快速倍运算
        '''
        newDHFE = DHPFE([],[])
        newDHFE.MD = (((1+self.MD**2)**l-(1-self.MD**2)**l)/((1+self.MD**2)**l+(1-self.MD**2)**l))**(1/2)
        newDHFE.NMD = ((2*(self.NMD**2)**l)/((2-self.NMD**2)**l+(self.NMD**2)**l))**(1/2)
        return newDHFE


########### Intersection and Union ###########
## Dual Hesitant Pythagorean Fuzzy Elements
## Intersection
def DHPFE_Intersection(dhfe1,dhfe2):
    
    assert dhfe1.qrung == dhfe2.qrung  and dhfe1.qrung==2 ,'ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!'          ## 判断是否是同一种对偶犹豫模糊集
    newDHFE = DHPFE([],[])
    
    md_min = min(max(dhfe1.MD),max(dhfe2.MD))
    nmd_max = max(min(dhfe1.NMD),min(dhfe2.NMD))
            
    md1 = dhfe1.MD[dhfe1.MD<=md_min] 
    md2 = dhfe2.MD[dhfe2.MD<=md_min]
    newDHFE.MD = np.unique(np.concatenate((md1,md2)))
    
    nmd1 = dhfe1.NMD[dhfe1.NMD>=nmd_max]
    nmd2 = dhfe2.NMD[dhfe2.NMD>=nmd_max]
    newDHFE.NMD = np.unique(np.concatenate((nmd1,nmd2)))
    
    return newDHFE

## Union
def DHPFE_Union(dhfe1,dhfe2):

    assert dhfe1.qrung == dhfe2.qrung  and dhfe1.qrung==2,'ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!'           ## 判断是否是同一种对偶犹豫模糊集
    newDHFE = DHPFE([],[])
    
    md_max = max(min(dhfe1.MD),min(dhfe2.MD))
    nmd_min = min(max(dhfe1.NMD),max(dhfe2.NMD))
    
    md1 = dhfe1.MD[dhfe1.MD>=md_max]
    md2 = dhfe2.MD[dhfe2.MD>=md_max]
    newDHFE.MD = np.unique(np.concatenate((md1,md2)))
    
    nmd1 = dhfe1.NMD[dhfe1.NMD<=nmd_min]
    nmd2 = dhfe2.NMD[dhfe2.NMD<=nmd_min]
    newDHFE.NMD = np.unique(np.concatenate((nmd1,nmd2)))
    
    return newDHFE


##########  Basic multiplication and addition operations
## Basic multiplication and addition operations
## Algebraic Basic Operations 
def DHPFE_Algebraic_Multiply(dhfe1,dhfe2):

    assert dhfe1.qrung == dhfe2.qrung  and dhfe1.qrung==2,'ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!'         ## 判断是否是同一种对偶犹豫模糊集

    dh_md1,dh_nmd1 = list(dhfe1.MD),list(dhfe1.NMD)
    dh_md2,dh_nmd2 = list(dhfe2.MD),list(dhfe2.NMD)
    
    MD = []
    NMD = []
    
    for md1 in dh_md1:
        for md2 in dh_md2:
            MD.append(pithy_algebraic_T(md1**2,md2**2)**(1/2))
    for nmd1 in dh_nmd1:
        for nmd2 in dh_nmd2:
            NMD.append(pithy_algebraic_S(nmd1**2,nmd2**2)**(1/2))

    newDHFE = DHPFE([],[])
    newDHFE.MD = np.asarray(MD)
    newDHFE.NMD = np.asarray(NMD)

    return newDHFE

def DHPFE_Algebraic_Plus(dhfe1,dhfe2):

    assert dhfe1.qrung == dhfe2.qrung  and dhfe1.qrung==2,'ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!'          ## 判断是否是同一种对偶犹豫模糊集

    dh_md1,dh_nmd1 = list(dhfe1.MD),list(dhfe1.NMD)
    dh_md2,dh_nmd2 = list(dhfe2.MD),list(dhfe2.NMD)
    
    MD = []
    NMD = []
    
    for md1 in dh_md1:
        for md2 in dh_md2:
            MD.append(pithy_algebraic_S(md1**2,md2**2)**(1/2))
    for nmd1 in dh_nmd1:
        for nmd2 in dh_nmd2:
            NMD.append(pithy_algebraic_T(nmd1**2,nmd2**2)**(1/2))
    
    newDHFE = DHPFE([],[])
    newDHFE.MD = np.asarray(MD)
    newDHFE.NMD = np.asarray(NMD)
    
    return newDHFE

## Einstein Basic Operations
def DHPFE_Einstein_Multiply(dhfe1,dhfe2):

    assert dhfe1.qrung == dhfe2.qrung  and dhfe1.qrung==2,'ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!'         ## 判断是否是同一种对偶犹豫模糊集
    
    dh_md1,dh_nmd1 = list(dhfe1.MD),list(dhfe1.NMD)
    dh_md2,dh_nmd2 = list(dhfe2.MD),list(dhfe2.NMD)
    
    MD = []
    NMD = []
    
    for md1 in dh_md1:
        for md2 in dh_md2:
            MD.append(pithy_einstein_T(md1**2,md2**2)**(1/2))
    for nmd1 in dh_nmd1:
        for nmd2 in dh_nmd2:
            NMD.append(pithy_einstein_S(nmd1**2,nmd2**2)**(1/2))

    newDHFE = DHPFE([],[])
    newDHFE.MD = np.asarray(MD)
    newDHFE.NMD = np.asarray(NMD)

    return newDHFE

def DHPFE_Einstein_Plus(dhfe1,dhfe2):

    assert dhfe1.qrung == dhfe2.qrung  and dhfe1.qrung==2,'ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!'          ## 判断是否是同一种对偶犹豫模糊集

    dh_md1,dh_nmd1 = list(dhfe1.MD),list(dhfe1.NMD)
    dh_md2,dh_nmd2 = list(dhfe2.MD),list(dhfe2.NMD)
    
    MD = []
    NMD = []
    
    for md1 in dh_md1:
        for md2 in dh_md2:
            MD.append(pithy_einstein_S(md1**2,md2**2)**(1/2))
    for nmd1 in dh_nmd1:
        for nmd2 in dh_nmd2:
            NMD.append(pithy_einstein_T(nmd1**2,nmd2**2)**(1/2))
    
    newDHFE = DHPFE([],[])
    newDHFE.MD = np.asarray(MD)
    newDHFE.NMD = np.asarray(NMD)
    
    return newDHFE
