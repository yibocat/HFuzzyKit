import numpy as np
import FNS.Archimedean
from FNS.Archimedean import *

class PFN:
    '''
        Intuitionistic Fuzzy Number
        MD denotes menbership degree, NMD denotes non-membership
    '''
    qrung=2
    def __init__(self,MD,NMD):
        MD = np.asarray(MD)
        NMD = np.asarray(NMD)
        assert ((MD.size==0 or MD.size==1) and (NMD.size==0 or NMD.size==1) and 0<=MD<=1 and 0<=NMD<=1) and 0<=MD**2+NMD**2<=1, \
        'ERROR: Both of MD and NMD and MD^2+NMD^2 must have be in the interval[0,1] and the number of MD or NMD must have be 1.'
        self.MD = MD
        self.NMD = NMD
        
    def __repr__(self):
        return 'PFN:('+str(np.around(self.MD,4))+','+str(np.around(self.NMD,4))+')'
    
    ## Algebraic Basic Operations 
    def Algebraic_Power(self,l):
        '''
            PFN Algebraic power operation
            l: parameter, denote the power of PFN
            return: new PFN
        '''
        newPFN =PFN(0,0)
        newPFN.MD = in_algebraic_tau(l*algebraic_tau(self.MD**2))**(1/2)
        newPFN.NMD = in_algebraic_s(l*algebraic_s(self.NMD**2))**(1/2)
        return newPFN
    def Algebraic_Times(self,l):
        '''
            PFN Algebraic times operation
            l: parameter, denote the times of PFN
            return: new PFN
        '''
        newPFN =PFN(0,0)
        newPFN.MD = in_algebraic_s(l*algebraic_s(self.MD**2))**(1/2)
        newPFN.NMD = in_algebraic_tau(l*algebraic_tau(self.NMD**2))**(1/2)
        return newPFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        '''
            PFN Einstein power operation
            l: parameter, denote the Einstein power of PFN
            return: new PFN
        '''
        newPFN =PFN(0,0)
        newPFN.MD = in_einstein_tau(l*einstein_tau(self.MD**2))**(1/2)
        newPFN.NMD = in_einstein_s(l*einstein_s(self.NMD**2))**(1/2)
        return newPFN
    def Einstein_Times(self,l):
        '''
            PFN Einstein times operation
            l: parameter, denote the Einstein times of PFN
            return: new PFN
        '''
        newPFN =PFN(0,0)
        newPFN.MD = in_einstein_s(l*einstein_s(self.MD**2))**(1/2)
        newPFN.NMD = in_einstein_tau(l*einstein_tau(self.NMD**2))**(1/2)
        return newPFN
    
    def Fast_Algebraic_Power(self,l):
        '''
           Fast PFN Algebraic power operation
        '''
        newPFN=PFN(0,0)
        newPFN.MD=self.MD**l
        newPFN.NMD=(1-(1-self.NMD**2)**l)**(1/2)
        return newPFN
    def Fast_Algebraic_Times(self,l):
        '''
            Fast PFN Algebraic times operation
        '''
        newPFN=PFN(0,0)
        newPFN.MD = (1-(1-self.MD**2)**l)**(1/2)
        newPFN.NMD = self.NMD**l
        return newPFN

    def Fast_Einstein_Power(self,l):
        '''
             Fast PFN Einstein power operation
        '''
        newPFN=PFN(0,0)
        newPFN.MD = ((2*(self.MD**2)**l)/((2-self.MD**2)**l+(self.MD**2)**l))**(1/2)
        newPFN.NMD = (((1+self.NMD**2)**l-(1-self.NMD**2)**l)/((1+self.NMD**2)**l+(1-self.NMD**2)**l))**(1/2)
        return newPFN

    def Fast_Einstein_Times(self,l):
        '''
            Fast PFN Einstein times operation
        '''
        newPFN=PFN(0,0)
        newPFN.MD = (((1+self.MD**2)**l-(1-self.MD**2)**l)/((1+self.MD**2)**l+(1-self.MD**2)**l))**(1/2)
        newPFN.NMD = ((2*(self.NMD**2)**l)/((2-self.NMD**2)**l+(self.NMD**2)**l))**(1/2)
        return newPFN
    
########### Intersection and Union ###########
## Intersection
def PFN_Intersection(ifn1,ifn2):
    '''
        毕达哥拉斯模糊数交运算
    '''
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR:the two PFNs are not the same PFN or they are not PFN!'
    newPFN = PFN(0,0)
    newPFN.MD = min(ifn1.MD, ifn2.MD)
    newPFN.NMD = max(ifn1.NMD, ifn2.NMD)
    return newPFN

def PFN_Union(ifn1,ifn2):
    '''
        毕达哥拉斯模糊数并运算
    '''
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR:the two PFNs are not the same PFN or they are not PFN!'
    newPFN = PFN(0,0)
    newPFN.MD = max(ifn1.MD, ifn2.MD)
    newPFN.NMD = min(ifn1.NMD, ifn2.NMD)
    return newPFN

##########  Basic multiplication and addition operations #################
def PFN_Algebraic_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR:the two PFNs are not the same PFN or they are not PFN!'
    newPFN = PFN(0,0)
    newPFN.MD = pithy_algebraic_T(ifn1.MD**2,ifn2.MD**2)**(1/2)
    newPFN.NMD = pithy_algebraic_S(ifn1.NMD**2,ifn2.NMD**2)**(1/2)
    return newPFN

def PFN_Algebraic_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR:the two PFNs are not the same PFN or they are not PFN!'
    
    newPFN = PFN(0,0)
    newPFN.MD = pithy_algebraic_S(ifn1.MD**2,ifn2.MD**2)**(1/2)
    newPFN.NMD = pithy_algebraic_T(ifn1.NMD**2,ifn2.NMD**2)**(1/2)
    return newPFN

## Einstein Basic Operations
def PFN_Einstein_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR:the two PFNs are not the same PFN or they are not PFN!'
    
    newPFN = PFN(0,0)
    newPFN.MD = pithy_einstein_T(ifn1.MD**2,ifn2.MD**2)**(1/2)
    newPFN.NMD = pithy_einstein_S(ifn1.NMD**2,ifn2.NMD**2)**(1/2)
    return newPFN
def PFN_Einstein_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR:the two PFNs are not the same PFN or they are not PFN!'
    
    newPFN = PFN(0,0)
    newPFN.MD = pithy_einstein_S(ifn1.MD**2,ifn2.MD**2)**(1/2)
    newPFN.NMD = pithy_einstein_T(ifn1.NMD**2,ifn2.NMD**2)**(1/2)
    return newPFN