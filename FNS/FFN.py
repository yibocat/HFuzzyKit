import numpy as np
import FNS.Archimedean
from FNS.Archimedean import *

class FFN:
    '''
        Intuitionistic Fuzzy Number
        MD denotes menbership degree, NMD denotes non-membership
    '''
    qrung=3
    def __init__(self,MD,NMD):
        MD = np.asarray(MD)
        NMD = np.asarray(NMD)
        assert ((MD.size==0 or MD.size==1) and (NMD.size==0 or NMD.size==1) and 0<=MD<=1 and 0<=NMD<=1) and 0<=MD**3+NMD**3<=1, \
        'ERROR: Both of MD and NMD and MD^3+NMD^3 must have be in the interval[0,1] and the number of MD or NMD must have be 1.'
        self.MD = MD
        self.NMD = NMD
        
    def __repr__(self):
        return 'FFN:('+str(np.around(self.MD,4))+','+str(np.around(self.NMD,4))+')'
    
    ## Algebraic Basic Operations 
    def Algebraic_Power(self,l):
        '''
            FFN Algebraic power operation
            l: parameter, denote the power of FFN
            return: new FFN
        '''
        newFFN =FFN(0,0)
        newFFN.MD = in_algebraic_tau(l*algebraic_tau(self.MD**3))**(1/3)
        newFFN.NMD = in_algebraic_s(l*algebraic_s(self.NMD**3))**(1/3)
        return newFFN
    def Algebraic_Times(self,l):
        '''
            FFN Algebraic times operation
            l: parameter, denote the times of FFN
            return: new FFN
        '''
        newFFN =FFN(0,0)
        newFFN.MD = in_algebraic_s(l*algebraic_s(self.MD**3))**(1/3)
        newFFN.NMD = in_algebraic_tau(l*algebraic_tau(self.NMD**3))**(1/3)
        return newFFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        '''
            FFN Einstein power operation
            l: parameter, denote the Einstein power of FFN
            return: new FFN
        '''
        newFFN =FFN(0,0)
        newFFN.MD = in_einstein_tau(l*einstein_tau(self.MD**3))**(1/3)
        newFFN.NMD = in_einstein_s(l*einstein_s(self.NMD**3))**(1/3)
        return newFFN
    def Einstein_Times(self,l):
        '''
            FFN Einstein times operation
            l: parameter, denote the Einstein times of FFN
            return: new FFN
        '''
        newFFN =FFN(0,0)
        newFFN.MD = in_einstein_s(l*einstein_s(self.MD**3))**(1/3)
        newFFN.NMD = in_einstein_tau(l*einstein_tau(self.NMD**3))**(1/3)
        return newFFN
    
    def Fast_Algebraic_Power(self,l):
        '''
           Fast FFN Algebraic power operation
        '''
        newFFN=FFN(0,0)
        newFFN.MD=self.MD**l
        newFFN.NMD=(1-(1-self.NMD**3)**l)**(1/3)
        return newFFN
    def Fast_Algebraic_Times(self,l):
        '''
            Fast FFN Algebraic times operation
        '''
        newFFN=FFN(0,0)
        newFFN.MD = (1-(1-self.MD**3)**l)**(1/3)
        newFFN.NMD = self.NMD**l
        return newFFN

    def Fast_Einstein_Power(self,l):
        '''
             Fast FFN Einstein power operation
        '''
        newFFN=FFN(0,0)
        newFFN.MD = ((2*(self.MD**3)**l)/((2-self.MD**3)**l+(self.MD**3)**l))**(1/3)
        newFFN.NMD = (((1+self.NMD**3)**l-(1-self.NMD**3)**l)/((1+self.NMD**3)**l+(1-self.NMD**3)**l))**(1/3)
        return newFFN

    def Fast_Einstein_Times(self,l):
        '''
            Fast FFN Einstein times operation
        '''
        newFFN=FFN(0,0)
        newFFN.MD = (((1+self.MD**3)**l-(1-self.MD**3)**l)/((1+self.MD**3)**l+(1-self.MD**3)**l))**(1/3)
        newFFN.NMD = ((2*(self.NMD**3)**l)/((2-self.NMD**3)**l+(self.NMD**3)**l))**(1/3)
        return newFFN
    
########### Intersection and Union ###########
## Intersection
def FFN_Intersection(ifn1,ifn2):
    '''
        毕达哥拉斯模糊数交运算
    '''
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR:the two FFNs are not the same FFN or they are not FFN!'
    newFFN = FFN(0,0)
    newFFN.MD = min(ifn1.MD, ifn2.MD)
    newFFN.NMD = max(ifn1.NMD, ifn2.NMD)
    return newFFN

def FFN_Union(ifn1,ifn2):
    '''
        毕达哥拉斯模糊数并运算
    '''
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR:the two FFNs are not the same FFN or they are not FFN!'
    newFFN = FFN(0,0)
    newFFN.MD = max(ifn1.MD, ifn2.MD)
    newFFN.NMD = min(ifn1.NMD, ifn2.NMD)
    return newFFN

##########  Basic multiplication and addition operations #################
def FFN_Algebraic_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR:the two FFNs are not the same FFN or they are not FFN!'
    newFFN = FFN(0,0)
    newFFN.MD = pithy_algebraic_T(ifn1.MD**3,ifn2.MD**3)**(1/3)
    newFFN.NMD = pithy_algebraic_S(ifn1.NMD**3,ifn2.NMD**3)**(1/3)
    return newFFN

def FFN_Algebraic_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR:the two FFNs are not the same FFN or they are not FFN!'
    
    newFFN = FFN(0,0)
    newFFN.MD = pithy_algebraic_S(ifn1.MD**3,ifn2.MD**3)**(1/3)
    newFFN.NMD = pithy_algebraic_T(ifn1.NMD**3,ifn2.NMD**3)**(1/3)
    return newFFN

## Einstein Basic Operations
def FFN_Einstein_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR:the two FFNs are not the same FFN or they are not FFN!'
    
    newFFN = FFN(0,0)
    newFFN.MD = pithy_einstein_T(ifn1.MD**3,ifn2.MD**3)**(1/3)
    newFFN.NMD = pithy_einstein_S(ifn1.NMD**3,ifn2.NMD**3)**(1/3)
    return newFFN
def FFN_Einstein_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR:the two FFNs are not the same FFN or they are not FFN!'
    
    newFFN = FFN(0,0)
    newFFN.MD = pithy_einstein_S(ifn1.MD**3,ifn2.MD**3)**(1/3)
    newFFN.NMD = pithy_einstein_T(ifn1.NMD**3,ifn2.NMD**3)**(1/3)
    return newFFN