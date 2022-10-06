import numpy as np
import FNS.Archimedean
from FNS.Archimedean import *

class IFN:
    '''
        Intuitionistic Fuzzy Number
        MD denotes menbership degree, NMD denotes non-membership
    '''
    qrung=1
    def __init__(self,MD,NMD):
        MD = np.asarray(MD)
        NMD = np.asarray(NMD)
        assert ((MD.size==0 or MD.size==1) and (NMD.size==0 or NMD.size==1) and 0<=MD<=1 and 0<=NMD<=1) and 0<=MD+NMD<=1, \
        'ERROR: Both of MD and NMD and MD+NMD must have be in the interval[0,1] and the number of MD or NMD must have be 1.'
        self.MD = MD
        self.NMD = NMD
        
    def __repr__(self):
        return 'IFN:('+str(np.around(self.MD,4))+','+str(np.around(self.NMD,4))+')'
    
    ## Algebraic Basic Operations 
    def Algebraic_Power(self,l):
        '''
            IFN Algebraic power operation
            l: parameter, denote the power of IFN
            return: new IFN
        '''
        newIFN =IFN(0,0)
        newIFN.MD = in_algebraic_tau(l*algebraic_tau(self.MD))
        newIFN.NMD = in_algebraic_s(l*algebraic_s(self.NMD))
        return newIFN
    def Algebraic_Times(self,l):
        '''
            IFN Algebraic times operation
            l: parameter, denote the times of IFN
            return: new IFN
        '''
        newIFN =IFN(0,0)
        newIFN.MD = in_algebraic_s(l*algebraic_s(self.MD))
        newIFN.NMD = in_algebraic_tau(l*algebraic_tau(self.NMD))
        return newIFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        '''
            IFN Einstein power operation
            l: parameter, denote the Einstein power of IFN
            return: new IFN
        '''
        newIFN =IFN(0,0)
        newIFN.MD = in_einstein_tau(l*einstein_tau(self.MD))
        newIFN.NMD = in_einstein_s(l*einstein_s(self.NMD))
        return newIFN
    def Einstein_Times(self,l):
        '''
            IFN Einstein times operation
            l: parameter, denote the Einstein times of IFN
            return: new IFN
        '''
        newIFN =IFN(0,0)
        newIFN.MD = in_einstein_s(l*einstein_s(self.MD))
        newIFN.NMD = in_einstein_tau(l*einstein_tau(self.NMD))
        return newIFN
    
    def Fast_Algebraic_Power(self,l):
        '''
           Fast IFN Algebraic power operation
        '''
        newIFN=IFN(0,0)
        newIFN.MD=self.MD**l
        newIFN.NMD=(1-(1-self.NMD)**l)
        return newIFN
    def Fast_Algebraic_Times(self,l):
        '''
            Fast IFN Algebraic times operation
        '''
        newIFN=IFN(0,0)
        newIFN.MD = (1-(1-self.MD)**l)
        newIFN.NMD = self.NMD**l
        return newIFN

    def Fast_Einstein_Power(self,l):
        '''
             Fast IFN Einstein power operation
        '''
        newIFN=IFN(0,0)
        newIFN.MD = ((2*(self.MD)**l)/((2-self.MD)**l+(self.MD)**l))
        newIFN.NMD = (((1+self.NMD)**l-(1-self.NMD)**l)/((1+self.NMD)**l+(1-self.NMD)**l))
        return newIFN

    def Fast_Einstein_Times(self,l):
        '''
            Fast IFN Einstein times operation
        '''
        newIFN=IFN(0,0)
        newIFN.MD = (((1+self.MD)**l-(1-self.MD)**l)/((1+self.MD)**l+(1-self.MD)**l))
        newIFN.NMD = ((2*(self.NMD)**l)/((2-self.NMD)**l+(self.NMD)**l))
        return newIFN
    
########### Intersection and Union ###########
## Intersection
def IFN_Intersection(ifn1,ifn2):
    '''
        直觉模糊数交运算
    '''
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR:the two IFNs are not the same IFN or they are not IFN!'
    newIFN = IFN(0,0)
    newIFN.MD = min(ifn1.MD, ifn2.MD)
    newIFN.NMD = max(ifn1.NMD, ifn2.NMD)
    return newIFN

def IFN_Union(ifn1,ifn2):
    '''
        直觉模糊并运算
    '''
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR:the two IFNs are not the same IFN or they are not IFN!'
    newIFN = IFN(0,0)
    newIFN.MD = max(ifn1.MD, ifn2.MD)
    newIFN.NMD = min(ifn1.NMD, ifn2.NMD)
    return newIFN

##########  Basic multiplication and addition operations #################
def IFN_Algebraic_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR:the two IFNs are not the same IFN or they are not IFN!'
    newIFN = IFN(0,0)
    newIFN.MD = pithy_algebraic_T(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_algebraic_S(ifn1.NMD,ifn2.NMD)
    return newIFN

def IFN_Algebraic_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR:the two IFNs are not the same IFN or they are not IFN!'
    
    newIFN = IFN(0,0)
    newIFN.MD = pithy_algebraic_S(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_algebraic_T(ifn1.NMD,ifn2.NMD)
    return newIFN

## Einstein Basic Operations
def IFN_Einstein_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR:the two IFNs are not the same IFN or they are not IFN!'
    
    newIFN = IFN(0,0)
    newIFN.MD = pithy_einstein_T(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_einstein_S(ifn1.NMD,ifn2.NMD)
    return newIFN
def IFN_Einstein_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR:the two IFNs are not the same IFN or they are not IFN!'
    
    newIFN = IFN(0,0)
    newIFN.MD = pithy_einstein_S(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_einstein_T(ifn1.NMD,ifn2.NMD)
    return newIFN