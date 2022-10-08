import numpy as np
from FNS.Archimedean import *

class IVIFN:
    '''
        Interval-Value Intuitionistic Fuzzy Number
        MD denotes menbership degree, NMD denotes non-membership
    '''
    qrung=1
    def __init__(self,MDL,MDU,NMDL,NMDU):
        assert MDL<=MDU and NMDL<=NMDU, 'ERROR: The upper limit of membership degree is greater than the lower limit.'
        assert 0<=MDL<=1 and 0<=MDU<=1 and 0<=NMDL<=1 and 0<=NMDU<=1, 'ERROR: The membership degree and non-membership degree must be in the interval[0,1]'
        assert 0<=MDL+NMDL<=1 and 0<=MDU+NMDU<=1, 'ERROR: The sum of membership degree and non-membership degree must be in the interval[0,1]'
        self.MDL = np.asarray(MDL)
        self.MDU = np.asarray(MDU)
        self.NMDL = np.asarray(NMDL)
        self.NMDU = np.asarray(NMDU)
        
    def __repr__(self):
        return 'IVIFN:(MD:['+str(np.around(self.MDL,4))+','+str(np.around(self.MDU,4))+'],NMD:['+str(np.around(self.NMDL,4))+','+str(np.around(self.NMDU,4)) +'])'
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newIVIFN =IVIFN(0,0,0,0)
        newIVIFN.MDL = in_algebraic_tau(l*algebraic_tau(self.MDL))
        newIVIFN.MDU = in_algebraic_tau(l*algebraic_tau(self.MDU))
        newIVIFN.NMDL = in_algebraic_s(l*algebraic_s(self.NMDL))
        newIVIFN.NMDU = in_algebraic_s(l*algebraic_s(self.NMDU))
        return newIVIFN
    
    def Algebraic_Times(self,l):
        newIVIFN =IVIFN(0,0,0,0)
        newIVIFN.MDL = in_algebraic_s(l*algebraic_s(self.MDL))
        newIVIFN.MDU = in_algebraic_s(l*algebraic_s(self.MDU))
        newIVIFN.NMDL = in_algebraic_tau(l*algebraic_tau(self.NMDL))
        newIVIFN.NMDU = in_algebraic_tau(l*algebraic_tau(self.NMDU))
        return newIVIFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newIVIFN =IVIFN(0,0,0,0)
        newIVIFN.MDL = in_einstein_tau(l*einstein_tau(self.MDL))
        newIVIFN.MDU = in_einstein_tau(l*einstein_tau(self.MDU))
        newIVIFN.NMDL = in_einstein_s(l*einstein_s(self.NMDL))
        newIVIFN.NMDU = in_einstein_s(l*einstein_s(self.NMDU))
        return newIVIFN
    
    def Einstein_Times(self,l):
        newIVIFN =IVIFN(0,0,0,0)
        newIVIFN.MDL = in_einstein_s(l*einstein_s(self.MDL))
        newIVIFN.MDU = in_einstein_s(l*einstein_s(self.MDU))
        newIVIFN.NMDL = in_einstein_tau(l*einstein_tau(self.NMDL))
        newIVIFN.NMDU = in_einstein_tau(l*einstein_tau(self.NMDU))
        return newIVIFN
        
# Basic multiplication and addition operations
## Interval-Value Intuitionistic Fuzzy Number
## Algebraic Basic Operations 
def IVIFN_Algebraic_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR: The two IVIFNs are not the same IVIFN!'
    newIVIFN = IVIFN(0,0,0,0)
    newIVIFN.MDL = pithy_algebraic_T(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_algebraic_T(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_algebraic_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_algebraic_S(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN
def IVIFN_Algebraic_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR: The two IVIFNs are not the same IVIFN!'
    newIVIFN = IVIFN(0,0,0,0)
    newIVIFN.MDL = pithy_algebraic_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_algebraic_S(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_algebraic_T(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_algebraic_T(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN

## Einstein Basic Operations
def IVIFN_Einstein_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR: The two IVIFNs are not the same IVIFN!'
    
    newIVIFN = IVIFN(0,0,0,0)
    newIVIFN.MDL = pithy_einstein_T(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_einstein_T(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_einstein_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_einstein_S(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN
def IVIFN_Einstein_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==1, 'ERROR: The two IVIFNs are not the same IVIFN!'
    
    newIVIFN = IVIFN(0,0,0,0)
    newIVIFN.MDL = pithy_einstein_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_einstein_S(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_einstein_T(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_einstein_T(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN