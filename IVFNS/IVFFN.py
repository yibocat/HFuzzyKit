import numpy as np
from FNS.Archimedean import *

class IVFFN:
    '''
        Interval-Value Intuitionistic Fuzzy Number
        MD denotes menbership degree, NMD denotes non-membership
    '''
    qrung=3
    def __init__(self,MDL,MDU,NMDL,NMDU):
        assert MDL<=MDU and NMDL<=NMDU, 'ERROR: The upper limit of membership degree is greater than the lower limit.'
        assert 0<=MDL<=1 and 0<=MDU<=1 and 0<=NMDL<=1 and 0<=NMDU<=1, 'ERROR: The membership degree and non-membership degree must be in the interval[0,1]'
        assert 0<=MDL**3+NMDL**3<=1 and 0<=MDU**3+NMDU**3<=1, 'ERROR: The sum of membership degree and non-membership degree must be in the interval[0,1]'
        self.MDL = np.asarray(MDL)
        self.MDU = np.asarray(MDU)
        self.NMDL = np.asarray(NMDL)
        self.NMDU = np.asarray(NMDU)
        
    def __repr__(self):
        return 'IVFFN:(MD:['+str(np.around(self.MDL,4))+','+str(np.around(self.MDU,4))+'],NMD:['+str(np.around(self.NMDL,4))+','+str(np.around(self.NMDU,4)) +'])'
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newIVFFN =IVFFN(0,0,0,0)
        newIVFFN.MDL = in_algebraic_tau(l*algebraic_tau(self.MDL**3))**(1/3)
        newIVFFN.MDU = in_algebraic_tau(l*algebraic_tau(self.MDU**3))**(1/3)
        newIVFFN.NMDL = in_algebraic_s(l*algebraic_s(self.NMDL**3))**(1/3)
        newIVFFN.NMDU = in_algebraic_s(l*algebraic_s(self.NMDU**3))**(1/3)
        return newIVFFN
    
    def Algebraic_Times(self,l):
        newIVFFN =IVFFN(0,0,0,0)
        newIVFFN.MDL = in_algebraic_s(l*algebraic_s(self.MDL**3))**(1/3)
        newIVFFN.MDU = in_algebraic_s(l*algebraic_s(self.MDU**3))**(1/3)
        newIVFFN.NMDL = in_algebraic_tau(l*algebraic_tau(self.NMDL**3))**(1/3)
        newIVFFN.NMDU = in_algebraic_tau(l*algebraic_tau(self.NMDU**3))**(1/3)
        return newIVFFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newIVFFN =IVFFN(0,0,0,0)
        newIVFFN.MDL = in_einstein_tau(l*einstein_tau(self.MDL**3))**(1/3)
        newIVFFN.MDU = in_einstein_tau(l*einstein_tau(self.MDU**3))**(1/3)
        newIVFFN.NMDL = in_einstein_s(l*einstein_s(self.NMDL**3))**(1/3)
        newIVFFN.NMDU = in_einstein_s(l*einstein_s(self.NMDU**3))**(1/3)
        return newIVFFN
    
    def Einstein_Times(self,l):
        newIVFFN =IVFFN(0,0,0,0)
        newIVFFN.MDL = in_einstein_s(l*einstein_s(self.MDL**3))**(1/3)
        newIVFFN.MDU = in_einstein_s(l*einstein_s(self.MDU**3))**(1/3)
        newIVFFN.NMDL = in_einstein_tau(l*einstein_tau(self.NMDL**3))**(1/3)
        newIVFFN.NMDU = in_einstein_tau(l*einstein_tau(self.NMDU**3))**(1/3)
        return newIVFFN
        
# Basic multiplication and addition operations
## Interval-Value Intuitionistic Fuzzy Number
## Algebraic Basic Operations 
def IVFFN_Algebraic_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR: The two IVFFNs are not the same IVFFN!'
    newIVFFN = IVFFN(0,0,0,0)
    newIVFFN.MDL = pithy_algebraic_T(ifn1.MDL**3,ifn2.MDL**3)**(1/3)
    newIVFFN.MDU = pithy_algebraic_T(ifn1.MDU**3,ifn2.MDU**3)**(1/3)
    newIVFFN.NMDL = pithy_algebraic_S(ifn1.NMDL**3,ifn2.NMDL**3)**(1/3)
    newIVFFN.NMDU = pithy_algebraic_S(ifn1.NMDU**3,ifn2.NMDU**3)**(1/3)
    return newIVFFN
def IVFFN_Algebraic_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR: The two IVFFNs are not the same IVFFN!'
    newIVFFN = IVFFN(0,0,0,0)
    newIVFFN.MDL = pithy_algebraic_S(ifn1.MDL**3,ifn2.MDL**3)**(1/3)
    newIVFFN.MDU = pithy_algebraic_S(ifn1.MDU**3,ifn2.MDU**3)**(1/3)
    newIVFFN.NMDL = pithy_algebraic_T(ifn1.NMDL**3,ifn2.NMDL**3)**(1/3)
    newIVFFN.NMDU = pithy_algebraic_T(ifn1.NMDU**3,ifn2.NMDU**3)**(1/3)
    return newIVFFN

## Einstein Basic Operations
def IVFFN_Einstein_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR: The two IVFFNs are not the same IVFFN!'
    
    newIVFFN = IVFFN(0,0,0,0)
    newIVFFN.MDL = pithy_einstein_T(ifn1.MDL**3,ifn2.MDL**3)**(1/3)
    newIVFFN.MDU = pithy_einstein_T(ifn1.MDU**3,ifn2.MDU**3)**(1/3)
    newIVFFN.NMDL = pithy_einstein_S(ifn1.NMDL**3,ifn2.NMDL**3)**(1/3)
    newIVFFN.NMDU = pithy_einstein_S(ifn1.NMDU**3,ifn2.NMDU**3)**(1/3)
    return newIVFFN
def IVFFN_Einstein_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==3, 'ERROR: The two IVFFNs are not the same IVFFN!'
    
    newIVFFN = IVFFN(0,0,0,0)
    newIVFFN.MDL = pithy_einstein_S(ifn1.MDL**3,ifn2.MDL**3)**(1/3)
    newIVFFN.MDU = pithy_einstein_S(ifn1.MDU**3,ifn2.MDU**3)**(1/3)
    newIVFFN.NMDL = pithy_einstein_T(ifn1.NMDL**3,ifn2.NMDL**3)**(1/3)
    newIVFFN.NMDU = pithy_einstein_T(ifn1.NMDU**3,ifn2.NMDU**3)**(1/3)
    return newIVFFN