import numpy as np
from FNS.Archimedean import *

class IVPFN:
    '''
        Interval-Value Intuitionistic Fuzzy Number
        MD denotes menbership degree, NMD denotes non-membership
    '''
    qrung=2
    def __init__(self,MDL,MDU,NMDL,NMDU):
        assert MDL<=MDU and NMDL<=NMDU, 'ERROR: The upper limit of membership degree is greater than the lower limit.'
        assert 0<=MDL<=1 and 0<=MDU<=1 and 0<=NMDL<=1 and 0<=NMDU<=1, 'ERROR: The membership degree and non-membership degree must be in the interval[0,1]'
        assert 0<=MDL**2+NMDL**2<=1 and 0<=MDU**2+NMDU**2<=1, 'ERROR: The sum of membership degree and non-membership degree must be in the interval[0,1]'
        self.MDL = np.asarray(MDL)
        self.MDU = np.asarray(MDU)
        self.NMDL = np.asarray(NMDL)
        self.NMDU = np.asarray(NMDU)
        
    def __repr__(self):
        return 'IVPFN:(MD:['+str(np.around(self.MDL,4))+','+str(np.around(self.MDU,4))+'],NMD:['+str(np.around(self.NMDL,4))+','+str(np.around(self.NMDU,4)) +'])'
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newIVPFN =IVPFN(0,0,0,0)
        newIVPFN.MDL = in_algebraic_tau(l*algebraic_tau(self.MDL**2))**(1/2)
        newIVPFN.MDU = in_algebraic_tau(l*algebraic_tau(self.MDU**2))**(1/2)
        newIVPFN.NMDL = in_algebraic_s(l*algebraic_s(self.NMDL**2))**(1/2)
        newIVPFN.NMDU = in_algebraic_s(l*algebraic_s(self.NMDU**2))**(1/2)
        return newIVPFN
    
    def Algebraic_Times(self,l):
        newIVPFN =IVPFN(0,0,0,0)
        newIVPFN.MDL = in_algebraic_s(l*algebraic_s(self.MDL**2))**(1/2)
        newIVPFN.MDU = in_algebraic_s(l*algebraic_s(self.MDU**2))**(1/2)
        newIVPFN.NMDL = in_algebraic_tau(l*algebraic_tau(self.NMDL**2))**(1/2)
        newIVPFN.NMDU = in_algebraic_tau(l*algebraic_tau(self.NMDU**2))**(1/2)
        return newIVPFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newIVPFN =IVPFN(0,0,0,0)
        newIVPFN.MDL = in_einstein_tau(l*einstein_tau(self.MDL**2))**(1/2)
        newIVPFN.MDU = in_einstein_tau(l*einstein_tau(self.MDU**2))**(1/2)
        newIVPFN.NMDL = in_einstein_s(l*einstein_s(self.NMDL**2))**(1/2)
        newIVPFN.NMDU = in_einstein_s(l*einstein_s(self.NMDU**2))**(1/2)
        return newIVPFN
    
    def Einstein_Times(self,l):
        newIVPFN =IVPFN(0,0,0,0)
        newIVPFN.MDL = in_einstein_s(l*einstein_s(self.MDL**2))**(1/2)
        newIVPFN.MDU = in_einstein_s(l*einstein_s(self.MDU**2))**(1/2)
        newIVPFN.NMDL = in_einstein_tau(l*einstein_tau(self.NMDL**2))**(1/2)
        newIVPFN.NMDU = in_einstein_tau(l*einstein_tau(self.NMDU**2))**(1/2)
        return newIVPFN
        
# Basic multiplication and addition operations
## Interval-Value Intuitionistic Fuzzy Number
## Algebraic Basic Operations 
def IVPFN_Algebraic_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR: The two IVPFNs are not the same IVPFN!'
    newIVPFN = IVPFN(0,0,0,0)
    newIVPFN.MDL = pithy_algebraic_T(ifn1.MDL**2,ifn2.MDL**2)**(1/2)
    newIVPFN.MDU = pithy_algebraic_T(ifn1.MDU**2,ifn2.MDU**2)**(1/2)
    newIVPFN.NMDL = pithy_algebraic_S(ifn1.NMDL**2,ifn2.NMDL**2)**(1/2)
    newIVPFN.NMDU = pithy_algebraic_S(ifn1.NMDU**2,ifn2.NMDU**2)**(1/2)
    return newIVPFN
def IVPFN_Algebraic_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR: The two IVPFNs are not the same IVPFN!'
    newIVPFN = IVPFN(0,0,0,0)
    newIVPFN.MDL = pithy_algebraic_S(ifn1.MDL**2,ifn2.MDL**2)**(1/2)
    newIVPFN.MDU = pithy_algebraic_S(ifn1.MDU**2,ifn2.MDU**2)**(1/2)
    newIVPFN.NMDL = pithy_algebraic_T(ifn1.NMDL**2,ifn2.NMDL**2)**(1/2)
    newIVPFN.NMDU = pithy_algebraic_T(ifn1.NMDU**2,ifn2.NMDU**2)**(1/2)
    return newIVPFN

## Einstein Basic Operations
def IVPFN_Einstein_Multiply(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR: The two IVPFNs are not the same IVPFN!'
    
    newIVPFN = IVPFN(0,0,0,0)
    newIVPFN.MDL = pithy_einstein_T(ifn1.MDL**2,ifn2.MDL**2)**(1/2)
    newIVPFN.MDU = pithy_einstein_T(ifn1.MDU**2,ifn2.MDU**2)**(1/2)
    newIVPFN.NMDL = pithy_einstein_S(ifn1.NMDL**2,ifn2.NMDL**2)**(1/2)
    newIVPFN.NMDU = pithy_einstein_S(ifn1.NMDU**2,ifn2.NMDU**2)**(1/2)
    return newIVPFN
def IVPFN_Einstein_Plus(ifn1,ifn2):
    assert ifn1.qrung == ifn2.qrung and ifn1.qrung==2, 'ERROR: The two IVPFNs are not the same IVPFN!'
    
    newIVPFN = IVPFN(0,0,0,0)
    newIVPFN.MDL = pithy_einstein_S(ifn1.MDL**2,ifn2.MDL**2)**(1/2)
    newIVPFN.MDU = pithy_einstein_S(ifn1.MDU**2,ifn2.MDU**2)**(1/2)
    newIVPFN.NMDL = pithy_einstein_T(ifn1.NMDL**2,ifn2.NMDL**2)**(1/2)
    newIVPFN.NMDU = pithy_einstein_T(ifn1.NMDU**2,ifn2.NMDU**2)**(1/2)
    return newIVPFN