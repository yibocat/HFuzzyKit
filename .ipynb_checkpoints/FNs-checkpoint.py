import math
from math import *
import Archimedean

## Intuitionistic Fuzzy Number
class IFN:
    ## MD denotes menbership degree, NMD denotes non-membership
    def __init__(self,MD,NMD):
        # try:
        self.MD = MD
        self.NMD = NMD
        #     assert 0 <= MD + NMD <= 1
        # except AssertionError as er:
        #     print("ERROE: IFN MD + NMD must be in [0,1]",er)
        
    def show(self):
        print((self.MD,self.NMD))
        
    ## Algebraic Basic Operations 
    def Algebraic_Power(self,l):
        # Power operation
        newIFN =IFN(None,None)
        newIFN.MD = in_algebraic_tau(l*algebraic_tau(self.MD))
        newIFN.NMD = in_algebraic_s(l*algebraic_s(self.NMD))
        return newIFN
    def Algebraic_Times(self,l):
        # Times operation
        newIFN =IFN(None,None)
        newIFN.MD = in_algebraic_s(l*algebraic_s(self.MD))
        newIFN.NMD = in_algebraic_tau(l*algebraic_tau(self.NMD))
        return newIFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newIFN =IFN(None,None)
        newIFN.MD = in_einstein_tau(l*einstein_tau(self.MD))
        newIFN.NMD = in_einstein_s(l*einstein_s(self.NMD))
        return newIFN
    def Einstein_Times(self,l):
        newIFN =IFN(None,None)
        newIFN.MD = in_einstein_s(l*einstein_s(self.MD))
        newIFN.NMD = in_einstein_tau(l*einstein_tau(self.NMD))
        return newIFN

## Interval-Value Intuitionistic Fuzzy Number
class IVIFN:
    ## MDL denotes lower limit of menbership degree, MDU denotes upper limit of menbership degree,
    ## NMDL denotes lower limit of non-membership degree, NMDU denotes lower limit of non-membership degree
    def __init__(self,MDL,MDU,NMDL,NMDU):
        # try:
        self.MDL = MDL
        self.MDU = MDU
        self.NMDL = NMDL
        self.NMDU = NMDU
        #     assert 0<=MDL+NMDL<=1 and 0<=MDU+NMDU<=1
        # except AssertionError as er:
        #     print("ERROE: IFN MDL + NMDL and MDU + NMDU must be in [0,1]",er)
        
    def show(self):
        print(([self.MDL,self.MDU],[self.NMDL,self.NMDU]))
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newIVIFN =IVIFN(None,None,None,None)
        newIVIFN.MDL = in_algebraic_tau(l*algebraic_tau(self.MDL))
        newIVIFN.MDU = in_algebraic_tau(l*algebraic_tau(self.MDU))
        newIVIFN.NMDL = in_algebraic_s(l*algebraic_s(self.NMDL))
        newIVIFN.NMDU = in_algebraic_s(l*algebraic_s(self.NMDU))
        return newIVIFN
    
    def Algebraic_Times(self,l):
        newIVIFN =IVIFN(None,None,None,None)
        newIVIFN.MDL = in_algebraic_s(l*algebraic_s(self.MDL))
        newIVIFN.MDU = in_algebraic_s(l*algebraic_s(self.MDU))
        newIVIFN.NMDL = in_algebraic_tau(l*algebraic_tau(self.NMDL))
        newIVIFN.NMDU = in_algebraic_tau(l*algebraic_tau(self.NMDU))
        return newIVIFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newIVIFN =IVIFN(None,None,None,None)
        newIVIFN.MDL = in_einstein_tau(l*einstein_tau(self.MDL))
        newIVIFN.MDU = in_einstein_tau(l*einstein_tau(self.MDU))
        newIVIFN.NMDL = in_einstein_s(l*einstein_s(self.NMDL))
        newIVIFN.NMDU = in_einstein_s(l*einstein_s(self.NMDU))
        return newIVIFN
    
    def Einstein_Times(self,l):
        newIVIFN =IVIFN(None,None,None,None)
        newIVIFN.MDL = in_einstein_s(l*einstein_s(self.MDL))
        newIVIFN.MDU = in_einstein_s(l*einstein_s(self.MDU))
        newIVIFN.NMDL = in_einstein_tau(l*einstein_tau(self.NMDL))
        newIVIFN.NMDU = in_einstein_tau(l*einstein_tau(self.NMDU))
        return newIVIFN

## Pythagorean Fuzzy Number
class PFN:
    def __init__(self,MD,NMD):
    # try:
        self.MD = MD
        self.NMD = NMD
        #     assert 0 <= pow(MD,2) + pow(NMD,2) <= 1
        # except AssertionError as er:
        #     print("ERROE: IFN MD^2 + NMD^2 must be in [0,1]",er)
    
    def show(self):
        print((self.MD,self.NMD))
        
    ## Algebraic Basic Operations 
    def Algebraic_Power(self,l):
        newPFN = PFN(None,None)
        newPFN.MD = pow(in_algebraic_tau(l*algebraic_tau(pow(self.MD),2)),1/2)
        newPFN.NMD = pow(in_algebraic_s(l*algebraic_s(pow(self.NMD),2)),1/2)
        return newPFN
    def Algebraic_Times(self,l):
        # Times operation
        newPFN = PFN(None,None)
        newPFN.MD = pow(in_algebraic_s(l*algebraic_s(pow(self.MD,2))),1/2)
        newPFN.NMD = pow(in_algebraic_tau(l*algebraic_tau(pow(self.NMD,2))),1/2)
        return newPFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newPFN = PFN(None,None)
        newPFN.MD = pow(in_einstein_tau(l*einstein_tau(pow(self.MD,2))),1/2)
        newPFN.NMD = pow(in_einstein_s(l*einstein_s(pow(self.NMD,2))),1/2)
        return newPFN
    
    def Einstein_Times(self,l):
        newPFN = PFN(None,None)
        newPFN.MD = pow(in_einstein_s(l*einstein_s(pow(self.MD,1/2))),1/2)
        newPFN.NMD = pow(in_einstein_tau(l*einstein_tau(pow(self.NMD,1/2))),1/2)
        return newPFN

## Interval-Value Pythagorean Fuzzy Number
class IVPFN:
    ## MDL denotes lower limit of menbership degree, MDU denotes upper limit of menbership degree,
    ## NMDL denotes lower limit of non-membership degree, NMDU denotes lower limit of non-membership degree
    def __init__(self,MDL,MDU,NMDL,NMDU):
        # try:
        self.MDL = MDL
        self.MDU = MDU
        self.NMDL = NMDL
        self.NMDU = NMDU
        #     assert 0<=pow(MDL,2)+pow(NMDL,2)<=1 and 0<=pow(MDU,2)+pow(NMDU,2)<=1
        # except AssertionError as er:
        #     print("ERROE: IFN MDL^2 + NMDL^2 and MDU^2 + NMDU^2 must be in [0,1]",er)
        
    def show(self):
        print(([self.MDL,self.MDU],[self.NMDL,self.NMDU]))
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newIVPFN =IVPFN(None,None,None,None)
        newIVPFN.MDL = pow(in_algebraic_tau(l*algebraic_tau(pow(self.MDL,2))),1/2)
        newIVPFN.MDU = pow(in_algebraic_tau(l*algebraic_tau(pow(self.MDU,2))),1/2)
        newIVPFN.NMDL = pow(in_algebraic_s(l*algebraic_s(pow(self.NMDL,2))),1/2)
        newIVPFN.NMDU = pow(in_algebraic_s(l*algebraic_s(pow(self.NMDU,2))),1/2)
        return newIVPFN
    
    def Algebraic_Times(self,l):
        newIVPFN =IVPFN(None,None,None,None)
        newIVPFN.MDL = pow(in_algebraic_s(l*algebraic_s(pow(self.MDL,2))),1/2)
        newIVPFN.MDU = pow(in_algebraic_s(l*algebraic_s(pow(self.MDU,2))),1/2)
        newIVPFN.NMDL = pow(in_algebraic_tau(l*algebraic_tau(pow(self.NMDL,2))),1/2)
        newIVPFN.NMDU = pow(in_algebraic_tau(l*algebraic_tau(pow(self.NMDU,2))),1/2)
        return newIVPFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newIVPFN =IVPFN(None,None,None,None)
        newIVPFN.MDL = pow(in_einstein_tau(l*einstein_tau(pow(self.MDL,2))),1/2)
        newIVPFN.MDU = pow(in_einstein_tau(l*einstein_tau(pow(self.MDU,2))),1/2)
        newIVPFN.NMDL = pow(in_einstein_s(l*einstein_s(pow(self.NMDL,2))),1/2)
        newIVPFN.NMDU = pow(in_einstein_s(l*einstein_s(pow(self.NMDU,2))),1/2)
        return newIVPFN
    
    def Einstein_Times(self,l):
        newIVPFN =IVPFN(None,None,None,None)
        newIVPFN.MDL = pow(in_einstein_s(l*einstein_s(pow(self.MDL,2))),1/2)
        newIVPFN.MDU = pow(in_einstein_s(l*einstein_s(pow(self.MDU,2))),1/2)
        newIVPFN.NMDL = pow(in_einstein_tau(l*einstein_tau(pow(self.NMDL,2))),1/2)
        newIVPFN.NMDU = pow(in_einstein_tau(l*einstein_tau(pow(self.NMDU,2))),1/2)
        return newIVPFN

## Fermatean Fuzzy Number
class FFN:
    def __init__(self,MD,NMD):
        # try:
        self.MD = MD
        self.NMD = NMD
        #     assert 0 <= pow(MD,3) + pow(NMD,3) <= 1
        # except AssertionError as er:
        #     print("ERROE: IFN MD^3 + NMD^3 must be in [0,1]",er)
    
    def show(self):
        print((self.MD,self.NMD))
        
    ## Algebraic Basic Operations 
    def Algebraic_Power(self,l):
        newFFN = FFN(None,None)
        newFFN.MD = pow(in_algebraic_tau(l*algebraic_tau(pow(self.MD),3)),1/3)
        newFFN.NMD = pow(in_algebraic_s(l*algebraic_s(pow(self.NMD),3)),1/3)
        return newFFN
    def Algebraic_Times(self,l):
        # Times operation
        newFFN = FFN(None,None)
        newFFN.MD = pow(in_algebraic_s(l*algebraic_s(pow(self.MD,3))),1/3)
        newFFN.NMD = pow(in_algebraic_tau(l*algebraic_tau(pow(self.NMD,3))),1/3)
        return newFFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newFFN = FFN(None,None)
        newFFN.MD = pow(in_einstein_tau(l*einstein_tau(pow(self.MD,3))),1/3)
        newFFN.NMD = pow(in_einstein_s(l*einstein_s(pow(self.NMD,3))),1/3)
        return newFFN
    
    def Einstein_Times(self,l):
        newFFN = FFN(None,None)
        newFFN.MD = pow(in_einstein_s(l*einstein_s(pow(self.MD,1/3))),1/3)
        newFFN.NMD = pow(in_einstein_tau(l*einstein_tau(pow(self.NMD,1/3))),1/3)
        return newFFN

## Interval-Value Fermatean Fuzzy Number
class IVFFN:
    ## MDL denotes lower limit of menbership degree, MDU denotes upper limit of menbership degree,
    ## NMDL denotes lower limit of non-membership degree, NMDU denotes lower limit of non-membership degree
    def __init__(self,MDL,MDU,NMDL,NMDU):
        # try:
        self.MDL = MDL
        self.MDU = MDU
        self.NMDL = NMDL
        self.NMDU = NMDU
        #     assert 0<=pow(MDL,2)+pow(NMDL,2)<=1 and 0<=pow(MDU,2)+pow(NMDU,2)<=1
        # except AssertionError as er:
        #     print("ERROE: IFN MDL^2 + NMDL^2 and MDU^2 + NMDU^2 must be in [0,1]",er)
        
    def show(self):
        print(([self.MDL,self.MDU],[self.NMDL,self.NMDU]))
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newIVFFN =IVFFN(None,None,None,None)
        newIVFFN.MDL = pow(in_algebraic_tau(l*algebraic_tau(pow(self.MDL,3))),1/3)
        newIVFFN.MDU = pow(in_algebraic_tau(l*algebraic_tau(pow(self.MDU,3))),1/3)
        newIVFFN.NMDL = pow(in_algebraic_s(l*algebraic_s(pow(self.NMDL,3))),1/3)
        newIVFFN.NMDU = pow(in_algebraic_s(l*algebraic_s(pow(self.NMDU,3))),1/3)
        return newIVFFN
    
    def Algebraic_Times(self,l):
        newIVFFN =IVFFN(None,None,None,None)
        newIVFFN.MDL = pow(in_algebraic_s(l*algebraic_s(pow(self.MDL,3))),1/3)
        newIVFFN.MDU = pow(in_algebraic_s(l*algebraic_s(pow(self.MDU,3))),1/3)
        newIVFFN.NMDL = pow(in_algebraic_tau(l*algebraic_tau(pow(self.NMDL,3))),1/3)
        newIVFFN.NMDU = pow(in_algebraic_tau(l*algebraic_tau(pow(self.NMDU,3))),1/3)
        return newIVFFN
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newIVFFN =IVFFN(None,None,None,None)
        newIVFFN.MDL = pow(in_einstein_tau(l*einstein_tau(pow(self.MDL,3))),1/3)
        newIVFFN.MDU = pow(in_einstein_tau(l*einstein_tau(pow(self.MDU,3))),1/3)
        newIVFFN.NMDL = pow(in_einstein_s(l*einstein_s(pow(self.NMDL,3))),1/3)
        newIVFFN.NMDU = pow(in_einstein_s(l*einstein_s(pow(self.NMDU,3))),1/3)
        return newIVFFN
    
    def Einstein_Times(self,l):
        newIVFFN =IVFFN(None,None,None,None)
        newIVFFN.MDL = pow(in_einstein_s(l*einstein_s(pow(self.MDL,3))),1/3)
        newIVFFN.MDU = pow(in_einstein_s(l*einstein_s(pow(self.MDU,3))),1/3)
        newIVFFN.NMDL = pow(in_einstein_tau(l*einstein_tau(pow(self.NMDL,2))),1/3)
        newIVFFN.NMDU = pow(in_einstein_tau(l*einstein_tau(pow(self.NMDU,2))),1/3)
        return newIVFFN


##########  Basic multiplication and addition operations

## Basic multiplication and addition operations
## Intuitionistic Fuzzy Number
## Algebraic Basic Operations 
def IFN_Algebraic_Multiply(ifn1,ifn2):
    newIFN = IFN(None,None)
    newIFN.MD = pithy_algebraic_T(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_algebraic_S(ifn1.NMD,ifn2.NMD)
    return newIFN
def IFN_Algebraic_Plus(ifn1,ifn2):
    newIFN = IFN(None,None)
    newIFN.MD = pithy_algebraic_S(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_algebraic_T(ifn1.NMD,ifn2.NMD)
    return newIFN

## Einstein Basic Operations
def IFN_Einstein_Multiply(ifn1,ifn2):
    newIFN = IFN(None,None)
    newIFN.MD = pithy_einstein_T(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_einstein_S(ifn1.NMD,ifn2.NMD)
    return newIFN
def IFN_Einstein_Plus(ifn1,ifn2):
    newIFN = IFN(None,None)
    newIFN.MD = pithy_einstein_S(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_einstein_T(ifn1.NMD,ifn2.NMD)
    return newIFN

## Basic multiplication and addition operations
## Interval-Value Intuitionistic Fuzzy Number
## Algebraic Basic Operations 
def IFN_Algebraic_Multiply(ifn1,ifn2):
    newIVIFN = IVIFN(None,None,None,None)
    newIVIFN.MDL = pithy_algebraic_T(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_algebraic_T(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_algebraic_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_algebraic_S(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN
def IFN_Algebraic_Plus(ifn1,ifn2):
    newIVIFN = IVIFN(None,None,None,None)
    newIVIFN.MDL = pithy_algebraic_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_algebraic_S(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_algebraic_T(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_algebraic_T(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN

## Einstein Basic Operations
def IFN_Einstein_Multiply(ifn1,ifn2):
    newIVIFN = IVIFN(None,None,None,None)
    newIVIFN.MDL = pithy_einstein_T(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_einstein_T(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_einstein_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_einstein_S(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN
def IFN_Einstein_Plus(ifn1,ifn2):
    newIVIFN = IVIFN(None,None,None,None)
    newIVIFN.MDL = pithy_einstein_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_einstein_S(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_einstein_T(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_einstein_T(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN

## Basic multiplication and addition operations
## Pythagorean Fuzzy Number
## Algebraic Basic Operations 
def IFN_Algebraic_Multiply(ifn1,ifn2):
    newPFN = PFN(None,None)
    newPFN.MD = pow(pithy_algebraic_T(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(pithy_algebraic_S(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN
def IFN_Algebraic_Plus(ifn1,ifn2):
    newPFN = PFN(None,None)
    newPFN.MD = pow(pithy_algebraic_S(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(pithy_algebraic_T(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN

## Einstein Basic Operations
def IFN_Einstein_Multiply(ifn1,ifn2):
    newPFN = PFN(None,None)
    newPFN.MD = pow(pithy_einstein_T(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(pithy_einstein_S(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN
def IFN_Einstein_Plus(ifn1,ifn2):
    newPFN = PFN(None,None)
    newPFN.MD = pow(pithy_einstein_S(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(pithy_einstein_T(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN

## Basic multiplication and addition operations
## Interval-Value Pythagorean Fuzzy Number
## Algebraic Basic Operations 
def IFN_Algebraic_Multiply(ifn1,ifn2):
    newIVPFN = IVPFN(None,None,None,None)
    newIVPFN.MDL = pow(pithy_algebraic_T(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(pithy_algebraic_T(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(pithy_algebraic_S(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(pithy_algebraic_S(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN
def IFN_Algebraic_Plus(ifn1,ifn2):
    newIVPFN = IVPFN(None,None,None,None)
    newIVPFN.MDL = pow(pithy_algebraic_S(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(pithy_algebraic_S(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(pithy_algebraic_T(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(pithy_algebraic_T(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN

## Einstein Basic Operations
def IFN_Einstein_Multiply(ifn1,ifn2):
    newIVPFN = IVPFN(None,None,None,None)
    newIVPFN.MDL = pow(pithy_einstein_T(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(pithy_einstein_T(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(pithy_einstein_S(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(pithy_einstein_S(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN
def IFN_Einstein_Plus(ifn1,ifn2):
    newIVPFN = IVPFN(None,None,None,None)
    newIVPFN.MDL = pow(pithy_einstein_S(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(pithy_einstein_S(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(pithy_einstein_T(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(pithy_einstein_T(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN

## Basic multiplication and addition operations
## Fermatean Fuzzy Number
## Algebraic Basic Operations 
def IFN_Algebraic_Multiply(ifn1,ifn2):
    newFFN = FFN(None,None)
    newFFN.MD = pow(pithy_algebraic_T(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(pithy_algebraic_S(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN
def IFN_Algebraic_Plus(ifn1,ifn2):
    newFFN = FFN(None,None)
    newFFN.MD = pow(pithy_algebraic_S(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(pithy_algebraic_T(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN

## Einstein Basic Operations
def IFN_Einstein_Multiply(ifn1,ifn2):
    newFFN = FFN(None,None)
    newFFN.MD = pow(pithy_einstein_T(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(pithy_einstein_S(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN
def IFN_Einstein_Plus(ifn1,ifn2):
    newFFN = FFN(None,None)
    newFFN.MD = pow(pithy_einstein_S(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(pithy_einstein_T(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN

## Basic multiplication and addition operations
## Interval-Value Fermatean Fuzzy Number
## Algebraic Basic Operations 
def IFN_Algebraic_Multiply(ifn1,ifn2):
    newIVFFN = IVFFN(None,None,None,None)
    newIVFFN.MDL = pow(pithy_algebraic_T(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(pithy_algebraic_T(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(pithy_algebraic_S(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(pithy_algebraic_S(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN
def IFN_Algebraic_Plus(ifn1,ifn2):
    newIVFFN = IVFFN(None,None,None,None)
    newIVFFN.MDL = pow(pithy_algebraic_S(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(pithy_algebraic_S(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(pithy_algebraic_T(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(pithy_algebraic_T(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN

## Einstein Basic Operations
def IFN_Einstein_Multiply(ifn1,ifn2):
    newIVFFN = IVFFN(None,None,None,None)
    newIVFFN.MDL = pow(pithy_einstein_T(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(pithy_einstein_T(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(pithy_einstein_S(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(pithy_einstein_S(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN
def IFN_Einstein_Plus(ifn1,ifn2):
    newIVFFN = IVFFN(None,None,None,None)
    newIVFFN.MDL = pow(pithy_einstein_S(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(pithy_einstein_S(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(pithy_einstein_T(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(pithy_einstein_T(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN





















