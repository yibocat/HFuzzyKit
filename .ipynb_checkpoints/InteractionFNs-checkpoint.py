import math
from math import *
from FNs import *
from Archimedean import *

class InteractionIFN(IFN):
    ## Interaction Algebraic Basic Operations
    def Interaction_Algebraic_Power(self,l):
        newIFN = IFN(None,None)
        newIFN.MD = in_algebraic_s(l*algebraic_s(self.MD+self.NMD))-in_algebraic_s(l*algebraic_s(self.NMD))
        newIFN.NMD = in_algebraic_s(l*algebraic_s(self.NMD))
        return newIFN
    def Interaction_Algebraic_Times(self,l):
        newIFN = IFN(None,None)
        newIFN.MD = in_algebraic_s(l*algebraic_s(self.MD))
        newIFN.NMD = in_algebraic_s(l*algebraic_s(self.MD+self.NMD))-in_algebraic_s(l*algebraic_s(self.MD))
        return newIFN
    
    ## Interaction Einstein Basic Operations
    def Interaction_Einstein_Power(self,l):
        newIFN = IFN(None,None)
        newIFN.MD = in_einstein_s(l*einstein_s(self.MD+self.NMD))-in_einstein_s(l*einstein_s(self.NMD))
        newIFN.NMD = in_einstein_s(l*einstein_s(self.NMD))
        return newIFN
    def Interaction_Einstein_Times(self,l):
        newIFN = IFN(None,None)
        newIFN.MD = in_einstein_s(l*einstein_s(self.MD))
        newIFN.NMD = in_einstein_s(l*einstein_s(self.MD+self.NMD))-in_einstein_s(l*einstein_s(self.MD))
        return newIFN
    
class InteractionIVIFN(IVFFN):
    ## Interaction Algebraic Basic Operations
    def Interaction_Algebraic_Power(self,l):
        newIVIFN = IVIFN(None,None,None,None)
        newIVIFN.MDL = in_algebraic_s(l*algebraic_s(self.MDL+self.NMDL))-in_algebraic_s(l*algebraic_s(self.NMDL))
        newIVIFN.MDU = in_algebraic_s(l*algebraic_s(self.MDU+self.NMDU))-in_algebraic_s(l*algebraic_s(self.NMDU))
        newIVIFN.NMDL = in_algebraic_s(l*algebraic_s(self.NMDL))
        newIVIFN.NMDU = in_algebraic_s(l*algebraic_s(self.NMDU))
        return newIVIFN
    def Interaction_Algebraic_Times(self,l):
        newIVIFN = IFN(None,None,None,None)
        newIVIFN.MDL = in_algebraic_s(l*algebraic_s(self.MDL))
        newIVIFN.MDU = in_algebraic_s(l*algebraic_s(self.MDU))
        newIVIFN.NMDL = in_algebraic_s(l*algebraic_s(self.MDL+self.NMDL))-in_algebraic_s(l*algebraic_s(self.MDL))
        newIVIFN.NMDU = in_algebraic_s(l*algebraic_s(self.MDU+self.NMDU))-in_algebraic_s(l*algebraic_s(self.MDU))
        return newIVIFN
    
    ## Interaction Einstein Basic Operations
    def Interaction_Einstein_Power(self,l):
        newIVIFN = IVIFN(None,None,None,None)
        newIVIFN.MDL = in_einstein_s(l*einstein_s(self.MDL+self.NMDL))-in_einstein_s(l*einstein_s(self.NMDL))
        newIVIFN.MDU = in_einstein_s(l*einstein_s(self.MDU+self.NMDU))-in_einstein_s(l*einstein_s(self.NMDU))
        newIVIFN.NMDL = in_einstein_s(l*einstein_s(self.NMDL))
        newIVIFN.NMDU = in_einstein_s(l*einstein_s(self.NMDU))
        return newIVIFN
    def Interaction_Einstein_Times(self,l):
        newIVIFN = IVIFN(None,None,None,None)
        newIVIFN.MDL = in_einstein_s(l*einstein_s(self.MDL))
        newIVIFN.MDU = in_einstein_s(l*einstein_s(self.MDU))
        newIVIFN.NMDL = in_einstein_s(l*einstein_s(self.MDL+self.NMDL))-in_einstein_s(l*einstein_s(self.MDL))
        newIVIFN.NMDU = in_einstein_s(l*einstein_s(self.MDU+self.NMDU))-in_einstein_s(l*einstein_s(self.MDU))
        return newIVIFN

class InteractionPFN(PFN):
    ## Interaction Algebraic Basic Operations
    def Interaction_Algebraic_Power(self,l):
        newPFN = PFN(None,None)
        newPFN.MD = pow(in_algebraic_s(l*algebraic_s(pow(self.MD,2)+pow(self.NMD,2)))-in_algebraic_s(l*algebraic_s(pow(self.NMD,2))),1/2)
        newPFN.NMD = pow(in_algebraic_s(l*algebraic_s(pow(self.NMD,2))),1/2)
        return newPFN
    def Interaction_Algebraic_Times(self,l):
        newPFN = PFN(None,None)
        newPFN.MD = pow(in_algebraic_s(l*algebraic_s(pow(self.MD,2))),1/2)
        newPFN.NMD = pow(in_algebraic_s(l*algebraic_s(pow(self.MD,2)+pow(self.NMD,2)))-in_algebraic_s(l*algebraic_s(pow(self.MD,2))),1/2)
        return newPFN
    
    ## Interaction Einstein Basic Operations
    def Interaction_Einstein_Power(self,l):
        newPFN = PFN(None,None)
        newPFN.MD = pow(in_einstein_s(l*einstein_s(pow(self.MD,2)+pow(self.NMD,2)))-in_einstein_s(l*einstein_s(pow(self.NMD,2))),1/2)
        newPFN.NMD = pow(in_einstein_s(l*einstein_s(pow(self.NMD,2))),1/2)
        return newPFN
    def Interaction_Einstein_Times(self,l):
        newPFN = PFN(None,None)
        newPFN.MD = pow(in_einstein_s(l*einstein_s(pow(self.MD,2))),1/2)
        newPFN.NMD = pow(in_einstein_s(l*einstein_s(pow(self.MD,2)+pow(self.NMD,2)))-in_einstein_s(l*einstein_s(pow(self.MD,2))),1/2)
        return newPFN
    
class InteractionIVPFN(IVPFN):
    ## Interaction Algebraic Basic Operations
    def Interaction_Algebraic_Power(self,l):
        newIVPFN = IVPFN(None,None,None,None)
        newIVPFN.MDL = pow(in_algebraic_s(l*algebraic_s(pow(self.MDL,2)+pow(self.NMDL,2)))-in_algebraic_s(l*algebraic_s(pow(self.NMDL,2))),1/2)
        newIVPFN.MDU = pow(in_algebraic_s(l*algebraic_s(pow(self.MDU,2)+pow(self.NMDU,2)))-in_algebraic_s(l*algebraic_s(pow(self.NMDU,2))),1/2)
        newIVPFN.NMDL = pow(in_algebraic_s(l*algebraic_s(pow(self.NMDL,2))),1/2)
        newIVPFN.NMDU = pow(in_algebraic_s(l*algebraic_s(pow(self.NMDU,2))),1/2)
        return newIVPFN
    def Interaction_Algebraic_Times(self,l):
        newIVPFN = IVPFN(None,None,None,None)
        newIVPFN.MDL = pow(in_algebraic_s(l*algebraic_s(pow(self.MDL,2))),1/2)
        newIVPFN.MDU = pow(in_algebraic_s(l*algebraic_s(pow(self.MDU,2))),1/2)
        newIVPFN.NMDL = pow(in_algebraic_s(l*algebraic_s(pow(self.MDL,2)+pow(self.NMDL,2)))-in_algebraic_s(l*algebraic_s(pow(self.MDL,2))),1/2)
        newIVPFN.NMDU = pow(in_algebraic_s(l*algebraic_s(pow(self.MDU,2)+pow(self.NMDU,2)))-in_algebraic_s(l*algebraic_s(pow(self.MDU,2))),1/2)
        return newIVPFN
    
    ## Interaction Einstein Basic Operations
    def Interaction_Einstein_Power(self,l):
        newIVPFN = IVPFN(None,None,None,None)
        newIVPFN.MDL = pow(in_einstein_s(l*einstein_s(pow(self.MDL,2)+pow(self.NMDL,2)))-in_einstein_s(l*einstein_s(pow(self.NMDL,2))),1/2)
        newIVPFN.MDU = pow(in_einstein_s(l*einstein_s(pow(self.MDU,2)+pow(self.NMDU,2)))-in_einstein_s(l*einstein_s(pow(self.NMDU,2))),1/2)
        newIVPFN.NMDL = pow(in_einstein_s(l*einstein_s(pow(self.NMDL,2))),1/2)
        newIVPFN.NMDU = pow(in_einstein_s(l*einstein_s(pow(self.NMDU,2))),1/2)
        return newIVPFN
    def Interaction_Einstein_Times(self,l):
        newIVPFN = IVPFN(None,None,None,None)
        newIVPFN.MDL = pow(in_einstein_s(l*einstein_s(pow(self.MDL,2))),1/2)
        newIVPFN.MDU = pow(in_einstein_s(l*einstein_s(pow(self.MDU,2))),1/2)
        newIVPFN.NMDL = pow(in_einstein_s(l*einstein_s(pow(self.MDL,2)+pow(self.NMDL,2)))-in_einstein_s(l*einstein_s(pow(self.MDL,2))),1/2)
        newIVPFN.NMDU = pow(in_einstein_s(l*einstein_s(pow(self.MDU,2)+pow(self.NMDU,2)))-in_einstein_s(l*einstein_s(pow(self.MDU,2))),1/2)
        return newIVPFN
    
class InteractionFFN(FFN):
    ## Interaction Algebraic Basic Operations
    def Interaction_Algebraic_Power(self,l):
        newFFN = FFN(None,None)
        newFFN.MD = pow(in_algebraic_s(l*algebraic_s(pow(self.MD,3)+pow(self.NMD,3)))-in_algebraic_s(l*algebraic_s(pow(self.NMD,3))),1/3)
        newFFN.NMD = pow(in_algebraic_s(l*algebraic_s(pow(self.NMD,3))),1/3)
        return newFFN
    def Interaction_Algebraic_Times(self,l):
        newFFN = FFN(None,None)
        newFFN.MD = pow(in_algebraic_s(l*algebraic_s(pow(self.MD,3))),1/3)
        newFFN.NMD = pow(in_algebraic_s(l*algebraic_s(pow(self.MD,3)+pow(self.NMD,3)))-in_algebraic_s(l*algebraic_s(pow(self.MD,3))),1/3)
        return newFFN
    
    ## Interaction Einstein Basic Operations
    def Interaction_Einstein_Power(self,l):
        newFFN = FFN(None,None)
        newFFN.MD = pow(in_einstein_s(l*einstein_s(pow(self.MD,3)+pow(self.NMD,3)))-in_einstein_s(l*einstein_s(pow(self.NMD,3))),1/3)
        newFFN.NMD = pow(in_einstein_s(l*einstein_s(pow(self.NMD,3))),1/3)
        return newFFN
    def Interaction_Einstein_Times(self,l):
        newFFN = FFN(None,None)
        newFFN.MD = pow(in_einstein_s(l*einstein_s(pow(self.MD,3))),1/3)
        newFFN.NMD = pow(in_einstein_s(l*einstein_s(pow(self.MD,3)+pow(self.NMD,3)))-in_einstein_s(l*einstein_s(pow(self.MD,3))),1/3)
        return newFFN
    
class InteractionIVFFN(IVFFN):
    ## Interaction Algebraic Basic Operations
    def Interaction_Algebraic_Power(self,l):
        newIVFFN = IVFFN(None,None,None,None)
        newIVFFN.MDL = pow(in_algebraic_s(l*algebraic_s(pow(self.MDL,3)+pow(self.NMDL,3)))-in_algebraic_s(l*algebraic_s(pow(self.NMDL,3))),1/3)
        newIVFFN.MDU = pow(in_algebraic_s(l*algebraic_s(pow(self.MDU,3)+pow(self.NMDU,3)))-in_algebraic_s(l*algebraic_s(pow(self.NMDU,3))),1/3)
        newIVFFN.NMDL = pow(in_algebraic_s(l*algebraic_s(pow(self.NMDL,3))),1/3)
        newIVFFN.NMDU = pow(in_algebraic_s(l*algebraic_s(pow(self.NMDU,3))),1/3)
        return newIVFFN
    def Interaction_Algebraic_Times(self,l):
        newIVFFN = IVFFN(None,None,None,None)
        newIVFFN.MDL = pow(in_algebraic_s(l*algebraic_s(pow(self.MDL,3))),1/3)
        newIVFFN.MDU = pow(in_algebraic_s(l*algebraic_s(pow(self.MDU,3))),1/3)
        newIVFFN.NMDL = pow(in_algebraic_s(l*algebraic_s(pow(self.MDL,3)+pow(self.NMDL,3)))-in_algebraic_s(l*algebraic_s(pow(self.MDL,3))),1/3)
        newIVFFN.NMDU = pow(in_algebraic_s(l*algebraic_s(pow(self.MDU,3)+pow(self.NMDU,3)))-in_algebraic_s(l*algebraic_s(pow(self.MDU,3))),1/3)
        return newIVFFN
    
    ## Interaction Einstein Basic Operations
    def Interaction_Einstein_Power(self,l):
        newIVFFN = IVFFN(None,None,None,None)
        newIVFFN.MDL = pow(in_einstein_s(l*einstein_s(pow(self.MDL,3)+pow(self.NMDL,3)))-in_einstein_s(l*einstein_s(pow(self.NMDL,3))),1/3)
        newIVFFN.MDU = pow(in_einstein_s(l*einstein_s(pow(self.MDU,3)+pow(self.NMDU,3)))-in_einstein_s(l*einstein_s(pow(self.NMDU,3))),1/3)
        newIVFFN.NMDL = pow(in_einstein_s(l*einstein_s(pow(self.NMDL,3))),1/3)
        newIVFFN.NMDU = pow(in_einstein_s(l*einstein_s(pow(self.NMDU,3))),1/3)
        return newIVIFN
    def Interaction_Einstein_Times(self,l):
        newIVFFN = IVFFN(None,None,None,None)
        newIVFFN.MDL = pow(in_einstein_s(l*einstein_s(pow(self.MDL,3))),1/3)
        newIVFFN.MDU = pow(in_einstein_s(l*einstein_s(pow(self.MDU,3))),1/3)
        newIVFFN.NMDL = pow(in_einstein_s(l*einstein_s(pow(self.MDL,3)+pow(self.NMDL,3)))-in_einstein_s(l*einstein_s(pow(self.MDL,3))),1/3)
        newIVFFN.NMDU = pow(in_einstein_s(l*einstein_s(pow(self.MDU,3)+pow(self.NMDU,3)))-in_einstein_s(l*einstein_s(pow(self.MDU,3))),1/3)
        return newIVFFN
    
##########  Basic multiplication and addition operations

## Intuitionistic Fuzzy Number
## Interaction Algebraic Basic Operations
def Interaction_IFN_Algebraic_Multiply(ifn1,ifn2):
    newIFN = IFN(None,None)
    newIFN.MD = algebraic_S(ifn1.MD+ifn1.NMD,ifn2.MD+ifn2.NMD)-algebraic_S(ifn1.NMD,ifn2.NMD)
    newIFN.NMD = algebraic_S(ifn1.NMD,ifn2.NMD)
    return newIFN
def Interaction_IFN_Algebraic_Plus(ifn1,ifn2):
    newIFN = IFN(None,None)
    newIFN.MD = algebraic_S(ifn1.MD,ifn2.MD)
    newIFN.NMD = algebraic_S(ifn1.MD+ifn1.NMD,ifn2.MD+ifn2.NMD)-algebraic_S(ifn1.MD,ifn2.MD)
    return newIFN

def Interaction_IFN_Einstein_Multiply(ifn1,ifn2):
    newIFN = IFN(None,None)
    newIFN.MD = einstein_S(ifn1.MD+ifn1.NMD,ifn2.MD+ifn2.NMD)-einstein_S(ifn1.NMD,ifn2.NMD)
    newIFN.NMD = einstein_S(ifn1.NMD,ifn2.NMD)
    return newIFN
def Interaction_IFN_Einstein_Plus(ifn1,ifn2):
    newIFN = IFN(None,None)
    newIFN.MD = einstein_S(ifn1.MD,ifn2.MD)
    newIFN.NMD = einstein_S(ifn1.MD+ifn1.NMD,ifn2.MD+ifn2.NMD)-einstein_S(ifn1.MD,ifn2.MD)
    return newIFN

## Interval-value Intuitionistic Fuzzy Number
## Interaction Algebraic Basic Operations
def Interaction_IVIFN_Algebraic_Multiply(ifn1,ifn2):
    newIVIFN = IVIFN(None,None,None,None)
    newIVIFN.MDL = algebraic_S(ifn1.MDL+ifn1.NMDL,ifn2.MDL+ifn2.NMDL)-algebraic_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.MDU = algebraic_S(ifn1.MDU+ifn1.NMDU,ifn2.MDU+ifn2.NMDU)-algebraic_S(ifn1.NMDU,ifn2.NMDU)
    newIVIFN.NMDL = algebraic_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = algebraic_S(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN
def Interaction_IVIFN_Algebraic_Plus(ifn1,ifn2):
    newIVIFN = IVIFN(None,None,None,None)
    newIVIFN.MDL = algebraic_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = algebraic_S(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = algebraic_S(ifn1.MDL+ifn1.NMDL,ifn2.MDL+ifn2.NMDL)-algebraic_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.NMDU = algebraic_S(ifn1.MDU+ifn1.NMDU,ifn2.MDU+ifn2.NMDU)-algebraic_S(ifn1.MDU,ifn2.MDU)
    return newIVIFN

def Interaction_IVIFN_Einstein_Multiply(ifn1,ifn2):
    newIVIFN = IVIFN(None,None,None,None)
    newIVIFN.MDL = einstein_S(ifn1.MDL+ifn1.NMDL,ifn2.MDL+ifn2.NMDL)-einstein_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.MDU = einstein_S(ifn1.MDU+ifn1.NMDU,ifn2.MDU+ifn2.NMDU)-einstein_S(ifn1.NMDU,ifn2.NMDU)
    newIVIFN.NMDL = einstein_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = einstein_S(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN
def Interaction_IVIFN_Einstein_Plus(ifn1,ifn2):
    newIVIFN = IVIFN(None,None,None,None)
    newIVIFN.MDL = einstein_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = einstein_S(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = einstein_S(ifn1.MDL+ifn1.NMDL,ifn2.MDL+ifn2.NMDL)-einstein_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.NMDU = einstein_S(ifn1.MDU+ifn1.NMDU,ifn2.MDU+ifn2.NMDU)-einstein_S(ifn1.MDU,ifn2.MDU)
    return newIVIFN

## Pythagorean Fuzzy Number
## Interaction Algebraic Basic Operations
def Interaction_PFN_Algebraic_Multiply(ifn1,ifn2):
    newPFN = PFN(None,None)
    newPFN.MD = pow(algebraic_S(pow(ifn1.MD,2)+pow(ifn1.NMD,2),pow(ifn2.MD,2)+pow(ifn2.NMD,2))-algebraic_S(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    newPFN.NMD = pow(algebraic_S(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN
def Interaction_PFN_Algebraic_Plus(ifn1,ifn2):
    newPFN = PFN(None,None)
    newPFN.MD = pow(algebraic_S(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(algebraic_S(pow(ifn1.MD,2)+pow(ifn1.NMD,2),pow(ifn2.MD,2)+pow(ifn2.NMD,2))-algebraic_S(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    return newPFN

def Interaction_PFN_Einstein_Multiply(ifn1,ifn2):
    newPFN = PFN(None,None)
    newPFN.MD = pow(einstein_S(pow(ifn1.MD,2)+pow(ifn1.NMD,2),pow(ifn2.MD,2)+pow(ifn2.NMD,2))-einstein_S(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    newPFN.NMD = pow(einstein_S(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN
def Interaction_PFN_Einstein_Plus(ifn1,ifn2):
    newPFN = PFN(None,None)
    newPFN.MD = pow(einstein_S(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(einstein_S(pow(ifn1.MD,2)+pow(ifn1.NMD,2),pow(ifn2.MD,2)+pow(ifn2.NMD,2))-einstein_S(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    return newPFN

## Interval-value Pythagorean Fuzzy Number
## Interaction Algebraic Basic Operations
def Interaction_IVPFN_Algebraic_Multiply(ifn1,ifn2):
    newIVPFN = IVPFN(None,None,None,None)
    newIVPFN.MDL = pow(algebraic_S(pow(ifn1.MDL,2)+pow(ifn1.NMDL,2),pow(ifn2.MDL,2)+pow(ifn2.NMDL,2))-algebraic_S(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.MDU = pow(algebraic_S(pow(ifn1.MDU,2)+pow(ifn1.NMDU,2),pow(ifn2.MDU,2)+pow(ifn2.NMDU,2))-algebraic_S(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    newIVPFN.NMDL = pow(algebraic_S(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(algebraic_S(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN
def Interaction_IVPFN_Algebraic_Plus(ifn1,ifn2):
    newIVPFN = IVPFN(None,None,None,None)
    newIVPFN.MDL = pow(algebraic_S(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(algebraic_S(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(algebraic_S(pow(ifn1.MDL,2)+pow(ifn1.NMDL,2),pow(ifn2.MDL,2)+pow(ifn2.NMDL,2))-algebraic_S(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.NMDU = pow(algebraic_S(pow(ifn1.MDU,2)+pow(ifn1.NMDU,2),pow(ifn2.MDU,2)+pow(ifn2.NMDU,2))-algebraic_S(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    return newIVPFN

def Interaction_IVPFN_Einstein_Multiply(ifn1,ifn2):
    newIVPFN = IVPFN(None,None,None,None)
    newIVPFN.MDL = pow(einstein_S(pow(ifn1.MDL,2)+pow(ifn1.NMDL,2),pow(ifn2.MDL,2)+pow(ifn2.NMDL,2))-einstein_S(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.MDU = pow(einstein_S(pow(ifn1.MDU,2)+pow(ifn1.NMDU,2),pow(ifn2.MDU,2)+pow(ifn2.NMDU,2))-einstein_S(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    newIVPFN.NMDL = pow(einstein_S(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(einstein_S(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN
def Interaction_IVPFN_Einstein_Plus(ifn1,ifn2):
    newIVPFN = IVPFN(None,None,None,None)
    newIVPFN.MDL = pow(einstein_S(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(einstein_S(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(einstein_S(pow(ifn1.MDL,2)+pow(ifn1.NMDL,2),pow(ifn2.MDL,2)+pow(ifn2.NMDL,2))-einstein_S(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.NMDU = pow(einstein_S(pow(ifn1.MDU,2)+pow(ifn1.NMDU,2),pow(ifn2.MDU,2)+pow(ifn2.NMDU,2))-einstein_S(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    return newIVPFN

## Fermatean Fuzzy Number
## Interaction Algebraic Basic Operations
def Interaction_FFN_Algebraic_Multiply(ifn1,ifn2):
    newFFN = FFN(None,None)
    newFFN.MD = pow(algebraic_S(pow(ifn1.MD,3)+pow(ifn1.NMD,3),pow(ifn2.MD,3)+pow(ifn2.NMD,3))-algebraic_S(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    newFFN.NMD = pow(algebraic_S(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN
def Interaction_FFN_Algebraic_Plus(ifn1,ifn2):
    newFFN = FFN(None,None)
    newFFN.MD = pow(algebraic_S(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(algebraic_S(pow(ifn1.MD,3)+pow(ifn1.NMD,3),pow(ifn2.MD,3)+pow(ifn2.NMD,3))-algebraic_S(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    return newFFN

def Interaction_FFN_Einstein_Multiply(ifn1,ifn2):
    newFFN = FFN(None,None)
    newFFN.MD = pow(einstein_S(pow(ifn1.MD,3)+pow(ifn1.NMD,3),pow(ifn2.MD,3)+pow(ifn2.NMD,3))-einstein_S(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    newFFN.NMD = pow(einstein_S(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN
def Interaction_FFN_Einstein_Plus(ifn1,ifn2):
    newFFN = FFN(None,None)
    newFFN.MD = pow(einstein_S(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(einstein_S(pow(ifn1.MD,3)+pow(ifn1.NMD,3),pow(ifn2.MD,3)+pow(ifn2.NMD,3))-einstein_S(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    return newPFN

## Interval-value Fermatean Fuzzy Number
## Interaction Algebraic Basic Operations
def Interaction_IVFFN_Algebraic_Multiply(ifn1,ifn2):
    newIVFFN = IVFFN(None,None,None,None)
    newIVFFN.MDL = pow(algebraic_S(pow(ifn1.MDL,3)+pow(ifn1.NMDL,3),pow(ifn2.MDL,3)+pow(ifn2.NMDL,3))-algebraic_S(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.MDU = pow(algebraic_S(pow(ifn1.MDU,3)+pow(ifn1.NMDU,3),pow(ifn2.MDU,3)+pow(ifn2.NMDU,3))-algebraic_S(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    newIVFFN.NMDL = pow(algebraic_S(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(algebraic_S(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN
def Interaction_IVFFN_Algebraic_Plus(ifn1,ifn2):
    newIVFFN = IVFFN(None,None,None,None)
    newIVFFN.MDL = pow(algebraic_S(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(algebraic_S(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(algebraic_S(pow(ifn1.MDL,3)+pow(ifn1.NMDL,3),pow(ifn2.MDL,3)+pow(ifn2.NMDL,3))-algebraic_S(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.NMDU = pow(algebraic_S(pow(ifn1.MDU,3)+pow(ifn1.NMDU,3),pow(ifn2.MDU,3)+pow(ifn2.NMDU,3))-algebraic_S(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    return newIVFFN

def Interaction_IVFFN_Einstein_Multiply(ifn1,ifn2):
    newIVFFN = IVFFN(None,None,None,None)
    newIVFFN.MDL = pow(einstein_S(pow(ifn1.MDL,3)+pow(ifn1.NMDL,3),pow(ifn2.MDL,3)+pow(ifn2.NMDL,3))-einstein_S(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.MDU = pow(einstein_S(pow(ifn1.MDU,3)+pow(ifn1.NMDU,3),pow(ifn2.MDU,3)+pow(ifn2.NMDU,3))-einstein_S(pow(ifn1.NMDU,3),pow(ifn2.NMDU,2)),1/3)
    newIVFFN.NMDL = pow(einstein_S(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(einstein_S(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN
def Interaction_IVFFN_Einstein_Plus(ifn1,ifn2):
    newIVFFN = IVFFN(None,None,None,None)
    newIVFFN.MDL = pow(einstein_S(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(einstein_S(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(einstein_S(pow(ifn1.MDL,3)+pow(ifn1.NMDL,3),pow(ifn2.MDL,3)+pow(ifn2.NMDL,3))-einstein_S(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.NMDU = pow(einstein_S(pow(ifn1.MDU,3)+pow(ifn1.NMDU,3),pow(ifn2.MDU,3)+pow(ifn2.NMDU,3))-einstein_S(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    return newIVFFN