from math import pow
from FNS.Archimedean import *
from FNS.FNs import *


##########  Basic multiplication and addition operations

## Basic multiplication and addition operations
## Intuitionistic Fuzzy Number
## Algebraic Basic Operations 
def IFN_Algebraic_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIFN = IFN(0,0)
    newIFN.MD = pithy_algebraic_T(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_algebraic_S(ifn1.NMD,ifn2.NMD)
    return newIFN
def IFN_Algebraic_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIFN = IFN(0,0)
    newIFN.MD = pithy_algebraic_S(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_algebraic_T(ifn1.NMD,ifn2.NMD)
    return newIFN

## Einstein Basic Operations
def IFN_Einstein_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIFN = IFN(0,0)
    newIFN.MD = pithy_einstein_T(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_einstein_S(ifn1.NMD,ifn2.NMD)
    return newIFN
def IFN_Einstein_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIFN = IFN(0,0)
    newIFN.MD = pithy_einstein_S(ifn1.MD,ifn2.MD)
    newIFN.NMD = pithy_einstein_T(ifn1.NMD,ifn2.NMD)
    return newIFN

## Basic multiplication and addition operations
## Interval-Value Intuitionistic Fuzzy Number
## Algebraic Basic Operations 
def IVIFN_Algebraic_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVIFN = IVIFN(0,0,0,0)
    newIVIFN.MDL = pithy_algebraic_T(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_algebraic_T(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_algebraic_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_algebraic_S(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN
def IVIFN_Algebraic_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVIFN = IVIFN(0,0,0,0)
    newIVIFN.MDL = pithy_algebraic_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_algebraic_S(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_algebraic_T(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_algebraic_T(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN

## Einstein Basic Operations
def IVIFN_Einstein_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVIFN = IVIFN(0,0,0,0)
    newIVIFN.MDL = pithy_einstein_T(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_einstein_T(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_einstein_S(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_einstein_S(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN
def IVIFN_Einstein_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVIFN = IVIFN(0,0,0,0)
    newIVIFN.MDL = pithy_einstein_S(ifn1.MDL,ifn2.MDL)
    newIVIFN.MDU = pithy_einstein_S(ifn1.MDU,ifn2.MDU)
    newIVIFN.NMDL = pithy_einstein_T(ifn1.NMDL,ifn2.NMDL)
    newIVIFN.NMDU = pithy_einstein_T(ifn1.NMDU,ifn2.NMDU)
    return newIVIFN

## Basic multiplication and addition operations
## Pythagorean Fuzzy Number
## Algebraic Basic Operations 
def PFN_Algebraic_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newPFN = PFN(0,0)
    newPFN.MD = pow(pithy_algebraic_T(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(pithy_algebraic_S(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN
def PFN_Algebraic_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newPFN = PFN(0,0)
    newPFN.MD = pow(pithy_algebraic_S(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(pithy_algebraic_T(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN

## Einstein Basic Operations
def PFN_Einstein_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newPFN = PFN(0,0)
    newPFN.MD = pow(pithy_einstein_T(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(pithy_einstein_S(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN
def PFN_Einstein_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newPFN = PFN(0,0)
    newPFN.MD = pow(pithy_einstein_S(pow(ifn1.MD,2),pow(ifn2.MD,2)),1/2)
    newPFN.NMD = pow(pithy_einstein_T(pow(ifn1.NMD,2),pow(ifn2.NMD,2)),1/2)
    return newPFN

## Basic multiplication and addition operations
## Interval-Value Pythagorean Fuzzy Number
## Algebraic Basic Operations 
def IVPFN_Algebraic_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVPFN = IVPFN(0,0,0,0)
    newIVPFN.MDL = pow(pithy_algebraic_T(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(pithy_algebraic_T(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(pithy_algebraic_S(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(pithy_algebraic_S(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN
def IVPFN_Algebraic_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVPFN = IVPFN(0,0,0,0)
    newIVPFN.MDL = pow(pithy_algebraic_S(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(pithy_algebraic_S(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(pithy_algebraic_T(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(pithy_algebraic_T(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN

## Einstein Basic Operations
def IVPFN_Einstein_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVPFN = IVPFN(0,0,0,0)
    newIVPFN.MDL = pow(pithy_einstein_T(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(pithy_einstein_T(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(pithy_einstein_S(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(pithy_einstein_S(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN
def IVPFN_Einstein_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVPFN = IVPFN(0,0,0,0)
    newIVPFN.MDL = pow(pithy_einstein_S(pow(ifn1.MDL,2),pow(ifn2.MDL,2)),1/2)
    newIVPFN.MDU = pow(pithy_einstein_S(pow(ifn1.MDU,2),pow(ifn2.MDU,2)),1/2)
    newIVPFN.NMDL = pow(pithy_einstein_T(pow(ifn1.NMDL,2),pow(ifn2.NMDL,2)),1/2)
    newIVPFN.NMDU = pow(pithy_einstein_T(pow(ifn1.NMDU,2),pow(ifn2.NMDU,2)),1/2)
    return newIVPFN

## Basic multiplication and addition operations
## Fermatean Fuzzy Number
## Algebraic Basic Operations 
def FFN_Algebraic_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newFFN = FFN(0,0)
    newFFN.MD = pow(pithy_algebraic_T(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(pithy_algebraic_S(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN
def FFN_Algebraic_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newFFN = FFN(0,0)
    newFFN.MD = pow(pithy_algebraic_S(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(pithy_algebraic_T(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN

## Einstein Basic Operations
def FFN_Einstein_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newFFN = FFN(0,0)
    newFFN.MD = pow(pithy_einstein_T(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(pithy_einstein_S(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN
def FFN_Einstein_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newFFN = FFN(0,0)
    newFFN.MD = pow(pithy_einstein_S(pow(ifn1.MD,3),pow(ifn2.MD,3)),1/3)
    newFFN.NMD = pow(pithy_einstein_T(pow(ifn1.NMD,3),pow(ifn2.NMD,3)),1/3)
    return newFFN

## Basic multiplication and addition operations
## Interval-Value Fermatean Fuzzy Number
## Algebraic Basic Operations 
def IVFFN_Algebraic_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVFFN = IVFFN(0,0,0,0)
    newIVFFN.MDL = pow(pithy_algebraic_T(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(pithy_algebraic_T(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(pithy_algebraic_S(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(pithy_algebraic_S(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN
def IVFFN_Algebraic_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVFFN = IVFFN(0,0,0,0)
    newIVFFN.MDL = pow(pithy_algebraic_S(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(pithy_algebraic_S(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(pithy_algebraic_T(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(pithy_algebraic_T(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN

## Einstein Basic Operations
def IVFFN_Einstein_Multiply(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVFFN = IVFFN(0,0,0,0)
    newIVFFN.MDL = pow(pithy_einstein_T(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(pithy_einstein_T(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(pithy_einstein_S(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(pithy_einstein_S(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN
def IVFFN_Einstein_Plus(ifn1,ifn2):
    try:
        assert ifn1.getQRung() == ifn2.getQRung()           ## 判断是否是同一种模糊集
    except AssertionError as es:
        print('ERROR! The two FNs are not the same FN!',es)
        return None
    
    newIVFFN = IVFFN(0,0,0,0)
    newIVFFN.MDL = pow(pithy_einstein_S(pow(ifn1.MDL,3),pow(ifn2.MDL,3)),1/3)
    newIVFFN.MDU = pow(pithy_einstein_S(pow(ifn1.MDU,3),pow(ifn2.MDU,3)),1/3)
    newIVFFN.NMDL = pow(pithy_einstein_T(pow(ifn1.NMDL,3),pow(ifn2.NMDL,3)),1/3)
    newIVFFN.NMDU = pow(pithy_einstein_T(pow(ifn1.NMDU,3),pow(ifn2.NMDU,3)),1/3)
    return newIVFFN