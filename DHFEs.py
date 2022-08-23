# from InteractionFNs import *
# from FNs import *
import math
from math import pow
import Archimedean
from Archimedean import *

## Dual Hesitant Intuitionistic Fuzzy Elements
class DHIFE:
    def __init__(self,MD,NMD):
        try:
            assert (not MD and not NMD) or (max(MD)<=1 and max(NMD)<=1 and min(MD)>=0 and min(NMD)>=0) and (0<=max(MD)+max(NMD)<=1 and 0<=min(MD)+min(NMD)<=1)
            self.MD = MD
            self.NMD = NMD
        except AssertionError as er:
            print("ERROE:Construction failed!\n"+
                  "DHIFE max(MD)+max(NMD) and min(MD)+min(NMD) must be in interval[0,1]!\n"+
                  "It will return a variable of type 'NoneType'!",er)
            self = None
    
    def show(self):
        print('DHIFE:{\n' + 'MD: '+str(self.MD) + ',\n' + 'NMD:' + str(self.NMD) +'}')
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(in_algebraic_tau(l*algebraic_tau(md)))
        for nmd in self.NMD:
            newDHIFE.NMD.append(in_algebraic_s(l*algebraic_s(nmd)))
        return newDHIFE
    def Algebraic_Times(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(in_algebraic_s(l*algebraic_s(md)))
        for nmd in self.NMD:
            newDHIFE.NMD.append(in_algebraic_tau(l*algebraic_tau(nmd)))
        return newDHIFE
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(in_einstein_tau(l*einstein_tau(md)))
        for nmd in self.NMD:
            newDHIFE.NMD.append(in_einstein_s(l*einstein_s(nmd)))
        return newDHIFE
    def Einstein_Times(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(in_einstein_s(l*einstein_s(md)))
        for nmd in self.NMD:
            newDHIFE.NMD.append(in_einstein_tau(l*einstein_tau(nmd)))
        return newDHIFE
    
    ## Fast Algebraic Basic Operations
    def Fast_Algebraic_Power(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(pow(md,l))
        for nmd in self.NMD:
            newDHIFE.NMD.append(1-pow(1-nmd,l))
        return newDHIFE

    def Fast_Algebraic_Times(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append(1-pow(1-md,l))
        for nmd in self.NMD:
            newDHIFE.NMD.append(pow(nmd,l))
        return newDHIFE
    
    ## Fast Einstein Basic Operations
    def Fast_Einstein_Power(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append((2*pow(md,l))/(pow(2-md,l)+pow(md,l)))
        for nmd in self.NMD:
            newDHIFE.NMD.append((pow(1+nmd,l)-pow(1-nmd,l))/(pow(1+nmd,l)+pow(1-nmd,l)))
        return newDHIFE
    def Fast_Einstein_Times(self,l):
        newDHIFE = DHIFE([],[])
        for md in self.MD:
            newDHIFE.MD.append((pow(1+md,l)-pow(1-md,l))/(pow(1+md,l)+pow(1-md,l)))
        for nmd in self.NMD:
            newDHIFE.NMD.append((2*pow(nmd,l))/(pow(2-nmd,l)+pow(nmd,l)))
        return newDHIFE
    

    
## Dual Hesitant Pythagorean Fuzzy Elements
class DHPFE:
    def __init__(self,MD,NMD):
        try:
            assert (not MD and not NMD) or (max(MD)<=1 and max(NMD)<=1 and min(MD)>=0 and min(NMD)>=0) and (0<=max(MD)**2+max(NMD)**2<=1 and 0<=min(MD)**2+min(NMD)**2<=1)
            self.MD = MD
            self.NMD = NMD
        except AssertionError as er:
            print("ERROE:Construction failed!\n"+
                  "DHPFE max(MD)^2+max(NMD)^2 and min(MD)^2+min(NMD)^2 must be in interval[0,1]!\n"+
                  "It will return a variable of type 'NoneType'!",er)
            self = None
    
    def show(self):
        print('DHPFE:{\n' + 'MD: '+str(self.MD) + ',\n' + 'NMD:' + str(self.NMD) +'}')
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(in_algebraic_tau(l*algebraic_tau(pow(md,2))),1/2))
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(in_algebraic_s(l*algebraic_s(pow(nmd,2))),1/2))
        return newDHPFE
    def Algebraic_Times(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(in_algebraic_s(l*algebraic_s(pow(md,2))),1/2))
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(in_algebraic_tau(l*algebraic_tau(pow(nmd,2))),1/2))
        return newDHPFE
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(in_einstein_tau(l*einstein_tau(pow(md,2))),1/2))
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(in_einstein_s(l*einstein_s(pow(nmd,2))),1/2))
        return newDHPFE
    def Einstein_Times(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(in_einstein_s(l*einstein_s(pow(md,2))),1/2))
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(in_einstein_tau(l*einstein_tau(pow(nmd,2))),1/2))
        return newDHPFE
    
    ## Fast Algebraic Basic Operations
    def Fast_Algebraic_Power(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(md,l))
        for nmd in self.NMD:
            newDHPFE.NMD.append(
                pow(1-pow(1-nmd**2,l),1/2)
            )
        return newDHPFE
    def Fast_Algebraic_Times(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(
                pow(1-pow(1-md**2,l),1/2)
            )
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(nmd,l))
        return newDHPFE
    
    ## Fast Einstein Basic Operations
    def Fast_Einstein_Power(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(
                pow((2*pow(md**2,l))/(pow(2-md**2,l)+pow(md**2,l)),1/2)
            )
        for nmd in self.NMD:
            newDHPFE.NMD.append(pow(
                (pow(1+nmd**2,l)-pow(1-nmd**2,l))/(pow(1+nmd**2,l)+pow(1-nmd**2,l)),1/2)
            )
        return newDHPFE
    def Fast_Einstein_Times(self,l):
        newDHPFE = DHPFE([],[])
        for md in self.MD:
            newDHPFE.MD.append(pow(
                (pow(1+md**2,l)-pow(1-md**2,l))/(pow(1+md**2,l)+pow(1-md**2,l)),1/2)
            )
        for nmd in self.NMD:
            newDHPFE.NMD.append(
                pow((2*pow(nmd**2,l))/(pow(2-nmd**2,l)+pow(nmd**2,l)),1/2)
            )
        return newDHPFE
    
## Dual Hesitant Fermatean Fuzzy Elements
class DHFFE:
    def __init__(self,MD,NMD):
        try:
            assert (not MD and not NMD) or (max(MD)<=1 and max(NMD)<=1 and min(MD)>=0 and min(NMD)>=0) and (0<=max(MD)**3+max(NMD)**3<=1 and 0<=min(MD)**3+min(NMD)**3<=1)
            self.MD = MD
            self.NMD = NMD
        except AssertionError as er:
            print("ERROE:Construction failed!\n"+
                  "DHFFE max(MD)^3+max(NMD)^3 and min(MD)^3+min(NMD)^3 must be in interval[0,1]!\n"+
                  "It will return a variable of type 'NoneType'!",er)
            self = None
    
    def show(self):
        print('DHFFE:{\n' + 'MD: '+str(self.MD) + ',\n' + 'NMD:' + str(self.NMD) +'}')
        
    ## Algebraic Basic Operations
    def Algebraic_Power(self,l):
        newDHFFE = DHFFE([],[])
        for md in self.MD:
            newDHFFE.MD.append(pow(in_algebraic_tau(l*algebraic_tau(pow(md,3))),1/3))
        for nmd in self.NMD:
            newDHFFE.NMD.append(pow(in_algebraic_s(l*algebraic_s(pow(nmd,3))),1/3))
        return newDHFFE
    def Algebraic_Times(self,l):
        newDHFFE = DHFFE([],[])
        for md in self.MD:
            newDHFFE.MD.append(pow(in_algebraic_s(l*algebraic_s(pow(md,3))),1/3))
        for nmd in self.NMD:
            newDHFFE.NMD.append(pow(in_algebraic_tau(l*algebraic_tau(pow(nmd,3))),1/3))
        return newDHFFE
    
    ## Einstein Basic Operations
    def Einstein_Power(self,l):
        newDHFFE = DHFFE([],[])
        for md in self.MD:
            newDHFFE.MD.append(pow(in_einstein_tau(l*einstein_tau(pow(md,3))),1/3))
        for nmd in self.NMD:
            newDHFFE.NMD.append(pow(in_einstein_s(l*einstein_s(pow(nmd,3))),1/3))
        return newDHFFE
    def Einstein_Times(self,l):
        newDHFFE = DHFFE([],[])
        for md in self.MD:
            newDHFFE.MD.append(pow(in_einstein_s(l*einstein_s(pow(md,3))),1/3))
        for nmd in self.NMD:
            newDHFFE.NMD.append(pow(in_einstein_tau(l*einstein_tau(pow(nmd,3))),1/3))
        return newDHFFE
    
    ## Fast Algebraic Basic Operations
    def Fast_Algebraic_Power(self,l):
        newDHFFE = DHFFE([],[])
        for md in self.MD:
            newDHFFE.MD.append(pow(md,l))
        for nmd in self.NMD:
            newDHFFE.NMD.append(
                pow(1-pow(1-nmd**3,l),1/3)
            )
        return newDHFFE
    def Fast_Algebraic_Times(self,l):
        newDHFFE = DHFFE([],[])
        for md in self.MD:
            newDHFFE.MD.append(
                pow(1-pow(1-md**3,l),1/3)
            )
        for nmd in self.NMD:
            newDHFFE.NMD.append(pow(nmd,l))
        return newDHFFE
    
    ## Fast Einstein Basic Operations
    def Fast_Einstein_Power(self,l):
        newDHFFE = DHFFE([],[])
        for md in self.MD:
            newDHFFE.MD.append(
                pow((2*pow(md**3,l))/(pow(2-md**3,l)+pow(md**3,l)),1/3)
            )
        for nmd in self.NMD:
            newDHFFE.NMD.append(pow(
                (pow(1+nmd**3,l)-pow(1-nmd**3,l))/(pow(1+nmd**3,l)+pow(1-nmd**3,l)),1/3)
            )
        return newDHFFE
    def Fast_Einstein_Times(self,l):
        newDHFFE = DHFFE([],[])
        for md in self.MD:
            newDHFFE.MD.append(pow(
                (pow(1+md**3,l)-pow(1-md**3,l))/(pow(1+md**3,l)+pow(1-md**3,l)),1/3)
            )
        for nmd in self.NMD:
            newDHFFE.NMD.append(pow(
                (2*pow(nmd**3,l))/(pow(2-nmd**3,l)+pow(nmd**3,l)),1/3)
            )
        return newDHFFE
    
##########  Basic multiplication and addition operations
## Basic multiplication and addition operations
## Dual Hesitant Intuitionistic Fuzzy Elements
## Algebraic Basic Operations 
def DHIFE_Algebraic_Multiply(dhfe1,dhfe2):
    newDHIFE = DHIFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHIFE.MD.append(pithy_algebraic_T(md1,md2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHIFE.NMD.append(pithy_algebraic_S(nmd1,nmd2))
    return newDHIFE
def DHIFE_Algebraic_Plus(dhfe1,dhfe2):
    newDHIFE = DHIFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHIFE.MD.append(pithy_algebraic_S(md1,md2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHIFE.NMD.append(pithy_algebraic_T(nmd1,nmd2))
    return newDHIFE

## Einstein Basic Operations
def DHIFE_Einstein_Multiply(dhfe1,dhfe2):
    newDHIFE = DHIFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHIFE.MD.append(pithy_einstein_T(md1,md2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHIFE.NMD.append(pithy_einstein_S(nmd1,nmd2))
    return newDHIFE
def DHIFE_Einstein_Plus(dhfe1,dhfe2):
    newDHIFE = DHIFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHIFE.MD.append(pithy_einstein_S(md1,md2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHIFE.NMD.append(pithy_einstein_T(nmd1,nmd2))
    return newDHIFE

## Basic multiplication and addition operations
## Dual Hesitant Pythagorean Fuzzy Elements
## Algebraic Basic Operations 
def DHPFE_Algebraic_Multiply(dhfe1,dhfe2):
    newDHPFE = DHPFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHPFE.MD.append(pow(pithy_algebraic_T(pow(md1,2),pow(md2,2)),1/2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHPFE.NMD.append(pow(pithy_algebraic_S(pow(nmd1,2),pow(nmd2,2)),1/2))
    return newDHPFE
def DHPFE_Algebraic_Plus(dhfe1,dhfe2):
    newDHPFE = DHPFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHPFE.MD.append(pow(pithy_algebraic_S(pow(md1,2),pow(md2,2)),1/2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHPFE.NMD.append(pow(pithy_algebraic_T(pow(nmd1,2),pow(nmd2,2)),1/2))
    return newDHPFE

## Einstein Basic Operations
def DHPFE_Einstein_Multiply(dhfe1,dhfe2):
    newDHPFE = DHPFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHPFE.MD.append(pow(pithy_einstein_T(pow(md1,2),pow(md2,2)),1/2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHPFE.NMD.append(pow(pithy_einstein_S(pow(nmd1,2),pow(nmd2,2)),1/2))
    return newDHPFE
def DHPFE_Einstein_Plus(dhfe1,dhfe2):
    newDHPFE = DHPFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHPFE.MD.append(pow(pithy_einstein_S(pow(md1,2),pow(md2,2)),1/2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHPFE.NMD.append(pow(pithy_einstein_T(pow(nmd1,2),pow(nmd2,2)),1/2))
    return newDHPFE

## Basic multiplication and addition operations
## Dual Hesitant Fermatean Fuzzy Elements
## Algebraic Basic Operations 
def DHFFE_Algebraic_Multiply(dhfe1,dhfe2):
    newDHFFE = DHFFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHFFE.MD.append(pow(pithy_algebraic_T(pow(md1,3),pow(md2,3)),1/3))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHFFE.NMD.append(pow(pithy_algebraic_S(pow(nmd1,3),pow(nmd2,3)),1/3))
    return newDHFFE
def DHFFE_Algebraic_Plus(dhfe1,dhfe2):
    newDHFFE = DHFFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHFFE.MD.append(pow(pithy_algebraic_S(pow(md1,3),pow(md2,3)),1/3))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHFFE.NMD.append(pow(pithy_algebraic_T(pow(nmd1,3),pow(nmd2,3)),1/3))
    return newDHFFE

## Einstein Basic Operations
def DHFFE_Einstein_Multiply(dhfe1,dhfe2):
    newDHFFE = DHFFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHFFE.MD.append(pow(pithy_einstein_T(pow(md1,3),pow(md2,3)),1/3))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHFFE.NMD.append(pow(pithy_einstein_S(pow(nmd1,3),pow(nmd2,3)),1/3))
    return newDHFFE
def DHFFE_Einstein_Plus(dhfe1,dhfe2):
    newDHFFE = DHFFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHFFE.MD.append(pow(pithy_einstein_S(pow(md1,3),pow(md2,3)),1/3))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHFFE.NMD.append(pow(pithy_einstein_T(pow(nmd1,3),pow(nmd2,3)),1/3))
    return newDHFFE