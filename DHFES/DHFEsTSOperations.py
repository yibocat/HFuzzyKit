import math
import DHFES.DHFEs
import Archimedean

from math import pow
from Archimedean import *
from DHFES.DHFEs import *

##########  Basic multiplication and addition operations
## Basic multiplication and addition operations
## Dual Hesitant Intuitionistic Fuzzy Elements
## Algebraic Basic Operations 
def DHIFE_Algebraic_Multiply(dhfe1,dhfe2):
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==1           ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHIFEs!',es)
        return None
    
    newDHIFE = DHIFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHIFE.MD.append(pithy_algebraic_T(md1,md2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHIFE.NMD.append(pithy_algebraic_S(nmd1,nmd2))
    return newDHIFE
def DHIFE_Algebraic_Plus(dhfe1,dhfe2):
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==1          ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHIFEs!',es)
        return None

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
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==1         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHIFEs!',es)
        return None
    
    newDHIFE = DHIFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHIFE.MD.append(pithy_einstein_T(md1,md2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHIFE.NMD.append(pithy_einstein_S(nmd1,nmd2))
    return newDHIFE
def DHIFE_Einstein_Plus(dhfe1,dhfe2):
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==1         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHIFEs!',es)
        return None
    
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
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==2         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!',es)
        return None
    
    newDHPFE = DHPFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHPFE.MD.append(pow(pithy_algebraic_T(pow(md1,2),pow(md2,2)),1/2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHPFE.NMD.append(pow(pithy_algebraic_S(pow(nmd1,2),pow(nmd2,2)),1/2))
    return newDHPFE
def DHPFE_Algebraic_Plus(dhfe1,dhfe2):
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==2         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!',es)
        return None
    
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
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==2         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!',es)
        return None
    
    newDHPFE = DHPFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHPFE.MD.append(pow(pithy_einstein_T(pow(md1,2),pow(md2,2)),1/2))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHPFE.NMD.append(pow(pithy_einstein_S(pow(nmd1,2),pow(nmd2,2)),1/2))
    return newDHPFE
def DHPFE_Einstein_Plus(dhfe1,dhfe2):
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==2         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHPFEs!',es)
        return None
    
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
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==3         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHFFEs!',es)
        return None
    
    newDHFFE = DHFFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHFFE.MD.append(pow(pithy_algebraic_T(pow(md1,3),pow(md2,3)),1/3))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHFFE.NMD.append(pow(pithy_algebraic_S(pow(nmd1,3),pow(nmd2,3)),1/3))
    return newDHFFE
def DHFFE_Algebraic_Plus(dhfe1,dhfe2):
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==3         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHFFEs!',es)
        return None
    
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
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==3         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHFFEs!',es)
        return None
    
    newDHFFE = DHFFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHFFE.MD.append(pow(pithy_einstein_T(pow(md1,3),pow(md2,3)),1/3))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHFFE.NMD.append(pow(pithy_einstein_S(pow(nmd1,3),pow(nmd2,3)),1/3))
    return newDHFFE
def DHFFE_Einstein_Plus(dhfe1,dhfe2):
    try:
        assert dhfe1.getQRung() == dhfe2.getQRung()  and dhfe1.getQRung()==3         ## 判断是否是同一种对偶犹豫模糊集
    except AssertionError as es:
        print('ERROR! The two DHFEs are not the same DHFE or they may be not DHFFEs!',es)
        return None
    
    newDHFFE = DHFFE([],[])
    for md1 in dhfe1.MD:
        for md2 in dhfe2.MD:
            newDHFFE.MD.append(pow(pithy_einstein_S(pow(md1,3),pow(md2,3)),1/3))
    for nmd1 in dhfe1.NMD:
        for nmd2 in dhfe2.NMD:
            newDHFFE.NMD.append(pow(pithy_einstein_T(pow(nmd1,3),pow(nmd2,3)),1/3))
    return newDHFFE
