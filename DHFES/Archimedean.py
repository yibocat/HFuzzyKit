import numpy as np
np.set_printoptions(suppress=True) ## 非科学计数法
np.set_printoptions(precision=4)

## Algebraic t-norm & t-conorm
def algebraic_tau(x):
    return -np.log2(x)
def in_algebraic_tau(x):
    return 1/(2**x)
def algebraic_s(x):
    return algebraic_tau(1-x)
def in_algebraic_s(x):
    return 1-in_algebraic_tau(x)

def algebraic_T(x,y):
    return in_algebraic_tau(algebraic_tau(x)+algebraic_tau(y))
def algebraic_S(x,y):
    return in_algebraic_s(algebraic_s(x)+algebraic_s(y))

def pithy_algebraic_T(x,y):
    # the pithy T function
    return x*y
def pithy_algebraic_S(x,y):
    # the pithy S function
    return x+y-x*y

## Einstein t-norm & t-conorm
def einstein_tau(x):
    return np.log2((2-x)/x)
def in_einstein_tau(x):
    return 2/((2**x)+1)
def einstein_s(x):
    return np.log2((1+x)/(1-x))
def in_einstein_s(x):
    return (2**x-1)/(2**x+1)

def einstein_T(x,y):
    return in_einstein_tau(einstein_tau(x)+einstein_tau(y))
def einstein_S(x,y):
    return in_einstein_s(einstein_s(x)+einstein_s(y))

def pithy_einstein_T(x,y):
    # the pithy einstein T function
    return (x*y)/(1+(1-x)*(1-y))
def pithy_einstein_S(x,y):
    # the pithy einstein S function
    return (x+y)/(1+x*y)


