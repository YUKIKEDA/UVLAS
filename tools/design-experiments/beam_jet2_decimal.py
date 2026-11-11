"""Decimal80 second-order AD; precision-preserving pose input, no float cast."""
from decimal import Decimal,localcontext,getcontext
import math
import numpy as np
from beam_decimal_reference import sincos
getcontext().prec=80
N=12
def number(x):return x if isinstance(x,Decimal) else Decimal(str(float(x)))
def zeros(shape):return np.full(shape,Decimal(0),dtype=object)
class Jet:
    def __init__(self,v,g=None,H=None):self.v=number(v);self.g=zeros(N) if g is None else np.asarray(g,dtype=object);self.H=zeros((N,N)) if H is None else np.asarray(H,dtype=object)
    @staticmethod
    def variable(v,i):
        g=zeros(N);g[i]=Decimal(1);return Jet(v,g)
    def __add__(self,x):
        x=asjet(x);return Jet(self.v+x.v,self.g+x.g,self.H+x.H)
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,-self.g,-self.H)
    def __sub__(self,x):return self+-asjet(x)
    def __rsub__(self,x):return asjet(x)+-self
    def __mul__(self,x):
        x=asjet(x);return Jet(self.v*x.v,self.g*x.v+self.v*x.g,self.H*x.v+self.v*x.H+np.outer(self.g,x.g)+np.outer(x.g,self.g))
    __rmul__=__mul__
    def __truediv__(self,x):return self*power(x,-1)
    def __rtruediv__(self,x):return asjet(x)*power(self,-1)
    def __pow__(self,p):return power(self,p)
def asjet(x):return x if isinstance(x,Jet) else Jet(x)
def value(x):return x.v if isinstance(x,Jet) else number(x)
def unary(x,f,d,dd):
    if not isinstance(x,Jet):return f(number(x))
    return Jet(f(x.v),d(x.v)*x.g,d(x.v)*x.H+dd(x.v)*np.outer(x.g,x.g))
def power(x,p):
    p=number(p);return unary(x,lambda t:t**p,lambda t:p*t**(p-1),lambda t:p*(p-1)*t**(p-2))
def sin(x):return unary(x,lambda t:sincos(t)[0],lambda t:sincos(t)[1],lambda t:-sincos(t)[0])
def cos(x):return unary(x,lambda t:sincos(t)[1],lambda t:-sincos(t)[0],lambda t:-sincos(t)[1])
def arc(t):
    x=number(math.acos(float(t)))
    for _ in range(8):
        s,c=sincos(x);x+=(c-t)/s
    return x
def acos(x):return unary(x,arc,lambda t:-1/(1-t*t).sqrt(),lambda t:-t/((1-t*t)**number('1.5')))
def sqrt(x):return power(x,.5)
def polynomial(x,coeff):
    out=Decimal(0)
    for c in coeff[::-1]:out=out*x+number(c)
    return out
