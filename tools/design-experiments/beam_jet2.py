"""Small second-order forward AD for an independently implemented beam probe."""
import math
import numpy as np
N=12
class Jet:
    def __init__(self,v,g=None,H=None):
        self.v=float(v);self.g=np.zeros(N) if g is None else np.asarray(g);self.H=np.zeros((N,N)) if H is None else np.asarray(H)
    @staticmethod
    def variable(v,i):
        g=np.zeros(N);g[i]=1;return Jet(v,g)
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
def value(x):return x.v if isinstance(x,Jet) else float(x)
def unary(x,f,d,dd):
    if not isinstance(x,Jet):return f(x)
    return Jet(f(x.v),d(x.v)*x.g,d(x.v)*x.H+dd(x.v)*np.outer(x.g,x.g))
def power(x,p):return unary(x,lambda t:t**p,lambda t:p*t**(p-1),lambda t:p*(p-1)*t**(p-2))
def sin(x):return unary(x,math.sin,math.cos,lambda t:-math.sin(t))
def cos(x):return unary(x,math.cos,lambda t:-math.sin(t),lambda t:-math.cos(t))
def acos(x):return unary(x,math.acos,lambda t:-1/math.sqrt(1-t*t),lambda t:-t/(1-t*t)**1.5)
def sqrt(x):return power(x,.5)
def polynomial(x,coeff):
    out=0.
    for c in coeff[::-1]:out=out*x+c
    return out
