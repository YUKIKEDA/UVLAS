"""First-order real forward differentiation for 13 local geometric variables."""
import math
import numpy as np
N=13
class Dual:
    def __init__(self,v,g=None):self.v=float(v);self.g=np.zeros(N) if g is None else np.asarray(g)
    @staticmethod
    def variable(v,index,scale=1.):
        g=np.zeros(N);g[index]=scale;return Dual(v,g)
    def __add__(self,x):
        x=dual(x);return Dual(self.v+x.v,self.g+x.g)
    __radd__=__add__
    def __neg__(self):return Dual(-self.v,-self.g)
    def __sub__(self,x):return self+-dual(x)
    def __rsub__(self,x):return dual(x)+-self
    def __mul__(self,x):
        x=dual(x);return Dual(self.v*x.v,self.g*x.v+self.v*x.g)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=dual(x);return Dual(self.v/x.v,(self.g-self.v/x.v*x.g)/x.v)
    def __rtruediv__(self,x):return dual(x)/self
def dual(x):return x if isinstance(x,Dual) else Dual(x)
def val(x):return x.v if isinstance(x,Dual) else float(x)
def sqrt(x):
    x=dual(x);r=math.sqrt(x.v);return Dual(r,x.g/(2*r))
def log1p(x):
    x=dual(x);return Dual(math.log1p(x.v),x.g/(1+x.v))
def atan2(y,x):
    y=dual(y);x=dual(x);return Dual(math.atan2(y.v,x.v),(x.v*y.g-y.v*x.g)/(x.v*x.v+y.v*y.v))
def absolute(x):return x if val(x)>=0 else -x
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Dual(0.))
def cross(a,b):return np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],object)
def norm(x):return sqrt(dot(x,x))
def scale(a,x):return np.array([v*x for v in a],object)
