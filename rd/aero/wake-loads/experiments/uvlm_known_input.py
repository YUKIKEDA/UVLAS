"""AE-07 design witness: immutable C1 input coefficients, derivatives, identity."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import hashlib,json,numpy as np
from scipy.interpolate import PchipInterpolator,PPoly

class KnownFlow:
    def __setattr__(self,key,value):
        if getattr(self,'_sealed',False):raise AttributeError('immutable_known_input')
        object.__setattr__(self,key,value)
    def __init__(self,times,values,wave=(0,0,0),direction=(0,0,1),omega=0.,phase=0.):
        t=np.array(times,float,copy=True);v=np.array(values,float,copy=True)
        if t.ndim!=1 or len(t)<2 or v.shape!=(len(t),4) or not np.isfinite(t).all() or not np.isfinite(v).all() or np.any(np.diff(t)<=0):raise ValueError('invalid_input_samples')
        self._initialize(t,PchipInterpolator(t,v,axis=0).c,wave,direction,omega,phase)
    def _initialize(self,t,c,wave,direction,omega,phase):
        self.times=np.array(t,float,copy=True);self.coefficients=np.array(c,float,copy=True);self.wave=np.array(wave,float,copy=True);self.direction=np.array(direction,float,copy=True);self.omega=float(omega);self.phase=float(phase)
        if self.wave.shape!=(3,) or self.direction.shape!=(3,) or not all(np.isfinite(a).all() for a in [self.times,self.coefficients,self.wave,self.direction,self.omega,self.phase]):raise ValueError('invalid_input_parameters')
        if abs(self.wave@self.direction)>64*np.finfo(float).eps*max(np.linalg.norm(self.wave)*np.linalg.norm(self.direction),1e-300):raise ValueError('non_solenoidal_gust')
        self.polynomial=PPoly(self.coefficients,self.times,extrapolate=False);self.derivative=self.polynomial.derivative()
        for a in [self.times,self.coefficients,self.wave,self.direction,self.polynomial.x,self.polynomial.c,self.derivative.x,self.derivative.c]:a.flags.writeable=False
        self.identity=hashlib.sha256(self.encode().encode()).hexdigest();self._sealed=True
    def encode(self):
        return json.dumps(dict(version='ae07-c1-pchip-coeff-v1',times=self.times.tolist(),coefficients=self.coefficients.tolist(),wave=self.wave.tolist(),direction=self.direction.tolist(),omega=self.omega,phase=self.phase,range_policy='reject'),sort_keys=True,separators=(',',':'))
    @classmethod
    def restore(cls,encoded):
        d=json.loads(encoded)
        if d.get('version')!='ae07-c1-pchip-coeff-v1' or d.get('range_policy')!='reject':raise ValueError('unsupported_input_policy')
        t=np.asarray(d['times']);c=np.asarray(d['coefficients'])
        if t.ndim!=1 or len(t)<2 or np.any(np.diff(t)<=0) or c.shape!=(4,len(t)-1,4):raise ValueError('invalid_saved_input')
        # A saved polynomial must also be C1 at every interior knot.
        for j in range(len(t)-2):
            h=t[j+1]-t[j];left=((c[0,j]*h+c[1,j])*h+c[2,j])*h+c[3,j];slope=(3*c[0,j]*h+2*c[1,j])*h+c[2,j]
            if not np.allclose(left,c[3,j+1],rtol=1e-12,atol=1e-12) or not np.allclose(slope,c[2,j+1],rtol=1e-12,atol=1e-12):raise ValueError('discontinuous_saved_input')
        obj=cls.__new__(cls);obj._initialize(t,c,d['wave'],d['direction'],d['omega'],d['phase']);return obj
    def evaluate(self,p,t):
        p=np.asarray(p,float)
        if p.shape[-1:]!=(3,) or not np.isfinite(p).all() or not np.isfinite(t) or t<self.times[0] or t>self.times[-1]:raise ValueError('input_outside_domain')
        v=self.polynomial(t);d=self.derivative(t);phase=p@self.wave-self.omega*t+self.phase;s=np.sin(phase);c=np.cos(phase);U=v[:3]+v[3]*s[...,None]*self.direction;J=v[3]*c[...,None,None]*np.outer(self.direction,self.wave);dt=d[:3]+(d[3]*s-self.omega*v[3]*c)[...,None]*self.direction;return U,J,dt
    def value(self,p,t):return self.evaluate(p,t)[0]
    def jvp(self,p,t,dp,dU=None):return np.einsum('...ij,...j->...i',self.evaluate(p,t)[1],dp)+(0 if dU is None else dU)
