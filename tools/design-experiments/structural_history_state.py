"""Numerical histories and transactional structural state design prototype.

This JSON round trip is a local binary64-qualified witness, not the product
checkpoint schema or the coupled time integrator's full restart state.
"""
import copy,json,hashlib
import numpy as np

class History:
    def __init__(self,time,value,slopes=None,role='force'):
        self.time=np.array(time,float,copy=True);self.value=np.array(value,float,copy=True);self.slopes=None if slopes is None else np.array(slopes,float,copy=True)
        if self.time.ndim!=1 or len(self.time)<2 or self.value.ndim<1 or self.value.shape[0]!=len(self.time) or not np.isfinite(self.time).all() or not np.isfinite(self.value).all() or np.any(np.diff(self.time)<=0):raise ValueError('invalid_history_samples')
        if self.slopes is not None and (self.slopes.shape!=self.value.shape or not np.isfinite(self.slopes).all()):raise ValueError('invalid_history_slopes')
        if role not in ['force','prescribed_motion']:raise ValueError('invalid_history_role')
        if role=='prescribed_motion' and self.slopes is None:raise ValueError('motion_requires_continuous_velocity')
        for x in [self.time,self.value,self.slopes]:
            if x is not None:x.flags.writeable=False
        self.events=tuple(self.time[1:-1].tolist())
    def evaluate(self,t,derivative=0,side=None):
        if not np.isfinite(t) or t<self.time[0] or t>self.time[-1]:raise ValueError('history_extrapolation')
        if derivative not in [0,1,2] or side not in [None,'left','right']:raise ValueError('invalid_history_query')
        knot=np.flatnonzero(self.time==t)
        if len(knot) and 0<knot[0]<len(self.time)-1 and derivative>(1 if self.slopes is not None else 0) and side is None:raise ValueError('one_sided_derivative_required')
        j=int(np.searchsorted(self.time,t,side='left' if side=='left' else 'right')-1);j=min(max(j,0),len(self.time)-2);h=self.time[j+1]-self.time[j];s=(t-self.time[j])/h;y0=self.value[j];y1=self.value[j+1]
        if self.slopes is None:
            return (1-s)*y0+s*y1 if derivative==0 else (y1-y0)/h if derivative==1 else np.zeros_like(y0)
        a=self.slopes[j];b=self.slopes[j+1]
        if derivative==0:return (2*s**3-3*s*s+1)*y0+(s**3-2*s*s+s)*h*a+(-2*s**3+3*s*s)*y1+(s**3-s*s)*h*b
        if derivative==1:return ((6*s*s-6*s)*y0+(3*s*s-4*s+1)*h*a+(-6*s*s+6*s)*y1+(3*s*s-2*s)*h*b)/h
        return ((12*s-6)*y0+(6*s-4)*h*a+(-12*s+6)*y1+(6*s-2)*h*b)/(h*h)

def plain(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)
def encode(state):return json.dumps(state,sort_keys=True,separators=(',',':'),allow_nan=False,default=plain)

class Store:
    def __init__(self,state,model_identity):
        self._state=copy.deepcopy(state);self.identity=model_identity;self.generation=0;encode(self._state)
    def fingerprint(self):return hashlib.sha256(encode(dict(state=self._state,generation=self.generation,identity=self.identity)).encode()).hexdigest()
    def begin(self):return dict(state=copy.deepcopy(self._state),generation=self.generation,identity=self.identity,qualified=False)
    def evaluate(self,trial,callback):
        if trial['generation']!=self.generation or trial['identity']!=self.identity:raise ValueError('stale_trial')
        trial['qualified']=False;trial.pop('qualified_fingerprint',None);out=callback(trial['state']);trial['qualified_fingerprint']=hashlib.sha256(encode(trial['state']).encode()).hexdigest();trial['qualified']=True;return out
    def accept(self,trial):
        if trial['generation']!=self.generation or trial['identity']!=self.identity or not trial['qualified']:raise ValueError('unqualified_trial')
        if hashlib.sha256(encode(trial['state']).encode()).hexdigest()!=trial['qualified_fingerprint']:raise ValueError('trial_mutated_after_evaluation')
        nextstate=copy.deepcopy(trial['state']);encode(nextstate);self._state=nextstate;self.generation+=1
    def checkpoint(self):return encode(dict(state=self._state,generation=self.generation,identity=self.identity,precision='binary64-qualified',scope='structural-physical-state-prototype'))
    @classmethod
    def restore(cls,payload,model_identity,callback):
        data=json.loads(payload)
        if data.get('identity')!=model_identity or data.get('precision')!='binary64-qualified' or data.get('scope')!='structural-physical-state-prototype':raise ValueError('incompatible_checkpoint')
        state={k:np.array(v) if isinstance(v,list) else v for k,v in data['state'].items()};callback(state);out=cls(state,model_identity);out.generation=data['generation'];return out
