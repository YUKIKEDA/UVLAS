"""Search-control oracle with a fully known two-mode spectrum; not aero physics."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import json,hashlib,copy,numpy as np
from scipy.linalg import expm,eig
from scipy.optimize import linear_sum_assignment
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)
class Search:
    def __init__(self):self.state=dict(version='two-mode-search-v1',h=.02,uncertainty=.0001,scan=[6.,9.,12.,15.,18.],completed={},queue=[],brackets=None,finished=False,actions=0,evaluations=0)
    def evaluate(self,U):
        key=float(U).hex();s=self.state
        if key in s['completed']:return s['completed'][key]
        theta=.03*U;R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]);growth=np.array([.03*(U-10),.05*(U-14)]);freq=np.array([1+.03*U,2-.07*U]);omega=2*np.pi*freq;A=np.block([[np.zeros((2,2)),np.eye(2)],[-R@np.diag(growth**2+omega**2)@R.T,R@np.diag(2*growth)@R.T]]);P=expm(s['h']*A);z,V=eig(P);ids=np.flatnonzero(z.imag>0);assert len(z)==4 and len(ids)==2;zz=z[ids];VV=V[:,ids];res=np.linalg.norm(P@VV-VV*zz,axis=0)/(np.linalg.norm(P)*np.linalg.norm(VV,axis=0));assert max(res)<1e-12;Q=VV[:2]/np.linalg.norm(VV[:2],axis=0)
        if s['completed']:
            near=min(s['completed'].values(),key=lambda p:(abs(p['U']-U),p['U']));raw=np.array(near['modes']);old=raw[...,0]+1j*raw[...,1];MAC=abs(old.conj().T@Q)**2;a,b=linear_sum_assignment(-MAC);order=b[np.argsort(a)];Q=Q[:,order];zz=zz[order]
            for j in range(2):phase=np.vdot(old[:,j],Q[:,j]);Q[:,j]*=phase.conjugate()/abs(phase)
        else:
            order=np.argsort(np.angle(zz));Q=Q[:,order];zz=zz[order]
            for j in range(2):k=np.argmax(abs(Q[:,j]));Q[:,j]*=Q[k,j].conjugate()/abs(Q[k,j])
        point=dict(U=float(U),growth=(np.log(abs(zz))/s['h']).tolist(),frequency=(np.angle(zz)/(2*np.pi*s['h'])).tolist(),modes=np.stack([Q.real,Q.imag],axis=-1).tolist(),original_eigen_residuals=res.tolist());s['completed'][key]=point;s['evaluations']+=1;return point
    def sign(self,p,j):
        g=p['growth'][j];u=self.state['uncertainty'];return -1 if g+u<0 else 1 if g-u>0 else 0
    def step(self):
        s=self.state
        if s['finished']:return
        if s['scan']:
            self.evaluate(s['scan'].pop(0))
            if not s['scan']:
                points=sorted(s['completed'].values(),key=lambda p:p['U']);brackets=[]
                for j in range(2):
                    pairs=[(a,b) for a,b in zip(points[:-1],points[1:]) if self.sign(a,j)<0 and self.sign(b,j)>0];assert len(pairs)==1;lo,hi=pairs[0];brackets.append(dict(mode=j,lo=lo['U'],hi=hi['U'],done=False))
                s['brackets']=brackets
        elif s['queue']:
            job=s['queue'].pop(0);j=job['mode'];b=s['brackets'][j];U=job['U']
            if b['lo']<U<b['hi']:
                point=self.evaluate(U);sign=self.sign(point,j)
                if sign<0:b['lo']=U
                elif sign>0:b['hi']=U
                else:s['queue'].extend([dict(mode=j,U=b['lo']+(b['hi']-b['lo'])/4),dict(mode=j,U=b['hi']-(b['hi']-b['lo'])/4)])
        else:
            for j,b in enumerate(s['brackets']):
                if b['done']:continue
                if b['hi']-b['lo']<=.04:
                    mid=(b['lo']+b['hi'])/2;point=self.evaluate(mid);b.update(done=True,estimate=mid,frequency=point['frequency'][j]);assert self.sign(self.evaluate(b['lo']),j)<0 and self.sign(self.evaluate(b['hi']),j)>0
                else:s['queue'].append(dict(mode=j,U=(b['lo']+b['hi'])/2))
                break
            s['finished']=all(b['done'] for b in s['brackets'])
        s['actions']+=1
        if s['actions']>500:raise ValueError('search_work_limit')
    def run(self):
        while not self.state['finished']:self.step()
        return self.state
    def encode(self):return canonical(dict(state=self.state,checksum=hashlib.sha256(canonical(self.state).encode()).hexdigest()))
    @classmethod
    def restore(cls,encoded):
        d=json.loads(encoded);s=d['state']
        if hashlib.sha256(canonical(s).encode()).hexdigest()!=d['checksum']:raise ValueError('search_checksum')
        required={'version','h','uncertainty','scan','completed','queue','brackets','finished','actions','evaluations'}
        if set(s)!=required or s['version']!='two-mode-search-v1' or s['h']!=.02 or s['uncertainty']!=.0001:raise ValueError('search_schema_or_policy')
        for key,p in s['completed'].items():
            if key!=float(p['U']).hex() or np.array(p['modes']).shape!=(2,2,2) or len(p['growth'])!=2 or len(p['frequency'])!=2 or max(p['original_eigen_residuals'])>=1e-12:raise ValueError('search_point_identity')
        if s['finished'] and (s['brackets'] is None or not all(b['done'] for b in s['brackets'])):raise ValueError('premature_search_completion')
        obj=cls();obj.state=copy.deepcopy(s);return obj
