"""Finite multi-surface P0-bound/W-P1-wake design implementation.

Fixed material position maps are independent of the trial q. Euler convection
and BDF1 are the declared time-step witness, not the product integrator.
"""
from dataclasses import dataclass
import hashlib,json,copy
import numpy as np
from types import MappingProxyType
import uvlm_analytic_induction_v2 as induction
import uvlm_dual13 as ad
import uvlm_geometry_guard as guard
_INHERIT=object()

def readonly(a,dtype=float):
    a=np.array(a,dtype=dtype,copy=True);a.flags.writeable=False;return a
def digest(*values):
    h=hashlib.sha256()
    for value in values:
        if isinstance(value,np.ndarray):h.update(str((value.shape,value.dtype)).encode());h.update(value.tobytes())
        else:h.update(repr(value).encode())
    return h.hexdigest()

@dataclass(frozen=True)
class Wake:
    sid:str
    vertices:np.ndarray
    mu:np.ndarray
    vertex_ids:tuple
    row_ids:tuple
    birth:tuple
    def __post_init__(self):
        object.__setattr__(self,'vertex_ids',tuple(tuple(row) for row in self.vertex_ids));object.__setattr__(self,'row_ids',tuple(self.row_ids));object.__setattr__(self,'birth',tuple(self.birth))
        object.__setattr__(self,'vertices',readonly(self.vertices));object.__setattr__(self,'mu',readonly(self.mu))
        if self.vertices.ndim!=3 or self.vertices.shape[2]!=3 or self.mu.shape!=(len(self.vertices),self.vertices.shape[1]-1) or len(self.vertex_ids)!=len(self.vertices) or any(len(row)!=self.vertices.shape[1] for row in self.vertex_ids) or len(self.row_ids)!=len(self.vertices) or len(self.birth)!=len(self.vertices):raise ValueError('invalid_wake_shape')
        if not np.isfinite(self.vertices).all() or not np.isfinite(self.mu).all():raise ValueError('nonfinite_wake')

@dataclass(frozen=True)
class State:
    identity:str
    epoch:int
    time:float
    q:np.ndarray
    v:np.ndarray
    gamma:np.ndarray
    wakes:tuple
    gamma_previous:np.ndarray=None
    last_h:float=None
    initialization:str='explicit_state'
    input_identity:str=''
    convection_mode:str='unselected'
    max_intervals:object=None
    def __post_init__(self):
        object.__setattr__(self,'wakes',tuple(self.wakes))
        for name in ['q','v','gamma']:object.__setattr__(self,name,readonly(getattr(self,name)))
        if self.gamma_previous is not None:
            object.__setattr__(self,'gamma_previous',readonly(self.gamma_previous))
            if self.gamma_previous.shape!=self.gamma.shape or not np.isfinite(self.gamma_previous).all() or self.last_h is None or not np.isfinite(self.last_h) or self.last_h<=0:raise ValueError('invalid_circulation_history')
        if not isinstance(self.epoch,(int,np.integer)) or self.epoch<0 or not np.isfinite(self.time) or not all(np.isfinite(getattr(self,n)).all() for n in ['q','v','gamma']):raise ValueError('invalid_state')
        if self.gamma_previous is None and self.last_h is not None:raise ValueError('missing_circulation_history')
        if self.convection_mode not in ['unselected','free','prescribed'] or self.epoch>0 and self.convection_mode=='unselected' or self.max_intervals is not None and (not isinstance(self.max_intervals,(int,np.integer)) or isinstance(self.max_intervals,bool) or self.max_intervals<1):raise ValueError('invalid_saved_wake_policy')
    def fingerprint(self):return digest(self.identity,self.input_identity,self.convection_mode,self.max_intervals,self.epoch,self.time,self.q,self.v,self.gamma,self.gamma_previous,self.last_h,self.initialization,*[(w.sid,w.vertex_ids,w.row_ids,w.birth,digest(w.vertices,w.mu)) for w in self.wakes])
    def gamma_dot(self):
        if self.gamma_previous is not None:return (self.gamma-self.gamma_previous)/self.last_h
        if self.initialization=='cold_zero_history':return np.zeros_like(self.gamma)
        raise ValueError('missing_circulation_derivative_history')

class Flow:
    def __setattr__(self,name,value):
        if getattr(self,'_sealed',False):raise AttributeError('immutable_flow_definition')
        object.__setattr__(self,name,value)
    def __init__(self,U,amplitude=None,wave=None,omega=0.):
        self.U=readonly(U);self.amplitude=readonly(np.zeros(3) if amplitude is None else amplitude);self.wave=readonly(np.zeros(3) if wave is None else wave);self.omega=float(omega)
        if self.U.shape!=(3,) or self.amplitude.shape!=(3,) or self.wave.shape!=(3,) or not all(np.isfinite(a).all() for a in [self.U,self.amplitude,self.wave,self.omega]):raise ValueError('invalid_flow')
        if abs(self.amplitude@self.wave)>1e-12*max(1.,np.linalg.norm(self.amplitude)*np.linalg.norm(self.wave)):raise ValueError('non_solenoidal_gust')
        self.identity=digest(self.U,self.amplitude,self.wave,self.omega);self._sealed=True
    def value(self,p,t):return self.U+self.amplitude*np.sin(np.asarray(p)@self.wave-self.omega*t)[...,None]
    def jvp(self,p,t,dp,dU=None):return self.amplitude*np.cos(np.asarray(p)@self.wave-self.omega*t)[...,None]*(np.asarray(dp)@self.wave)[...,None]+(0 if dU is None else dU)

class Model:
    def __setattr__(self,name,value):
        if getattr(self,'_sealed',False):raise AttributeError('immutable_prepared_aerodynamic_model')
        object.__setattr__(self,name,value)
    def __init__(self,bindings,eps,rho=1.2):
        if not np.isfinite(eps) or eps<0 or not np.isfinite(rho) or rho<=0:raise ValueError('invalid_core_or_density')
        self.eps=float(eps);self.rho=float(rho);self.bindings=tuple((b['sid'],b['m'],b['n']) for b in bindings)
        if len(set(s for s,m,n in self.bindings))!=len(bindings):raise ValueError('duplicate_surface_id')
        self.ns=sum(b['Nv'].shape[1] for b in bindings);nv=sum(len(b['Nv']) for b in bindings);nb=sum(len(b['Nc']) for b in bindings);self.Nv=np.zeros((nv,self.ns));self.Nc=np.zeros((nb,self.ns));self.Dx=self.Nc.copy();self.Dy=self.Nc.copy();faces=[];vids=[];self.surface_info=[];iv=ib=iq=it=0
        for b in bindings:
            sid,m,n=b['sid'],b['m'],b['n'];ns=b['Nv'].shape[1];nvert=(m+1)*(n+1);ncell=m*n
            if m<1 or n<1 or b['Nv'].shape!=(nvert,ns) or any(b[key].shape!=(ncell,ns) for key in ['Nc','Dx','Dy']):raise ValueError('invalid_surface_binding')
            if any(not np.isfinite(b[key]).all() for key in ['Nv','Nc','Dx','Dy']):raise ValueError('nonfinite_surface_binding')
            self.Nv[iv:iv+nvert,iq:iq+ns]=b['Nv'];self.Nc[ib:ib+ncell,iq:iq+ns]=b['Nc'];self.Dx[ib:ib+ncell,iq:iq+ns]=b['Dx'];self.Dy[ib:ib+ncell,iq:iq+ns]=b['Dy'];ids=[(sid,'bound',i,j) for i in range(m+1) for j in range(n+1)];vids+=ids
            for i in range(m):
                for j in range(n):a=iv+i*(n+1)+j;faces.append([a,a+1,a+n+2,a+n+1])
            te=np.arange(iv+m*(n+1),iv+(m+1)*(n+1));self.surface_info.append(dict(sid=sid,m=m,n=n,te=te,mu_slice=slice(it,it+n),gamma_te=np.arange(ib+(m-1)*n,ib+m*n)));iv+=nvert;ib+=ncell;iq+=ns;it+=n
        self.vertex_ids=tuple(vids);self.faces=np.array(faces,int);self.nb=nb;self.nte=it;self.T=np.zeros((it,nb))
        for info in self.surface_info:self.T[info['mu_slice'],info['gamma_te']]=np.eye(info['n'])
        self.edges=sorted({tuple(sorted((f[i],f[(i+1)%4]))) for f in faces for i in range(4)});index={e:i for i,e in enumerate(self.edges)};self.D2=np.zeros((len(self.edges),nb),int)
        for j,f in enumerate(faces):
            for i in range(4):pair=(f[i],f[(i+1)%4]);e=tuple(sorted(pair));self.D2[index[e],j]+=1 if pair==e else -1
        self.D1=np.zeros((nv,len(self.edges)),int)
        for i,(a,b) in enumerate(self.edges):self.D1[a,i]=-1;self.D1[b,i]=1
        if np.any(self.D1@self.D2):raise ValueError('invalid_integer_boundary')
        te_edges={tuple(sorted((int(i),int(j)))) for info in self.surface_info for i,j in zip(info['te'][:-1],info['te'][1:])};self.force_edges=tuple(i for i,e in enumerate(self.edges) if e not in te_edges);self.formulation='dual-quarter-P0/WP1-RM-fixed-global-core/Euler-BDF1-witness-v1';self.identity=digest(self.formulation,self.bindings,self.Nv,self.Nc,self.Dx,self.Dy,self.faces,self.D2,self.eps,self.rho)
        for name in ['Nv','Nc','Dx','Dy','T','D1','D2','faces']:getattr(self,name).flags.writeable=False
        for info in self.surface_info:
            for name in ['te','gamma_te']:info[name].flags.writeable=False
        self.surface_info=tuple(MappingProxyType(info) for info in self.surface_info);self.edges=tuple(self.edges);self._sealed=True
    def geometry(self,q,v=None,dq=None,dv=None):
        q=np.asarray(q);v=np.zeros_like(q) if v is None else np.asarray(v);dq=np.zeros_like(q) if dq is None else np.asarray(dq);dv=np.zeros_like(q) if dv is None else np.asarray(dv)
        if any(a.shape!=(self.ns,3) or not np.isfinite(a).all() for a in [q,v,dq,dv]):raise ValueError('invalid_bound_state_shape')
        a=self.Dx@q;b=self.Dy@q;cross=np.cross(a,b);length=np.linalg.norm(cross,axis=1)
        if np.any(length<1e-12*np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)) or np.any(length==0):raise ValueError('degenerate_collocation_frame')
        normals=cross/length[:,None];dc=np.cross(self.Dx@dq,b)+np.cross(a,self.Dy@dq);dn=(dc-normals*np.sum(normals*dc,axis=1)[:,None])/length[:,None]
        return dict(X=self.Nv@q,V=self.Nv@v,P=self.Nc@q,CV=self.Nc@v,n=normals,dX=self.Nv@dq,dV=self.Nv@dv,dP=self.Nc@dq,dCV=self.Nc@dv,dn=dn)
    def wake_faces(self,wakes):
        for iw,w in enumerate(wakes):
            n=w.mu.shape[1]
            for r in range(len(w.vertices)-1):
                for j in range(n):
                    xyz=np.array([w.vertices[r,j],w.vertices[r,j+1],w.vertices[r+1,j+1],w.vertices[r+1,j]]);ids=[w.vertex_ids[r][j],w.vertex_ids[r][j+1],w.vertex_ids[r+1][j+1],w.vertex_ids[r+1][j]];idx=[r*n+j,r*n+j,(r+1)*n+j,(r+1)*n+j]
                    for tri in [[0,2,1],[0,3,2]]:yield iw,np.array(idx)[tri],tuple(ids[k] for k in tri),xyz[tri],(w.sid,w.row_ids[r],w.row_ids[r+1],j,tuple(tri)),r,j,tri
    def validate_geometry(self,g,wakes):
        facets=[]
        for k,f in enumerate(self.faces):
            for j,tri in enumerate([[0,2,1],[0,3,2]]):
                ids=f[tri];x=g['X'][ids];normal=np.cross(x[1]-x[0],x[2]-x[0])
                if normal@g['n'][k]<=0:raise ValueError('inverted_bound_facet')
                facets.append((('bound',k,j),tuple(self.vertex_ids[i] for i in ids),x))
        for iw,idx,ids,x,name,r,j,tri in self.wake_faces(wakes):facets.append((('wake',name),ids,x))
        return guard.validate_facets(facets)
    def validate_state(self,state):
        if state.identity!=self.identity or state.gamma.shape!=(self.nb,) or state.q.shape!=(self.ns,3) or len(state.wakes)!=len(self.surface_info):raise ValueError('state_model_mismatch')
        g=self.geometry(state.q,state.v)
        if state.initialization=='cold_zero_history' and (state.epoch!=0 or np.any(state.gamma!=0) or any(len(w.vertices)!=1 or np.any(w.mu!=0) for w in state.wakes)):raise ValueError('false_cold_history')
        for w,info in zip(state.wakes,self.surface_info):
            if w.sid!=info['sid'] or w.mu.shape[1]!=info['n'] or w.vertex_ids[0]!=tuple(self.vertex_ids[i] for i in info['te']) or np.linalg.norm(w.vertices[0]-g['X'][info['te']])>1e-10 or np.linalg.norm(w.mu[0]-state.gamma[info['gamma_te']])>1e-10:raise ValueError('invalid_attached_wake_state')
            if len(set(w.row_ids))!=len(w.row_ids) or len(set(x for row in w.vertex_ids for x in row))!=w.vertices.shape[0]*w.vertices.shape[1]:raise ValueError('duplicate_material_id')
            if not np.isfinite(w.birth).all() or abs(w.birth[0]-state.time)>1e-12*max(1.,abs(state.time)) or np.any(np.diff(w.birth)>1e-12):raise ValueError('invalid_wake_age_history')
        return self.validate_geometry(g,state.wakes)

def segment_basis(p,a,b,eps,incident=False,dp=None,da=None,db=None,fraction=None):
    length=np.linalg.norm(b-a)
    if length<=0:raise ValueError('zero_edge')
    contact=guard.segment_distance(p,a,b)<=1e-11*length
    if contact and not incident:raise ValueError('nonincident_edge_contact')
    if incident and not contact:raise ValueError('edge_incidence_geometry_mismatch')
    if contact and incident:
        if fraction is None:raise ValueError('missing_incident_material_fraction')
        scale=max(np.linalg.norm(a),np.linalg.norm(b),length)
        if np.linalg.norm(p-((1-fraction)*a+fraction*b))>64*np.finfo(float).eps*scale:raise ValueError('incident_point_position_precision')
        if dp is not None and np.linalg.norm(dp-((1-fraction)*da+fraction*db))>1e-10*max(1.,np.linalg.norm(dp),np.linalg.norm(da),np.linalg.norm(db)):raise ValueError('incident_tangent_mismatch')
        # Exact straight self-line value and derivative ALONG the declared
        # material point map. Other edges/sheets remain in the field sum.
        return np.zeros(3),np.zeros(3)
    pp=np.array([ad.Dual.variable((p[j]-a[j])/length,j,1/length) for j in range(3)],object);aa=np.array([ad.Dual.variable(0.,3+j,1/length) for j in range(3)],object);bb=np.array([ad.Dual.variable((b[j]-a[j])/length,6+j,1/length) for j in range(3)],object);ee=ad.Dual.variable(eps/length,9,1/length);v=induction.segment(pp,aa,bb,ee,1.,1.);value=np.array([x.v for x in v])/(4*np.pi*length)
    if dp is None:return value,np.zeros(3)
    direction=np.r_[dp,da,db,0.,np.zeros(3)];j=np.array([x.g for x in v])/(4*np.pi*length);return value,j@direction

def field_basis(model,g,wakes,p,vertex_id=None,edge_id=None,dp=None,dwakes=None):
    B=np.zeros((3,model.nb));dB=B.copy();sizes=[w.mu.size for w in wakes];offset=np.cumsum([0]+sizes);W=np.zeros((3,sum(sizes)));dW=W.copy();dp=np.zeros(3) if dp is None else dp;dwakes=[np.zeros_like(w.vertices) for w in wakes] if dwakes is None else dwakes
    for i,(a,b) in enumerate(model.edges):
        ids=(model.vertex_ids[a],model.vertex_ids[b]);inc=vertex_id in ids or edge_id==i;fraction=.5 if edge_id==i else 0. if vertex_id==ids[0] else 1. if vertex_id==ids[1] else None;v,dv=segment_basis(p,g['X'][a],g['X'][b],model.eps,inc,dp,g['dX'][a],g['dX'][b],fraction);B+=v[:,None]*model.D2[i];dB+=dv[:,None]*model.D2[i]
    for iw,idx,ids,x,name,r,j,tri in model.wake_faces(wakes):
        v,jac,_=induction.triangle(p,x,model.eps,vertex_id in ids);der=np.array([dwakes[iw][r,j],dwakes[iw][r,j+1],dwakes[iw][r+1,j+1],dwakes[iw][r+1,j]])[tri];direction=np.r_[dp,der.ravel(),0.];dv=jac@direction
        for k,index in enumerate(idx):W[:,offset[iw]+index]+=v[k];dW[:,offset[iw]+index]+=dv[k]
    return B,W,dB,dW

def field(model,g,wakes,gamma,p,vertex_id=None,edge_id=None,dp=None,dwakes=None,dgamma=None,dmu=None):
    B,W,dB,dW=field_basis(model,g,wakes,p,vertex_id,edge_id,dp,dwakes);mu=np.concatenate([w.mu.ravel() for w in wakes]);dgamma=np.zeros_like(gamma) if dgamma is None else dgamma;dmu=np.zeros_like(mu) if dmu is None else dmu;return B@gamma+W@mu,dB@gamma+dW@mu+B@dgamma+W@dmu

class Workspace:
    def __init__(self):self.key=None;self.A=None;self.W=None;self.hits=0;self.misses=0
    def assemble(self,model,g,wakes):
        key=digest(model.identity,g['X'],g['P'],g['n'],*[(w.vertex_ids,digest(w.vertices)) for w in wakes])
        if key==self.key:self.hits+=1;return self.A,self.W
        Wrows=[];Brows=[]
        for p,n in zip(g['P'],g['n']):B,W,_,_=field_basis(model,g,wakes,p);Brows.append(n@B);Wrows.append(n@W)
        B=np.array(Brows);W=np.array(Wrows);A=B.copy();offset=0
        for w,info in zip(wakes,model.surface_info):A+=W[:,offset:offset+info['n']]@model.T[info['mu_slice']];offset+=w.mu.size
        self.key=key;self.A=readonly(A);self.W=readonly(W);self.misses+=1;return self.A,self.W

def cold_state(model,q,v=None,time=0.):
    g=model.geometry(q,v);gamma=np.zeros(model.nb);wakes=[]
    for info in model.surface_info:
        ids=tuple(model.vertex_ids[i] for i in info['te']);wakes.append(Wake(info['sid'],g['X'][info['te']][None],np.zeros((1,info['n'])),(ids,),((info['sid'],'attached'),),(float(time),)))
    return State(model.identity,0,float(time),q,np.zeros_like(q) if v is None else v,gamma,tuple(wakes),initialization='cold_zero_history')

def convection(model,state,flow,mode=None):
    if mode is None:mode='free' if state.convection_mode=='unselected' else state.convection_mode
    if mode not in ['free','prescribed']:raise ValueError('unknown_convection_mode')
    if mode=='free' and model.eps<=0:raise ValueError('free_self_requires_positive_core')
    model.validate_state(state);g=model.geometry(state.q,state.v);out=[]
    for w in state.wakes:
        vel=flow.value(w.vertices,state.time).copy()
        if mode=='free':
            for r in range(len(w.vertices)):
                for j in range(w.vertices.shape[1]):vel[r,j]+=field(model,g,state.wakes,state.gamma,w.vertices[r,j],w.vertex_ids[r][j])[0]
        out.append(vel)
    return tuple(out)

def begin_trial(model,parent,q,v,h,flow,mode=None,max_intervals=_INHERIT):
    if mode is None:mode='free' if parent.convection_mode=='unselected' else parent.convection_mode
    if max_intervals is _INHERIT:max_intervals=parent.max_intervals
    if not np.isfinite(h) or h<=0 or max_intervals is not None and (not isinstance(max_intervals,(int,np.integer)) or isinstance(max_intervals,bool) or max_intervals<1):raise ValueError('invalid_step_or_tail')
    model.validate_state(parent);g=model.geometry(q,v);vel=convection(model,parent,flow,mode);wakes=[];retired=[]
    for w,info,V in zip(parent.wakes,model.surface_info,vel):
        release=(w.sid,'released',parent.epoch+1,float(parent.time),float(h));ids=tuple((release,j) for j in range(info['n']+1));xyz=np.concatenate([g['X'][info['te']][None],w.vertices+h*V]);mu=np.vstack([np.zeros(info['n']),w.mu]);vertex_ids=(w.vertex_ids[0],ids)+w.vertex_ids[1:];row_ids=(w.row_ids[0],release)+w.row_ids[1:];birth=(parent.time+h,parent.time)+w.birth[1:]
        if max_intervals is not None and len(xyz)>max_intervals+1:
            retired.extend(row_ids[max_intervals+1:]);xyz=xyz[:max_intervals+1];mu=mu[:max_intervals+1];vertex_ids=vertex_ids[:max_intervals+1];row_ids=row_ids[:max_intervals+1];birth=birth[:max_intervals+1]
        wakes.append(Wake(w.sid,xyz,mu,vertex_ids,row_ids,birth))
    model.validate_geometry(g,wakes);return dict(identity=model.identity,flow_identity=flow.identity,parent_epoch=parent.epoch,parent_fingerprint=parent.fingerprint(),time=parent.time+h,h=float(h),q=np.array(q,copy=True),v=np.array(v,copy=True),wakes=tuple(wakes),mode=mode,retired=tuple(retired),tail='growing' if max_intervals is None else 'finite_closed_window',max_intervals=max_intervals,qualified=False)

def solve_trial(model,parent,trial,flow,workspace=None):
    trial['qualified']=False
    if trial['identity']!=model.identity or trial['flow_identity']!=flow.identity or trial['parent_epoch']!=parent.epoch or trial['parent_fingerprint']!=parent.fingerprint():raise ValueError('stale_aerodynamic_trial')
    workspace=Workspace() if workspace is None else workspace;g=model.geometry(trial['q'],trial['v']);model.validate_geometry(g,trial['wakes']);A,W=workspace.assemble(model,g,trial['wakes']);mu=np.concatenate([w.mu.ravel() for w in trial['wakes']]);offset=0
    for w,info in zip(trial['wakes'],model.surface_info):mu[offset:offset+info['n']]=0.;offset+=w.mu.size
    rhs=np.sum(g['n']*(g['CV']-flow.value(g['P'],trial['time'])),axis=1)-W@mu;gamma=np.linalg.solve(A,rhs);res=A@gamma-rhs;error=np.linalg.norm(res)/max(np.linalg.norm(rhs),1e-30)
    if not np.isfinite(gamma).all() or error>1e-11:raise ValueError('nonpenetration_solve_unqualified')
    wakes=[]
    for w,info in zip(trial['wakes'],model.surface_info):mu=w.mu.copy();mu[0]=gamma[info['gamma_te']];wakes.append(Wake(w.sid,w.vertices,mu,w.vertex_ids,w.row_ids,w.birth))
    trial.update(gamma=gamma,wakes=tuple(wakes),gamma_dot=(gamma-parent.gamma)/trial['h'],residual=float(error));trial['state']=State(model.identity,parent.epoch+1,trial['time'],trial['q'],trial['v'],gamma,tuple(wakes),gamma_previous=parent.gamma,last_h=trial['h'],initialization='continued',input_identity=flow.identity,convection_mode=trial['mode'],max_intervals=trial['max_intervals']);model.validate_state(trial['state']);trial['loads']=loads(model,g,tuple(wakes),gamma,trial['gamma_dot'],flow,trial['time']);trial['qualified_fingerprint']=trial_fingerprint(trial);trial['qualified']=True;return trial

def trial_fingerprint(t):return digest(t['identity'],t['flow_identity'],t['parent_epoch'],t['parent_fingerprint'],t['time'],t['h'],t['mode'],t['tail'],t['max_intervals'],t['retired'],t['q'],t['v'],t.get('gamma'),t.get('gamma_dot'),None if 'state' not in t else t['state'].fingerprint(),None if 'loads' not in t else digest(*[t['loads'][name] for name in ['Q','Fv','positions','forces','F','M']]),*[(w.vertex_ids,w.row_ids,w.birth,digest(w.vertices,w.mu)) for w in t['wakes']])

def accept(model,parent,trial):
    if not trial['qualified'] or trial['identity']!=model.identity or trial['parent_epoch']!=parent.epoch or trial['parent_fingerprint']!=parent.fingerprint() or trial_fingerprint(trial)!=trial['qualified_fingerprint']:raise ValueError('unqualified_or_changed_trial')
    model.validate_state(trial['state']);return trial['state']

def loads(model,g,wakes,gamma,gamma_dot,flow,time,dgamma=None,dgamma_dot=None,dwakes=None,dmu=None,dU=None):
    dgamma=np.zeros_like(gamma) if dgamma is None else dgamma;dgamma_dot=np.zeros_like(gamma_dot) if dgamma_dot is None else dgamma_dot;Fv=np.zeros_like(g['X']);dFv=Fv.copy();positions=[];forces=[];dforces=[];kappa=model.D2@gamma;dkappa=model.D2@dgamma
    for edge in model.force_edges:
        a,b=model.edges[edge];p=(g['X'][a]+g['X'][b])/2;dp=(g['dX'][a]+g['dX'][b])/2;length=g['X'][b]-g['X'][a];dl=g['dX'][b]-g['dX'][a];ind,dind=field(model,g,wakes,gamma,p,edge_id=edge,dp=dp,dwakes=dwakes,dgamma=dgamma,dmu=dmu);rel=flow.value(p,time)+ind-(g['V'][a]+g['V'][b])/2;drel=flow.jvp(p,time,dp,dU)+dind-(g['dV'][a]+g['dV'][b])/2;force=model.rho*kappa[edge]*np.cross(rel,length);der=model.rho*(dkappa[edge]*np.cross(rel,length)+kappa[edge]*(np.cross(drel,length)+np.cross(rel,dl)))
        Fv[a]+=force/2;Fv[b]+=force/2;dFv[a]+=der/2;dFv[b]+=der/2;positions.append(p);forces.append(force);dforces.append(der)
    circulatory_count=len(forces)
    for panel,f in enumerate(model.faces):
        for tri in [[0,2,1],[0,3,2]]:
            ids=f[tri];a,b,c=g['X'][ids];da,db,dc=g['dX'][ids];area=np.cross(b-a,c-a)/2;darea=(np.cross(db-da,c-a)+np.cross(b-a,dc-da))/2;force=model.rho*gamma_dot[panel]*area;der=model.rho*(dgamma_dot[panel]*area+gamma_dot[panel]*darea);Fv[ids]+=force/3;dFv[ids]+=der/3;positions.append((a+b+c)/3);forces.append(force);dforces.append(der)
    positions=np.array(positions);forces=np.array(forces);out=dict(Q=model.Nv.T@Fv,dQ=model.Nv.T@dFv,Fv=Fv,dFv=dFv,positions=positions,forces=forces,dforces=np.array(dforces),circulatory_count=circulatory_count,F=forces.sum(axis=0),M=np.cross(positions,forces).sum(axis=0))
    if any(not np.isfinite(v).all() for v in out.values() if isinstance(v,np.ndarray)):raise ValueError('nonfinite_aerodynamic_load')
    return out

def trial_jvp(model,parent,trial,flow,dq,dv=None,dU=None,workspace=None):
    if not trial['qualified'] or trial_fingerprint(trial)!=trial['qualified_fingerprint'] or trial['flow_identity']!=flow.identity or trial['parent_fingerprint']!=parent.fingerprint():raise ValueError('changed_trial_linearization')
    g=model.geometry(trial['q'],trial['v'],dq,dv);workspace=Workspace() if workspace is None else workspace;A,_=workspace.assemble(model,g,trial['wakes']);dw=[np.zeros_like(w.vertices) for w in trial['wakes']]
    for d,info in zip(dw,model.surface_info):d[0]=g['dX'][info['te']]
    if dU is not None:
        for d in dw:d[1:]=trial['h']*np.asarray(dU)
    residual=[]
    for i,p in enumerate(g['P']):
        ind,der=field(model,g,trial['wakes'],trial['gamma'],p,dp=g['dP'][i],dwakes=dw);relative=ind+flow.value(p,trial['time'])-g['CV'][i];residual.append(g['dn'][i]@relative+g['n'][i]@(der+flow.jvp(p,trial['time'],g['dP'][i],dU)-g['dCV'][i]))
    dgamma=np.linalg.solve(A,-np.array(residual));dmu=[]
    for w,info in zip(trial['wakes'],model.surface_info):a=np.zeros_like(w.mu);a[0]=dgamma[info['gamma_te']];dmu.extend(a.ravel())
    out=loads(model,g,trial['wakes'],trial['gamma'],trial['gamma_dot'],flow,trial['time'],dgamma,dgamma/trial['h'],dw,np.array(dmu),dU);out['dgamma']=dgamma;return out

def convection_jvp(model,state,flow,dq,dgamma,dwakes,dmu,mode=None):
    if mode is None:mode='free' if state.convection_mode=='unselected' else state.convection_mode
    if mode not in ['free','prescribed']:raise ValueError('unknown_convection_mode')
    model.validate_state(state);g=model.geometry(state.q,state.v,dq);out=[]
    for iw,w in enumerate(state.wakes):
        dv=flow.jvp(w.vertices,state.time,dwakes[iw])
        if mode=='free':
            for r in range(len(w.vertices)):
                for j in range(w.vertices.shape[1]):dv[r,j]+=field(model,g,state.wakes,state.gamma,w.vertices[r,j],w.vertex_ids[r][j],dp=dwakes[iw][r,j],dwakes=dwakes,dgamma=dgamma,dmu=dmu)[1]
        out.append(dv)
    return tuple(out)

def checkpoint(model,state):
    model.validate_state(state);data=dict(format='wp1-physical-state-design-v1',identity=state.identity,input_identity=state.input_identity,convection_mode=state.convection_mode,max_intervals=state.max_intervals,epoch=state.epoch,time=state.time,q=state.q.tolist(),v=state.v.tolist(),gamma=state.gamma.tolist(),gamma_previous=None if state.gamma_previous is None else state.gamma_previous.tolist(),last_h=state.last_h,initialization=state.initialization,wakes=[dict(sid=w.sid,vertices=w.vertices.tolist(),mu=w.mu.tolist(),vertex_ids=w.vertex_ids,row_ids=w.row_ids,birth=w.birth) for w in state.wakes]);payload=json.dumps(data,sort_keys=True,allow_nan=False,separators=(',',':'));return json.dumps(dict(payload=payload,checksum=hashlib.sha256(payload.encode()).hexdigest()),sort_keys=True)

def restore(model,encoded,flow=None):
    wrap=json.loads(encoded)
    if hashlib.sha256(wrap['payload'].encode()).hexdigest()!=wrap['checksum']:raise ValueError('checkpoint_integrity')
    data=json.loads(wrap['payload'])
    if data.get('format')!='wp1-physical-state-design-v1' or data['identity']!=model.identity:raise ValueError('incompatible_checkpoint')
    if data.get('input_identity') and (flow is None or flow.identity!=data['input_identity']):raise ValueError('checkpoint_input_mismatch')
    def tuples(x):return tuple(tuples(y) for y in x) if isinstance(x,list) else x
    wakes=tuple(Wake(w['sid'],w['vertices'],w['mu'],tuples(w['vertex_ids']),tuples(w['row_ids']),tuples(w['birth'])) for w in data['wakes']);state=State(data['identity'],data['epoch'],data['time'],data['q'],data['v'],data['gamma'],wakes,data['gamma_previous'],data['last_h'],data['initialization'],data.get('input_identity',''),data['convection_mode'],data['max_intervals']);model.validate_state(state);return state

class Store:
    def __init__(self,model,state):model.validate_state(state);self.model=model;self._state=state
    @property
    def state(self):return self._state
    def begin(self,q,v,h,flow,**kwargs):return begin_trial(self.model,self._state,q,v,h,flow,**kwargs)
    def solve(self,trial,flow,workspace=None):return solve_trial(self.model,self._state,trial,flow,workspace)
    def accept(self,trial):
        qualified=accept(self.model,self._state,trial);self._state=qualified;return qualified
