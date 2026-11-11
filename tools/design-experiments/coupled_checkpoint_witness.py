import copy,json,hashlib,numpy as np
from dataclasses import replace
import uvlm_wake_assembly_v2 as aero
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)
def seal(payload):return canonical(dict(payload=payload,checksum=hashlib.sha256(canonical(payload).encode()).hexdigest()))
def model_id(p):return aero.digest('J1-BDF2-physical-velocity-v1',p.model.identity,p.data['C'],p.data['e0'],p.data['weights'],p.data['G'],p.B,p.mass,p.fixed,p.flow.identity,p.h)
def physical_next(p,parent,s):
    state=p.decode(s);old=parent.wakes[0];release=('pg','released',parent.epoch+1,parent.time,p.h);w=aero.Wake('pg',state.wakes[0].vertices,state.wakes[0].mu,(old.vertex_ids[0],tuple((release,j) for j in range(2)),old.vertex_ids[1]),(old.row_ids[0],release,old.row_ids[1]),(parent.time+p.h,parent.time,old.birth[1]));out=replace(state,epoch=parent.epoch+1,time=parent.time+p.h,wakes=(w,),initialization='continued',gamma_previous=parent.gamma);p.model.validate_state(out);return out
def packet(p,s,parent_s,physical):
    nb=p.nb;gamma=[s[8:8+nb],parent_s[8:8+nb],parent_s[8+nb:8+2*nb]];vel=[np.einsum('nci,i->nc',p.B,u) for u in [s[2:4],parent_s[2:4],parent_s[6:8]]];gd=(3*gamma[0]-4*gamma[1]+gamma[2])/(2*p.h);av=(3*vel[0]-4*vel[1]+vel[2])/(2*p.h);g=p.model.geometry(physical.q,physical.v);loads=aero.loads(p.model,g,physical.wakes,physical.gamma,gd,p.flow,physical.time);return dict(version='joint-checkpoint-witness-v1',model_id=model_id(p),input_id=p.flow.identity,scheme='BDF2-physical-velocity-v1',h=p.h,epoch=physical.epoch,time=physical.time,state=s.tolist(),physical=aero.checkpoint(p.model,physical),stencil=dict(gamma=[x.tolist() for x in gamma],velocity=[x.tolist() for x in vel]),gamma_dot=gd.tolist(),acceleration_algorithmic=av.tolist(),loads={key:loads[key].tolist() for key in ['Q','F','M']})
def restore(p,encoded):
    d=json.loads(encoded)
    if set(d)!=set(['payload','checksum']) or hashlib.sha256(canonical(d['payload']).encode()).hexdigest()!=d['checksum']:raise ValueError('joint_checksum')
    a=d['payload'];required={'version','model_id','input_id','scheme','h','epoch','time','state','physical','stencil','gamma_dot','acceleration_algorithmic','loads'}
    if set(a)!=required or a['version']!='joint-checkpoint-witness-v1' or a['model_id']!=model_id(p) or a['input_id']!=p.flow.identity or a['scheme']!='BDF2-physical-velocity-v1' or a['h']!=p.h:raise ValueError('joint_schema_or_policy')
    s=np.array(a['state']);physical=aero.restore(p.model,a['physical'],p.flow)
    if s.shape!=(p.dimension,) or len(a['stencil']['gamma'])!=3 or len(a['stencil']['velocity'])!=3:raise ValueError('missing_joint_history')
    G=np.array(a['stencil']['gamma']);V=np.array(a['stencil']['velocity']);gd=(3*G[0]-4*G[1]+G[2])/(2*p.h);av=(3*V[0]-4*V[1]+V[2])/(2*p.h)
    if not np.array_equal(gd,a['gamma_dot']) or not np.array_equal(av,a['acceleration_algorithmic']):raise ValueError('joint_derivative_stencil')
    if not np.array_equal(G[0],s[8:8+p.nb]) or not np.array_equal(G[1],s[8+p.nb:8+2*p.nb]) or not np.array_equal(physical.gamma,G[0]) or not np.array_equal(physical.gamma_previous,G[1]) or physical.time!=a['time'] or physical.epoch!=a['epoch']:raise ValueError('joint_physical_history_mismatch')
    g=p.model.geometry(physical.q,physical.v);loads=aero.loads(p.model,g,physical.wakes,physical.gamma,gd,p.flow,physical.time)
    for key in ['Q','F','M']:
        if not np.array_equal(loads[key],np.array(a['loads'][key])):raise ValueError('joint_load_packet_mismatch')
    return s,physical,a
class Store:
    def __init__(self,encoded):self.encoded=encoded
    def qualify(self,p,payload):
        encoded=seal(payload);restore(p,encoded);return dict(parent_hash=hashlib.sha256(self.encoded.encode()).hexdigest(),encoded=encoded,qualification=hashlib.sha256(encoded.encode()).hexdigest())
    def commit(self,p,trial):
        if trial['parent_hash']!=hashlib.sha256(self.encoded.encode()).hexdigest() or trial['qualification']!=hashlib.sha256(trial['encoded'].encode()).hexdigest():raise ValueError('stale_or_modified_joint_trial')
        _,old,_=restore(p,self.encoded);_,new,_=restore(p,trial['encoded'])
        if new.epoch!=old.epoch+1 or abs(new.time-old.time-p.h)>1e-12:raise ValueError('joint_generation')
        self.encoded=trial['encoded']
