"""All shared q01 constraints with compatible global geometric kernels retained."""
import numpy as np
from scipy.linalg import qr
import mixed_mean_stress_bending_v2 as t
import mixed_mean_stress_geometric_kernel as known
def assemble(cells,n,kappa):
    ny=sum(c['ny'] for c in cells);registry={};offset=0;trace_count=0;fixed=[];equations=[]
    for index,cell in enumerate(cells):
        cell['yoffset']=offset;offset+=cell['ny'];nodes=[];transforms=[]
        for face,key in zip(cell['faces'],t.fit.edge_keys(index,n)):
            Mc=face['Mcoeff'];nd=3*face['n'];sl=slice(3*face['offset'],3*(face['offset']+face['n']))
            if key not in registry:
                ids=np.arange(trace_count,trace_count+face['n']);trace_count+=face['n'];registry[key]=dict(M=Mc.copy(),nodes=ids)
                if key[0]==0 and key[1]==0:fixed.extend(ids)
            r=registry[key];T=Mc.T@r['M'];Tc=np.kron(T,np.eye(3));nodes.extend(r['nodes']);transforms.append(Tc)
            if np.linalg.norm(T.T@T-np.eye(face['n']))>1e-10:raise ValueError('MMS1_bending_face_transform')
            local=np.vstack((Tc.T@cell['Z'][sl],Tc.T@cell['Z'][cell['nc']+sl.start:cell['nc']+sl.stop]))
            if 'seen' in r:
                prior,block=r['seen'];G=np.zeros((2*nd,ny));G[:,cell['yoffset']:offset]=local;G[:,prior:prior+block.shape[1]]-=block;equations.append(G)
            else:r['seen']=(cell['yoffset'],local)
            if key[0]==0 and key[1]==0:
                G=np.zeros((2*nd,ny));G[:,cell['yoffset']:offset]=local;equations.append(G)
        cell['trace_ids']=np.array(nodes);cell['trace_transform']=np.zeros((cell['nb'],cell['nb']));p=0
        for T in transforms:cell['trace_transform'][p:p+len(T),p:p+len(T)]=T;p+=len(T)
    G=np.vstack(equations);Yknown=[];representation=0.;names=None
    for cell in cells:
        s=t.f.prepare(cell['q'],.01);field,names=known.native_global_fields(cell['q'],kappa,s);nc=cell['nc'];n0=cell['n0'];y0=cell['Z'][:nc,:n0].T@field[:nc];res=field[nc:]-cell['Z'][nc:,:n0]@y0;yh=cell['Z'][nc:,n0:].T@res;y=np.vstack((y0,yh));representation=max(representation,float(np.linalg.norm(cell['Z']@y-field)/max(1.,np.linalg.norm(field))));Yknown.append(y)
    Yknown=np.vstack(Yknown);scale=np.linalg.norm(Yknown,axis=0);Q,R=qr(Yknown/scale,mode='full',check_finite=False);nk=len(scale);Q0=Q[:,:nk];Qp=Q[:,nk:];Zc,info=t.nullspace(G@Qp);Z=np.column_stack((Q0,Qp@Zc));residual=float(np.linalg.norm(G@Z)/max(1.,np.linalg.norm(G)));known_residual=float(np.linalg.norm(G@Yknown)/max(1.,np.linalg.norm(G)*np.linalg.norm(Yknown)))
    if max(representation,residual,known_residual)>1e-9:raise ValueError('MMS1_known_shared_kernel_residual:'+str((representation,residual,known_residual)))
    info.update(complement_columns=ny-nk,columns=ny,nullity=Z.shape[1],residual=residual,known_compatible_fields=names,known_kernel_residual=known_residual,known_local_representation=representation)
    free=np.setdiff1d(np.arange(trace_count),fixed);free=free[free%3==2];return Z,info,free,trace_count
