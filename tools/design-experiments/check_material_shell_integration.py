"""Material input, effective laws, section work, and adopted full PG space."""
from pathlib import Path
from math import comb
from types import SimpleNamespace
import copy,gc,hashlib,json,time,traceback
import numpy as np
import material_equivalent_input as inp
import mass_paired_ancf as shell
import exact_paired_finite as engine
import pg_ritz_mechanics as kernel
ROOT=Path(__file__).resolve().parents[2]
def rel(a,b):return float(np.linalg.norm(a)/max(np.linalg.norm(b),1e-300))
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,allow_nan=False).encode()).hexdigest()
def requests():
    common=dict(inp.TAGS,e0=[.001,-.0003,.0002,.0001,-.00005,.00008])
    iso=dict(common,kind='isotropic',E=12.,nu=.29,rho=1.2)
    orth=dict(common,kind='orthotropic',E1=90.,E2=9.,E3=6.,nu12=.27,nu13=.23,nu23=.29,G23=2.5,G13=3.5,G12=4.,rho=.8)
    v=np.array([1.,.5,-.3,.2,.1,.4]);C=np.diag([30.,18.,11.,5.,7.,8.])+3*np.outer(v,v)
    eff=dict(common,kind='effective_3d',C=C.tolist(),rho=1.6,e0=[-.0004,.0006,-.0002,.00003,.0001,-.0002])
    return [dict(id=f'L{i}',lower=lo,upper=hi,fiber_global=[np.cos(theta),np.sin(theta),.1],material=m)
        for i,(lo,hi,theta,m) in enumerate(zip(['-1','-.37','.22'],['-.37','.22','1'],[.23,-.41,.67],[iso,orth,eff]))]

def independent_Ceng(request):
    if request['kind']=='isotropic':
        E=request['E'];nu=request['nu'];mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu));C=np.diag([2*mu]*3+[mu]*3);C[:3,:3]+=lam;return C
    if request['kind']=='orthotropic':
        E1,E2,E3=[request[k] for k in ['E1','E2','E3']]
        H=np.diag([1/E1,1/E2,1/E3,1/request['G23'],1/request['G13'],1/request['G12']])
        H[0,1]=H[1,0]=-request['nu12']/E1;H[0,2]=H[2,0]=-request['nu13']/E1;H[1,2]=H[2,1]=-request['nu23']/E2
        return np.linalg.inv(H)
    C=np.array(request['C']);D=np.diag([1,1,1,1/np.sqrt(2),1/np.sqrt(2),1/np.sqrt(2)])
    return D@C@D if request['components'].startswith('mandel') else C

def direct_section(req,h,q,v):
    energy=0.;force=np.zeros((3,6));massenergy=0.;nodes,weights=np.polynomial.legendre.leggauss(5)
    for l in req:
        m=l['material'];lo=float(l['lower'])*h/2;hi=float(l['upper'])*h/2
        theta=np.arctan2(l['fiber_global'][1],l['fiber_global'][0]);c=np.cos(theta);s=np.sin(theta);A=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
        e0=inp.tensor(np.array(m['e0']));C=independent_Ceng(m)
        for z,w in zip((lo+hi)/2+(hi-lo)*nodes/2,weights*(hi-lo)/2):
            powers=np.array([1.,z,z*z]);E=A.T@inp.tensor(powers@q)@A-e0
            ee=np.array([E[0,0],E[1,1],E[2,2],2*E[1,2],2*E[0,2],2*E[0,1]]);ss=C@ee
            Sm=np.array([[ss[0],ss[5],ss[4]],[ss[5],ss[1],ss[3]],[ss[4],ss[3],ss[2]]]);S=A@Sm@A.T
            sm=np.einsum('kij,ij->k',inp.MB,S)
            energy+=w*(ee@ss)/2;force+=w*powers[:,None]*sm;massenergy+=w*m['rho']*np.sum((powers@v)**2)/2
    return energy,force,massenergy

def data_for(model,stack,order=45,zorder=7):
    plane,wp=shell.geometry.quadrature(model['n'],order);z,wz=np.polynomial.legendre.leggauss(zorder)
    blocks={key:[] for key in ['p','weights','N','G','C','e0','rho','X','A','T','layer']}
    for li,l in enumerate(stack):
        lo=float(l.lower);hi=float(l.upper);zz=(lo+hi)/2+(hi-lo)*z/2
        for begin in range(0,len(plane),128):
            xy=plane[begin:begin+128];F=shell.fields(model,xy)
            p,N,G,det,X=engine.kinematics_from_fields(model,xy,zz,model['h'],model['variable'],F)
            J=shell.geometry.reference(p[:,:2],p[:,2],model['h'],model['parameters'],model['variable'])[1]
            C,e0,A,T=inp.point_fields(l,J);ww=(wp[begin:begin+128,None]*wz*(hi-lo)/2).ravel()*det
            values=dict(p=p,N=N,G=G,C=C,e0=e0,A=A,T=T,weights=ww,rho=np.full(len(p),l.material.rho),X=X,layer=np.full(len(p),li))
            for key in blocks:blocks[key].append(values[key])
    data={key:np.concatenate(value) for key,value in blocks.items()};data['M']=data['N'].T@((data['weights']*data['rho'])[:,None]*data['N'])
    return data

def independent_point(data,u,req):
    F=np.eye(3)+np.einsum('na,pnb->pab',u,data['G']);E=(F.transpose(0,2,1)@F-np.eye(3))/2;P=np.empty_like(F);W=np.zeros(len(F))
    for li,l in enumerate(req):
        mask=data['layer']==li;A=data['A'][mask];e=A.transpose(0,2,1)@E[mask]@A-inp.tensor(np.array(l['material']['e0']))
        eng=np.stack([e[:,0,0],e[:,1,1],e[:,2,2],2*e[:,1,2],2*e[:,0,2],2*e[:,0,1]],axis=1)
        stress=eng@independent_Ceng(l['material']).T;S=np.zeros((sum(mask),3,3));S[:,0,0]=stress[:,0];S[:,1,1]=stress[:,1];S[:,2,2]=stress[:,2]
        for k,(i,j) in enumerate([(1,2),(0,2),(0,1)],3):S[:,i,j]=S[:,j,i]=stress[:,k]
        P[mask]=F[mask]@A@S@A.transpose(0,2,1);W[mask]=np.sum(eng*stress,axis=1)/2
    U=np.sum(data['weights']*W);force=np.einsum('p,pnb,pab->na',data['weights'],data['G'],P,optimize=True)
    return U,force

def main():
    dest=ROOT/'docs/verification/material-shell-integration-result.json';assert not dest.exists();records=[];failures=[];start=time.monotonic();artifact=None
    req=requests();original=digest(req);stack=inp.prepare_stack(req);h=.05;rng=np.random.default_rng(20261105)
    try:
        equivalent=copy.deepcopy(req)
        for l,m in zip(equivalent,stack):l['material']=dict(inp.TAGS,kind='effective_3d',C=m.material.C.tolist(),rho=m.material.rho,e0=m.material.e0.tolist())
        equivalent_stack=inp.prepare_stack(equivalent);assert inp.identity(stack)==inp.identity(equivalent_stack)
        engineering=copy.deepcopy(equivalent);D=np.diag([1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2)])
        for l in engineering:
            m=l['material'];m['C']=(np.linalg.solve(D,np.array(m['C']))@np.linalg.inv(D)).tolist();m['e0']=(D@np.array(m['e0'])).tolist();m['components']='engineering_xx_yy_zz_yz_xz_xy'
        eng_stack=inp.prepare_stack(engineering);bridge=max(rel(a.material.C-b.material.C,a.material.C) for a,b in zip(stack,eng_stack));assert bridge<1e-12
        records.append(dict(case='explicit_3d_input_and_engineering_bridge',canonical_effective_identity=True,bridge_error=bridge))
        C,g,rho,w0=inp.section_moments(stack,h);K=np.block([[C[i+j] for j in range(3)] for i in range(3)]);load=g.ravel();mass=np.kron(np.array([[rho[i+j] for j in range(3)] for i in range(3)]),np.eye(3))
        q=rng.normal(size=(3,6))*np.array([.003,.04,.2])[:,None];v=rng.normal(size=(3,3));U,f,ke=direct_section(req,h,q,v)
        metrics=dict(energy=rel(np.array(U-(q.ravel()@K@q.ravel()/2-q.ravel()@load+w0)),np.array(U)),force=rel(f.ravel()-(K@q.ravel()-load),f),mass=rel(np.array(ke-v.ravel()@mass@v.ravel()/2),np.array(ke)))
        assert max(metrics.values())<=1e-10;records.append(dict(case='quadratic_strain_velocity_section_work',metrics=metrics))
        cert=dict(reference='mid_surface',units='SI',h=h,C_moments=C.tolist(),initial_moments=g.tolist(),mass_moments=rho.tolist(),initial_constant=w0)
        assert inp.verify_certificate(stack,h,cert)
        shift=.007;Cp=[sum(comb(k,j)*(-shift)**(k-j)*C[j] for j in range(k+1)) for k in range(5)];gp=[sum(comb(k,j)*(-shift)**(k-j)*g[j] for j in range(k+1)) for k in range(3)]
        qp=np.array([q[0]+shift*q[1]+shift**2*q[2],q[1]+2*shift*q[2],q[2]])
        Kp=np.block([[Cp[i+j] for j in range(3)] for i in range(3)]);reference=qp.ravel()@Kp@qp.ravel()/2-qp.ravel()@np.array(gp).ravel()+w0
        error=abs(reference-U)/abs(U);assert error<1e-10;records.append(dict(case='reference_plane_shift',energy_error=error))
        # Same ABD and transverse shear, different normal elasticity.
        C1=np.diag([10.,12.,3.,4.,5.,6.]);C1[0,1]=C1[1,0]=1.;C2=C1.copy();C2[2,2]=30.;ia=[0,1,5]
        assert np.array_equal(C1[np.ix_(ia,ia)],C2[np.ix_(ia,ia)]) and np.array_equal(C1[3:5,3:5],C2[3:5,3:5])
        e=np.array([0.,0.,.01,0.,0.,0.]);difference=float(e@(C2-C1)@e/2);assert difference>1e-4
        # Same areal mass and mass centroid do not fix the second mass moment.
        base=copy.deepcopy(req[0]);base['lower']='-1';base['upper']='1';base['material']['rho']=1.;uniform=inp.prepare_stack([base]);variable=[]
        for i,(lo,hi,density) in enumerate([('-1','-1/3',1.25),('-1/3','1/3',.5),('1/3','1',1.25)]):
            l=copy.deepcopy(base);l.update(id=f'm{i}',lower=lo,upper=hi);l['material']['rho']=density;variable.append(l)
        ma=inp.section_moments(uniform,h)[2];mb=inp.section_moments(inp.prepare_stack(variable),h)[2]
        assert abs(ma[0]-mb[0])<1e-14 and abs(ma[1]-mb[1])<1e-14 and abs(ma[2]-mb[2])>1e-7
        records.append(dict(case='nonunique_ABD_and_areal_mass',same_ABD_shear=True,normal_energy_difference=difference,second_mass_moment_difference=float(mb[2]-ma[2])))
        # A normal reversal changes the material coordinate basis, too.
        A=np.array([[np.cos(.37),-np.sin(.37),0.],[np.sin(.37),np.cos(.37),0.],[0.,0.,1.]]);Q=np.diag([1.,-1.,-1.]);T=inp.transform(A);Tq=inp.transform(Q);m=stack[2].material
        transformed=inp.transform(A@Q);corrected=transformed@(Tq.T@m.C@Tq)@transformed.T;unmodified=transformed@m.C@transformed.T
        error=rel(corrected-T@m.C@T.T,T@m.C@T.T);negative=rel(unmodified-T@m.C@T.T,T@m.C@T.T);assert error<1e-12 and negative>1e-3
        records.append(dict(case='normal_and_material_basis_reversal',tensor_error=error,forgotten_component_transform_defect=negative))
        # Explicit negative inputs, each with unchanged caller data.
        negatives=[]
        def reject(name,fn):
            try:fn()
            except (ValueError,KeyError):negatives.append(name)
            else:raise AssertionError('accepted '+name)
        for key,value in [('units','mm_N'),('components','unknown'),('measure','small_strain'),('rho',0.),('rho',-1.),('e0',[0,0,-1,0,0,0]),('e0',[float('nan')]*6),('C',np.eye(5).tolist())]:
            m=copy.deepcopy(equivalent[2]['material']);m[key]=value;reject(key+'_'+str(value)[:20],lambda m=m:inp.prepare_material(m))
        for name,Cbad in [('nonsymmetric',np.eye(6)),('indefinite',np.diag([1,1,-1,1,1,1])),('near_null',np.ones((6,6)))]:
            if name=='nonsymmetric':Cbad[0,1]=.1
            m=copy.deepcopy(equivalent[2]['material']);m['C']=Cbad.tolist();reject(name,lambda m=m:inp.prepare_material(m))
        for kind in ['abd','abd_shear']:
            m=dict(inp.TAGS,kind=kind,rho=1.,e0=[0.]*6);reject(kind,lambda m=m:inp.prepare_material(m))
        m=copy.deepcopy(req[1]['material']);del m['E3'];reject('missing_E3',lambda:inp.prepare_material(m))
        m=copy.deepcopy(req[0]['material']);m['nu']=.5;reject('incompressible_input',lambda:inp.prepare_material(m))
        for name,index,key,value in [('gap',1,'lower','-.36'),('overlap',1,'lower','-.38'),('reverse',1,'upper','-.8'),('incomplete',2,'upper','.9'),('duplicate_id',1,'id','L0'),('zero_fiber',0,'fiber_global',[0,0,0])]:
            bad=copy.deepcopy(req);bad[index][key]=value;before=repr(bad);reject(name,lambda bad=bad:inp.prepare_stack(bad));assert repr(bad)==before
        reject('parallel_fiber',lambda:inp.axes(np.eye(3)[None],np.array([0.,0.,1.])))
        for name,key,value in [('certificate_h','h',.051),('certificate_reference','reference','offset'),('certificate_units','units','mm_N'),('certificate_moment','initial_constant',w0*1.1)]:
            bad=copy.deepcopy(cert);bad[key]=value;reject(name,lambda bad=bad:inp.verify_certificate(stack,h,bad))
        reject('outside_layer',lambda:inp.layer_at(stack,1.1));reject('side_without_interface',lambda:inp.layer_at(stack,0.,'upper'))
        assert len(negatives)==28;records.append(dict(case='invalid_inputs',rejected=negatives,count=len(negatives)))
        model=shell.prepare(1,(.4,.2,.07,.1,-.06,.025),h,True,'extended');q0=model['reference_coordinates'];data=data_for(model,stack)
        # Explicit effective input must give the same physical-point law, including E0/rho/axes.
        defect=0.
        for li,l in enumerate(equivalent_stack):
            mask=data['layer']==li;p=data['p'][mask][::103];J=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],True)[1]
            Ce,e0e,Ae,Te=inp.point_fields(l,J)
            defect=max(defect,rel(Ce-data['C'][mask][::103],Ce),rel(e0e-data['e0'][mask][::103],e0e));assert l.material.rho==stack[li].material.rho
        assert defect<1e-12
        direction=rng.normal(size=q0.shape);direction*=.003/np.max(abs(np.einsum('na,pnb->pab',direction,data['G'],optimize=True)))
        u=q0@(kernel.FSTAR-np.eye(3)).T+direction;U,f,K,info=kernel.evaluate(data,u,True);eps=1e-28;Uc,fc,_,_=kernel.evaluate(data,u+1j*eps*direction)
        Ui,fi=independent_point(data,u,req)
        metrics=dict(effective_point_law=defect,independent_tensor_energy=abs(Ui-U)/abs(Ui),independent_tensor_force=rel(fi-f,fi),
            energy_direction=abs(Uc.imag/eps-np.sum(f*direction))/abs(np.sum(f*direction)),tangent_direction=rel(fc.imag/eps-(K@direction.ravel()).reshape(u.shape),fc.imag/eps))
        assert max(metrics.values())<=1e-8;records.append(dict(case='three_input_kinds_same_adopted_S',metrics=metrics,ns=model['ns'],min_det_F=info['min_det']))
        print('SHELL_LAW',metrics,flush=True)
        # Independent physical inertia: different thickness quadrature, analytic X, no N@q0.
        N=data['N'];ww=data['weights']*data['rho'];translation=np.linalg.solve(data['M'],N.T@ww);rigid=np.zeros((model['ns'],3,6));rigid[:,:,:3]=translation[:,None,None]*np.eye(3)[None]
        for j in range(3):rigid[:,:,3+j]=np.cross(np.eye(3)[j],q0)
        fullM=np.kron(data['M'],np.eye(3));Qrigid=rigid.reshape(3*model['ns'],6);direct=np.zeros((6,6));plane,wp=shell.geometry.quadrature(1,45);z,wz=np.polynomial.legendre.leggauss(9)
        for l in stack:
            lo=float(l.lower);hi=float(l.upper);p=engine.points(plane,(lo+hi)/2+(hi-lo)*z/2)
            X,J,det,H,Ha=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],True);w=(wp[:,None]*wz*(hi-lo)/2).ravel()*det*l.material.rho
            V=np.zeros((len(p),3,6));V[:,:,:3]=np.eye(3)[None]
            for j in range(3):V[:,:,3+j]=np.cross(np.eye(3)[j],X)
            direct+=np.einsum('p,pai,paj->ij',w,V,V)
        inertia=rel(Qrigid.T@fullM@Qrigid-direct,direct);assert inertia<1e-9
        records.append(dict(case='equivalent_density_and_physical_spatial_inertia',inertia_error=inertia,mass_eigen_bounds=np.linalg.eigvalsh(data['M'])[[0,-1]].tolist()))
        recovery=[]
        for zeta in ['-.37','.22']:
            lower=inp.layer_at(stack,zeta,'lower');upper=inp.layer_at(stack,zeta,'upper');assert inp.layer_at(stack,zeta).id==upper.id
            p=np.array([[.31,.47,float(zeta)],[.67,.73,float(zeta)]]);N,G,B,det,X=shell.kinematics(model,p,h,True);J=shell.geometry.reference(p[:,:2],p[:,2],h,model['parameters'],True)[1]
            F=np.eye(3)+np.einsum('na,pnb->pab',u,G);E=np.einsum('kij,pij->pk',inp.MB,(F.transpose(0,2,1)@F-np.eye(3))/2);values=[]
            for l in [lower,upper]:
                Cg,e0,A,T=inp.point_fields(l,J);values.append(np.einsum('pij,pj->pi',Cg,E-e0))
            difference=float(np.linalg.norm(values[0]-values[1]));assert difference>1e-6
            recovery.append(dict(zeta=zeta,lower=lower.id,upper=upper.id,stress_difference=difference))
        records.append(dict(case='unequal_layers_two_sided_recovery',interfaces=recovery))
        assert digest(req)==original and inp.identity(stack)==inp.identity(equivalent_stack)
        try:stack[0].material.C[0,0]=1.
        except ValueError:pass
        else:raise AssertionError('writable prepared material')
        records.append(dict(case='immutable_input_and_canonical_identity',input_unchanged=True,prepared_matrix_readonly=True,physics_identity_equal=True))
        out=ROOT/'build/analysis'/('material-shell-integration-'+time.strftime('%Y%m%d-%H%M%S')+'.npz')
        np.savez_compressed(out,section_C=C,section_g=g,section_mass=rho,section_K=K,fullM=fullM,q0=q0,u=u,force=f,spatial_inertia=direct)
        artifact=dict(artifact=out.relative_to(ROOT).as_posix(),artifact_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
    except Exception as exc:failures.append(dict(reason=repr(exc),traceback=traceback.format_exc()));print('FAIL',failures[-1],flush=True)
    paths=[Path(__file__),Path(inp.__file__),Path(shell.__file__),Path(engine.__file__),Path(kernel.__file__),ROOT/'docs/verification/material-shell-integration-protocol.md']
    result=dict(passed=len(records)==10 and not failures,records=records,failures=failures,artifact=artifact,seconds=time.monotonic()-start,
        material_input=req,real_material_validation=False,element_plane_stress_adopted=False,
        sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print('MATERIAL_INTEGRATION_RESULT',result['passed'],len(records),flush=True)

if __name__=='__main__':main()
