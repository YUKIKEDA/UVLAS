"""Orthogonal constitutive/thickness and element-mean error diagnostics."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import q9_eas21 as method
import q9_developable_bending_approximation as base
import check_q9_plate_convergence as memory
ROOT=Path(__file__).resolve().parents[2]
LABELS=['membrane','bending','higher_thickness','transverse_normal','transverse_shear']
MU=method.common.mu;LAM=method.common.lam
PS=method.C[np.ix_([0,1,5],[0,1,5])]-np.outer([LAM,LAM,0.],[LAM,LAM,0.])/(LAM+2*MU)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def eq(p):return np.einsum('...i,ij,...j->...',p,PS,p)
def norm2(e):return np.einsum('...i,ij,...j->...',e,method.C,e)
def split(error,weight,x,h,points,order,kappa):
    D=1+kappa*kappa*x*x;frame=np.zeros((len(x),3,3));frame[:,0,0]=frame[:,2,2]=1/np.sqrt(D);frame[:,1,1]=1;frame[:,2,0]=kappa*x/np.sqrt(D);frame[:,0,2]=-kappa*x/np.sqrt(D)
    tensor=np.einsum('kij,pk->pij',method.MB,error);local=frame.transpose(0,2,1)@tensor@frame
    p=np.stack((local[:,0,0],local[:,1,1],np.sqrt(2.)*local[:,0,1]),axis=1).reshape(order,order,order,3)
    w=weight.reshape(order,order,order);z=(h*points[:,2]/2).reshape(order,order,order);w0=w.sum(axis=2);z0=(w*z).sum(axis=2)/w0;zc=z-z0[:,:,None];w2=(w*zc*zc).sum(axis=2)
    p0=np.einsum('ijk,ijka->ija',w,p)/w0[:,:,None];p1=np.einsum('ijk,ijk,ijka->ija',w,zc,p)/w2[:,:,None];rest=p-p0[:,:,None,:]-zc[:,:,:,None]*p1[:,:,None,:]
    normal=local[:,2,2]+LAM/(LAM+2*MU)*(local[:,0,0]+local[:,1,1]);shear=2*MU*(2*local[:,0,2]**2+2*local[:,1,2]**2)
    parts=np.array([np.sum(w0*eq(p0)),np.sum(w2*eq(p1)),np.sum(w*eq(rest)),weight@((LAM+2*MU)*normal**2),weight@shear])
    mean=(weight@error)/weight.sum();mean_energy=float(weight.sum()*norm2(mean));fluct=float(weight@norm2(error-mean));total=float(weight@norm2(error))
    # The frame varies across a cell: this is a separate local-component statistic.
    pm=np.einsum('ij,ija->a',w0,p0)/w0.sum();membrane_mean=float(w0.sum()*eq(pm));membrane_fluct=float(np.sum(w0*eq(p0-pm)))
    return dict(total=total,parts=parts,mean=mean_energy,fluctuation=fluct,membrane_mean=membrane_mean,membrane_fluctuation=membrane_fluct,reconstruction=float(abs(parts.sum()-total)/max(total,1e-300)),mean_reconstruction=float(abs(mean_energy+fluct-total)/max(total,1e-300)),membrane_reconstruction=float(abs(membrane_mean+membrane_fluct-parts[0])/max(total,1e-300)))
def observe(n,h,u,kappa=.4):
    xy,cells,coords,interp,_=base.mesh(kappa,n,.2);fields={'saved_raw':u,'saved_eas':u,'interpolation_raw':interp,'interpolation_eas':interp};totals={};cell_logs=[];reference={9:0.,11:0.};maxmean=0.;checks=True
    for ci,(_,_,ids,_) in enumerate(cells):
        state=method.prepare(coords[ids].reshape(81),h);observed={}
        for order in [9,11]:
            points,w,N,_=base.table(order);J,B,R,_=method.evaluate(state,points);weight=w*np.linalg.det(J);assert np.all(weight>0);x=(N@xy[ids])[:,0];star=base.target(kappa,x,h*points[:,2]/2);ref=float(weight@norm2(star));reference[order]+=ref
            for name,field in fields.items():
                operator=R if name.endswith('eas') else B;error=operator@field[ids].reshape(81)-star;log=split(error,weight,x,h,points,order,kappa);observed[name,order]=log
                key=(name,order)
                if key not in totals:totals[key]=dict(total=0.,parts=np.zeros(5),mean=0.,fluctuation=0.,membrane_mean=0.,membrane_fluctuation=0.)
                for item in totals[key]:totals[key][item]+=log[item]
                # Near-exact flat controls have negligible energies, so use reference scaling.
                checks=checks and abs(log['parts'].sum()-log['total'])<=1e-10*max(log['total'],ref*1e-20)
                checks=checks and abs(log['mean']+log['fluctuation']-log['total'])<=1e-10*max(log['total'],ref*1e-20)
            for field in [u,interp]:
                diff=(R-B)@field[ids].reshape(81);mean=weight@diff/weight.sum();maxmean=max(maxmean,float(np.sqrt(weight.sum()*norm2(mean)/ref)))
        cell_logs.append(dict(cell=ci,**{name:{'normalized_parts':(observed[name,11]['parts']/ref).tolist(),'mean_fraction':observed[name,11]['mean']/max(observed[name,11]['total'],1e-300)} for name in fields}))
    logs={};maxquad=0.
    for name in fields:
        l=totals[name,11];q=totals[name,9];ref=reference[11];maxquad=max(maxquad,float(np.max(abs(l['parts']/ref-q['parts']/reference[9]))))
        logs[name]=dict(error=float(np.sqrt(l['total']/ref)),energy_parts={label:float(v/ref) for label,v in zip(LABELS,l['parts'])},fractions={label:float(v/max(l['total'],1e-300)) for label,v in zip(LABELS,l['parts'])},mean_fraction=float(l['mean']/max(l['total'],1e-300)),normalized_mean=float(l['mean']/ref),normalized_fluctuation=float(l['fluctuation']/ref),membrane_mean_fraction=float(l['membrane_mean']/max(l['parts'][0],1e-300)),normalized_membrane_mean=float(l['membrane_mean']/ref),normalized_membrane_fluctuation=float(l['membrane_fluctuation']/ref))
    return dict(n=n,h=h,kappa=kappa,routes=logs,quadrature_difference=maxquad,mean_preservation_error=maxmean,reconstruction_pass=bool(checks),cells=cell_logs)
def main():
    dest=ROOT/'docs/verification/q9-eas21-residual-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parent=ROOT/'docs/verification/q9-eas21-shared-approximation-result.json';data=json.loads(parent.read_text());records=[]
    for path,sha in data['sources'].items():assert digest(ROOT/path)==sha
    for record in data['records']:
        path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];u=np.load(path)['u'];log=observe(record['n'],record['h'],u);log['parent_error_difference']=abs(log['routes']['saved_eas']['error']-record['observation']['error']);log['checks']=dict(reconstruction=log['reconstruction_pass'],quadrature=log['quadrature_difference']<=1e-7,mean_preservation=log['mean_preservation_error']<=1e-8,parent=log['parent_error_difference']<=1e-8);assert all(log['checks'].values());records.append(log);print(log['n'],log['h'],log['routes']['saved_eas'],flush=True)
    _,_,_,u,_=base.mesh(0.,2,.2);control=observe(2,.0001,u,0.);assert control['routes']['saved_eas']['error']<=1e-8
    sources=[Path(__file__),Path(method.__file__),Path(base.__file__),parent,ROOT/'docs/verification/q9-eas21-residual-experiment.md']
    out=dict(records=records,control=control,all_checks_pass=True,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start,peak_bytes=memory.memory_bytes());dest.write_text(json.dumps(out,indent=2),encoding='utf-8')
if __name__=='__main__':main()
