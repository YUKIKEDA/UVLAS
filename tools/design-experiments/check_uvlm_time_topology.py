"""Limited UVLM interface experiment; protocol in docs/verification/uvlm-time-topology-experiment.md."""
import copy
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import check_q9_plate_convergence as memory

ROOT=Path(__file__).resolve().parents[2]
START=time.monotonic()
DIMS=[(2,2),(1,3)]
OFF=[0,4]
RHO=1.2
FLOW=np.array([8.,0.,0.])
Z64,W64=leggauss(64)


def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def clone(value):return json.loads(json.dumps(value))


def gamma(t):
    p=np.arange(7);return 1+.1*p+(.3+.02*p)*t+(.2+.01*p)*t*t


def gamma_dot(t):return .3+.02*np.arange(7)+2*(.2+.01*np.arange(7))*t


def bound(s,t):
    m,n=DIMS[s];x=np.linspace(0,1,m+1)+(.2 if s else 0)
    y=np.linspace(-.5,.5,n+1)+(2. if s else 0)
    points=np.array([[[a,b,.5*s] for b in y] for a in x])
    c=np.cos(.05*t);d=np.sin(.05*t);rotation=np.array([[c,0,d],[0,1,0],[-d,0,c]])
    return points@rotation.T+np.array([0,0,.02*t])


def te(g,s):
    m,n=DIMS[s];return np.array(g)[OFF[s]+(m-1)*n:OFF[s]+m*n]


def initial():
    surfaces=[]
    for s in range(2):
        nodes=[bound(s,0)[-1]]+[bound(s,t)[-1]-t*FLOW for t in [-.02,-.04]]
        surfaces.append(dict(nodes=np.array(nodes).tolist(),circulation=[te(gamma(t),s).tolist() for t in [-.02,-.04]],births=[-1,-2],
            intervals=[[-.02,0.],[-.04,-.02]]))
    return dict(t=0.,generation=0,gamma=gamma(0).tolist(),gamma_previous=gamma(-.02).tolist(),h_previous=.02,wake=surfaces)


def segment(p,a,b):
    r1=p-a;r2=p-b;cross=np.cross(r1,r2);den=cross@cross
    if den<1e-26 or min(np.linalg.norm(r1),np.linalg.norm(r2))<1e-12:raise ValueError('singular kernel outside fixture')
    return cross*((b-a)@(r1/np.linalg.norm(r1)-r2/np.linalg.norm(r2)))/(4*np.pi*den)


def quad_segment(p,a,b):
    t=(Z64+1)/2;dl=b-a;r=p[None,:]-(a[None,:]+t[:,None]*dl)
    return np.sum(W64[:,None]*np.cross(dl,r)/np.linalg.norm(r,axis=1)[:,None]**3,axis=0)/(8*np.pi)


def ring_velocity(p,points,fn=segment):return sum((fn(p,points[i],points[(i+1)%4]) for i in range(4)),np.zeros(3))


def ring_nodes(grid,i,j):return np.array([grid[i,j],grid[i,j+1],grid[i+1,j+1],grid[i+1,j]])


def bdf(current,old,previous,h,k,order=2):
    if order==1:return (current-old)/h
    if previous is None or k is None or k<=0:raise ValueError('missing BDF2 history')
    return (2*h+k)/(h*(h+k))*current-(h+k)/(h*k)*old+h/(k*(h+k))*previous


def topology(state):
    coords={};faces=[];strength=[];bound_faces=[];te_edges=set()
    for s,(m,n) in enumerate(DIMS):
        grid=bound(s,state['t']);wake=np.array(state['wake'][s]['nodes']);births=state['wake'][s]['births']
        for i in range(m+1):
            for j in range(n+1):coords[(s,'b',i,j)]=grid[i,j]
        for r,birth in enumerate(births):
            for j in range(n+1):coords[(s,'w',birth,j)]=wake[r+1,j]
        for i in range(m):
            for j in range(n):
                face=[(s,'b',i,j),(s,'b',i,j+1),(s,'b',i+1,j+1),(s,'b',i+1,j)]
                g=state['gamma'][OFF[s]+i*n+j];faces.append(face);strength.append(g);bound_faces.append((face,g))
        for j in range(n):te_edges.add(((s,'b',m,j),(s,'b',m,j+1)))
        for r,birth in enumerate(births):
            front=[(s,'b',m,j) if r==0 else (s,'w',births[r-1],j) for j in range(n+1)]
            back=[(s,'w',birth,j) for j in range(n+1)]
            for j in range(n):
                faces.append([front[j],front[j+1],back[j+1],back[j]]);strength.append(state['wake'][s]['circulation'][r][j])
    edges=sorted({tuple(sorted((face[i],face[(i+1)%4]))) for face in faces for i in range(4)})
    vertices=sorted(coords);ei={e:i for i,e in enumerate(edges)};vi={v:i for i,v in enumerate(vertices)}
    d1=np.zeros((len(vertices),len(edges)),dtype=int);d2=np.zeros((len(edges),len(faces)),dtype=int)
    for i,(a,b) in enumerate(edges):d1[vi[a],i]=-1;d1[vi[b],i]=1
    for j,face in enumerate(faces):
        for i in range(4):
            a,b=face[i],face[(i+1)%4];edge=tuple(sorted((a,b)));d2[ei[edge],j]+=1 if (a,b)==edge else -1
    return coords,bound_faces,edges,d1,d2,np.array(strength),te_edges


def load_checks(state,dot):
    coords,faces,edges,d1,d2,g,te_edges=topology(state)
    nodes=sorted(coords);fnode={key:np.zeros(3) for key in nodes};direct_f=np.zeros(3);direct_m=np.zeros(3)
    def relative_u(x):return np.array([8+.1*x[0],.2*x[1],.05*x[2]])
    kappas={edge:0. for edge in edges}
    for face,circ in faces:
        for i in range(4):
            a,b=face[i],face[(i+1)%4];edge=tuple(sorted((a,b)))
            if edge in te_edges:continue
            kappas[edge]+=circ*(1 if (a,b)==edge else -1)
            x=.5*(coords[a]+coords[b]);f=RHO*circ*np.cross(relative_u(x),coords[b]-coords[a])
            direct_f+=f;direct_m+=np.cross(x,f)
    for (a,b),circ in kappas.items():
        if circ==0:continue
        x=.5*(coords[a]+coords[b]);f=RHO*circ*np.cross(relative_u(x),coords[b]-coords[a])
        fnode[a]+=.5*f;fnode[b]+=.5*f
    reference_area_force=np.zeros(3);omitted=np.zeros(3)
    for p,(face,circ) in enumerate(faces):
        points=np.array([coords[key] for key in face]);area=-.5*sum((np.cross(points[i],points[(i+1)%4]) for i in range(4)),np.zeros(3))
        force=RHO*dot[p]*area;center=points.mean(axis=0)
        reference_area_force+=force;direct_f+=force;direct_m+=np.cross(center,force)
        for key in face:
            fnode[key]+=.25*force
            if key[1]=='b' and key[2]==DIMS[key[0]][0]:omitted+=.25*force
    total=sum(fnode.values(),np.zeros(3));moment=sum((np.cross(coords[key],fnode[key]) for key in nodes),np.zeros(3))
    translation=np.array([.1,-.2,.3]);rotation=np.array([.04,.05,-.06])
    work=sum((translation+np.cross(rotation,coords[key]))@fnode[key] for key in nodes)
    expected=translation@direct_f+rotation@direct_m
    return dict(boundary_squared_zero=bool(np.all(d1@d2==0)),force_error=float(np.linalg.norm(total-direct_f)),
        moment_error=float(np.linalg.norm(moment-direct_m)),virtual_work_error=float(abs(work-expected)),
        omitted_trailing_vertices_force_difference=float(np.linalg.norm(omitted)),
        reversed_unsteady_normal_difference=float(2*np.linalg.norm(reference_area_force)))


def build_trial(old,h):
    t=old['t']+h;n=old['generation'];new=clone(old);new.update(t=t,generation=n+1)
    for s in range(2):
        w=old['wake'][s];nodes=np.array(w['nodes'])
        new['wake'][s]=dict(nodes=np.concatenate([bound(s,t)[-1][None,:,:],nodes+h*FLOW]).tolist(),
            circulation=[te(old['gamma'],s).tolist()]+clone(w['circulation']),births=[n]+w['births'],
            intervals=[[old['t'],t]]+clone(w['intervals']))
    bound_rings=[];normals=[];centers=[];wake_rings=[];wake_gamma=[];wake_owner=[]
    for s,(m,ns) in enumerate(DIMS):
        grid=bound(s,t)
        for i in range(m):
            for j in range(ns):
                ring=ring_nodes(grid,i,j);bound_rings.append(ring);centers.append(ring.mean(axis=0))
                normal=np.cross(ring[3]-ring[0],ring[1]-ring[0]);normals.append(normal/np.linalg.norm(normal))
        grid=np.array(new['wake'][s]['nodes'])
        for i in range(len(grid)-1):
            for j in range(ns):
                wake_rings.append(ring_nodes(grid,i,j));wake_gamma.append(new['wake'][s]['circulation'][i][j]);wake_owner.append(s)
    ab=np.array([[ring_velocity(c,ring)@normal for ring in bound_rings] for c,normal in zip(centers,normals)])
    aw=np.array([[ring_velocity(c,ring)@normal for ring in wake_rings] for c,normal in zip(centers,normals)])
    ab_ref=np.array([[ring_velocity(c,ring,quad_segment)@normal for ring in bound_rings] for c,normal in zip(centers,normals)])
    aw_ref=np.array([[ring_velocity(c,ring,quad_segment)@normal for ring in wake_rings] for c,normal in zip(centers,normals)])
    target=gamma(t);wg=np.array(wake_gamma);rhs=-aw_ref@wg-ab_ref@target
    solution=np.linalg.solve(ab,-aw@wg-rhs)
    new['gamma']=solution.tolist();new['gamma_previous']=old['gamma'];new['h_previous']=h
    dot=bdf(solution,np.array(old['gamma']),np.array(old['gamma_previous']),h,old['h_previous'])
    metrics=load_checks(new,dot)
    owner=np.array([0]*4+[1]*3)
    omitted_ab=ab*(owner[:,None]==owner[None,:]);omitted_aw=aw*(owner[:,None]==np.array(wake_owner)[None,:])
    omitted_gamma=np.linalg.solve(omitted_ab,-omitted_aw@wg-rhs)
    metrics.update(aic_quadrature_error=float(max(np.max(np.abs(ab-ab_ref)),np.max(np.abs(aw-aw_ref)))),
        omitted_inter_surface_gamma_error=float(np.max(np.abs(omitted_gamma-target))))
    metrics.update(circulation_error=float(np.max(np.abs(solution-target))),
        nonpenetration_relative=float(np.linalg.norm(ab@solution+aw@wg+rhs)/max(1.,np.linalg.norm(rhs))),
        bdf2_error=float(np.max(np.abs(dot-gamma_dot(t)))),
        bdf1_known_error=float(np.max(np.abs((solution-np.array(old['gamma']))/h-gamma_dot(t)+(.2+.01*np.arange(7))*h))),
        fixed_step_bdf2_mismatch=float(np.max(np.abs((1.5*solution-2*np.array(old['gamma'])+.5*np.array(old['gamma_previous']))/h-gamma_dot(t)))),
        newest_wrong_time_mismatch=float(max(np.max(np.abs(te(solution,s)-te(old['gamma'],s))) for s in range(2))))
    coords,faces,edges,d1,d2,g,te_edges=topology(new);kappa=d2@g;index={edge:i for i,edge in enumerate(edges)}
    te_error=0.
    for s,(m,ns) in enumerate(DIMS):
        for j in range(ns):
            edge=((s,'b',m,j),(s,'b',m,j+1))
            expected=new['wake'][s]['circulation'][0][j]-te(solution,s)[j]
            te_error=max(te_error,abs(kappa[index[edge]]-expected))
    metrics['trailing_induction_error']=float(te_error)
    geometry_error=0.;transport_error=0.
    for s,wake in enumerate(new['wake']):
        geometry_error=max(geometry_error,float(np.max(np.abs(np.array(wake['nodes'])[0]-bound(s,t)[-1]))))
        for row,interval in enumerate(wake['intervals']):
            born=interval[0]
            expected=bound(s,born)[-1]+(t-born)*FLOW
            geometry_error=max(geometry_error,float(np.max(np.abs(np.array(wake['nodes'])[row+1]-expected))))
            transport_error=max(transport_error,float(np.max(np.abs(np.array(wake['circulation'][row])-te(gamma(born),s)))))
    metrics.update(material_geometry_error=geometry_error,material_circulation_error=transport_error)
    steady=clone(new);steady['gamma']=[1.]*7
    for wake in steady['wake']:wake['circulation']=[[1.]*len(row) for row in wake['circulation']]
    _,_,se,_,sd,sg,ste=topology(steady);sk=sd@sg;si={edge:i for i,edge in enumerate(se)}
    metrics['steady_trailing_induction_error']=float(max(abs(sk[si[edge]]) for edge in ste))
    return dict(parent=digest(old),candidate=new,metrics=metrics)


def commit(old,trial):
    if trial['parent']!=digest(old):raise ValueError('stale_parent')
    return clone(trial['candidate'])


def main():
    known=segment(np.array([1.,0,0]),np.array([0.,-1,0]),np.array([0.,1,0]))*2*np.pi
    reference_error=float(np.max(np.abs(known-np.array([0,0,-1/np.sqrt(2)]))))
    kernel_error=0.
    for s in range(2):
        ring=ring_nodes(bound(s,0),0,0)
        for point in [np.array([2.,.3,1.]),np.array([-1.,-.2,2.]),np.array([3.,2.,-1.])]:
            kernel_error=max(kernel_error,float(np.max(np.abs(ring_velocity(point,ring)-ring_velocity(point,ring,quad_segment)))))
    state=initial();cases=[];unchanged=True;repeat_equal=True;double_refused=False
    saved=None
    for h in [.02,.013,.027]:
        parent=digest(state);trials=[build_trial(state,h) for _ in range(5)]
        unchanged=unchanged and digest(state)==parent
        repeat_equal=repeat_equal and len({digest(t['candidate']) for t in trials})==1
        # Discard all preceding candidates and accept exactly the final candidate.
        state=commit(state,trials[-1]);cases.append(trials[-1]['metrics'])
        try:commit(state,trials[-1])
        except ValueError:double_refused=True
        if state['generation']==2:saved=clone(state)
    restored=json.loads(json.dumps(saved));resumed=commit(restored,build_trial(restored,.027))
    history_refused=False
    try:bdf(gamma(.02),gamma(0),None,.02,None)
    except ValueError:history_refused=True
    linear_error=float(np.max(np.abs(bdf(np.array([1.04]),np.array([1.]),None,.02,None,1)-2)))
    old_nodes=np.array(saved['wake'][0]['nodes']);wrong_nodes=old_nodes.copy()
    for iteration in range(5):wrong_nodes+=.027*FLOW
    repeated_drift=float(np.max(np.abs(wrong_nodes-(old_nodes+.027*FLOW))))
    checks=dict(kernel_reference=reference_error<=1e-12,kernel_quadrature=kernel_error<=1e-10,
        topology=all(c['boundary_squared_zero'] and c['trailing_induction_error']<=1e-12 and c['steady_trailing_induction_error']<=1e-12 for c in cases),
        material_transport=all(c['material_geometry_error']<=1e-12 and c['material_circulation_error']<=1e-12 for c in cases),
        circulation=all(c['circulation_error']<=1e-12 and c['nonpenetration_relative']<=1e-12 for c in cases),
        aic_quadrature=all(c['aic_quadrature_error']<=1e-10 for c in cases),
        negative_omitted_inter_surface=all(c['omitted_inter_surface_gamma_error']>1e-12 for c in cases),
        derivatives=all(c['bdf2_error']<=1e-12 and c['bdf1_known_error']<=1e-12 for c in cases) and linear_error<=1e-12,
        force_work=all(max(c['force_error'],c['moment_error'],c['virtual_work_error'])<=1e-10 for c in cases),
        parent_unchanged=unchanged,repeat_trial_identical=repeat_equal,double_commit_refused=double_refused,
        checkpoint_equal=digest(state)==digest(resumed),final_rows=all(len(w['circulation'])==5 and len(set(w['births']))==5 for w in state['wake']),
        negative_time_transfer=all(c['newest_wrong_time_mismatch']>1e-12 for c in cases),
        negative_fixed_step_coefficients=cases[1]['fixed_step_bdf2_mismatch']>1e-12,
        negative_normal=all(c['reversed_unsteady_normal_difference']>1e-10 for c in cases),
        negative_omitted_vertices=all(c['omitted_trailing_vertices_force_difference']>1e-10 for c in cases),
        negative_repeated_convection=repeated_drift>1e-12,missing_history_refused=history_refused)
    source=Path(__file__)
    result=dict(protocol='docs/verification/uvlm-time-topology-experiment.md',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        numpy=np.__version__,memory_helper_sha256=hashlib.sha256(Path(memory.__file__).read_bytes()).hexdigest(),
        kernel_reference_error=reference_error,kernel_quadrature_error=kernel_error,cases=cases,cold_start_linear_error=linear_error,
        repeated_convection_drift_m=repeated_drift,checks=checks,all_criteria_pass=all(checks.values()),final_state=state,
        seconds=time.monotonic()-START,peak_process_working_bytes=memory.memory_bytes(),
        sources=['https://spiral.imperial.ac.uk/bitstreams/68d9db23-1616-4cc0-ab88-4f24970aff85/download',
            'https://orbi.uliege.be/bitstream/2268/200767/1/AIAAjournal_InducedDrag_accepted.pdf'],
        limitations='prescribed planar motion, manufactured normal RHS, uniform wake convection, prescribed edge load velocity; no physical aero-load or ANCF coupled solve, no core/near-field adoption')
    if result['seconds']>120 or result['peak_process_working_bytes']>1024**3:result['all_criteria_pass']=False
    (ROOT/'docs/verification/uvlm-time-topology-result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,seconds=result['seconds'],cases=cases)),flush=True)
    if not result['all_criteria_pass']:raise SystemExit(1)


if __name__=='__main__':main()
