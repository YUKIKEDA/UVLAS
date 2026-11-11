from pathlib import Path
import json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parents[2]
def initial():return dict(point=[1.,0.],tangent=(np.array([-1.,2.])/np.sqrt(5)).tolist(),steps=0,records=[])
def step(state):
    old=np.array(state['point']);t=np.array(state['tangent']);pred=old+.04*t;x=pred.copy()
    for it in range(10):
        z,U=x;R=np.array([z**3-z+U,t@(x-pred)]);J=np.array([[3*z*z-1,1],t]);dx=np.linalg.solve(J,-R);x+=dx
        if np.linalg.norm(dx)<1e-13:break
    else:raise ValueError('arc_Newton_unresolved')
    z,U=x;newt=np.array([-1.,3*z*z-1]);newt/=np.linalg.norm(newt)
    if newt@t<0:newt=-newt
    residual=abs(z**3-z+U);assert residual<1e-12 and newt@t>0;return dict(point=x.tolist(),tangent=newt.tolist(),steps=state['steps']+1,records=state['records']+[dict(z=float(z),U=float(U),stiffness=float(3*z*z-1),residual=float(residual))])
def main():
    s=initial()
    for _ in range(40):s=step(s)
    points=np.array([[r['z'],r['U']] for r in s['records']]);assert min(np.diff(points[:,1]))<0 and max(np.diff(points[:,1]))>0 and min(np.diff(points[:,0]))<0;assert any(r['stiffness']<0 for r in s['records']) and any(r['stiffness']>0 for r in s['records'])
    for stop in [7,17,29]:
        r=initial()
        for _ in range(stop):r=step(r)
        r=json.loads(json.dumps(r))
        for _ in range(stop,40):r=step(r)
        assert r==s
    foldz=1/np.sqrt(3);need_arc=abs(3*foldz*foldz-1)<1e-10;assert need_arc;paths=[Path(__file__),ROOT/'docs/j1-static-continuation-error-contract.md'];result=dict(passed=True,scope='continuation control on a known static fold; not a physical wing branch test',steps=40,max_equilibrium_residual=max(r['residual'] for r in s['records']),known_fold_U=float(foldz-foldz**3),direction_continued_through_fold=True,static_not_flutter=True,restart_bitwise=True,singular_fixed_speed_Newton_rejected=bool(need_arc),sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths});dest=ROOT/'docs/verification/static-branch-continuation-result.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
