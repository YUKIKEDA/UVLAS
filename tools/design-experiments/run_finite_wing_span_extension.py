"""One added span level and separately controlled approximations."""
from pathlib import Path
import argparse,json,time,gc
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/finite-wing-2-3'
def main():
    p=argparse.ArgumentParser();p.add_argument('side');args=p.parse_args()
    if args.side=='reference':from probe_aerosandbox_reference import run
    else:from uvlm_steady_reference import run
    for d in [0.,5.]:
        tasks=[('span384',16,384,1.114285714285714e-5,1024.,1.,1e-8)]
        if args.side=='candidate':
            tasks += [(f'length{L}',16,192,1.114285714285714e-5,L,1.,1e-8) for L in [256.,512.]]
            tasks += [(f'core{c}',16,192,c*1.114285714285714,1024.,1.,1e-8) for c in [1e-3,1e-4]]
        else:
            tasks += [('core1e-7',16,192,0,0,1.,1e-7)]
            tasks += [(f'angle{angle}',4,12,0,0,angle,1e-8) for angle in [-2.,-1.,2.]]
        for label,m,n,eps,L,angle,core in tasks:
            name=f'{args.side}-d{int(d)}-{label}';p=OUT/(name+'.json')
            if p.exists():continue
            start=time.monotonic();print('START',name,flush=True)
            if args.side=='candidate':
                a=run(m,n,d,eps,L);r=dict(F=a['F'].tolist(),M=a['M'].tolist(),bands=a['bands'].tolist(),residual=float(a['residual']),eps=eps,length=L);np.savez_compressed(p.with_suffix('.npz'),gamma=a['gamma'],horse=a['horse'],left=a['g']['left'],right=a['g']['right'],collocation=a['g']['P'],F=a['F'],M=a['M'],bands=a['bands'])
            else:
                a,r0,U,elapsed=run(m,n,d,angle,core);horse=a.vortex_strengths/np.sin(np.deg2rad(angle));F=1.2*horse[:,None]*np.cross([20.,0,0],a.right_vortex_vertices-a.left_vortex_vertices);M=np.cross((a.left_vortex_vertices+a.right_vortex_vertices)/2-[.375,0,0],F);bands=[]
                for lo in np.linspace(-4,3.6,20):
                    fraction=np.maximum(0,np.minimum(a.right_vortex_vertices[:,1],lo+.4)-np.maximum(a.left_vortex_vertices[:,1],lo))/(a.right_vortex_vertices[:,1]-a.left_vortex_vertices[:,1]);bands.append(float(fraction@F[:,2]))
                r=dict(F=F.sum(axis=0).tolist(),M=M.sum(axis=0).tolist(),bands=bands,core=core,alpha=angle,raw_F=np.asarray(r0['F_g']).tolist(),raw_M=np.asarray(r0['M_g']).tolist());np.savez_compressed(p.with_suffix('.npz'),horse=horse,left=a.left_vortex_vertices,right=a.right_vortex_vertices,collocation=a.collocation_points,F=F.sum(axis=0),M=M.sum(axis=0),bands=bands,raw_F=a.forces_geometry,raw_M=a.moments_geometry)
            r.update(side=args.side,dihedral=d,m=m,n=n,seconds=time.monotonic()-start);p.write_text(json.dumps(r,indent=2),encoding='utf-8');print('DONE',name,flush=True);del a;gc.collect()
if __name__=='__main__':main()
