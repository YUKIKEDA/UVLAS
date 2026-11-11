"""Reference-only physical expansion, with raw/derived provenance kept separate."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,shutil
import numpy as np

R=np.array([1.,-1.,1.])
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expand(raw,m,physical_n,out):
    out.mkdir(parents=True,exist_ok=False)
    rows=np.loadtxt(raw/'uv-rows.txt',skiprows=2);forces=np.loadtxt(raw/'uv-forces.txt',skiprows=2)
    half=physical_n//2;size=m*half
    assert rows.shape==(size,21) and forces.shape==(3*size,9)
    order=np.array([(half-1-s)*m+j for s in range(half) for j in range(m)])
    negative=rows[order].copy()
    for start in [1,4,14,17]:negative[:,start:start+3]*=R
    negative[:,7:10]=rows[order,10:13]*R;negative[:,10:13]=rows[order,7:10]*R
    expanded=np.vstack([negative,rows]);expanded[:,0]=np.arange(1,2*size+1)
    pos=forces.copy();pos[:,0]+=size;pos[:,1]+=half
    neg=forces.copy();source=forces[:,0].astype(int)-1
    neg[:,0]=(half-1-source//m)*m+source%m+1
    neg[:,1]=half-forces[:,1]+1
    neg[:,2]=np.where(forces[:,2]==0,0,3-forces[:,2])
    neg[:,3:6]*=R;neg[:,6:9]*=R
    ff=np.vstack([neg,pos]);ff=ff[np.lexsort((ff[:,2],ff[:,0]))]
    np.savetxt(out/'uv-rows.txt',expanded,fmt=['%d']+['%.17e']*20,header='UVLAS_AVL_ROWS_SYMMETRY_EXPANDED_V1\ni rc n rv1 rv2 gamma_over_U known induced res',comments='')
    np.savetxt(out/'uv-forces.txt',ff,fmt=['%d']*3+['%.17e']*6,header='UVLAS_AVL_FORCE_SYMMETRY_EXPANDED_V1\ni strip kind position[3] force_over_q[3]',comments='')
    strips=np.array([[s+1,expanded[s*m,8],expanded[s*m,11],m,s*m+1] for s in range(physical_n)])
    np.savetxt(out/'uv-strips.txt',strips,fmt=['%d','%.17e','%.17e','%d','%d'],header='UVLAS_AVL_STRIPS_SYMMETRY_EXPANDED_V1\nstrip y_left y_right count first_vortex',comments='')
    shutil.copyfile(raw/'total.txt',out/'total.txt')
    info=dict(kind='derived-by-symmetry',raw_kind='native-half',native_unknowns=size,physical_panels=2*size,
        reflection=[1,-1,1],swapped_endpoints=True,swapped_legs=True,
        total_file_role='unchanged native half-mode FT, already containing physical full-wing totals',
        parent_sha256={name:digest(raw/name) for name in ['uv-rows.txt','uv-forces.txt','uv-strips.txt','total.txt','uh-meta.txt']})
    (out/'derivation.json').write_text(json.dumps(info,indent=2)+'\n')
    return expanded,ff,info

def canonical_forces(f):return f[np.lexsort((f[:,2],f[:,0]))]

def synthetic():
    p=np.array([.3,.7,.2]);f=np.array([1.,2.,3.]);o=np.array([.375,0.,0.])
    moment=np.cross(p-o,f);reflected=np.cross(p*R-o,f*R)
    return dict(axial_moment_error=float(np.max(abs(reflected+moment*R))),
        paired_force_error=float(np.max(abs(f+f*R-[2*f[0],0,2*f[2]]))),
        paired_moment_error=float(np.max(abs(moment+reflected-[0,2*moment[1],0]))))
