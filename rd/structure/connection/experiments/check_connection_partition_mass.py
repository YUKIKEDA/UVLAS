"""No mass duplication under actual PG thickness/footprint partitioning."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import numpy as np
import json,hashlib,time
import connection_shell as pg
from check_connection_fold import Q
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    start=time.monotonic();dest=ROOT/'rd/structure/connection/records/connection-partition-mass-result.json';assert not dest.exists();rows=[]
    for h in [.1,.05]:
        for layers in [1,2,4]:
            boxes=[]
            for part in [0,1,2]:
                for j in range(layers):
                    if part<2:
                        length=1-h if part==0 else h;o=np.array([0 if part==0 else 1-h,0.,-h/2+(j+.5)*h/layers]);rot=np.eye(3)
                    else:length=.5;o=np.array([1-h+(j+.5)*h/layers,0.,h/2]);rot=Q
                    boxes.append((length,o,rot))
            v=np.array([.3,-.2,.7]);w=np.array([.4,.8,-.1]);energy=0.;volume=0.;minimum=1e300
            for length,o,rot in boxes:
                m=pg.flat_model(1,h/layers,length,.2);d=pg.data(m,origin=o,Q=rot,order=10);velocity=d['translation'][:,None]*v+np.cross(w,d['q0']);energy+=float(np.sum(velocity*(d['M']@velocity)));volume+=float(sum(d['weights']));minimum=min(minimum,float(np.linalg.eigvalsh(d['M']).min()))
            reference=0.;totalmass=0.
            # Independent moments of the TWO original non-overlapping boxes.
            for dims,o,rot in [(np.array([1.,.2,h]),np.zeros(3),np.eye(3)),(np.array([.5,.2,h]),np.array([1-h/2,0.,h/2]),Q)]:
                center=o+rot@np.array([dims[0]/2,dims[1]/2,0.]);mass=2*np.prod(dims);cov=rot@np.diag(dims*dims/12)@rot.T;J=mass*(np.trace(cov)*np.eye(3)-cov);reference+=mass*np.dot(v+np.cross(w,center),v+np.cross(w,center))+w@J@w;totalmass+=mass
            error=abs(energy/reference-1);volume_error=abs(volume/(1.5*.2*h)-1);assert error<1e-9 and volume_error<1e-10 and minimum>0
            r=dict(h=h,layers=layers,physical_mass=totalmass,child_patch_count=len(boxes),rigid_energy_error=error,volume_error=volume_error,min_child_mass_eigenvalue=minimum);rows.append(r);print('PARTITION_MASS',r,flush=True)
    paths=[Path(__file__),Path(pg.__file__),ROOT/'rd/structure/connection/experiments/check_connection_fold.py'];dest.write_text(json.dumps(dict(passed=True,records=rows,sources={p.relative_to(ROOT).as_posix():sha(p) for p in paths},seconds=time.monotonic()-start),indent=2),encoding='utf-8');print('CONNECTION_PARTITION_MASS_RESULT',True)
if __name__=='__main__':main()
