"""Export explicit displacement/constant input axes without changing witnesses."""
from pathlib import Path
import hashlib,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    dest=ROOT/'docs/verification/mq3-zero-first-map-export.json'
    if dest.exists():raise FileExistsError(dest)
    parent=ROOT/'docs/verification/mq3-zero-first-map-result.json';data=json.loads(parent.read_text());assert data['all_checks_pass'];out=ROOT/'build/analysis'/('mq3-zero-first-export-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir();records=[]
    for record in data['records']:
        path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];saved=np.load(path);B0=saved['B0'];B1=saved['B1'];assert B0.shape==B1.shape==(11,11,6,150);b0=B0[...,:144].reshape(121,6,144);b1=B1[...,:144].reshape(121,6,144);c0=B0[...,144:].reshape(121,6,6);c1=B1[...,144:].reshape(121,6,6);assert np.max(abs(c0-np.eye(6)))<=1e-9 and np.max(abs(c1))<=1e-9;file=out/Path(record['artifact']).name;np.savez_compressed(file,q=saved['q'],T=saved['T'],rows=saved['rows'],U=saved['U'],X0=saved['X0'],Xd=saved['Xd'],B0_displacement=b0,B1_displacement=b1,B0_constant=c0,B1_constant=c1,plane_points=np.array([(s,t) for s in np.polynomial.legendre.leggauss(11)[0] for t in np.polynomial.legendre.leggauss(11)[0]]));records.append(dict(kappa=record['kappa'],n=record['n'],cell=record['cell'],original_sha256=record['artifact_sha256'],artifact=file.relative_to(ROOT).as_posix(),artifact_sha256=digest(file)))
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=len(records)==10,displacement_shape=[121,6,144],constant_shape=[121,6,6],meaning='11x11 observation witnesses, not interpolation tables for arbitrary points',sources={Path(__file__).relative_to(ROOT).as_posix():digest(__file__),parent.relative_to(ROOT).as_posix():digest(parent)}),indent=2),encoding='utf-8');print('Exported 10 explicit 144+6 input maps; original witnesses unchanged.')
if __name__=='__main__':main()
