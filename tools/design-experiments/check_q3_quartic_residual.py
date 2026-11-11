"""Shared Q3 residual split with a second tensor/thickness-QR calculation."""
from pathlib import Path
import hashlib,json,time,itertools
import numpy as np
import q3_curved_element as element
import developable_general_reference as reference
import check_q9_eas21_residual as partition
import audit_q9_eas21_residual_v2 as tensor_partition
ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    dest=ROOT/'docs/verification/q3-quartic-residual-result.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.monotonic();parent=ROOT/'docs/verification/q3-quartic-shared-approximation-result.json';data=json.loads(parent.read_text());assert data['all_numerical_checks_pass'];assert np.max(abs(partition.method.MB-element.MB))<1e-15 and np.max(abs(partition.method.C-element.C))<1e-15;records=[];points,w=element.quadrature(11)
    for record in data['records']:
        path=ROOT/record['artifact'];assert digest(path)==record['artifact_sha256'];saved=np.load(path);parts=np.zeros(5);audit_parts=np.zeros(5);mean=total=refnorm=memmean=0.;maxrecon=0.
        for ids in saved['ids']:
            q=saved['coords'][ids].reshape(144);_,J,B=element.evaluate(q,record['h'],points);weight=w*np.linalg.det(J);xy,_=element.bilinear(points,saved['xy'][ids][[0,3,15,12]]);x=xy[:,0];z=record['h']*points[:,2]/2;star=np.einsum('kij,pij->pk',element.MB,reference.strain('quartic',.4,x,z));error=B@saved['u'][ids].reshape(144)-star;log=partition.split(error,weight,x,record['h'],points,11,.4);tensor=np.einsum('kij,pk->pij',element.MB,error);other,_,_,_=tensor_partition.local(tensor,weight,x,z,11);parts+=log['parts'];audit_parts+=other;mean+=log['mean'];total+=log['total'];memmean+=log['membrane_mean'];refnorm+=weight@partition.norm2(star);maxrecon=max(maxrecon,log['reconstruction'])
        err=float(np.sqrt(total/refnorm));auditdiff=float(np.max(abs(parts-audit_parts))/refnorm);checks=dict(reconstruction=maxrecon<=1e-10,partition_audit=auditdiff<=1e-8,parent=abs(err-record['observation']['error'])<=1e-8);assert all(checks.values());records.append(dict(n=record['n'],h=record['h'],error=err,normalized_parts={label:float(v/refnorm) for label,v in zip(partition.LABELS,parts)},fractions={label:float(v/total) for label,v in zip(partition.LABELS,parts)},element_mean_fraction=float(mean/total),membrane_internal_variation_fraction=float(1-memmean/parts[0]),reconstruction=maxrecon,audit_difference=auditdiff,checks=checks));print(record['n'],record['h'],records[-1]['fractions'],flush=True)
    sources=[Path(__file__),Path(element.__file__),Path(reference.__file__),Path(partition.__file__),Path(tensor_partition.__file__),parent,ROOT/'docs/verification/q3-quartic-residual-experiment.md']
    dest.write_text(json.dumps(dict(records=records,all_checks_pass=True,physical_adoption=False,sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},seconds=time.monotonic()-start),indent=2),encoding='utf-8')
if __name__=='__main__':main()
