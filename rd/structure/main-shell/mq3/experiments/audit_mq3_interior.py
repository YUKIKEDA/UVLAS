"""Separate raw-constraint construction and direct-complex saved-field evaluation."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import itertools
import hashlib
import json
import time
import numpy as np
from scipy.linalg import solve
import mq3_factor_family as family

ROOT = rd_paths.ROOT


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def raw(q, plane):
    N, d = family.base.geometry.shape(plane)
    grad = np.zeros((len(plane), 48, 3))
    grad[:, ::3, :2] = d
    grad[:, 1::3, 2] = N
    J = np.einsum('na,pnb->pab', q.reshape(48, 3), grad)
    gradient = grad @ np.linalg.inv(J)
    return J, family.base.strain_matrix(gradient)


def cov_component(J, B, i, j):
    return np.einsum('pa,kab,pb,pkn->pn', J[:, :, i], family.base.MB, J[:, :, j], B)


def local(q):
    p, w = np.polynomial.legendre.leggauss(11)
    plane = np.array(list(itertools.product(p, p)))
    weights = np.array([a*b for a, b in itertools.product(w, w)])
    J, B = raw(q, plane)
    volume = weights*np.linalg.det(J)
    volume /= volume.sum()
    out = dict(N=cov_component(J, B, 2, 2)*np.sqrt(weights)[:, None],
        W=np.einsum('p,pkn->kn', volume, B))
    membrane, shear = [], []
    for axis in range(2):
        for sign in [-1., 1.]:
            edge = np.zeros((11, 2))
            edge[:, axis] = p
            edge[:, 1-axis] = sign
            Je, Be = raw(q, edge)
            m = cov_component(Je, Be, axis, axis)
            s = cov_component(Je, Be, axis, 2)
            for i in range(3):
                membrane.append((w*p**i) @ m)
            for i in range(3):
                # Analytic thickness integral at h=0; no side-face raw evaluation.
                for j in range(3):
                    shear.append(((w*p**i) @ s) * (0. if j%2 else 2./(j+1)))
    out.update(M=np.array(membrane), S=np.array(shear))
    return out


def constraints(saved, level):
    index = {int(node): i for i, node in enumerate(saved['free'])}
    groups = {name: [] for name in level}
    for ids in saved['ids']:
        blocks = local(saved['coords'][ids].reshape(144))
        for name in level:
            block = blocks[name]
            H = np.zeros((len(block), 6*len(index)))
            for a, node in enumerate(ids):
                if int(node) in index:
                    j = index[int(node)]
                    H[:, 6*j:6*j+6] = block[:, 9*a:9*a+6]
            groups[name].append(H)
    return np.vstack([np.vstack(groups[name]) for name in level])


def energy(E):
    return family.base.LAM*np.trace(E, axis1=1, axis2=2)**2 + 2*family.base.MU*np.einsum('pij,pij->p', E, E)



def interior_constraints(saved, record, data):
    Hbase = constraints(saved, record['base_level'])
    index = {int(node): i for i,node in enumerate(saved['free'])}
    groups = {name: [] for name in record['level']}
    points, w = family.base.quadrature(11)
    for ids, key in zip(saved['ids'], record['map_keys']):
        item = data['maps'][key]
        split_item = data['splits'][key]
        assert digest(ROOT/item['artifact']) == item['artifact_sha256']
        assert digest(ROOT/split_item['artifact']) == split_item['artifact_sha256']
        m = np.load(ROOT/item['artifact'])
        z = np.load(ROOT/split_item['artifact'])
        q, T = m['q'], m['T']
        target = family.data(q, 0., points)
        J, B = raw(q, points[:,:2])
        center, _ = raw(q, np.zeros((1,2)))
        length = np.linalg.norm(center[0],axis=0);length[2]=1.
        Jh = J/length[None,None,:]
        cw = family.base.covariant(Jh, target['W']@T@z['K'])
        cr = family.base.covariant(Jh, B)
        l = np.einsum('p,pabi,pabn->in',w,cw,cr,optimize=True)
        for name in groups:
            block = z['Z'+name+'11'].T@l
            global_block = np.zeros((len(block),6*len(index)))
            for local,node in enumerate(ids):
                if int(node) in index:
                    j=index[int(node)]
                    global_block[:,6*j:6*j+6]=block[:,9*local:9*local+6]
            groups[name].append(global_block)
    return np.vstack([Hbase]+[np.vstack(groups[name]) for name in record['level']])

def main():
    dest = ROOT/'rd/structure/main-shell/mq3/records/mq3-interior-audit.json'
    if dest.exists():
        raise FileExistsError(dest)
    start = time.monotonic()
    parent_path = ROOT/'rd/structure/main-shell/mq3/records/mq3-interior-result.json'
    data = json.loads(parent_path.read_text())
    for path, sha in data['sources'].items():
        assert digest(ROOT/path) == sha
    p, w = np.polynomial.legendre.leggauss(11)
    plane = np.array(list(itertools.product(p, p)))
    weights = np.array([a*b for a, b in itertools.product(w, w)])
    points = np.column_stack((plane, np.ones(len(plane))))
    N, g = family.base.geometry.shape(plane)
    records, traces, cache, stop = [], [], {}, None
    try:
        for record in data['records']:
            if not all(record['checks'].values()):
                continue
            family.configure('covariant', record['work'])
            assert digest(ROOT/record['artifact']) == record['artifact_sha256']
            saved = np.load(ROOT/record['artifact'])
            H = interior_constraints(saved, record, data)
            hd = float(np.linalg.norm(H-saved['H11'])/max(1., np.linalg.norm(H)))
            kres = float(np.linalg.norm(H@saved['kernel'])/(np.linalg.norm(H)*np.linalg.norm(saved['kernel'])))
            totals = {kind: dict(error=0., reference=0., h1=0.) for kind in ['quadratic', 'quartic']}
            area = 0.
            for ids, key in zip(saved['ids'], record['map_keys']):
                if key not in cache:
                    item = data['maps'][key]
                    assert digest(ROOT/item['artifact']) == item['artifact_sha256']
                    m = np.load(ROOT/item['artifact'])
                    state = {name: m[name] for name in ['q', 'T', 'rows', 'U']}
                    Sc, bc, _, _ = family.system(state, 1e-25j)
                    Xc = solve(Sc, bc, assume_a='gen', check_finite=False)
                    B1 = 2*family.evaluate(state, 1e-25j, Xc, points).imag[:, :, :144]/1e-25
                    J, _ = raw(state['q'], plane)
                    cache[key] = B1, J
                B1, J = cache[key]
                weight = weights*np.linalg.det(J)
                x = (N@saved['coords'][ids, 0])[:, 0]
                vec = np.column_stack((np.ones_like(x), np.zeros_like(x), .4*x))
                D = 1+.16*x*x
                Js = np.einsum('na,pnb->pab', saved['coords'][ids, 0], g)
                G = Js.transpose(0, 2, 1) @ Js
                aw = weights*np.sqrt(np.linalg.det(G))
                area += aw.sum()
                for kind, total in totals.items():
                    u = saved['u_'+kind][ids]
                    curvature = np.ones_like(x) if kind=='quadratic' else x*x
                    target = -curvature[:, None, None]*vec[:, :, None]*vec[:, None, :]/D[:, None, None]
                    E = np.einsum('kij,pk->pij', family.base.MB, B1@u.reshape(144))
                    total['error'] += weight@energy(E-target)
                    total['reference'] += weight@energy(target)
                    val = N@u[:, 3:9]
                    du = np.einsum('na,pnb->pab', u[:, 3:9], g)
                    total['h1'] += aw@(np.sum(val*val, axis=1)+np.einsum('pai,pij,paj->p', du, np.linalg.inv(G), du))
            checks = dict(constraint_matrix=hd<=1e-10, kernel=kres<=1e-9)
            observations = {}
            for kind, total in totals.items():
                error = float(np.sqrt(total['error']/total['reference']))
                h1 = float(np.sqrt(total['h1']/area))
                ed = abs(error-record['solutions'][kind]['observed_error'])
                hdiff = abs(h1/record['solutions'][kind]['observed_h1']-1)
                rd = saved['u_'+kind][saved['free'], :6].reshape(-1)
                residual = float(np.linalg.norm(H@rd)/(np.linalg.norm(H)*np.linalg.norm(rd)))
                checks[kind] = ed<=1e-8 and hdiff<=1e-8 and residual<=1e-9
                observations[kind] = dict(error=error, error_difference=ed,
                    h1_difference=hdiff, necessary_residual=residual)
            records.append(dict(amplitude=record['amplitude'], work=record['work'],
                level=record['level'], constraint_matrix_difference=hd,
                kernel_residual=kres, observations=observations,
                checks={k: bool(v) for k, v in checks.items()}))
            print(record['amplitude'], record['work'], record['level'], checks, flush=True)
            if not all(checks.values()):
                raise RuntimeError('necessary_audit_gate_failed')
    except Exception as exc:
        stop = repr(exc)
    complete = sum(all(r['checks'].values()) for r in data['records'])
    sources = [Path(__file__), Path(family.__file__), parent_path,
        ROOT/'rd/structure/main-shell/mq3/records/mq3-interior-audit-protocol.md']
    dest.write_text(json.dumps(dict(records=records,
        validated_primary_problems=complete, primary_complete=data['all_checks_pass'],
        all_checks_pass=len(records)==complete and stop is None,
        stop_reason=stop, independent_global_resolve=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},
        seconds=time.monotonic()-start), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
