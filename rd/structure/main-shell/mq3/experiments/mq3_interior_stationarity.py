"""Field-defined pure membrane/shear decomposition of constrained stationarity."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
from scipy.linalg import cholesky, solve_triangular, svd
import mq3_factor_space as method


def null(A):
    norms = np.linalg.norm(A, axis=1)
    mask = norms > 1e-13*max(float(norms.max()), 1e-300)
    B = A[mask]/norms[mask, None]
    _, s, V = svd(B, full_matrices=True, check_finite=False)
    ranks = [int(np.count_nonzero(s > t*s[0])) for t in [1e-9, 1e-10, 1e-11]]
    return V[ranks[1]:].T, ranks, s


def data(state, order):
    f = method.family
    points, w = f.base.quadrature(order)
    d = f.data(state['q'], 0., points)
    center, _ = f.raw(state['q'], 0., np.zeros((1, 3)))
    lengths = np.linalg.norm(center[0], axis=0)
    lengths[2] = 1.
    Jhat = d['J']/lengths[None, None, :]
    W = d['W']@state['T']
    cw = f.base.covariant(Jhat, W)
    cf = f.base.covariant(Jhat, d['F'])
    weight = w*np.linalg.det(d['J'])
    weight /= weight.sum()
    F, g = f.moments(state['q'], 0., d, weight)
    D = state['U']@(F[state['rows']]@state['T'])
    H = np.einsum('p,pabi,pabj->ij', w, cw, cw, optimize=True)
    l = np.einsum('p,pabi,pabn->in', w, cw, cf, optimize=True)
    return dict(H=H, l=l, D=D, cw=cw, weight=w, points=points)


def split(state, order, K=None):
    d = data(state, order)
    H = d['H']
    if K is None:
        C = cholesky(H, lower=True, check_finite=False)
        K = solve_triangular(C.T, np.eye(len(C)), check_finite=False)
    cw = d['cw']@K
    sw = np.sqrt(d['weight'])
    M = (cw[:, :2, :2]*sw[:, None, None, None]).reshape(-1, len(K))
    S = (cw[:, :2, 2]*np.sqrt(2.)*sw[:, None, None]).reshape(-1, len(K))
    D = d['D']@K
    Z, rankD, sD = null(D)
    PM, rankS, sS = null(S)
    PS, rankM, sM = null(M)
    Qm, rankDm, sDm = null(D@PM)
    Qs, rankDs, sDs = null(D@PS)
    Zm, Zs = PM@Qm, PS@Qs
    pure = np.column_stack((Zm, Zs))
    # The complement lives within ker D; pure field spaces cannot intersect.
    Qx, rankPure, sPure = null(pure.T@Z)
    Zx = Z@Qx
    allZ = np.column_stack((Zm, Zs, Zx))
    orth = float(np.linalg.norm(allZ.T@allZ-np.eye(allZ.shape[1])))
    reconstruct = float(np.linalg.norm(allZ@allZ.T-Z@Z.T))
    gram = float(np.linalg.norm(K.T@H@K-np.eye(len(K))))
    field_split = float(np.linalg.norm(M.T@M+S.T@S-K.T@H@K))
    whitened_l = K.T@d['l']
    rows = {name: basis.T@whitened_l[:, :144] for name, basis in [('m', Zm), ('s', Zs), ('x', Zx)]}
    ranks = dict(D=rankD, membrane_field=rankM, shear_field=rankS,
        D_on_membrane=rankDm, D_on_shear=rankDs, pure_in_kernel=rankPure)
    metrics = dict(dimensions=dict(total=Z.shape[1], membrane=Zm.shape[1],
        shear=Zs.shape[1], mixed=Zx.shape[1]), ranks=ranks,
        orthogonality=orth, reconstruction=reconstruct, whitened_gram=gram,
        field_gram_split=field_split, constraint_residual=float(np.linalg.norm(D@allZ)/max(1., np.linalg.norm(D))))
    checks = dict(rank=all(v==[v[1]]*3 for v in ranks.values()),
        decomposition=allZ.shape[1]==Z.shape[1] and max(orth, reconstruct, gram, field_split)<=1e-10,
        constraints=metrics['constraint_residual']<=1e-10)
    return dict(K=K, Z=Z, Zm=Zm, Zs=Zs, Zx=Zx, rows=rows,
        metrics=metrics, checks={k: bool(v) for k,v in checks.items()},
        spectra=dict(D=sD, M=sM, S=sS, Dm=sDm, Ds=sDs, pure=sPure))


def assemble(a, cache, order, splits, split_paths, out):
    groups = {name: [] for name in ['m', 's', 'x']}
    free_index = {int(node): i for i,node in enumerate(a['free'])}
    for (_, _, ids), key in zip(a['cells'], a['map_keys']):
        state = cache[key]['state']
        if key not in splits:
            first = split(state, 9)
            second = split(state, 11, first['K'])
            stable = first['metrics']['ranks']==second['metrics']['ranks'] and first['metrics']['dimensions']==second['metrics']['dimensions']
            group_differences = {}
            for name in ['m', 's', 'x']:
                z, zo = first['Z'+name], second['Z'+name]
                group_differences[name] = float(np.linalg.norm(z@z.T-zo@zo.T))
            passed = all(first['checks'].values()) and all(second['checks'].values()) and stable and max(group_differences.values())<=1e-10
            path = out/(key+'.npz')
            saved = dict(K=first['K'])
            for tag, value in [('9', first), ('11', second)]:
                saved.update({n+tag:value[n] for n in ['Z', 'Zm', 'Zs', 'Zx']})
                saved.update({n+tag:v for n,v in value['rows'].items()})
                saved.update({'s_'+n+tag:v for n,v in value['spectra'].items()})
            np.savez_compressed(path, **saved)
            split_paths[key] = dict(artifact=path.relative_to(method.ROOT).as_posix(),
                artifact_sha256=method.digest(path), metrics9=first['metrics'], metrics11=second['metrics'],
                component_projector_differences=group_differences, passed=bool(passed))
            splits[key] = {9:first, 11:second}
            if not passed:
                raise ValueError('local_interior_decomposition_failed:'+key)
        for name in groups:
            block = splits[key][order]['rows'][name]
            global_block = np.zeros((len(block), 6*len(a['free'])))
            for local, node in enumerate(ids):
                if int(node) in free_index:
                    j = free_index[int(node)]
                    global_block[:, 6*j:6*j+6] = block[:, 9*local:9*local+6]
            groups[name].append(global_block)
    return {name:np.vstack(value) for name,value in groups.items()}
