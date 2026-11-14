"""Decompose the saved 9^3 operator; verify rather than redefine it at 11^3."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import numpy as np
import mq3_interior_stationarity as original

method = original.method


def assemble(a, cache, order, splits, split_paths, out):
    groups = {name: [] for name in ['m', 's', 'x']}
    index = {int(node): i for i, node in enumerate(a['free'])}
    for (_, _, ids), key in zip(a['cells'], a['map_keys']):
        if key not in splits:
            state = cache[key]['state']
            fixed = original.split(state, 9)
            K = fixed['K']
            d9, d11 = original.data(state, 9), original.data(state, 11)
            cw = d11['cw']@K
            sw = np.sqrt(d11['weight'])
            M = (cw[:, :2, :2]*sw[:, None, None, None]).reshape(-1, len(K))
            S = (cw[:, :2, 2]*np.sqrt(2.)*sw[:, None, None]).reshape(-1, len(K))
            PM11, rS11, _ = original.null(S)
            PS11, rM11, _ = original.null(M)
            c9 = d9['cw']@K
            w9 = np.sqrt(d9['weight'])
            M9 = (c9[:, :2, :2]*w9[:, None, None, None]).reshape(-1, len(K))
            S9 = (c9[:, :2, 2]*np.sqrt(2.)*w9[:, None, None]).reshape(-1, len(K))
            PM9, _, _ = original.null(S9)
            PS9, _, _ = original.null(M9)
            checks = dict(fixed['checks'])
            metrics = dict(
                H_change=float(np.linalg.norm(K.T@(d11['H']-d9['H'])@K)),
                D_change=float(np.linalg.norm((d11['D']-d9['D'])@K)/max(1., np.linalg.norm(d9['D']@K))),
                l_change=float(np.linalg.norm(K.T@(d11['l']-d9['l']))/max(1., np.linalg.norm(K.T@d9['l']))),
                membrane_target_projector=float(np.linalg.norm(PM9@PM9.T-PM11@PM11.T)),
                shear_target_projector=float(np.linalg.norm(PS9@PS9.T-PS11@PS11.T)),
                membrane_purity=float(np.linalg.norm(S@fixed['Zm'])/max(1., np.linalg.norm(S))),
                shear_purity=float(np.linalg.norm(M@fixed['Zs'])/max(1., np.linalg.norm(M))))
            checks['independent_quadrature'] = max(metrics.values())<=1e-10
            checks['independent_field_ranks'] = rS11==fixed['metrics']['ranks']['shear_field'] and rM11==fixed['metrics']['ranks']['membrane_field'] and rS11==[rS11[1]]*3 and rM11==[rM11[1]]*3
            path = out/(key+'.npz')
            saved = {name:fixed[name] for name in ['K', 'Z', 'Zm', 'Zs', 'Zx']}
            # Internal rows are frozen; the suffix records evaluation order only.
            for tag in ['9', '11']:
                saved.update({name+tag:fixed[name] for name in ['Z', 'Zm', 'Zs', 'Zx']})
                saved.update({name+tag:value for name, value in fixed['rows'].items()})
            saved.update({'s_'+name:value for name, value in fixed['spectra'].items()})
            np.savez_compressed(path, **saved)
            passed = all(checks.values())
            split_paths[key] = dict(artifact=path.relative_to(method.ROOT).as_posix(),
                artifact_sha256=method.digest(path), fixed_metrics=fixed['metrics'],
                independent_quadrature_metrics=metrics, checks={k:bool(v) for k,v in checks.items()},
                operator_frozen=True, passed=bool(passed))
            splits[key] = fixed
            if not passed:
                raise ValueError('fixed_operator_decomposition_failed:'+key)
        for name in groups:
            block = splits[key]['rows'][name]
            global_block = np.zeros((len(block), 6*len(a['free'])))
            for local, node in enumerate(ids):
                if int(node) in index:
                    j = index[int(node)]
                    global_block[:, 6*j:6*j+6] = block[:, 9*local:9*local+6]
            groups[name].append(global_block)
    return {name:np.vstack(value) for name,value in groups.items()}
