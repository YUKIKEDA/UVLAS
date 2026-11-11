"""Nested necessary-kernel relaxations with fixed saved M-Q3 first maps."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import time
import numpy as np
import mq3_necessary_constraints as necessary
import mq3_factor_space as method
import check_q9_plate_convergence as memory

ROOT = rd_paths.ROOT


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_parent():
    path = ROOT/'rd/structure/main-shell/mq3/records/mq3-factor-v2-result.json'
    parent = json.loads(path.read_text())
    assert parent['all_checks_pass']
    cache = {}
    for item in parent['maps'].values():
        saved_path = ROOT/item['artifact']
        assert digest(saved_path) == item['artifact_sha256']
        with np.load(saved_path) as m:
            state = {name: m[name] for name in ['q', 'T', 'rows', 'U']}
            cache[next(key for key, val in parent['maps'].items() if val is item)] = dict(
                state=state, X0=m['X0'], Xd=m['Xd'], validation=item['validation'],
                existing=True, metrics=item['metrics'])
    return parent, cache


def main():
    dest = ROOT/'rd/structure/main-shell/mq3/records/mq3-necessary-result.json'
    if dest.exists():
        raise FileExistsError(dest)
    start = time.monotonic()
    parent, cache = load_parent()
    out = ROOT/'build/analysis'/('mq3-necessary-'+time.strftime('%Y%m%d-%H%M%S'))
    out.mkdir()
    (out/'maps').mkdir()
    records, traces, paths, stop = [], [], {}, None
    try:
        for amp in [0., .2]:
            for work in [False, True]:
                method.family.configure('covariant', work)
                ref = next(r for r in parent['records'] if
                    (r['metric'], r['work'], r['kappa'], r['amplitude'], r['n']) ==
                    ('covariant', work, .4, amp, 2))
                assert digest(ROOT/ref['artifact']) == ref['artifact_sha256']
                previous = np.load(ROOT/ref['artifact'])
                a = method.assemble(.4, 2, amp, 9, cache, paths, out/'maps')
                b = method.assemble(.4, 2, amp, 11, cache, paths, out/'maps')
                groups9, groups11 = necessary.assemble(a, 9), necessary.assemble(b, 11)
                scale = previous['common_scale']
                last_error = 0.
                levels = ['N', 'NM', 'NMS'] + (['NMSW'] if work else [])
                for level in levels:
                    tick = time.monotonic()
                    H9 = np.vstack([groups9[c] for c in level])
                    H11 = np.vstack([groups11[c] for c in level])
                    N, info9, s9 = necessary.kernel(H9, scale, a)
                    No, info11, s11 = necessary.kernel(H11, scale, b)
                    V, Vo = N/scale[:, None], No/scale[:, None]
                    projection = float(np.linalg.norm(V@V.T-Vo@Vo.T))
                    inclusion = necessary.relative_constraint(H11, previous['kernel'])
                    checks = dict(rank=info9['ranks']==info11['ranks']==[info9['ranks'][1]]*3,
                        kernel=max(info9['residual'], info11['residual'])<=1e-10,
                        anchor=max(info9['anchor_shift'], info11['anchor_shift'])<=1e-7 and
                            max(info9['orthogonality'], info11['orthogonality'])<=1e-10,
                        projector=projection<=1e-7, inclusion=inclusion<=1e-9)
                    solutions = {}
                    saved = dict(xy=a['xy'], coords=a['coords'], ids=previous['ids'],
                        free=a['free'], kernel=N, common_scale=scale, H9=H9, H11=H11,
                        singular9=s9, singular11=s11)
                    if all(checks.values()):
                        T = np.column_stack((a['A1']@N, a['D1']))
                        for kind in ['quadratic', 'quartic']:
                            rh = method.common.reference_h1(kind, .4, 40)
                            q, fit = method.solve_semidefinite_ball(T, a['targets'][kind], N, a['M'], 2*rh)
                            u = method.common.full_field(a, N, q)
                            rd = u[a['free'], :6].reshape(-1)
                            c = u[a['free'], 6:].reshape(-1)
                            error = float(np.linalg.norm(b['A1']@rd+b['D1']@c-b['targets'][kind])/
                                np.linalg.norm(b['targets'][kind]))
                            h1 = method.common.h1tool.h1(u, b['free'], b['M'])
                            residual = necessary.relative_constraint(H11, rd)
                            sc = dict(kkt=fit['kkt']<=1e-8, sphere=h1<=2*rh*(1+1e-8),
                                quadrature=abs(error-fit['error'])<=1e-8,
                                H1=abs(h1/fit['h1']-1)<=1e-8, necessary=residual<=1e-9,
                                positive=kind!='quadratic' or error<=1e-8)
                            if kind == 'quartic':
                                sc['relaxation'] = error<=ref['solutions'][kind]['observed_error']+1e-8
                                sc['monotone'] = error+1e-8>=last_error
                                last_error = error
                            checks[kind] = all(sc.values())
                            solutions[kind] = dict(fit=fit, observed_error=error, observed_h1=h1,
                                reference_h1=rh, necessary_residual=residual,
                                full_zero_residual=necessary.relative_constraint(b['A0'], rd),
                                checks={key: bool(val) for key, val in sc.items()})
                            saved['u_'+kind] = u
                    artifact = out/f'a{amp:g}-w{int(work)}-{level}.npz'
                    np.savez_compressed(artifact, **saved)
                    checks['resources'] = time.monotonic()-tick<=180 and time.monotonic()-start<=900 and memory.memory_bytes()<=2*1024**3
                    record = dict(kappa=.4, amplitude=amp, n=2, metric='covariant', work=work,
                        level=level, nullity=N.shape[1], rank9=info9, rank11=info11,
                        projector_difference=projection, inclusion_residual=inclusion,
                        parent_error=ref['solutions']['quartic']['observed_error'],
                        solutions=solutions, map_keys=a['map_keys'],
                        checks={key: bool(val) for key, val in checks.items()},
                        artifact=artifact.relative_to(ROOT).as_posix(), artifact_sha256=digest(artifact),
                        seconds=time.monotonic()-tick)
                    records.append(record)
                    print(amp, work, level, 'null', N.shape[1],
                        'errors', {k: v['observed_error'] for k,v in solutions.items()}, checks, flush=True)
                    if not all(checks.values()):
                        raise RuntimeError('required_necessary_gate_failed')
                if not work:
                    u = previous['u_quartic']
                    cell_traces = []
                    for _, _, ids in a['cells']:
                        raw = necessary.local(a['coords'][ids].reshape(144), 11)
                        cell_traces.append((raw['W']@u[ids].reshape(144)).tolist())
                    traces.append(dict(amplitude=amp, source_artifact=ref['artifact'],
                        mean_raw_zero_mandel=cell_traces,
                        global_W_action=(groups11['W']@u[a['free'], :6].reshape(-1)).tolist()))
    except Exception as exc:
        stop = repr(exc)
    sources = [Path(__file__), Path(necessary.__file__), Path(method.__file__),
        Path(method.family.__file__), ROOT/'rd/structure/main-shell/mq3/mq3-necessary-constraint-diagnostic.md',
        ROOT/'rd/structure/main-shell/mq3/records/mq3-necessary-preparation.json',
        ROOT/'rd/structure/main-shell/mq3/records/mq3-factor-v2-result.json']
    dest.write_text(json.dumps(dict(records=records, maps=paths, work_traces=traces,
        all_checks_pass=len(records)==14 and stop is None, stop_reason=stop,
        physical_adoption=False, finite_thickness_lower_bound=False,
        sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources},
        seconds=time.monotonic()-start, peak_bytes=memory.memory_bytes()), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
