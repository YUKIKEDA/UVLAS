"""M-Q3I1: hard work/shear conditions and low-degree membrane moment fit."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
import itertools
import numpy as np
from scipy.linalg import qr, solve, solve_triangular, svd, svdvals
import mq3_factor_family as family
import mq3_interior_stationarity as interior


class PreparationFailure(ValueError):
    def __init__(self, report, state, A, b, C):
        super().__init__('coefficient_space_not_unisolvent:'+str(report))
        self.report = report
        self.witness = dict(q=state['q'], T=state['T'], Qs=state['Qs'], A=A, b=b, C=C,
            geometry_projection=state['geometry_projection'], geometry_vectors=state['geometry_vectors'])


def ranks(s, scale=None):
    if scale is None:
        scale = s[0]
    return [int(np.count_nonzero(s > t*scale)) for t in [1e-9, 1e-10, 1e-11]]


def plane(order):
    p, w = np.polynomial.legendre.leggauss(order)
    return np.array(list(itertools.product(p, p))), np.array([a*b for a,b in itertools.product(w, w)])


def polynomial(points, full=False):
    s, t = points[:, :2].T
    ls = np.polynomial.legendre.legvander(s, 3)*np.sqrt((2*np.arange(4)+1)/2)
    lt = np.polynomial.legendre.legvander(t, 3)*np.sqrt((2*np.arange(4)+1)/2)
    shapes = [(0,2,3), (1,3,2), (2,2,2)] if full else [(0,2,1), (1,1,2), (2,2,2)]
    rows = []
    for component, imax, jmax in shapes:
        for i,j in itertools.product(range(imax+1), range(jmax+1)):
            v = np.zeros((len(points), 3))
            v[:, component] = ls[:, i]*lt[:, j]
            rows.append(v)
    return np.stack(rows, axis=2)


def normalized_membrane(value, length):
    return value*np.array([1/length[0]**2, 1/length[1]**2,
        np.sqrt(2.)/(length[0]*length[1])])[None, :, None]


def membrane_data(state, points):
    pts = np.column_stack((points[:, :2], np.zeros(len(points))))
    d = family.data(state['q'], 0., pts)
    return normalized_membrane(d['M'], state['length']), normalized_membrane(d['m'], state['length'])


def tests(state, points):
    base = polynomial(points, True)
    canonical = polynomial(points)
    if not state['geometry_vectors'].shape[1]:
        return canonical
    _, raw = membrane_data(state, points)
    residual = raw[:, :, 144:]-base@state['geometry_projection']
    extra = residual@state['geometry_vectors']
    return np.concatenate((canonical, extra), axis=2)


def moment_rows(state, order=9):
    points, w = plane(order)
    target, raw = membrane_data(state, points)
    phi = tests(state, points)
    A = np.einsum('p,pci,pcj->ij', w, phi, target@state['T'], optimize=True)
    b = np.einsum('p,pci,pcn->in', w, phi, raw, optimize=True)
    return A, b


def gram(state, h, order=9):
    points, d, weight = family.volume(state['q'], h, order)
    J0, _ = family.raw(state['q'], 0., points)
    Jh = J0/state['length'][None, None, :]
    W = family.base.covariant(Jh, d['W']@state['T'])
    F = family.base.covariant(Jh, d['F'])
    _, w = family.base.quadrature(order)
    H = np.einsum('p,pabi,pabj->ij', w, W, W, optimize=True)
    l = np.einsum('p,pabi,pabn->in', w, W, F, optimize=True)
    full, rhs = family.moments(state['q'], h, d, weight)
    D = state['U']@(full[state['rows']]@state['T'])
    g = state['U']@rhs[state['rows']]
    C = np.vstack((D, state['Qs'].T@H))
    crhs = np.vstack((g, state['Qs'].T@l))
    return C, crhs, full, rhs, H, l


def system(state, h):
    C, rhs, full, g, _, _ = gram(state, h)
    C = state['hard_U']@C
    rhs = state['hard_U']@rhs
    A, b = state['A'], state['b']
    n, m = A.shape[1], len(C)
    S = np.zeros((n+m, n+m), dtype=np.result_type(h, float))
    S[:n, :n] = A.T@A
    S[:n, n:] = C.T
    S[n:, :n] = C
    return S, np.vstack((A.T@b, rhs)), full, g


def coefficients(state, h):
    S, rhs, F, g = system(state, h)
    X = solve(S, rhs, assume_a='sym', check_finite=False)
    return X, S, rhs, F, g


def choose(q):
    family.configure('covariant', True)
    state = family.choose(q)
    split = interior.split(state, 9)
    state['T'] = state['T']@split['K']
    state['Qs'] = split['Zs']
    Jc, _ = family.raw(q, 0., np.zeros((1,3)))
    length = np.linalg.norm(Jc[0], axis=0)
    length[2] = 1.
    state['length'] = length
    points, w = plane(9)
    base = polynomial(points, True)
    _, raw = membrane_data(state, points)
    const = raw[:, :, 144:]
    projection = np.einsum('p,pci,pcn->in', w, base, const)
    residual = const-base@projection
    U, s, V = svd((residual*np.sqrt(w)[:,None,None]).reshape(-1,6), full_matrices=False, check_finite=False)
    scale = svdvals((const*np.sqrt(w)[:,None,None]).reshape(-1,6), check_finite=False)[0]
    granks = ranks(s, scale)
    state['geometry_projection'] = projection
    state['geometry_vectors'] = V[:granks[1]].T/s[:granks[1]]
    A, b = moment_rows(state)
    C, _, _, _, _, _ = gram(state, 0.)
    Q, R = qr(C.T, mode='full', check_finite=False)
    c = len(C)
    cranks = ranks(svdvals(C, check_finite=False))
    free = Q[:, c:]
    AF = A@free
    aranks = ranks(svdvals(AF, check_finite=False))
    report = dict(geometry_ranks=granks, C_ranks=cranks, A_kernel_ranks=aranks,
        coefficient_dimension=C.shape[1], constraints=c, kernel_dimension=free.shape[1], moments=len(A))
    if granks!=[granks[1]]*3 or cranks!=[c]*3 or aranks!=[free.shape[1]]*3 or not all(split['checks'].values()):
        raise PreparationFailure(report,state,A,b,C)
    _, Ra = qr(AF, mode='economic', check_finite=False)
    # Whiten the objective on the constraint kernel, without a penalty.
    P = np.column_stack((Q[:, :c], free@solve_triangular(Ra, np.eye(len(Ra)), check_finite=False)))
    state['T'] = state['T']@P
    state['Qs'] = solve(P, state['Qs'], check_finite=False)
    state['hard_U'] = solve_triangular(R[:c, :c].T, np.eye(c), lower=True, check_finite=False)
    state['A'], state['b'] = moment_rows(state)
    state['preparation'] = dict(space_ranks=state['space_ranks'], base_constraint_ranks=state['constraint_ranks'],
        geometry_ranks=granks, C_ranks=cranks, A_kernel_ranks=aranks, kernel_dimension=free.shape[1],
        pure_shear_directions=split['Zs'].shape[1], moment_count=len(A),
        checks=dict(interior=all(split['checks'].values()), space=state['space_ranks']==[len(state['T'][0])]*3,
            geometry=granks==[granks[1]]*3, C=cranks==[c]*3, A_kernel=aranks==[free.shape[1]]*3))
    return state


def evaluate(state, h, X, points):
    return family.evaluate(state, h, X, points)
