"""Necessary zero-strain constraints from raw fields, without KKT row choices."""
import itertools
import numpy as np
from scipy.linalg import svd
import mq3_factor_space as method
from check_mq3_factor_v2 import anchor


def local(q, order):
    family = method.family
    p, w = np.polynomial.legendre.leggauss(order)
    plane = np.array(list(itertools.product(p, p)))
    weights = np.array([a*b for a, b in itertools.product(w, w)])
    pts = np.column_stack((plane, np.zeros(len(plane))))
    J, E = family.raw(q, 0., pts)
    cov = family.base.covariant(J, E[:, :, :144])
    normal = cov[:, 2, 2] * np.sqrt(weights)[:, None]
    volume_weights = weights * np.linalg.det(J)
    volume_weights /= volume_weights.sum()
    work = np.einsum('p,pkn->kn', volume_weights, E[:, :, :144])
    membrane, shear = [], []
    for tangent in range(2):
        for sign in [-1., 1.]:
            edge = np.zeros((order, 3))
            edge[:, tangent] = p
            edge[:, 1-tangent] = sign
            Je, Ee = family.raw(q, 0., edge)
            ce = family.base.covariant(Je, Ee[:, :, :144])
            for degree in range(3):
                membrane.append((w*p**degree) @ ce[:, tangent, tangent])
            # Retain every thickness moment, including dependent/zero odd rows.
            face = np.zeros((order**2, 3))
            face[:, tangent] = np.repeat(p, order)
            face[:, 1-tangent] = sign
            face[:, 2] = np.tile(p, order)
            Js, Es = family.raw(q, 0., face)
            cs = family.base.covariant(Js, Es[:, :, :144])
            for i in range(3):
                for j in range(3):
                    shear.append(np.outer(w*p**i, w*p**j).reshape(-1) @ cs[:, tangent, 2])
    return dict(N=normal, M=np.array(membrane), S=np.array(shear), W=work)


def assemble(a, order):
    free_index = {int(node): i for i, node in enumerate(a['free'])}
    groups = {name: [] for name in ['N', 'M', 'S', 'W']}
    for _, _, ids in a['cells']:
        blocks = local(a['coords'][ids].reshape(144), order)
        for name, block in blocks.items():
            global_block = np.zeros((len(block), 6*len(a['free'])))
            for node, point in enumerate(ids):
                if int(point) in free_index:
                    j = free_index[int(point)]
                    global_block[:, 6*j:6*j+6] = block[:, 9*node:9*node+6]
            groups[name].append(global_block)
    return {name: np.vstack(rows) for name, rows in groups.items()}


def normalized(H, scale):
    value = H * scale
    norm = np.linalg.norm(value, axis=1)
    active = norm > 1e-13 * norm.max()
    return value[active] / norm[active, None], int(np.count_nonzero(~active))


def kernel(H, scale, a):
    A, dropped = normalized(H, scale)
    # Full right vectors are needed when constraints have fewer rows than columns.
    _, s, V = svd(A, full_matrices=True, check_finite=False)
    ranks = [int(np.count_nonzero(s > t*s[0])) for t in [1e-9, 1e-10, 1e-11]]
    basis, shift, orth = anchor(a, V[ranks[1]:].T, scale, .4)
    residual = float(np.linalg.norm(A @ basis) / (np.linalg.norm(A)*np.linalg.norm(basis)))
    return scale[:, None]*basis, dict(ranks=ranks, dropped_rows=dropped,
        nullity=basis.shape[1], anchor_shift=shift, orthogonality=orth,
        residual=residual), s


def relative_constraint(H, x):
    return float(np.linalg.norm(H @ x) / max(1e-300, np.linalg.norm(H)*np.linalg.norm(x)))
