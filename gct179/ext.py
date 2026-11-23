"""Extend a fixed core subspace (rows of W0, k x 20) by e extra vectors to a (k+e)-dim solution.
usage: ext.py CORE E SEED NREST"""
import numpy as np, sys, time
from core import check_basis, NE
from sym import RowProb

z5 = np.exp(2j * np.pi / 5)
w3 = np.exp(2j * np.pi / 3)


def M(d):
    t = np.zeros((4, 5), complex)
    for (i, j), v in d.items():
        t[i, j] += v
    return t.ravel()


def core(name):
    if name == 'rowcol':
        return np.array([M({(0, 1): 1, (0, 2): -1}), M({(0, 1): 1, (0, 3): -1}), M({(0, 1): 1, (0, 4): -1}),
                         M({(0, 0): 1, (1, 0): -1}), M({(1, 0): 1, (2, 0): -1}), M({(1, 0): 1, (3, 0): -1})])
    if name == 'rowcol2':  # B=0, A+x=0 variant
        return np.array([M({(1, 0): 1, (2, 0): -1}), M({(1, 0): 1, (3, 0): -1}),
                         M({(0, 0): 1, (0, 1): -1}), M({(0, 1): 1, (0, 2): -1}), M({(0, 1): 1, (0, 3): -1}), M({(0, 1): 1, (0, 4): -1})])
    if name.startswith('fourier'):
        m = int(name[7:]) if len(name) > 7 else 1
        return np.array([M({(i, j): z5 ** (m * j) for j in range(5)}) for i in range(4)])
    if name == 'g180':
        # rows0..2 = a+b, w a + w^2 b, w^2 a + w b; a=p0(1,-1,0,0,0), b=p1(0,0,1,-1,0); row3 sum-zero
        a = np.array([1, -1, 0, 0, 0]); b = np.array([0, 0, 1, -1, 0])
        def stack(x, y):
            t = np.zeros((4, 5), complex)
            t[0] = x + y; t[1] = w3 * x + w3 ** 2 * y; t[2] = w3 ** 2 * x + w3 * y
            return t.ravel()
        vs = [stack(a, 0 * a), stack(0 * a, b)]
        for j in range(4):
            vs.append(M({(3, j): 1, (3, 4): -1}))
        return np.array(vs)
    raise ValueError(name)


if __name__ == '__main__':
    name, e, seed, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    W0 = core(name)
    print('core check', check_basis(W0))
    k = W0.shape[0]
    rng = np.random.default_rng(seed)
    # complement of core
    Q, _ = np.linalg.qr(W0.T, mode='complete')
    comp = Q[:, k:].T  # (20-k) x 20
    best = []
    t0 = time.time()
    for rs in range(n):
        Qr, _ = np.linalg.qr(rng.normal(size=(20 - k, 20 - k)) + 1j * rng.normal(size=(20 - k, 20 - k)))
        U = Qr.conj().T @ comp
        base = [W0[i] for i in range(k)] + [U[i] for i in range(e)]
        Vs = [np.zeros((0, NE))] * k + [U[e:] for _ in range(e)]
        P = RowProb(base, Vs, rng)
        th = rng.choice([0.3, 1.0, 3.0]) * (rng.normal(size=P.np_) + 1j * rng.normal(size=P.np_))
        th, f, it = P.solve(th, iters=300)
        if f < 1e-4:
            th, f, it = P.solve(th, iters=3000)
        best.append(f)
        if f < 1e-8:
            print('SOLUTION?', rs, f, check_basis(P.basis(th)), flush=True)
            np.save('ext_%s_%d_%d.npy' % (name, seed, rs), P.basis(th))
    print(name, 'e=%d' % e, 'n=%d' % n, 'best', ['%.2e' % x for x in sorted(best)[:8]], '%.0fs' % (time.time() - t0))
