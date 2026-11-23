"""Grow named seed families (from attackers A/C) to dim 7: LM over the full d=7 chart started at
span(seed, random extra vectors) with small perturbation. usage: seeds.py D SEED SECONDS OUTFILE"""
import numpy as np, sys, time
from core import Prob, check_basis

z5 = np.exp(2j * np.pi / 5)


def M(d):
    t = np.zeros((4, 5), complex)
    for (i, j), v in d.items():
        t[i, j] += v
    return t.ravel()


def Fb(b):
    a = b / (1 - b); c = 3 * b - 1
    v3 = np.zeros((4, 5), complex); v3[:, 2] = [-1 / 3, -1 / 3, 1, -1 / 3]
    v4 = np.zeros((4, 5), complex); v4[0] = [0, 1, 0, -a, a - 1]
    v5 = np.zeros((4, 5), complex); v5[2] = [1, -b, 0, b - 1, 0]
    v6 = np.zeros((4, 5), complex); v6[2] = [0, -c, 0, c - 1, 1]
    return np.array([M({(0, 0): 1, (1, 0): -1}), M({(3, 0): 1, (1, 0): -1}), v3.ravel(), v4.ravel(), v5.ravel(), v6.ravel()])


def W6(c=0, rows=(0, 1, 2), w=(0, -1, 1, -1j, 1j)):
    vs = []
    for i in range(3):
        t = np.zeros((4, 5), complex); t[i, c] = 1; t[3, c] = -1; vs.append(t.ravel())
    for r in rows:
        t = np.zeros((4, 5), complex); t[r] = w; vs.append(t.ravel())
    return np.array(vs)


def Vz(m=1):
    return np.array([M({(i, j): z5 ** (m * j) for j in range(5)}) for i in range(4)])


def seed_list(rng):
    out = []
    bs = [0.25, 0.5, 2.0, -1.0, 1 / 3 + 0.0j, 0.5 + 0.5j, np.exp(2j * np.pi / 3), np.exp(1j * np.pi / 3), 2 + 1j]
    for b in bs:
        out.append(('Fb(%s)' % np.round(b, 3), Fb(b)))
    for k in range(3):
        b = complex(rng.normal(), rng.normal())
        out.append(('Fb(rand %s)' % np.round(b, 3), Fb(b)))
    out.append(('W6', W6()))
    out.append(('Vz1', Vz(1)))
    out.append(('Vz2', Vz(2)))
    return out


def grow(H, d, rng, noise, iters=300):
    k = H.shape[0]
    e = d - k
    B0 = np.vstack([H, rng.normal(size=(e, 20)) + 1j * rng.normal(size=(e, 20))])
    Q, _ = np.linalg.qr(B0.T, mode='complete')
    Qr, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    P = Prob(d, rng=rng)
    P.U = np.vstack([(Q[:, :d] @ Qr).T, Q[:, d:].T])
    X = noise * (rng.normal(size=(d, 20 - d)) + 1j * rng.normal(size=(d, 20 - d)))
    X, f, it = P.solve(X, iters=iters)
    if f < 1e-3:
        X, f, it = P.solve(X, iters=3000, tol=1e-28)
    return P, X, f


if __name__ == '__main__':
    d, seed, secs, outf = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    rng = np.random.default_rng(seed)
    fo = open(outf, 'a')
    fo.write('# seeds d=%d seed=%d\n' % (d, seed))
    seeds = seed_list(rng)
    for name, H in seeds:
        fo.write('# check %s dim=%d chk=%s\n' % (name, np.linalg.matrix_rank(H), ' '.join('%.1e' % x for x in check_basis(H))))
    fo.flush()
    t0 = time.time(); best = {}; cnt = {}
    rs = 0
    while time.time() - t0 < secs:
        name, H = seeds[rs % len(seeds)]
        if name.startswith('Fb(rand'):
            b = complex(rng.normal(), rng.normal()); H = Fb(b); name = 'Fb(rand)'
        H = np.linalg.qr(H.T)[0].T
        noise = [0.0, 0.03, 0.1, 0.3][(rs // len(seeds)) % 4]
        P, X, f = grow(H, d, rng, noise)
        best[name] = min(best.get(name, 1e9), f); cnt[name] = cnt.get(name, 0) + 1
        line = '%s noise=%.2f f=%.3e t=%.0f' % (name, noise, f, time.time() - t0)
        if f < 1e-8:
            B = P.basis(X)
            line += ' chk=' + ' '.join('%.1e' % x for x in check_basis(B))
            line += '\nSOLUTION ' + ' '.join('%.17g,%.17g' % (z.real, z.imag) for z in B.ravel())
        fo.write(line + '\n'); fo.flush()
        rs += 1
    for k in best:
        fo.write('# best %s n=%d f=%.3e\n' % (k, cnt[k], best[k]))
    fo.close()
