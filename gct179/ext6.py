"""Seeded extension search: find a random (d-1)-dim solution H by LM, then try to grow it
to a d-dim solution with LM over the full chart started at span(H, v) for random v.
usage: ext6.py D SEED SECONDS OUTFILE"""
import numpy as np, sys, time
from core import Prob, check_basis


def lm_restart(d, rng, iters=200):
    P = Prob(d, rng=rng)
    th = rng.choice([0.3, 1.0, 3.0]) * (rng.normal(size=(d, P.m - d)) + 1j * rng.normal(size=(d, P.m - d)))
    th, f, it = P.solve(th, iters=iters)
    if f < 1e-3:
        th, f, it = P.solve(th, iters=2000, tol=1e-28)
    return P, th, f


def grow(H, rng, noise=0.0, iters=300):
    """H: (d-1) x 20. returns (basis, f) after LM at dim d started at span(H, v)."""
    k = H.shape[0]
    d = k + 1
    v = rng.normal(size=20) + 1j * rng.normal(size=20)
    B0 = np.vstack([H, v])
    Q, _ = np.linalg.qr(B0.T, mode='complete')  # 20x20 unitary, first d cols span B0
    Qr, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    top = (Q[:, :d] @ Qr).T
    P = Prob(d, rng=rng)
    P.U = np.vstack([top, Q[:, d:].T])
    X = noise * (rng.normal(size=(d, 20 - d)) + 1j * rng.normal(size=(d, 20 - d)))
    X, f, it = P.solve(X, iters=iters)
    if f < 1e-3:
        X, f, it2 = P.solve(X, iters=3000, tol=1e-28)
    return P, X, f


if __name__ == '__main__':
    d, seed, secs, outf = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    tries = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    rng = np.random.default_rng(seed)
    fo = open(outf, 'a')
    fo.write('# ext6 d=%d seed=%d\n' % (d, seed)); fo.flush()
    t0 = time.time(); hist = []; ncore = 0; nfail = 0
    while time.time() - t0 < secs:
        P, th, f = lm_restart(d - 1, rng)
        if f > 1e-9:
            nfail += 1
            continue
        ncore += 1
        H = np.linalg.qr(P.basis(th).T)[0].T
        fs = []
        for k in range(tries):
            noise = [0.0, 0.0, 0.1, 0.3][k % 4]
            P2, X, f2 = grow(H, rng, noise)
            fs.append(f2)
            hist.append(f2)
            if f2 < 1e-8:
                B = P2.basis(X)
                fo.write('SOLUTION f=%.3e chk=%s\n' % (f2, ' '.join('%.1e' % x for x in check_basis(B))))
                fo.write('BASIS ' + ' '.join('%.17g,%.17g' % (z.real, z.imag) for z in B.ravel()) + '\n')
                fo.flush()
        fo.write('ext core=%d fail=%d best=%.3e all=%s t=%.0f\n' % (ncore, nfail, min(fs), ' '.join('%.1e' % x for x in fs), time.time() - t0))
        fo.flush()
    h = np.sort(np.array(hist))
    fo.write('# done cores=%d fails=%d grows=%d best=%s\n' % (ncore, nfail, len(h), ' '.join('%.2e' % x for x in h[:10])))
