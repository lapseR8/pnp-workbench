"""Search for g-invariant d-dim subspaces (g in S4 x S5), L = sum of subspaces of eigenspaces.
usage: sym.py D SEED SECONDS OUTFILE"""
import numpy as np, sys, time, itertools
from math import comb
from core import forms, check_basis, R, C, NE

SIG = {'e': [0, 1, 2, 3], '(01)': [1, 0, 2, 3], '(01)(23)': [1, 0, 3, 2], '(012)': [1, 2, 0, 3], '(0123)': [1, 2, 3, 0]}
TAU = {'e': [0, 1, 2, 3, 4], '(01)': [1, 0, 2, 3, 4], '(01)(23)': [1, 0, 3, 2, 4], '(012)': [1, 2, 0, 3, 4],
       '(012)(34)': [1, 2, 0, 4, 3], '(0123)': [1, 2, 3, 0, 4], '(01234)': [1, 2, 3, 4, 0]}


def order(p):
    q = list(range(len(p))); k = 0
    while True:
        q = [p[x] for x in q]; k += 1
        if q == list(range(len(p))):
            return k


def eigenspaces(s, t):
    P = np.zeros((NE, NE))
    for i in range(R):
        for j in range(C):
            P[s[i] * C + t[j], i * C + j] = 1
    n = int(np.lcm(order(s), order(t)))
    spaces = {}
    for m in range(n):
        lam = np.exp(2j * np.pi * m / n)
        Pr = np.zeros((NE, NE), complex)
        Pk = np.eye(NE)
        for k in range(n):
            Pr += lam ** (-k) * Pk
            Pk = P @ Pk
        Pr /= n
        u, sv, vh = np.linalg.svd(Pr)
        r = int(np.sum(sv > 1e-8))
        if r:
            spaces[m] = u[:, :r].T  # r x 20 orthonormal rows
    return n, spaces


class RowProb:
    """basis row r: base[r] + th_r @ V[r] (V[r]: n_r x 20)."""

    def __init__(self, base, Vs, rng, npts=None, extra_lin=None):
        self.base = np.array(base, complex)
        self.Vs = Vs
        self.d = len(base)
        self.sizes = [V.shape[0] for V in Vs]
        self.off = np.concatenate([[0], np.cumsum(self.sizes)])
        self.np_ = int(self.off[-1])
        n4 = comb(self.d + 3, 4)
        self.npts = npts or int(n4 * 1.3) + 12
        c = rng.normal(size=(self.npts, self.d)) + 1j * rng.normal(size=(self.npts, self.d))
        self.c = c / np.linalg.norm(c, axis=1, keepdims=True)

    def basis(self, th):
        B = self.base.copy()
        for r in range(self.d):
            if self.sizes[r]:
                B[r] += th[self.off[r]:self.off[r + 1]] @ self.Vs[r]
        return B

    def resid(self, th, jac=False):
        B = self.basis(th)
        t = (self.c @ B).reshape(-1, R, C)
        if not jac:
            return forms(t).ravel()
        E, G = forms(t, True)
        J = np.zeros((self.npts, 4, self.np_), complex)
        for r in range(self.d):
            if self.sizes[r]:
                J[:, :, self.off[r]:self.off[r + 1]] = self.c[:, r][:, None, None] * (G @ self.Vs[r].T)
        return E.ravel(), J.reshape(-1, self.np_)

    def solve(self, th, iters=300, tol=1e-28):
        lam = 1e-2
        r, J = self.resid(th, True)
        f = np.vdot(r, r).real
        it = 0
        for it in range(iters):
            A = J.conj().T @ J
            g = J.conj().T @ r
            dA = np.diag(A).real + 1e-12
            ok = False
            while lam < 1e12:
                try:
                    dx = -np.linalg.solve(A + lam * np.diag(dA), g)
                except np.linalg.LinAlgError:
                    lam *= 10; continue
                rn = self.resid(th + dx)
                fn = np.vdot(rn, rn).real
                if fn < f:
                    th = th + dx; lam = max(lam / 3, 1e-12); ok = True
                    break
                lam *= 4
            if not ok:
                break
            r, J = self.resid(th, True)
            f = np.vdot(r, r).real
            if f < tol:
                break
        return th, f, it


def chart_block(E, k, rng):
    """random chart for k-dim subspace of rowspace(E) (m x 20): rows base_i + x_i @ V."""
    m = E.shape[0]
    Q, _ = np.linalg.qr(rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m)))
    U = Q.conj().T @ E
    return [U[i] for i in range(k)], [U[k:] for _ in range(k)]


def random_dimvec(dims, d, rng):
    keys = list(dims)
    for _ in range(1000):
        v = {k: 0 for k in keys}
        left = d
        order_ = list(rng.permutation(keys))
        for k in order_:
            v[k] = int(rng.integers(0, min(dims[k], left) + 1))
            left -= v[k]
        if left == 0:
            return v
    return None


if __name__ == '__main__':
    d, seed, secs, outf = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    only = sys.argv[5] if len(sys.argv) > 5 else None
    rng = np.random.default_rng(seed)
    classes = [(a, b) for a in SIG for b in TAU if (a, b) != ('e', 'e')]
    if only:
        classes = [tuple(only.split('|'))]
    fo = open(outf, 'a')
    t0 = time.time(); rs = 0
    while time.time() - t0 < secs:
        a, b = classes[rng.integers(len(classes))]
        n, sp = eigenspaces(SIG[a], TAU[b])
        dims = {m: E.shape[0] for m, E in sp.items()}
        v = random_dimvec(dims, d, rng)
        if v is None:
            continue
        base, Vs = [], []
        for m, k in v.items():
            if k:
                b_, V_ = chart_block(sp[m], k, rng)
                base += b_; Vs += V_
        P = RowProb(base, Vs, rng)
        th = rng.choice([0.3, 1.0, 2.0]) * (rng.normal(size=P.np_) + 1j * rng.normal(size=P.np_))
        if P.np_ == 0:
            r0 = P.resid(th); f = float(np.vdot(r0, r0).real); it = 0
        else:
            th, f, it = P.solve(th, iters=150)
        if f < 1e-3 and P.np_ > 0:
            th, f, it2 = P.solve(th, iters=2000)
        vs = ','.join('%d:%d' % (m, k) for m, k in sorted(v.items()) if k)
        line = '%s|%s n=%d dims=%s f=%.3e t=%.0f' % (a, b, n, vs, f, time.time() - t0)
        if f < 1e-8:
            B = P.basis(th)
            line += ' chk=' + ' '.join('%.1e' % x for x in check_basis(B))
            line += '\nSOLUTION ' + ' '.join('%.17g,%.17g' % (z.real, z.imag) for z in B.ravel())
        fo.write(line + '\n'); fo.flush()
        rs += 1
    fo.write('# done %d\n' % rs)
