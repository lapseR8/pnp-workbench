"""LM search for d-dim subspaces inside {S=0} that are Q-isotropic by construction
(Q = 2*e2 restricted), so only e3, e4 are residuals."""
import numpy as np, sys, time, os
from math import comb
from core import forms, check_basis, R, C, NE

# Q matrix: Q(t) = sum t^2 - sum R^2 - sum C^2 + S^2  (= 2 e2)
def qmat():
    Q = np.eye(NE)
    for i in range(R):
        v = np.zeros(NE); v[i * C:(i + 1) * C] = 1
        Q -= np.outer(v, v)
    for j in range(C):
        v = np.zeros(NE); v[j::C] = 1
        Q -= np.outer(v, v)
    Q += np.ones((NE, NE))
    return Q

QM = qmat()


def split_basis(rng):
    """random basis E(9), F(9), G(1) of {S=0} with Q(E_i,F_j)=d_ij, others 0, Q(G,G)=1."""
    # orthonormal basis of S=0 hyperplane wrt Euclidean, then Q-orthonormalize (complex)
    one = np.ones(NE) / np.sqrt(NE)
    M = rng.normal(size=(NE, NE)) + 1j * rng.normal(size=(NE, NE))
    M = M - np.outer(one, one @ M)
    # Q-Gram-Schmidt on 19 vectors
    vs = []
    k = 0
    while len(vs) < NE - 1:
        v = M[:, k]; k += 1
        for w in vs:
            v = v - (w @ QM @ v) * w
        n2 = v @ QM @ v
        if abs(n2) < 1e-8:
            continue
        vs.append(v / np.sqrt(n2))
    V = np.array(vs)  # 19 x 20, Q-orthonormal
    E = (V[0:18:2] + 1j * V[1:18:2]) / np.sqrt(2)
    F = (V[0:18:2] - 1j * V[1:18:2]) / np.sqrt(2)
    G = V[18:19]
    return E, F, G


class IsoProb:
    def __init__(self, d, rng, npts=None):
        self.d = d
        self.rng = rng
        self.E, self.F, self.G = split_basis(rng)
        n4 = comb(d + 3, 4)
        self.npts = npts or int(n4 * 1.3) + 12
        c = rng.normal(size=(self.npts, d)) + 1j * rng.normal(size=(self.npts, d))
        self.c = c / np.linalg.norm(c, axis=1, keepdims=True)
        self.iu = np.triu_indices(9, 1)
        self.np_ = (9 - d) * d + 36 + 9

    def unpack(self, th):
        d = self.d
        Cm = th[:(9 - d) * d].reshape(9 - d, d)
        s = th[(9 - d) * d:(9 - d) * d + 36]
        b = th[(9 - d) * d + 36:]
        S = np.zeros((9, 9), complex)
        S[self.iu] = s
        S = S - S.T
        return Cm, S, b

    def basis(self, th):
        Cm, S, b = self.unpack(th)
        M = np.vstack([np.eye(self.d), Cm])  # 9 x d
        A = S - 0.5 * np.outer(b, b)
        return M.T @ self.E + (A @ M).T @ self.F + np.outer(b @ M, self.G[0])

    def dbasis(self, th):
        n = len(th)
        out = np.zeros((n, self.d, NE), complex)
        for k in range(n):
            e = np.zeros(n, complex); e[k] = 1.0
            out[k] = (self.basis(th + e) - self.basis(th - e)) / 2
        return out

    def resid(self, th, jac=False):
        B = self.basis(th)
        t = (self.c @ B).reshape(-1, R, C)
        if not jac:
            return forms(t)[:, 2:].ravel()
        Ef, G = forms(t, True)
        dB = self.dbasis(th)  # n x d x 20
        # dt/dth_k = c @ dB_k : N x n x 20
        dt = np.einsum('nm,kmp->nkp', self.c, dB)
        J = np.einsum('nqp,nkp->nqk', G[:, 2:, :], dt)
        return Ef[:, 2:].ravel(), J.reshape(-1, len(th))

    def solve(self, th, iters=300, tol=1e-28, verbose=False):
        lam = 1e-2
        r, J = self.resid(th, True)
        f = np.vdot(r, r).real
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
                thn = th + dx
                rn = self.resid(thn)
                fn = np.vdot(rn, rn).real
                if fn < f:
                    th = thn; lam = max(lam / 3, 1e-12); ok = True
                    break
                lam *= 4
            if not ok:
                break
            r, J = self.resid(th, True)
            f = np.vdot(r, r).real
            if verbose and it % 10 == 0:
                print(it, f, lam, np.abs(th).max(), flush=True)
            if f < tol:
                break
        return th, f, it


if __name__ == '__main__':
    d = int(sys.argv[1]); n = int(sys.argv[2]); seed = int(sys.argv[3])
    out = sys.argv[4] if len(sys.argv) > 4 else None
    rng = np.random.default_rng(seed)
    t0 = time.time()
    fs = []
    for rs in range(n):
        P = IsoProb(d, rng)
        sc = rng.choice([0.2, 0.5, 1.0])
        th = sc * (rng.normal(size=P.np_) + 1j * rng.normal(size=P.np_))
        th, f, it = P.solve(th, iters=200)
        if f < 1e-4:
            th, f, it2 = P.solve(th, iters=2000)
            it += it2
        fs.append(f)
        B = P.basis(th)
        chk = check_basis(B) if f < 1e-8 else None
        print(rs, '%.3e' % f, it, '%.1f' % np.abs(th).max(), '%.0fs' % (time.time() - t0), chk, flush=True)
        if f < 1e-8 and out:
            np.save('%s_d%d_s%d_%d.npy' % (out, d, seed, rs), B)
    print('sorted', sorted(fs)[:10])
