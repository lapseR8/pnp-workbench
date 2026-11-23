import numpy as np, itertools
from math import comb

R, C = 4, 5
NE = R * C
INJ = np.array(list(itertools.permutations(range(C), R)))  # 120 x 4
# for each row i and column j: indices of injections with f(i)=j
SEL = [[np.where(INJ[:, i] == j)[0] for j in range(C)] for i in range(R)]
SM = [np.array([[1.0 if INJ[f, i] == j else 0.0 for j in range(C)] for f in range(len(INJ))]) for i in range(R)]


def esym(vals):
    m = vals.shape[-1]
    e = np.zeros(vals.shape[:-1] + (m + 1,), dtype=complex)
    e[..., 0] = 1
    for i in range(m):
        v = vals[..., i:i + 1]
        e[..., 1:] = e[..., 1:] + v * e[..., :-1]
    return e


def _es(cols):
    """cols: list of arrays (same shape). returns list [e1..em]."""
    e = [np.ones_like(cols[0])]
    for v in cols:
        ne = [e[0]]
        for k in range(1, len(e)):
            ne.append(e[k] + v * e[k - 1])
        ne.append(v * e[-1])
        e = ne
    return e[1:]


def forms(t, grad=False):
    """t: N x 4 x 5. E: N x 4 level forms (sum over injections of e_k), G: N x 4 x 20."""
    N = t.shape[0]
    V = [np.ascontiguousarray(t[:, i, INJ[:, i]]) for i in range(R)]  # each N x 120
    e = _es(V)
    E = np.stack([x.sum(axis=1) for x in e], axis=1)
    if not grad:
        return E
    G = np.zeros((N, R, NE), dtype=complex)
    for i in range(R):
        eo = _es([V[k] for k in range(R) if k != i])  # e1,e2,e3 of others
        eo = np.stack([np.ones_like(V[0])] + eo, axis=1)  # N x 4 x 120
        G[:, :, i * C:(i + 1) * C] = eo @ SM[i]
    return E, G


class Prob:
    """Find d-dim subspace inside span(rows of W) (W: m x 20) killing all level forms.
    Unknown X: d x (m-d); basis = [I X] @ U, U = random unitary mix of W rows."""

    def __init__(self, d, W=None, rng=None, npts=None, wts=(1., 1., 1., 1.)):
        self.d = d
        self.rng = rng if rng is not None else np.random.default_rng()
        if W is None:
            W = np.eye(NE)
        W = np.asarray(W, dtype=complex)
        Wq, _ = np.linalg.qr(W.T)
        m = Wq.shape[1]
        self.m = m
        Q, _ = np.linalg.qr(self.rng.normal(size=(m, m)) + 1j * self.rng.normal(size=(m, m)))
        self.U = Q.conj().T @ Wq.T  # m x 20
        n4 = comb(d + 3, 4)
        self.npts = npts or int(n4 * 1.3) + 12
        c = self.rng.normal(size=(self.npts, d)) + 1j * self.rng.normal(size=(self.npts, d))
        self.c = c / np.linalg.norm(c, axis=1, keepdims=True)
        self.w = np.array(wts)

    def basis(self, X):
        return np.hstack([np.eye(self.d), X]) @ self.U

    def resid(self, X, jac=False):
        B = self.basis(X)
        t = (self.c @ B).reshape(-1, R, C)
        if jac:
            E, G = forms(t, True)
            E = E * self.w
            G = G * self.w[None, :, None]
            GU = G @ self.U[self.d:].T  # N x 4 x (m-d)
            J = np.einsum('nm,nka->nkma', self.c, GU)
            return E.ravel(), J.reshape(E.size, -1)
        return (forms(t) * self.w).ravel()

    def solve(self, X0, iters=300, tol=1e-26, verbose=False):
        X = X0.copy()
        lam = 1e-2
        r, J = self.resid(X, True)
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
                    lam *= 10
                    continue
                Xn = X + dx.reshape(X.shape)
                rn = self.resid(Xn)
                fn = np.vdot(rn, rn).real
                if fn < f:
                    X = Xn
                    lam = max(lam / 3, 1e-12)
                    ok = True
                    break
                lam *= 4
            if not ok:
                break
            r, J = self.resid(X, True)
            f = np.vdot(r, r).real
            if verbose and it % 10 == 0:
                print(it, f, lam, np.abs(X).max(), flush=True)
            if f < tol:
                break
        return X, f, it


def check_basis(B, ntest=400, rng=None):
    """max |level form| on random unit combos of orthonormalized B."""
    rng = rng if rng is not None else np.random.default_rng(1)
    Q, _ = np.linalg.qr(B.T)
    Bo = Q.T
    d = Bo.shape[0]
    c = rng.normal(size=(ntest, d)) + 1j * rng.normal(size=(ntest, d))
    c /= np.linalg.norm(c, axis=1, keepdims=True)
    E = forms((c @ Bo).reshape(-1, R, C))
    return np.abs(E).max(axis=0)
