"""Time-bounded random-restart search. usage: cloud.py MODE D SEED SECONDS OUTFILE
MODE: gen (chart [I X]) | iso (isotropic param) | mix (alternate) | cont (iso, e3 first then ramp e4)."""
import numpy as np, sys, time
from core import Prob, check_basis
from iso import IsoProb

mode, d, seed, secs, outf = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
rng = np.random.default_rng(seed)
t0 = time.time()
hist = []
fo = open(outf, 'a')
fo.write('# mode=%s d=%d seed=%d\n' % (mode, d, seed)); fo.flush()
rs = 0
while time.time() - t0 < secs:
    m = mode if mode != 'mix' else ('iso' if rs % 2 == 0 else 'gen')
    extra = ''
    if m in ('iso', 'cont'):
        P = IsoProb(d, rng)
        th = rng.choice([0.2, 0.5, 1.0]) * (rng.normal(size=P.np_) + 1j * rng.normal(size=P.np_))
    else:
        P = Prob(d, rng=rng)
        th = rng.choice([0.3, 1.0, 3.0]) * (rng.normal(size=(d, P.m - d)) + 1j * rng.normal(size=(d, P.m - d)))
    if m == 'cont':
        P.w = np.array([1.0, 0.0])
        th, f3, it = P.solve(th, iters=300, tol=1e-24)
        extra = ' f3=%.2e' % f3
        if f3 < 1e-8:
            for w4 in [0.01, 0.03, 0.1, 0.3]:
                P.w = np.array([1.0, w4])
                th, f, it = P.solve(th, iters=200, tol=1e-28)
                extra += ' w%.2f:%.2e' % (w4, f)
        P.w = np.array([1.0, 1.0])
    th, f, it = P.solve(th, iters=200)
    if f < 1e-3:
        th, f, it2 = P.solve(th, iters=3000, tol=1e-28)
        it += it2
    hist.append(f)
    line = '%s %d f=%.3e it=%d t=%.0f%s' % (m, rs, f, it, time.time() - t0, extra)
    if f < 1e-5:
        B = P.basis(th)
        chk = check_basis(B)
        line += ' chk=' + ' '.join('%.1e' % x for x in chk)
        line += '\n%s ' % ('SOLUTION' if f < 1e-8 else 'CANDIDATE') + ' '.join('%.17g,%.17g' % (z.real, z.imag) for z in B.ravel())
    fo.write(line + '\n'); fo.flush()
    rs += 1
h = np.sort(np.array(hist))
fo.write('# done restarts=%d best=%s\n' % (rs, ' '.join('%.2e' % x for x in h[:10])))
fo.close()
