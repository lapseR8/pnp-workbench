"""Time-bounded random-restart search. usage: cloud.py MODE D SEED SECONDS OUTFILE
MODE: gen (chart [I X]) | iso (isotropic param) | mix (alternate)."""
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
    if m == 'iso':
        P = IsoProb(d, rng)
        th = rng.choice([0.2, 0.5, 1.0]) * (rng.normal(size=P.np_) + 1j * rng.normal(size=P.np_))
    else:
        P = Prob(d, rng=rng)
        th = rng.choice([0.3, 1.0, 3.0]) * (rng.normal(size=(d, P.m - d)) + 1j * rng.normal(size=(d, P.m - d)))
    th, f, it = P.solve(th, iters=200)
    if f < 1e-3:
        th, f, it2 = P.solve(th, iters=3000, tol=1e-28)
        it += it2
    hist.append(f)
    line = '%s %d f=%.3e it=%d t=%.0f' % (m, rs, f, it, time.time() - t0)
    if f < 1e-8:
        B = P.basis(th)
        chk = check_basis(B)
        line += ' chk=' + ' '.join('%.1e' % x for x in chk)
        line += '\nSOLUTION ' + ' '.join('%.17g,%.17g' % (z.real, z.imag) for z in B.ravel())
    fo.write(line + '\n'); fo.flush()
    rs += 1
h = np.sort(np.array(hist))
fo.write('# done restarts=%d best=%s\n' % (rs, ' '.join('%.2e' % x for x in h[:10])))
fo.close()
