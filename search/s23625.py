import json, time, os
src = open('upgen.py').read().replace("for N in range(1, K + 1):", "for N in range(K, K + 1):")
ns = {}; exec(compile(src, 'upgen_k', 'exec'), ns)
x = json.load(open('stmt_23625.json')); os.makedirs('prep', exist_ok=True); os.makedirs('results', exist_ok=True)
t = time.time(); txt = ns['solve_stmt'](x['statement'], x['module'])
print('witness' if txt else 'none', round(time.time() - t), flush=True)
if txt: open('prep/L-meta-23625.lean', 'w').write(txt)
