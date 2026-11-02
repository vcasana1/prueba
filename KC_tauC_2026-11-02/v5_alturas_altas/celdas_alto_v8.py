# celdas_alto_v8.py -- prueba de celdas en una ventana alta: N(tauc_k) in {k-1,k}, margenes y ternas apretadas
# usage: python3 celdas_alto_v8.py tag
import numpy as np, sys, mpmath as mp
tag = sys.argv[1]
d = np.load(f'z_{tag}_zeros.npz'); T0 = float(mp.mpf(str(d['T0']))); r = d['roots']; K0 = int(d['K0']); n0 = K0 + int(d['n0rel'])
tc = np.loadtxt(f'tc_{tag}.txt'); k = tc[:, 0].astype(np.int64); x = tc[:, 1]          # x = tauc - T0
# ceros con indice absoluto: zero n0+1+i esta en r[i]
idx = k - (n0 + 1)                                                                     # posicion de gamma_k en r
ok = (idx >= 1) & (idx < r.size - 1); k, x, idx = k[ok], x[ok], idx[ok]
g = r[idx]; gl = r[idx - 1]; gr = r[idx + 1]
sp = 2 * np.pi / np.log((T0 + g) / (2 * np.pi))
Nt = n0 + np.searchsorted(r, x, side='right')                                         # N(tauc_k)
off = Nt - k
u, c = np.unique(off, return_counts=True)
print(f"[{tag}] T0={T0:.6e}  taucs comprobados={k.size}  (indices {k[0]}..{k[-1]})")
print("  N(tauc_k)-k :", dict(zip(u.tolist(), c.tolist())), "   (fuera de {-1,0}:", int(((off != 0) & (off != -1)).sum()), ")")
err = (x - g) / sp; print(f"  RMS(tauc-gamma) = {np.sqrt((err**2).mean()):.5f} esp ; max |.| = {np.abs(err).max():.4f}")
m = np.minimum(x - gl, gr - x) / sp; o = np.argsort(m)[:8]
print("  margen minimo en la celda (gamma_{k-1},gamma_{k+1}):", f"{m.min():+.4f} esp;", [(int(k[i]), round(float(m[i]), 4)) for i in o])
spr = 2 * np.pi / np.log((T0 + r) / (2 * np.pi)); w = (r[2:] - r[:-2]) / spr[1:-1]; o = np.argsort(w)[:8]
pos = {int(kk): ii for ii, kk in enumerate(k)}
print("  ternas mas apretadas (k del central, anchura, margen del tauc central):")
for i in o:
    kc = n0 + 1 + i + 1
    mm = m[pos[kc]] if kc in pos else float('nan')
    print(f"     k={kc}  anchura={w[i]:.4f} esp  margen={mm:+.4f}")
bad = np.nonzero((off != 0) & (off != -1))[0]
if bad.size:
    print("  FUERA DE CELDA:"); [print(f"     k={int(k[i])} N(tauc)-k={int(off[i])} err={err[i]:+.4f}") for i in bad[:30]]
