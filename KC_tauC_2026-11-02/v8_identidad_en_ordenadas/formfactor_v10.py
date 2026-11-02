# formfactor_v10.py -- structure factor (pair-correlation form factor) of zeros and of taucs, alpha in (0, 2]
#   unfolded points x_k = N0(p_k), N0(t) = t/2pi log(t/2pi e) + 7/8 (mean spacing 1)
#   S(alpha) = |sum_k e^{2 pi i alpha x_k}|^2 / n, averaged over alpha-bins of width DA  (GUE form factor: min(alpha,1);
#   equispaced points: 0 on (0,1), Bragg peak at 1)
#   coherence C(alpha) = |(1/n) sum_k e^{2 pi i alpha (x_k^zero - x_k^set)}| (point-by-point alignment with the zeros)
#   continuous part (optional): Sc(alpha) = |int e^{2 pi i alpha N0(t)} dF(t)|^2 / n, F = N0 + clip(Phi), from a phi grid;
#   S_set - Sc = rearrangement + discreteness (saw-tooth) part
# usage: python3 formfactor_v10.py zeross.bin k1 k2 DA name1=file1 [name2=file2 ...] [grid=phiprefix]
import numpy as np, finufft, sys, math
g = np.fromfile(sys.argv[1]); k1, k2, DA = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
sets = [a.split('=') for a in sys.argv[5:] if not a.startswith('grid=')]
grid = [a[5:] for a in sys.argv[5:] if a.startswith('grid=')]
N0 = lambda t: t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 0.875
n = k2 - k1; AMAX = 2.0
xz = N0(g[k1:k2]); xc = 0.5 * (xz[0] + xz[-1]); P = 1.05 * (xz[-1] - xz[0]) + 10   # period (unfolded units)
Mmodes = 2 * int(math.ceil(AMAX * P)) + 2                                          # modes m <-> alpha = m/P
def ft(x, w=None):
    xp = 2 * math.pi * (x - xc) / P
    c = np.ones(x.size, complex) if w is None else w.astype(complex)
    f = finufft.nufft1d1(xp, c, Mmodes, isign=1, eps=1e-12)                       # modes -M/2..M/2-1
    m = np.arange(-Mmodes // 2, Mmodes // 2); keep = m > 0
    return m[keep] / P, f[keep]
def binned(al, v):
    edges = np.arange(0, AMAX + DA / 2, DA); idx = np.digitize(al, edges) - 1
    s = np.bincount(idx, weights=v, minlength=edges.size)[:edges.size - 1]; c = np.bincount(idx, minlength=edges.size)[:edges.size - 1]
    return 0.5 * (edges[:-1] + edges[1:]), s / np.maximum(c, 1)
al, fz = ft(xz); a_mid, Sz = binned(al, np.abs(fz)**2 / n)
cols = [('zeros', Sz)]; coh = []
for name, fn in sets:
    p = np.fromfile(fn)[k1:k2]; x = N0(p)
    _, f = ft(x); cols.append((name, binned(al, np.abs(f)**2 / n)[1]))
    d = xz - x
    coh.append((name, np.array([abs(np.exp(2j * math.pi * a * d).mean()) for a in a_mid])))
if grid:
    K, gmax, HINV = [int(v) for v in open(grid[0] + '_meta.txt').read().split()[:3]]; h = 1.0 / HINV
    phi = np.memmap(grid[0] + '_phi.bin', dtype=np.float64, mode='r'); dphi = np.memmap(grid[0] + '_dphi.bin', dtype=np.float64, mode='r')
    ga, gb = int(math.ceil(g[k1] * HINV)), int(math.floor(g[k2 - 1] * HINV))
    t = (np.arange(ga, gb) + 0.5) * h; ph = 0.5 * (phi[ga:gb] + phi[ga + 1:gb + 1]); dph = np.array(dphi[ga:gb])
    hc = np.log(np.log(t)); cl = np.clip(ph, -hc, hc); dcl = np.where(np.abs(ph) > hc, np.sign(ph) / (t * np.log(t)), dph)
    w = (np.log(t / (2 * math.pi)) / (2 * math.pi) + dcl) * h                    # dF = F' dt
    _, f = ft(N0(t), w); cols.append(('cont_dF', binned(al, np.abs(f)**2 / n)[1]))
print(f"# zeros k in [{k1},{k2}) n={n}  heights [{g[k1]:.0f},{g[k2-1]:.0f}]  bin DA={DA}")
print("# alpha   " + "  ".join(f"{c[0]:>12s}" for c in cols) + "   GUE  | " + "  ".join(f"coh_{c[0]:>9s}" for c in coh))
for i, a in enumerate(a_mid):
    print(f"{a:7.3f}  " + "  ".join(f"{c[1][i]:12.4f}" for c in cols) + f"   {min(a,1):.3f} | " + "  ".join(f"{c[1][i]:13.4f}" for c in coh))
