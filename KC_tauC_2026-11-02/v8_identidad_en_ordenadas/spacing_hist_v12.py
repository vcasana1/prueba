# spacing_hist_v12.py -- nearest-neighbour spacing distribution of zeros and taucs (unfolded by N0), compared with GUE
#   GUE reference: Wigner surmise p(s) = 32/pi^2 s^2 exp(-4 s^2/pi)  (beta = 2), and its CDF.
#   reports: histogram (bin 0.1), KS distance to the surmise and to the zeros, small-spacing fractions,
#   fitted repulsion exponent b in P(s < x) ~ x^{b+1} for x in [0.05, 0.25], and spacing of the taucs at the 10 closest zero pairs.
# usage: python3 spacing_hist_v12.py zeross.bin k1 k2 name=file ...
import numpy as np, sys, math
g = np.fromfile(sys.argv[1]); k1, k2 = int(sys.argv[2]), int(sys.argv[3])
N0 = lambda t: t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 0.875
sets = [('zeros', g)] + [(a.split('=')[0], np.fromfile(a.split('=')[1])) for a in sys.argv[4:]]
from scipy.special import erf
cdf_w = lambda s: erf(2 * s / math.sqrt(math.pi)) - (4 * s / math.pi) * np.exp(-4 * s * s / math.pi)
pdf_w = lambda s: 32 / math.pi**2 * s * s * np.exp(-4 * s * s / math.pi)
sp = {}
for name, p in sets:
    x = N0(p[k1:k2]); s = np.diff(x); sp[name] = s
edges = np.arange(0, 3.01, 0.1); mid = 0.5 * (edges[1:] + edges[:-1])
print(f"# k in [{k1},{k2})  n = {k2-k1-1} spacings (unfolded, mean 1)")
print("# s       GUE(surmise) " + " ".join(f"{n:>10s}" for n, _ in sets))
H = {n: np.histogram(s, edges, density=True)[0] for n, s in sp.items()}
for i, m in enumerate(mid):
    print(f"{m:5.2f}   {np.mean(pdf_w(np.linspace(edges[i], edges[i+1], 11))):9.4f}   " + " ".join(f"{H[n][i]:10.4f}" for n, _ in sets))
print("# name        mean    var     KS vs GUE   KS vs zeros   P(s<0.1)   P(s<0.2)   P(s<0.05)  min s    exponent b (P(s<x)~x^(b+1))")
zs = np.sort(sp['zeros'])
for n, _ in sets:
    s = np.sort(sp[n]); F = np.arange(1, s.size + 1) / s.size
    ks_g = np.max(np.abs(F - cdf_w(s)))
    Fz = np.searchsorted(zs, s, side='right') / zs.size; ks_z = np.max(np.abs(F - Fz))
    xs = np.array([0.05, 0.1, 0.15, 0.2, 0.25]); P = np.array([(s < v).mean() for v in xs])
    b = np.polyfit(np.log(xs), np.log(np.maximum(P, 1e-12)), 1)[0] - 1 if (P > 0).all() else float('nan')
    print(f"{n:10s}  {s.mean():.4f}  {s.var():.4f}   {ks_g:.4f}      {ks_z:.4f}       {(s<0.1).mean():.5f}    {(s<0.2).mean():.5f}    {(s<0.05).mean():.6f}  {s.min():.4f}   {b:.2f}")
print(f"# GUE surmise: var = {3*math.pi/8 - 1:.4f};  P(s<0.1) = {cdf_w(0.1):.5f}, P(s<0.2) = {cdf_w(0.2):.5f}, P(s<0.05) = {cdf_w(0.05):.6f}; exponent b = 2")
o = np.argsort(sp['zeros'])[:10]
print("# 10 closest zero pairs (spacing of zeros, then of each taucs set at the same k):")
for i in o:
    print(f"   k={k1+i+1:8d}  " + "  ".join(f"{n}={sp[n][i]:.4f}" for n, _ in sets))
