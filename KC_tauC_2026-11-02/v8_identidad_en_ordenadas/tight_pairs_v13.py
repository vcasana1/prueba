# tight_pairs_v13.py -- what the taucs know about tight zero pairs (gap g < GMAX mean spacings), several prime rules
#   for each tight pair (k, k+1): taucs gap s, taucs midpoint error |m_tau - m_zero|, individual errors |tau-gamma| (unfolded units)
#   Spearman correlation between s and g  (information on g carried by the taucs below their resolution floor)
#   "GUE mean" replacement: s -> E_GUE[g | g < s_floor] = 3/4 s_floor  (p(g) ~ g^2): shows it carries no information on g
# usage: python3 tight_pairs_v13.py zeross.bin k1 k2 GMAX name=file ...
import numpy as np, sys, math
from scipy.stats import spearmanr
g = np.fromfile(sys.argv[1]); k1, k2, GMAX = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
N0 = lambda t: t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 0.875
xz = N0(g[k1:k2]); dz = np.diff(xz); idx = np.nonzero(dz < GMAX)[0]
print(f"# k in [{k1},{k2}): {idx.size} zero pairs with gap < {GMAX}  (median gap {np.median(dz[idx]):.4f})")
print("# rule        floor(min s)  median s(pairs)  Spearman(s,g)   median|mid err|  median|indiv err|  RMS indiv err (all k)")
for a in sys.argv[5:]:
    name, fn = a.split('='); x = N0(np.fromfile(fn)[k1:k2]); s = np.diff(x)
    mid_t = 0.5 * (x[idx] + x[idx + 1]); mid_z = 0.5 * (xz[idx] + xz[idx + 1])
    ind = np.abs(np.r_[x[idx] - xz[idx], x[idx + 1] - xz[idx + 1]])
    rho = spearmanr(s[idx], dz[idx]).correlation
    print(f"{name:12s}  {s.min():.4f}        {np.median(s[idx]):.4f}           {rho:+.3f}          {np.median(np.abs(mid_t-mid_z)):.5f}          {np.median(ind):.4f}            {np.sqrt(((x-xz)**2).mean()):.4f}")
