# margin_triangle_v14.py -- criterion  M_k = (tau_{k+1}-tau_k) - |tau_k-gamma_k| - |tau_{k+1}-gamma_{k+1}| > 0  for all k  =>  gamma_k < gamma_{k+1}
#   (triangle inequality: a doubled ordinate gamma_k = gamma_{k+1} forces M_k <= 0 for ANY points tau_k < tau_{k+1})
#   prints min M_k, count M_k <= 0, and the relation M_k vs zero gap g_k on the tight pairs (unfolded units)
# usage: python3 margin_triangle_v14.py zeross.bin k1 k2 name=file ...
import numpy as np, sys, math
g = np.fromfile(sys.argv[1]); k1, k2 = int(sys.argv[2]), int(sys.argv[3])
N0 = lambda t: t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 0.875
xz = N0(g[k1:k2]); gz = np.diff(xz)
for a in sys.argv[4:]:
    name, fn = a.split('='); x = N0(np.fromfile(fn)[k1:k2]); e = np.abs(x - xz)
    M = np.diff(x) - e[:-1] - e[1:]
    t = gz < 0.05
    print(f"{name:10s}: min M = {M.min():+.4f}   #M<=0: {(M<=0).sum():6d} of {M.size}   "
          f"tight pairs (g<0.05): median M = {np.median(M[t]):+.4f}, median g = {np.median(gz[t]):.4f}, corr(M,g) = {np.corrcoef(M[t], gz[t])[0,1]:+.3f}")
