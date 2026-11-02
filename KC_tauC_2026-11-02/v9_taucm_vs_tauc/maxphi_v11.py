# maxphi_v11.py -- max |Phi^o(t)| / loglog t and / H(t) on the grid, by dyadic block (does any clip ever act?)
import numpy as np, sys, math
pp = sys.argv[1]; K, gmax, HINV = [int(v) for v in open(pp + '_meta.txt').read().split()[:3]]; h = 1.0 / HINV
phi = np.memmap(pp + '_phi.bin', dtype=np.float64, mode='r')
print(f"== {pp}:  block   max|Phi|   max |Phi|/loglog t   max |Phi|/H(t)")
for j in range(4, 21):
    a, b = int(2**j * HINV), min(int(2**(j + 1) * HINV), gmax)
    if a >= gmax: break
    t = np.arange(a, b) * h; p = np.abs(np.asarray(phi[a:b]))
    lt = np.log(np.log(np.maximum(t, math.e**math.e))); H = 0.08382 * np.log(np.maximum(t, math.e)) + 5.5
    print(f"   2^{j:2d}   {p.max():7.4f}   {(p/lt).max():8.4f}            {(p/H).max():8.4f}")
