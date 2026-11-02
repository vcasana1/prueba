# maxS_v11.py -- max |S(t)| (S = N - N0 at the jumps, both sides) against loglog t and H(t), by dyadic block, from zeross.bin
import numpy as np, math
g = np.fromfile('zeross.bin'); k = np.arange(1, g.size + 1)
N0 = g / (2*math.pi) * np.log(g / (2*math.pi*math.e)) + 0.875
Sp, Sm = k - N0, k - 1 - N0                        # S just after / just before gamma_k (extreme values occur there)
lt = np.log(np.log(np.maximum(g, math.e**math.e))); H = 0.08382 * np.log(g) + 5.5
print("block   max|S|   max |S|/loglog t   max |S|/H   #zeros with |S|>loglog t")
for j in range(4, 21):
    m = (g >= 2**j) & (g < 2**(j+1))
    if not m.any(): continue
    a = np.maximum(np.abs(Sp[m]), np.abs(Sm[m]))
    print(f" 2^{j:2d}   {a.max():6.3f}   {(a/lt[m]).max():8.3f}          {(a/H[m]).max():6.3f}     {(a > lt[m]).sum()}")
