# budget_j1_v9.py -- explicit budget for an unconditional j=1 proportion theorem, block (2^j, 2^{j+1}]
#   B/DN <= sqrt(mean D^4) sqrt(mean Phi_c'^2) / mean N0'      (Cauchy-Schwarz; D = N - F)
#   prints mean D^2, mean D^4, mean Phi_c'^2 / mean N0'^2, CS bound for B/DN, and the needed A/DN for 41.7% and >0
# usage: python3 budget_j1_v9.py zeross.bin phiprefix j
import numpy as np, sys, math
g = np.fromfile(sys.argv[1]); pp = sys.argv[2]; j = int(sys.argv[3])
K, gmax, HINV = [int(x) for x in open(pp + '_meta.txt').read().split()[:3]]; h = 1.0 / HINV
phi = np.memmap(pp + '_phi.bin', dtype=np.float64, mode='r'); dphi = np.memmap(pp + '_dphi.bin', dtype=np.float64, mode='r')
a, b = 2.0**j, min(2.0**(j + 1), g[-1]); ga, gb = int(math.ceil(a * HINV)), int(math.floor(b * HINV))
s2 = s4 = sp2 = sn = sn2 = sB = 0.0; L = 0
for c0 in range(ga, gb, 1 << 22):
    c1 = min(c0 + (1 << 22), gb); t = (np.arange(c0, c1) + 0.5) * h
    ph = 0.5 * (phi[c0:c1] + phi[c0 + 1:c1 + 1]); dph = dphi[c0:c1]
    hc = np.log(np.log(t)); cl = np.clip(ph, -hc, hc); dcl = np.where(np.abs(ph) > hc, np.sign(ph) / (t * np.log(t)), dph)
    N = np.searchsorted(g, t, side='right'); N0 = t / (2*math.pi) * np.log(t / (2*math.pi*math.e)) + 0.875
    N0p = np.log(t / (2*math.pi)) / (2*math.pi); D = N - N0 - cl
    s2 += (D**2).sum(); s4 += (D**4).sum(); sp2 += (dcl**2).sum(); sn += N0p.sum(); sn2 += (N0p**2).sum(); sB += (D**2 * dcl).sum(); L += c1 - c0
m2, m4, mp2, mn = s2 / L, s4 / L, sp2 / L, sn / L
print(f"block 2^{j}:  mean D^2 = {m2:.5f}   mean D^4 = {m4:.5f}  (Gaussian would be {3*m2*m2:.5f})")
print(f"  RMS Phi_c' / mean N0' = {math.sqrt(mp2)/mn:.4f}")
print(f"  B/DN measured = {sB/L/mn:+.5f}    Cauchy-Schwarz bound = {math.sqrt(m4*mp2)/mn:.5f}")
print(f"  A/DN + B/DN measured = {m2 + sB/L/mn:.5f}  (= 1/12 + eta)")
for name, eta in [("record 41.7%", (1 - 0.417) / 4), ("positive proportion", 0.25)]:
    print(f"  {name}: need A/DN + B/DN < {1/12 + eta:.4f};  with CS for B: need mean D^2 < {1/12 + eta - math.sqrt(m4*mp2)/mn:.4f}")
