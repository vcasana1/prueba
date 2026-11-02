# decomp_j1_v9.py -- decomposition of the j=1 identity on dyadic blocks, F = N0 + clip(Phi^o, +-hc):
#   sum_k u_k^2 = int (N-F)^2 F' dt - DeltaN/12 - bdry,   F' = N0' + (clip Phi)'
#   A = int (N-F)^2 N0' dt   (continuous mean square, the Selberg-type term)
#   B = int (N-F)^2 (clip Phi)' dt   (the term coming from the slope of the prime phase)
#   reports per block: mean_t (N-F)^2, A/DN, B/DN, (A+B)/DN - 1/12, eta measured from u_k (same F), GUE prediction of mean_t
# usage: python3 decomp_j1_v9.py zeross.bin phiprefix ukprefix W TH
import numpy as np, sys, math
from scipy.integrate import quad
zf, pp, up, prof, TH = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5])
g = np.fromfile(zf); u = np.fromfile(up + '_uk.bin')
K, gmax, HINV = [int(x) for x in open(pp + '_meta.txt').read().split()[:3]]; h = 1.0 / HINV
phi = np.memmap(pp + '_phi.bin', dtype=np.float64, mode='r'); dphi = np.memmap(pp + '_dphi.bin', dtype=np.float64, mode='r')
def Wf(x):
    v = min(max(x - 1.0, 0.0), 1.0)
    if x >= 2: return 0.0
    return (1 - v) if prof == 'selberg' else (1 - v + 0.25 * v * (1 - v))
# GUE / Montgomery heuristic: mean_t (S - Phi)^2 ~ (1/2pi^2) int F(a) (1 - W(a/TH))^2 / a^2 da, F(a) = min(a,1)  (primes only)
pred = (quad(lambda a: (1 - Wf(a / TH))**2 / a, TH, 1)[0] + quad(lambda a: (1 - Wf(a / TH))**2 / a**2, 1, 2 * TH)[0] * (TH > 0.5)
        + quad(lambda a: 1 / a**2, max(1, 2 * TH), np.inf)[0]) / (2 * math.pi**2)
print(f"== phi={pp}  u={up}   GUE prediction of mean_t (S-Phi)^2 (prime powers p only): {pred:.5f}   1/12 = {1/12:.5f}")
print("  block     DN      mean_t(N-F)^2   A/DN     B/DN    (A+B)/DN-1/12   eta(u_k)   diff")
for j in range(10, 21):
    a, b = 2.0**j, min(2.0**(j + 1), g[-1])
    ga, gb = int(math.ceil(a * HINV)), int(math.floor(b * HINV))
    if gb <= ga: continue
    A = B = C = 0.0; L = 0
    for c0 in range(ga, gb, 1 << 22):
        c1 = min(c0 + (1 << 22), gb)
        t = (np.arange(c0, c1) + 0.5) * h                      # midpoints; phi linear-interpolated with Hermite slope
        p0 = phi[c0:c1]; p1 = phi[c0 + 1:c1 + 1]; d0 = dphi[c0:c1]
        ph = 0.5 * (p0 + p1); dph = d0
        hc = np.log(np.log(np.maximum(t, math.e**math.e))); dhc = 1.0 / (t * np.log(t))
        cl = np.clip(ph, -hc, hc); dcl = np.where(ph > hc, dhc, np.where(ph < -hc, -dhc, dph))
        N = np.searchsorted(g, t, side='right')
        N0 = t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 0.875; N0p = np.log(t / (2 * math.pi)) / (2 * math.pi)
        D2 = (N - N0 - cl)**2
        A += (D2 * N0p).sum() * h; B += (D2 * dcl).sum() * h; C += D2.sum() * h; L += c1 - c0
    k1, k2 = np.searchsorted(g, a, side='right'), np.searchsorted(g, b, side='right')
    DN = k2 - k1; eta = (u[k1:k2]**2).mean()
    tot = (A + B) / DN - 1 / 12
    print(f"  2^{j:2d}  {DN:8d}    {C/(L*h):.5f}      {A/DN:.5f}  {B/DN:+.5f}    {tot:+.5f}      {eta:.5f}   {tot-eta:+.5f}")
