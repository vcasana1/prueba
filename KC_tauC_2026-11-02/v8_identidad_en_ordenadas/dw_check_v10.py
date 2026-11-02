# dw_check_v10.py -- Debye-Waller (Gaussian) prediction of the form factor of the taucs from the covariance of the prime phase
#   taucs (F monotone, unfolded):  x_k = k - 1/2 - phi_k,   phi = clip(Phi^o) evaluated at t_k = N0^{-1}(k - 1/2)
#   Gaussian stationary phi with covariance R(m) (m = lag in spacings), V = R(0):
#      S_DW(alpha) = sum_m cos(2 pi alpha m) exp(-4 pi^2 alpha^2 (V - R(m)))
#                  = Bragg comb * exp(-4pi^2 alpha^2 V)  +  diffuse,   diffuse = sum_m cos(2pi alpha m) e^{-4pi^2a^2V}(e^{4pi^2a^2R(m)} - 1)
#   compared with: measured taucs, the synthetic displaced lattice y_k = k - 1/2 - phi_k (tests the linearisation x ~ k-1/2-phi),
#   and the zeros.  R(m) is the empirical autocovariance of the pure-prime sequence phi_k.
# usage: python3 dw_check_v10.py zeross.bin rvm_tau.bin taucs.bin phiprefix k1 k2 DA
import numpy as np, finufft, sys, math
g = np.fromfile(sys.argv[1]); tau = np.fromfile(sys.argv[2]); tc = np.fromfile(sys.argv[3]); pp = sys.argv[4]
k1, k2, DA = int(sys.argv[5]), int(sys.argv[6]), float(sys.argv[7]); n = k2 - k1
K, gmax, HINV = [int(v) for v in open(pp + '_meta.txt').read().split()[:3]]; h = 1.0 / HINV
phi = np.memmap(pp + '_phi.bin', dtype=np.float64, mode='r')
N0 = lambda t: t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 0.875
t = tau[k1:k2]; gi = t * HINV; i0 = np.floor(gi).astype(np.int64); fr = gi - i0
ph = (1 - fr) * phi[i0] + fr * phi[i0 + 1]; hc = np.log(np.log(t)); ph = np.clip(ph, -hc, hc)
ph = ph - ph.mean()
k = np.arange(k1 + 1, k2 + 1); y = k - 0.5 - ph
# empirical autocovariance R(m), m = 0..MM
MM = 20000; L = 1 << int(math.ceil(math.log2(2 * n)))
f = np.fft.rfft(ph, L); ac = np.fft.irfft(np.abs(f)**2, L)[:MM + 1] / (n - np.arange(MM + 1))
V = ac[0]
AMAX = 2.0; xc = 0.5 * (k1 + k2); P = 1.05 * n + 10; Mmodes = 2 * int(math.ceil(AMAX * P)) + 2
def S_of(x):
    fx = finufft.nufft1d1(2 * math.pi * (x - xc) / P, np.ones(x.size, complex), Mmodes, isign=1, eps=1e-12)
    m = np.arange(-Mmodes // 2, Mmodes // 2); keep = m > 0
    al = m[keep] / P; v = np.abs(fx[keep])**2 / n
    edges = np.arange(0, AMAX + DA / 2, DA); idx = np.digitize(al, edges) - 1
    s = np.bincount(idx, weights=v, minlength=edges.size)[:edges.size - 1]; c = np.bincount(idx, minlength=edges.size)[:edges.size - 1]
    return 0.5 * (edges[:-1] + edges[1:]), s / np.maximum(c, 1)
am, Sz = S_of(N0(g[k1:k2])); _, St = S_of(N0(tc[k1:k2])); _, Sy = S_of(y)
mm = np.arange(1, MM + 1)
def diffuse(a):
    e = math.exp(-4 * math.pi**2 * a * a * V)
    return e * (math.exp(4 * math.pi**2 * a * a * V) - 1) + 2 * e * (np.cos(2 * math.pi * a * mm) * (np.exp(4 * math.pi**2 * a * a * ac[1:]) - 1)).sum()
def dw_bin(a0):
    return np.mean([diffuse(a) for a in np.linspace(a0 - DA / 2, a0 + DA / 2, 41)[1:-1]])
print(f"# k in [{k1},{k2}) n={n}; phi = clip Phi^o at tau_k (RvM); V = Var phi = {V:.4f};  Debye-Waller factor exp(-4pi^2 V) = {math.exp(-4*math.pi**2*V):.2e}")
print(f"# R(m)/V for m = 1,2,3,5,10,20,50,100,1000: " + " ".join(f"{ac[m]/V:.3f}" for m in [1, 2, 3, 5, 10, 20, 50, 100, 1000]))
print("# alpha     zeros     taucs   lattice-phi   DW-gauss(diffuse)")
for i, a in enumerate(am):
    print(f"{a:6.3f}   {Sz[i]:7.4f}   {St[i]:7.4f}   {Sy[i]:9.4f}   {dw_bin(a):9.4f}")
