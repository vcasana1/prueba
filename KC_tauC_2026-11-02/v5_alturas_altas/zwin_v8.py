# zwin_v8.py -- ceros de Z(t) en [T0, T0+W] a gran altura, sin tablas externas.
#   Z(t) = 2 Re[ e^{i theta(t)} sum_{n<=N} n^{-1/2} e^{-i t log n} ] + (-1)^{N-1} (2pi/t)^{1/4} C0(p),
#   N = floor(sqrt(t/2pi)), p = sqrt(t/2pi) - N, C0(p) = cos(2pi(p^2-p-1/16))/cos(2pi p)   (Riemann-Siegel, error O(t^{-3/4}))
#   La suma se evalua en la malla t = (g0+m) h por NUFFT de tipo 1 (FINUFFT); fases con long double.
#   theta(T0) en mpmath (40 cifras); theta(T0+u) = theta(T0) + u theta'(T0) + u^2 theta''(T0)/2 + u^3 theta'''(T0)/6.
#   Ceros: cambios de signo + raiz de la cubica de Hermite (Z, Z') en cada celda.
#   Indice (Turing): N(t) = n0 + #{ceros hallados <= t};  n0 - K0 = round(mean(N0(t) - K0 - count(t))),
#   N0(t) = (t/2pi) log(t/2pi e) + 7/8,  K0 = floor(N0(T0)).
# usage: python3 zwin_v8.py T0 W h_inv outprefix
import numpy as np, finufft, sys, math, mpmath as mp
mp.mp.dps = 40
T0s, W, HI, pref = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
T0 = mp.mpf(T0s); h = 1.0 / HI
g0 = int(mp.floor(T0 * HI)); T0 = mp.mpf(g0) / HI          # T0 exactamente en la malla
M = int(math.ceil(W * HI)) + 1
tm = np.arange(M) * h                                       # u = t - T0
th0 = mp.siegeltheta(T0); th1 = mp.log(T0 / (2 * mp.pi)) / 2; th2 = 1 / (2 * T0); th3 = -1 / (2 * T0 ** 2)
th0m = float(mp.fmod(th0, 2 * mp.pi))
theta = th0m + float(th1) * tm + float(th2) * tm ** 2 / 2 + float(th3) * tm ** 3 / 6
dtheta = float(th1) + float(th2) * tm
tf = float(T0) + tm
Nt = np.floor(np.sqrt(tf / (2 * np.pi))).astype(np.int64); Nmin = int(Nt.min()); Nmax = int(Nt.max())
n = np.arange(1, Nmin + 1, dtype=np.float64)
lnL = np.log(np.arange(1, Nmin + 1, dtype=np.longdouble))
c = n ** -0.5
# NUFFT tipo 1 con isign=-1:  f_k = sum_j c_j e^{-i k x_j},  x_j = om_j h ;  t = (g0 + Mc + k) h
M2 = 1 << (M - 1).bit_length(); Mc = M2 // 2
TWOPI = np.longdouble('6.2831853071795864769252867665590058')
ph = np.mod(lnL * np.longdouble(g0 + Mc) * np.longdouble(h), TWOPI).astype(np.float64)
c0 = c * np.exp(-1j * ph); om = lnL.astype(np.float64)
f = finufft.nufft1d1(om * h, np.stack([c0, c0 * (-1j) * om]), M2, isign=-1, eps=1e-13)
# modo k corresponde a m = k + Mc  (k = -Mc..Mc-1)  -> m = 0..M2-1
F = f[0][:M]; dF = f[1][:M]
for nn in range(Nmin + 1, Nmax + 1):                          # terminos extra donde N(t) > Nmin
    msk = Nt >= nn; ln = math.log(nn)
    e = np.exp(-1j * np.mod(np.longdouble(ln) * (np.longdouble(g0) + np.arange(M)) * np.longdouble(h), TWOPI).astype(np.float64)) / math.sqrt(nn)
    F[msk] += e[msk]; dF[msk] += (-1j * ln) * e[msk]
E = np.exp(1j * theta)
p = np.sqrt(tf / (2 * np.pi)) - Nt
C0 = np.cos(2 * np.pi * (p * p - p - 1 / 16)) / np.cos(2 * np.pi * p)
R = np.where(Nt % 2 == 1, 1.0, -1.0) * (2 * np.pi / tf) ** 0.25 * C0
Z = 2 * np.real(E * F) + R
dZ = 2 * np.real(1j * dtheta * E * F + E * dF)
# ceros
s = np.sign(Z); idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
roots = []
for i in idx:
    z0, z1, d0, d1 = Z[i], Z[i + 1], dZ[i] * h, dZ[i + 1] * h
    a, b = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (a + b)
        val = z0 * (2 * mid**3 - 3 * mid**2 + 1) + d0 * (mid**3 - 2 * mid**2 + mid) + z1 * (-2 * mid**3 + 3 * mid**2) + d1 * (mid**3 - mid**2)
        if np.sign(val) == np.sign(z0): a = mid
        else: b = mid
    roots.append(tm[i] + 0.5 * (a + b) * h)
roots = np.array(roots)
# indice (Turing): media de N0 - count sobre la malla
N0T0 = (T0 / (2 * mp.pi)) * mp.log(T0 / (2 * mp.pi * mp.e)) + mp.mpf(7) / 8
K0 = int(mp.floor(N0T0)); fr0 = float(N0T0 - K0)
dN0 = (tm * float(mp.log(T0 / (2 * mp.pi * mp.e))) + tf * np.log1p(tm / float(T0))) / (2 * np.pi)
cnt = np.searchsorted(roots, tm, side='right')
q = fr0 + dN0 - cnt
off = q.mean(); n0rel = int(round(off))
print(f"T0={mp.nstr(T0,20)} W={W} h=1/{HI} Nterms={Nmin}..{Nmax} ceros={roots.size}  K0={K0}  mean(N0-K0-count)={off:.4f} -> n0-K0={n0rel}  (resto {off-n0rel:+.4f})")
print(f"  S medio estimado en la ventana = {off-n0rel:+.4f} ; min/max de S: {q.min()-n0rel:+.3f} / {q.max()-n0rel:+.3f}")
np.savez(pref + '_zeros.npz', roots=roots, T0=mp.nstr(T0, 30), K0=K0, n0rel=n0rel, fr0=fr0, h=h)
