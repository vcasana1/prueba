# uk_at_zeros_v9.py -- u_k = k - 1/2 - F(gamma_k) at every zero of zeross.bin, for the taucs function F of a run.
#   F = N0 + clip(Phi^o, -hc, hc),  hc = loglog max(t, e^e)   (same F as tauC_v1 / tauC_fromphi_v5)
#   Phi^o on [2^j, 2^{j+1}] = (2-u) Phi^(j) + (u-1) Phi^(j+1),  u = t/2^j
#   Phi^(L) = -(1/pi) sum_{n=p^r} c_n W_L(n) sin(t log n),  c_n = 1/(r p^{r/2})
#   schemes:  count C      : first floor(C 2^{L/2}) primes, W = 1          (IJNT draft rule for C = 1)
#             wsize PROF TH: X_L = 2^{TH L}, W = PROF(log n / log X_L), n <= X_L^2   (PROF = selberg | w14)
#   Evaluation at the zeros by type-3 NUFFT (FINUFFT), phases centred per block in long double.
# usage: python3 uk_at_zeros_v9.py zeross.bin outprefix count 1
#        python3 uk_at_zeros_v9.py zeross.bin outprefix wsize w14 0.7
# output: outprefix_uk.bin (float64, u_k for k = 1..K), outprefix_meta.txt
import numpy as np, finufft, sys, math, time
zf, pref, sch = sys.argv[1], sys.argv[2], sys.argv[3]
if sch == 'count': prof, par = 'none', float(sys.argv[4])
else: prof, par = sys.argv[4], float(sys.argv[5])
g = np.fromfile(zf); K = g.size
jmax = int(math.floor(math.log2(g[-1])))
def logX(L): return L * par * math.log(2.0) if L >= 1 else 0.0
def Wf(u):
    v = np.clip(u - 1.0, 0.0, 1.0)
    w = (1.0 - v) if prof == 'selberg' else (1.0 - v + 0.25 * v * (1.0 - v))
    return np.where(u >= 2.0, 0.0, w)
def nL(L):
    if sch == 'count': return int(math.floor(par * 2 ** (L / 2) + 1e-12)) if L >= 0 else 0
    lx = logX(L)
    return 0 if lx <= 0 else int(np.searchsorted(PR, math.exp(2 * lx), side='right'))
if sch == 'count': LIM = 2000000
else: LIM = int(math.exp(2 * logX(jmax + 1))) + 10
s = np.ones(LIM + 1, bool); s[:2] = False
for i in range(2, int(LIM ** .5) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0].astype(np.float64); del s
def points(i0, i1, L):
    p = PR[i0:i1]; lp = np.log(p.astype(np.longdouble)); om = []; cf = []; r = 1
    lx = logX(L) if sch != 'count' else 0.0
    while True:
        c = 1.0 / (r * np.exp((r / 2.0) * np.log(p))); m = c > 1e-15
        if not m.any(): break
        rl = r * lp[m]
        w = np.ones(int(m.sum())) if sch == 'count' else Wf((rl / lx).astype(np.float64))
        om.append(rl); cf.append(c[m] * w); r += 1
    return np.concatenate(om), np.concatenate(cf)
TWOPI = np.longdouble('6.2831853071795864769252867665590058')
def phiL(L, t):
    """Phi^(L)(t) at the points t (float64 array) -- type-3 NUFFT"""
    n = nL(L); out = np.zeros(t.size)
    if n == 0 or t.size == 0: return out
    tc = 0.5 * (t[0] + t[-1]); tau = t - tc
    for b0 in range(0, n, 3000000):
        om, cf = points(b0, min(b0 + 3000000, n), L)
        ph = np.mod(om * np.longdouble(tc), TWOPI).astype(np.float64)
        c = (cf * np.exp(1j * ph)).astype(np.complex128)
        f = finufft.nufft1d3(om.astype(np.float64), c, tau, isign=1, eps=1e-13)
        out += -f.imag / math.pi
    return out
uk = np.empty(K); tic = time.time()
for j in range(2, jmax + 1):
    lo, hi = 2.0 ** j, 2.0 ** (j + 1)
    i0, i1 = np.searchsorted(g, lo), np.searchsorted(g, hi)
    if i1 <= i0: continue
    t = g[i0:i1]; u = t / lo
    ph = (2 - u) * phiL(j, t) + (u - 1) * phiL(j + 1, t)
    hc = np.log(np.log(np.maximum(t, math.e ** math.e)))
    F = t / (2 * math.pi) * np.log(t / (2 * math.pi * math.e)) + 0.875 + np.clip(ph, -hc, hc)
    k = np.arange(i0 + 1, i1 + 1)
    uk[i0:i1] = k - 0.5 - F
    print(f"  j={j:2d} zeros {i1-i0:8d} primes {nL(j)}->{nL(j+1)}  max|u|={np.abs(uk[i0:i1]).max():.4f}  t={time.time()-tic:.0f}s", flush=True)
uk.tofile(pref + '_uk.bin')
with open(pref + '_meta.txt', 'w') as fo: fo.write(f"{K} {sch} {prof} {par}\n")
print("done", time.time() - tic)
