# nufft_phi_v6.py -- (v6: primes processed in batches of PB to bound memory; terms with coefficient <= 1e-15 dropped)
# nufft_phi_v5.py -- Phi^o(t) and Phi^o'(t) on the grid t = g*h (h = 1/HINV), g = 0..gmax, via type-1 NUFFT (FINUFFT).
#   Phi^(L)(t) = -(1/pi) sum_{p in P_L} sum_{r>=1} sin(r t log p)/(r p^{r/2})          (exact Euler phases, all prime powers)
#   P_L = first n_L primes;  scheme 'count c': n_L = floor(c 2^{L/2});  scheme 'size th': n_L = pi(2^{L th})   (L<0: n_L = 0)
#   on [2^j, 2^{j+1}]:  Phi^o = (2-u) Phi^(j) + (u-1) Phi^(j+1),  u = t/2^j      (IJNT paper eq. F_def for 'count 1')
# usage: python3 nufft_phi_v5.py zeross.bin outprefix SCHEME PARAM [HINV=32]
# output: outprefix_phi.bin, outprefix_dphi.bin (float64, g = 0..gmax; dphi = derivative of the block containing [g, g+1]),
#         outprefix_bnd.txt  (g at block ends and the left derivative there), outprefix_meta.txt
import numpy as np, finufft, sys, math, time
zf, pref, sch, par = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
HINV = int(sys.argv[5]) if len(sys.argv) > 5 else 32; h = 1.0 / HINV
K = np.fromfile(zf).size
# Tmax: N0(Tmax) = K + 12  (same as tauC_v1)
lam = K + 12.0; T = 2 * math.pi * lam / math.log(lam) + 50
for _ in range(100):
    x = T / (2 * math.pi); T -= (x * (math.log(x) - 1) + 0.875 - lam) / (math.log(x) / (2 * math.pi))
gmax = int(math.ceil(T / h)); jmax = int(math.floor(math.log2(gmax * h)))
def nL(L):
    if L < 0: return 0
    if sch == 'count': return int(math.floor(par * 2 ** (L / 2) + 1e-12))
    X = 2 ** (L * par); return int(np.searchsorted(PR, X, side='right'))
# primes
LIM = int(2 ** ((jmax + 1) * par)) + 10 if sch == 'size' else int(1.5 * (par * 2 ** ((jmax + 1) / 2)) * 20 + 1000)
s = np.ones(LIM + 1, bool); s[:2] = False
for i in range(2, int(LIM ** .5) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0].astype(np.float64); del s
print(f"K={K} Tmax={T:.6f} gmax={gmax} jmax={jmax} scheme={sch} {par} primes sieved={PR.size} n_jmax+1={nL(jmax+1)}", flush=True)
# frequency/coefficient list for primes with index in [i0, i1): all prime powers with coefficient > 1e-17
def points(i0, i1):
    p = PR[i0:i1]; lp = np.log(p.astype(np.longdouble)); om = []; cf = []
    r = 1
    while True:
        c = 1.0 / (r * np.exp((r / 2.0) * np.log(p))); m = c > 1e-15
        if not m.any(): break
        om.append(r * lp[m]); cf.append(c[m]); r += 1
    return (np.concatenate(om) if om else np.zeros(0, np.longdouble)), (np.concatenate(cf) if cf else np.zeros(0))
EPS = 1e-13
def evalsum(om, cf, g0, M):
    """ returns A(t), A'(t) at t=(g0+m)h, m=0..M-1, where A = -(1/pi) sum cf sin(om t) """
    if om.size == 0: return np.zeros(M), np.zeros(M)
    gc = g0 + M // 2
    TWOPI = np.longdouble('6.2831853071795864769252867665590058')
    ph = np.mod(om * np.longdouble(gc) * np.longdouble(h), TWOPI).astype(np.float64)
    o64 = om.astype(np.float64)
    c0 = cf * np.exp(1j * ph); cv = np.stack([c0, c0 * o64])
    f = finufft.nufft1d1(o64 * h, cv, M, isign=1, eps=EPS)          # modes k = -M/2 .. M/2-1  <->  m = k + M/2
    return -f[0].imag / math.pi, -f[1].real / math.pi
phi = np.zeros(gmax + 1); dphi = np.zeros(gmax + 1); bnd = []
MB = 1 << 21
tic = time.time()
PB = 3000000
for j in range(-1, jmax + 1):
    ga = int(round(2.0 ** j * HINV)); gb = min(2 * ga, gmax)
    if ga >= gmax: break
    n1, n2 = nL(j), nL(j + 1)
    pj = 2.0 ** j; L = gb + 1 - ga
    A = np.zeros(L); dA = np.zeros(L); B = np.zeros(L); dB = np.zeros(L); npts = 0
    for (i0, i1, VA, VD) in [(0, n1, A, dA), (n1, n2, B, dB)]:
        for b0 in range(i0, i1, PB):
            om, cf = points(b0, min(b0 + PB, i1)); npts += om.size
            for g0 in range(ga, gb + 1, MB):
                M = min(MB, gb + 1 - g0); M2 = max(64, 1 << (M - 1).bit_length())
                v, d = evalsum(om, cf, g0, M2)
                VA[g0 - ga:g0 - ga + M] += v[:M]; VD[g0 - ga:g0 - ga + M] += d[:M]
            del om, cf
    t = (ga + np.arange(L)) * h; u = t / pj
    val = A + (u - 1) * B; der = B / pj + dA + (u - 1) * dB
    phi[ga:gb + 1] = val
    if gb == 2 * ga:
        bnd.append((gb, der[-1])); dphi[ga:gb] = der[:-1]
    else:
        dphi[ga:gb + 1] = der
    del A, dA, B, dB, val, der, t, u
    print(f"  block j={j:3d}  primes {n1}->{n2}  points {npts}  t={time.time() - tic:.0f}s", flush=True)
phi.tofile(pref + '_phi.bin'); dphi.tofile(pref + '_dphi.bin')
with open(pref + '_bnd.txt', 'w') as fo:
    for gbb, d in bnd: fo.write(f"{gbb} {d:.17g}\n")
with open(pref + '_meta.txt', 'w') as fo: fo.write(f"{K} {gmax} {HINV} {sch} {par}\n")
print("done", time.time() - tic)
