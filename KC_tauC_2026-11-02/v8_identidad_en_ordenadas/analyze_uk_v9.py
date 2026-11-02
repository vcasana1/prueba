# analyze_uk_v9.py -- coincidence test from the exact defect identity (IJNT v55, Thm identity) with F = taucs function.
#   u_k = k - 1/2 - F(gamma_k).  Psi_j(u) = ((1/2+u)^{2j+1} + (1/2-u)^{2j+1})/(2j+1).
#   E_j(window) = sum_k [Psi_j(u_k) - Psi_j(0)],  c_j = 2(1-4^{-j})/(2j+1).
#   Every coincidence gamma_k = gamma_{k+1} inside the window forces E_j >= c_j, for every j.
#   So E_j < c_j for some j certifies the window (no off-line zero, no multiple zero).
#   j = 1: #zeros not simple on the line <= 4 sum u_k^2   (proportion bound 1 - 4 eta, eta = mean u^2).
#   j*   : least j with E_j < c_j (exists iff max|u| < 1/2).
# usage: python3 analyze_uk_v9.py zeross.bin prefix [W=1000]
import numpy as np, sys, math
g = np.fromfile(sys.argv[1]); pref = sys.argv[2]; W = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
u = np.fromfile(pref + '_uk.bin'); K = u.size; k = np.arange(1, K + 1)
lam = np.log(g / (2 * math.pi)) / (2 * math.pi)
a = np.abs(u)
print(f"== {pref}:  K={K}  RMS(u)={np.sqrt((u**2).mean()):.5f}  mean u={u.mean():+.5f}  #|u|>=1/2: {(a>=0.5).sum()}  max|u|={a.max():.4f}")
# ---- j = 1 on dyadic blocks
print(" dyadic blocks (T, 2T]:  n   eta=mean u^2   1-4eta   max|u|   #|u|>=1/2")
for j in range(10, 21):
    m = (g > 2.0**j) & (g <= 2.0**(j + 1))
    if m.sum() == 0: continue
    eta = (u[m]**2).mean()
    print(f"   2^{j:2d}  {m.sum():8d}  {eta:.5f}   {1-4*eta:.4f}   {a[m].max():.4f}   {(a[m]>=0.5).sum()}")
# ---- j* per window
def logsum_terms(uw, j):
    """log of sum_k [(1/2+u)^{2j+1} + (1/2-u)^{2j+1}] for |u|<1/2"""
    e = 2 * j + 1
    l1 = e * np.log(0.5 + np.abs(uw)); l2 = e * np.log(0.5 - np.abs(uw))
    L = np.concatenate([l1, l2]); M = L.max()
    return M + math.log(np.exp(L - M).sum())
def certified(uw, j):
    n = uw.size; e = 2 * j + 1
    ls = logsum_terms(uw, j)                     # log sum
    base = math.log(n) - 2 * j * math.log(2.0)   # log(n 2^{-2j})
    lhs = math.exp(ls) - math.exp(base) if ls < 700 else float('inf')
    return lhs < 2 * (1 - 4.0**(-j))
def jstar(uw):
    if np.abs(uw).max() >= 0.5: return None
    if certified(uw, 1): return 1
    lo, hi = 1, 2
    while not certified(uw, hi): lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if certified(uw, mid): hi = mid
        else: lo = mid
    return hi
nw = K // W; js = []; dl = []; bad = 0
for w in range(nw):
    uw = u[w*W:(w+1)*W]; d = 0.5 - np.abs(uw).max(); dl.append(d)
    jj = jstar(uw)
    if jj is None: bad += 1; js.append(np.nan)
    else: js.append(jj)
js = np.array(js, float); dl = np.array(dl)
ok = ~np.isnan(js)
print(f" windows of {W} zeros: {nw}   certified (max|u|<1/2): {ok.sum()}  ({100*ok.mean():.2f}%)   not certifiable: {bad}")
if ok.any():
    q = np.percentile(js[ok], [0, 50, 90, 99, 100]); qd = np.percentile(dl[ok], [0, 1, 10, 50])
    print(f"   j*: min {q[0]:.0f}  median {q[1]:.0f}  p90 {q[2]:.0f}  p99 {q[3]:.0f}  max {q[4]:.0f}")
    print(f"   delta = 1/2 - max|u|: min {qd[0]:.4f}  p1 {qd[1]:.4f}  p10 {qd[2]:.4f}  median {qd[3]:.4f}")
    r = js[ok] / (math.log(W) / (2 * dl[ok]))
    print(f"   j* / (log W / (2 delta)):  median {np.median(r):.3f}  range [{r.min():.3f}, {r.max():.3f}]")
# by height: fraction of certifiable windows per dyadic block
print(" certifiable windows by height:")
wc = g[np.arange(nw) * W + W // 2]
for j in range(12, 21):
    m = (wc > 2.0**j) & (wc <= 2.0**(j + 1))
    if m.sum(): print(f"   2^{j:2d}: {m.sum():5d} windows, certified {100*ok[m].mean():6.2f}%,  median j* {np.nanmedian(js[m]):.0f}")
# worst points
o = np.argsort(-a)[:15]
print(" largest |u_k| (k, gamma_k, u_k, left gap, right gap in mean spacings):")
for i in o:
    gl = (g[i] - g[i-1]) * lam[i] if i > 0 else np.nan; gr = (g[i+1] - g[i]) * lam[i] if i + 1 < K else np.nan
    print(f"   k={i+1:8d}  gamma={g[i]:14.6f}  u={u[i]:+.4f}  gaps {gl:.3f} {gr:.3f}")
if (a >= 0.5).sum():
    i = np.nonzero(a >= 0.5)[0]
    gl = (g[i] - g[i-1]) * lam[i]; gr = (g[np.minimum(i+1, K-1)] - g[i]) * lam[i]
    mg = np.minimum(gl, gr); allg = np.minimum(np.diff(g)[:-1], np.diff(g)[1:]) * lam[1:-1]
    print(f" at |u|>=1/2: median of min neighbour gap {np.median(mg):.3f} (all zeros: {np.median(allg):.3f});  "
          f"sign of u: +{(u[i]>0).sum()} / -{(u[i]<0).sum()}")
np.savetxt(pref + f'_windows_W{W}.txt', np.column_stack([np.arange(nw) * W + 1, wc, dl, js]),
           fmt=['%d', '%.6f', '%.6f', '%.0f'], header='k_first  gamma_center  delta=1/2-max|u|  jstar(nan=not certifiable)')
