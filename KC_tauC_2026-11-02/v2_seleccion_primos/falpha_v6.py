# falpha_v6.py -- variance of S(t) - Phi_X(t) with a SHARP prime cutoff X = Tc^theta (all prime powers of p<=X),
# on a uniform grid t in [T1, T1+W) (step 1/32), S(t) = N(t) - N0(t) counted on zeross.bin.
# Heuristic (derived in the reply, not proved):  V(theta) = E[(S-Phi_X)^2] ~ (1/2pi^2) int_theta^inf F(a) a^-2 da,
# so  F_est(theta) = -2 pi^2 theta^2 V'(theta)  should reproduce Montgomery's F(alpha) (= alpha for alpha<=1).
# usage: python3 falpha_v6.py T1 W th1,th2,...
import numpy as np, finufft, sys, math, time
T1 = float(sys.argv[1]); W = float(sys.argv[2]); TH = [float(v) for v in sys.argv[3].split(',')]
h = 1 / 32; g0 = int(T1 * 32); M = 1 << int(math.ceil(math.log2(W * 32))); Tc = T1 + W / 2
g = np.fromfile('zeross.bin')
t = (g0 + np.arange(M)) * h; x = t / (2 * np.pi)
S = np.searchsorted(g, t, side='right') - (x * (np.log(x) - 1) + 0.875)
Xmax = Tc ** max(TH); LIM = int(Xmax) + 10
s = np.ones(LIM + 1, bool); s[:2] = False
for i in range(2, int(LIM ** .5) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0].astype(np.float64); del s
TWOPI = np.longdouble('6.2831853071795864769252867665590058')
def add(i0, i1, acc):
    PB = 3000000
    for b0 in range(i0, i1, PB):
        p = PR[b0:min(b0 + PB, i1)]; lp = np.log(p.astype(np.longdouble)); om = []; cf = []; r = 1
        while True:
            c = 1.0 / (r * np.exp((r / 2.0) * np.log(p))); m = c > 1e-15
            if not m.any(): break
            om.append(r * lp[m]); cf.append(c[m]); r += 1
        om = np.concatenate(om); cf = np.concatenate(cf)
        gc = g0 + M // 2
        ph = np.mod(om * np.longdouble(gc) * np.longdouble(h), TWOPI).astype(np.float64)
        o64 = om.astype(np.float64)
        f = finufft.nufft1d1(o64 * h, cf * np.exp(1j * ph), M, isign=1, eps=1e-12)
        acc += -f.imag / math.pi
phi = np.zeros(M); prev = 0; out = []
tic = time.time()
print(f"# window [{T1:g},{T1+W:g}) grid {M} pts, Tc={Tc:g}, primes up to {Xmax:.3g} ({PR.size})")
print("# theta    X          #primes     E[(S-Phi)^2]   E[(S-Phi)^4]")
for th in TH:
    n = int(np.searchsorted(PR, Tc ** th, side='right'))
    add(prev, n, phi); prev = n
    e = S - phi; v2 = float((e * e).mean()); v4 = float((e ** 4).mean()); out.append((th, v2, v4))
    print(f"{th:6.3f} {Tc**th:12.4g} {n:10d}   {v2:.6f}      {v4:.6f}   ({time.time()-tic:.0f}s)", flush=True)
print("# F_est(theta_mid) = -2 pi^2 theta_mid^2 dV/dtheta   vs  min(theta,1)")
for (a, va, _), (b, vb, _) in zip(out[:-1], out[1:]):
    tm = 0.5 * (a + b); print(f"  theta={tm:.3f}   F_est={-2*math.pi**2*tm*tm*(vb-va)/(b-a):.3f}   min(theta,1)={min(tm,1):.3f}")
