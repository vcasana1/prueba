# falpha_w_v7.py -- linea GUE con perfiles: V_W(theta) = E[(S - Phi^W_X)^2], X = Tc^theta, en t uniforme,
# frente a la prediccion heuristica  V_W = (1/2pi^2) int_0^inf (1 - W(a/theta))^2 F(a) a^-2 da
# con F(a) = min(a,1) (Montgomery) y con F = a (a<=1), Fp (a>1)  (meseta medida).
# usage: python3 falpha_w_v7.py T1 W profile th1,th2,... [Fp=0.87]
import numpy as np, finufft, sys, math
from scipy.integrate import quad
T1 = float(sys.argv[1]); WW = float(sys.argv[2]); prof = sys.argv[3]; TH = [float(v) for v in sys.argv[4].split(',')]
Fp = float(sys.argv[5]) if len(sys.argv) > 5 else 0.87
h = 1 / 32; g0 = int(T1 * 32); M = 1 << int(math.ceil(math.log2(WW * 32))); Tc = T1 + WW / 2
g = np.fromfile('zeross.bin')
t = (g0 + np.arange(M)) * h; x = t / (2 * np.pi)
S = np.searchsorted(g, t, side='right') - (x * (np.log(x) - 1) + 0.875)
def Wf(u):
    u = np.asarray(u, float)
    if prof == 'sharp': return (u <= 1).astype(float)
    v = np.clip(u - 1, 0, 1); w = 1 - v if prof == 'selberg' else 1 - v + 0.25 * v * (1 - v)
    return np.where(u >= 2, 0.0, w)
umax = 1.0 if prof == 'sharp' else 2.0
LIM = int(Tc ** (umax * max(TH))) + 10
s = np.ones(LIM + 1, bool); s[:2] = False
for i in range(2, int(LIM ** .5) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0].astype(np.float64); del s
TWOPI = np.longdouble('6.2831853071795864769252867665590058')
def phiW(th):
    lx = th * math.log(Tc); n = int(np.searchsorted(PR, math.exp(umax * lx), side='right')); acc = np.zeros(M)
    for b0 in range(0, n, 3000000):
        p = PR[b0:min(b0 + 3000000, n)]; lp = np.log(p.astype(np.longdouble)); om = []; cf = []; r = 1
        while True:
            c = 1.0 / (r * np.exp((r / 2.0) * np.log(p))); m = c > 1e-15
            if not m.any(): break
            rl = r * lp[m]; om.append(rl); cf.append(c[m] * Wf((rl / lx).astype(float))); r += 1
        om = np.concatenate(om); cf = np.concatenate(cf); keep = cf != 0; om = om[keep]; cf = cf[keep]
        ph = np.mod(om * np.longdouble(g0 + M // 2) * np.longdouble(h), TWOPI).astype(np.float64)
        f = finufft.nufft1d1(om.astype(np.float64) * h, cf * np.exp(1j * ph), M, isign=1, eps=1e-12)
        acc += -f.imag / math.pi
    return acc
def pred(th, Fhi):
    F = lambda a: a if a <= 1 else Fhi
    f = lambda a: (1 - float(Wf(a / th))) ** 2 * F(a) / a ** 2
    pts = sorted(set([th, umax * th, 1.0]))
    tot = 0; lo = 1e-9
    for b in pts + [50.0]:
        if b > lo: tot += quad(f, lo, b, limit=400)[0]; lo = b
    return (tot + Fhi / 50.0) / (2 * math.pi ** 2)
print(f"# {prof}  window [{T1:g},{T1+WW:g})  Tc={Tc:g}")
print("# theta   V_med      V_pred(F=Montgomery)  V_pred(F meseta %.2f)" % Fp)
for th in TH:
    e = S - phiW(th); v = float((e * e).mean())
    print(f"{th:6.3f}  {v:.6f}   {pred(th,1.0):.6f}             {pred(th,Fp):.6f}", flush=True)
