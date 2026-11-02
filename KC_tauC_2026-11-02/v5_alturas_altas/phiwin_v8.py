# phiwin_v8.py -- Phi^o(t), Phi^o'(t) del draft (m_j = floor(2^{j/2}), interpolacion diadica) en t = T0 + j h, j=0..M-1,
# por NUFFT tipo 1. Requiere que la ventana este dentro de un bloque diadico [2^J, 2^{J+1}].
# usage: python3 phiwin_v8.py T0 W h_inv outprefix
import numpy as np, finufft, sys, math
T0 = float(sys.argv[1]); W = float(sys.argv[2]); HI = int(sys.argv[3]); pref = sys.argv[4]; h = 1.0 / HI
g0 = int(round(T0 * HI)); assert abs(g0 / HI - T0) < 1e-9
M = int(math.ceil(W * HI)) + 1
J = int(math.floor(math.log2(T0))); assert T0 + W <= 2.0 ** (J + 1), "ventana cruza un bloque diadico"
m1 = int(math.floor(2 ** (J / 2) + 1e-12)); m2 = int(math.floor(2 ** ((J + 1) / 2) + 1e-12))
LIM = int(m2 * (math.log(m2) + math.log(math.log(m2)))) + 1000
s = np.ones(LIM + 1, bool); s[:2] = False
for i in range(2, int(LIM ** .5) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0].astype(np.float64)[:m2]; del s
assert PR.size == m2
TWOPI = np.longdouble('6.2831853071795864769252867665590058')
M2 = 1 << (M - 1).bit_length(); Mc = M2 // 2
def phase_sum(i0, i1):
    p = PR[i0:i1]; lp = np.log(p.astype(np.longdouble)); om = []; cf = []; r = 1
    while True:
        c = 1.0 / (r * np.exp((r / 2.0) * np.log(p))); m = c > 1e-15
        if not m.any(): break
        om.append(r * lp[m]); cf.append(c[m]); r += 1
    om = np.concatenate(om); cf = np.concatenate(cf)
    ph = np.mod(om * np.longdouble(g0 + Mc) * np.longdouble(h), TWOPI).astype(np.float64)
    o = om.astype(np.float64); c0 = cf * np.exp(1j * ph)
    f = finufft.nufft1d1(o * h, np.stack([c0, c0 * 1j * o]), M2, isign=1, eps=1e-13)
    return -f[0][:M].imag / math.pi, -f[1][:M].imag / math.pi   # A = -(1/pi) sum c sin(om t), A' = -(1/pi) sum c om cos(om t)
A, dA = phase_sum(0, m1); B, dB = phase_sum(m1, m2)                 # B = Phi_{m2} - Phi_{m1}
t = T0 + np.arange(M) * h; pj = 2.0 ** J; u = t / pj
phi = A + (u - 1) * B; dphi = B / pj + dA + (u - 1) * dB
phi.tofile(pref + '_phi.bin'); dphi.tofile(pref + '_dphi.bin')
print(f"T0={T0} W={W} M={M} bloque J={J} primos m_J={m1} m_J+1={m2}  |Phi| max={np.abs(phi).max():.3f}  loglog T0={math.log(math.log(T0)):.3f}")
