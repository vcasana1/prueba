# mv_dphi_v9.py -- Montgomery-Vaughan (unconditional) value of mean_t Phi'^2 on (2^j, 2^{j+1}] for W14, X_L = 2^{TH L}
#   Phi^o' ~ -(1/pi) sum_n a_n cos(t log n),  a_n = c_n log n * Wmix(n)   (c_n log n = Lambda(n)/sqrt n for n = p^r)
#   MV: int_T^{2T} |sum b_n n^{-it}|^2 = sum |b_n|^2 (T + theta 3pi n/... ) -> diagonal + explicit error (Preissmann constant)
import numpy as np, math, sys
TH, j = float(sys.argv[1]), int(sys.argv[2])
lx0, lx1 = j * TH * math.log(2), (j + 1) * TH * math.log(2); LIM = int(math.exp(2 * lx1)) + 10
s = np.ones(LIM + 1, bool); s[:2] = False
for i in range(2, int(LIM**.5) + 1):
    if s[i]: s[i*i::i] = False
p = np.nonzero(s)[0].astype(float); del s
def W(u):
    v = np.clip(u - 1, 0, 1); return np.where(u >= 2, 0, 1 - v + 0.25 * v * (1 - v))
D = 0.0; E = 0.0
for r in range(1, 40):
    n = p**r; m = n <= LIM
    if not m.any(): break
    ln = np.log(n[m]); a = np.log(p[m]) / np.sqrt(n[m]) * 0.5 * (W(ln / lx0) + W(ln / lx1)) / math.pi
    D += (a**2).sum(); E += (n[m] * a**2).sum()
T = 2.0**j; T2 = 2 * T
N0p = (T2 * math.log(T2 / (2*math.pi*math.e)) - T * math.log(T / (2*math.pi*math.e))) / (2 * math.pi) / T   # mean N0' on block
diag = 0.5 * D                                       # mean of cos^2 = 1/2
err = 0.5 * math.sqrt(1 + 2/3 * math.sqrt(6/5)) * math.pi * E / T   # Preissmann form of the MV constant, divided by block length
print(f"TH={TH} block 2^{j}: MV diagonal mean Phi'^2 = {diag:.4f};  explicit MV error <= {err:.4f};  mean N0' = {N0p:.4f}")
print(f"  RMS Phi'/N0': diagonal {math.sqrt(diag)/N0p:.4f},  upper bound (diag+err) {math.sqrt(diag+err)/N0p:.4f}")
