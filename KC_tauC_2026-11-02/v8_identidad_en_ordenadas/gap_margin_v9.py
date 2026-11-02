# gap_margin_v9.py -- margin delta_k = 1/2 - |u_k| against the gap to the nearer neighbour (mean spacings), close pairs only
# usage: python3 gap_margin_v9.py zeross.bin prefix [...]
import numpy as np, sys, math
g = np.fromfile(sys.argv[1]); lam = np.log(g / (2 * math.pi)) / (2 * math.pi)
gl = np.r_[np.inf, np.diff(g)] * lam; gr = np.r_[np.diff(g), np.inf] * lam; mg = np.minimum(gl, gr)
for pref in sys.argv[2:]:
    u = np.fromfile(pref + '_uk.bin'); d = 0.5 - np.abs(u)
    print(f"== {pref}")
    for lo, hi in [(0, .01), (.01, .02), (.02, .04), (.04, .08), (.08, .16), (.16, .32), (.32, 10)]:
        m = (mg >= lo) & (mg < hi)
        print(f"   gap in [{lo:.2f},{hi:.2f}): n={m.sum():8d}  delta: min {d[m].min():+.4f}  median {np.median(d[m]):.4f}   #delta<=0: {(d[m]<=0).sum()}")
    m = mg < 0.05; A = np.column_stack([mg[m], np.ones(m.sum())]); c, *_ = np.linalg.lstsq(A, d[m], rcond=None)
    q = np.quantile(d[m] - c[0] * mg[m], 0.01)
    print(f"   close pairs (gap<0.05, n={m.sum()}): delta ~ {c[0]:.2f}*gap {c[1]:+.4f};  1% lower envelope: delta >= {c[0]:.2f}*gap {q:+.4f}")
