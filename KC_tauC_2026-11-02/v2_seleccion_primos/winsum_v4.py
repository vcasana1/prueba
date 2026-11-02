# increment of B - B^C over each window: sum_{k in window} (e^{i gamma_k} - e^{i x_k}), and smooth-window version
import numpy as np, glob
g=np.fromfile('zeross.bin')
print(f"{'window_scheme':24s} {'|sum(e^ig - e^itC)|':>20s} {'|smooth sum|':>13s}   (tau: plain / smooth)")
for f in sorted(glob.glob('w1*_win.txt')):
    w=np.loadtxt(f); k=w[:,0].astype(int); tc=w[:,1]; t0=w[:,2]; gk=g[k-1]
    x=np.linspace(0,1,len(k)+2)[1:-1]; phi=np.exp(4-1/(x*(1-x)))
    a=np.exp(1j*gk)-np.exp(1j*tc); b=np.exp(1j*gk)-np.exp(1j*t0)
    print(f"{f[:-8]:24s} {abs(a.sum()):20.4f} {abs((phi*a).sum()):13.4f}   {abs(b.sum()):.4f} / {abs((phi*b).sum()):.4f}")
