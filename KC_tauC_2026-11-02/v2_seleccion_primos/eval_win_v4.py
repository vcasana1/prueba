# eval_win_v4.py : RMS(tau^C - gamma) in windows, per scheme, in t-units and in mean spacings
import numpy as np, glob, re
g=np.fromfile('zeross.bin')
rows=[]
for f in sorted(glob.glob('w*_win.txt')):
    w=np.loadtxt(f); k=w[:,0].astype(int); tc=w[:,1]; t0=w[:,2]; gk=g[k-1]
    sp=2*np.pi/np.log(gk/(2*np.pi))
    dC=(tc-gk)/sp; d0=(t0-gk)/sp
    rows.append((f[:-8], len(k), np.sqrt((dC**2).mean()), np.sqrt((d0**2).mean()), np.abs(dC).max(), (np.abs(dC)>0.5).mean(), (dC**4).mean()))
print(f"{'window_scheme':24s} {'#k':>5s} {'rms(tC-g)/sp':>12s} {'rms(tau-g)/sp':>13s} {'max|.|/sp':>9s} {'P(|.|>1/2)':>10s} {'E[d^4]':>8s}")
for r in rows: print(f"{r[0]:24s} {r[1]:5d} {r[2]:12.5f} {r[3]:13.5f} {r[4]:9.4f} {r[5]:10.5f} {r[6]:8.5f}")
