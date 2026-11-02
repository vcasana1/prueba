# error_por_gap_v8.py : RMS(tauc-gamma) en espaciados, por intervalos del hueco minimo vecino (en espaciados)
import numpy as np, mpmath as mp
bins=[0,0.2,0.4,0.6,0.8,1.0,1.5,10]
def tabla(nombre,g,x,sp):
    gl=np.r_[np.inf,np.diff(g)]; gr=np.r_[np.diff(g),np.inf]; gm=np.minimum(gl,gr)/sp; e=(x-g)/sp
    gm,e=gm[1:-1],e[1:-1]
    print(f"{nombre}: N={e.size}  RMS total={np.sqrt((e**2).mean()):.4f}")
    for a,b in zip(bins[:-1],bins[1:]):
        m=(gm>=a)&(gm<b)
        if m.sum(): print(f"   hueco min en [{a:.1f},{b:.1f}): {m.sum():7d} ({100*m.mean():5.1f} %)  RMS={np.sqrt((e[m]**2).mean()):.4f}")
# altura ~1e6: k en [1.5e6, 2e6] de tus ceros
g=np.fromfile('/home/claude/KC_v2/zeross.bin'); x=np.fromfile('/home/claude/KC_v1/zTaucs2e6.bin')
s=slice(1500000,2001052); gg=g[s]; sp=2*np.pi/np.log(gg/(2*np.pi)); tabla("t ~ 1e6 (k 1.5e6-2e6)",gg,x[s],sp)
# altura 2.7e11
d=np.load('z_h1e12_zeros.npz'); T0=float(mp.mpf(str(d['T0']))); r=d['roots']; n0=int(d['K0'])+int(d['n0rel'])
tc=np.loadtxt('tc_h1e12.txt'); k=tc[:,0].astype(np.int64); xx=tc[:,1]; i=k-(n0+1); ok=(i>=0)&(i<r.size); i,xx=i[ok],xx[ok]
gg=r[i]; sp=2*np.pi/np.log((T0+gg)/(2*np.pi)); tabla("t ~ 2.7e11",gg,xx,sp)
