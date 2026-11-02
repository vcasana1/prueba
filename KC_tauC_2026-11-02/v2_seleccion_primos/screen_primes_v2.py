# screen_primes_v2.py -- criba de reglas de seleccion/peso de primos por el error e = Phi_scheme(t) - (N(t)-N0(t))
# en puntos t uniformes de una ventana [t0, t0+W]; N(t) contado sobre zeross.bin (varpi/pi = O(1/t) despreciado).
import numpy as np, sys, math, time
t0=float(sys.argv[1]); W=float(sys.argv[2]); n=int(sys.argv[3]); PMAX=int(sys.argv[4]); seed=int(sys.argv[5]) if len(sys.argv)>5 else 1
g=np.fromfile('zeross.bin')
rng=np.random.default_rng(seed); t=np.sort(t0+W*rng.random(n))
N=np.searchsorted(g,t,side='right').astype(float)
x=t/(2*np.pi); R=N-(x*(np.log(x)-1)+0.875)          # = S(t) + O(1/t)
def primes(L):
    s=np.ones(L+1,bool); s[:2]=False
    for i in range(2,int(L**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
P=primes(PMAX); lt=np.log(t)
j=np.floor(np.log2(t)); u=t/2**j
mj=np.floor(2**(j/2)+1e-12); mj1=np.floor(2**((j+1)/2)+1e-12)
schemes={'none':None,'paper':None}
import os
CNT=[0.5,1,1.5,2,3,4,8]; SIZE=[float(v) for v in os.environ.get("SIZES","0.3,0.4,0.5,0.6,0.7,0.8,0.9").split(",")]; SEL=[0.2,0.25,0.3,0.35,0.4]
acc={k:np.zeros(n) for k in ['paper']+[f'cnt{c}' for c in CNT]+[f'size{q}' for q in SIZE]+[f'selb{q}' for q in SEL]+[f'selbd{q}' for q in SEL]}
nb=60; cb=np.zeros(nb); gb=np.zeros(nb); cpl=[]; gpl=[]
xs={q:np.exp(q*lt) for q in SEL}
tic=time.time()
for i,p in enumerate(P):
    lg=math.log(p); a=p**-0.5; th=t*lg
    s1=np.sin(th); c1=np.cos(th)
    phi=-np.arctan2(a*s1,1-a*c1)/np.pi                   # exact Euler phase of p
    # paper: weight 1 if i<mj ; (u-1) if mj<=i<mj1 ; 0 otherwise
    wpap=np.where(i<mj,1.0,np.where(i<mj1,u-1,0.0)); acc['paper']+=wpap*phi
    for c in CNT: acc[f'cnt{c}']+=np.where(i<np.floor(c*np.sqrt(t)),phi,0.0)
    for q in SIZE: acc[f'size{q}']+=np.where(lg<=q*lt,phi,0.0)
    # Selberg-type taper on prime powers n=p^r : w=1 (n<=x), log(x^2/n)/log x (x<n<=x^2), 0 beyond
    if p<=xs[SEL[-1]].max()**2:
        rmax=int(math.floor(2*SEL[-1]*lt.max()/lg))+1
        sr_prev=np.zeros(n); sr=s1; cr=c1
        for r in range(1,rmax+1):
            term=-(a**r)*sr/(r*np.pi); ln=r*lg
            for q in SEL:
                lx=q*lt; w=np.clip(2-ln/lx,0,1)
                acc[f'selb{q}']+=w*term
                acc[f'selbd{q}']+=w*np.exp(-ln/lx)*term
            sr,cr=sr*c1+cr*s1, cr*c1-sr*s1               # sin((r+1)th), cos((r+1)th)
    # empirical least-squares weight of the exact phase, binned in u=log p/log t0
    cpv=float((R*phi).sum()); gpv=float((phi*phi).sum()); cpl.append(cpv); gpl.append(gpv)
    b=min(int(lg/math.log(t0)/0.02),nb-1); cb[b]+=cpv; gb[b]+=gpv
print(f"# t0={t0:g} W={W:g} n={n} PMAX={PMAX} primes={len(P)}  time={time.time()-tic:.0f}s")
def st(e): return float((e*e).mean()), float((e**4).mean())
v2,v4=st(R); print(f"{'scheme':10s} {'E[e^2]':>9s} {'E[e^4]':>9s}  (ratio to none)")
print(f"{'none':10s} {v2:9.5f} {v4:9.5f}")
for k,v in acc.items():
    a2,a4=st(v-R); print(f"{k:10s} {a2:9.5f} {a4:9.5f}   {a2/v2:6.3f} {a4/v4:6.3f}")
cpl=np.array(cpl); gpl=np.array(gpl)
print("LS in-sample floor (diag. approx):  E[e^2] ~",(v2*n-(cpl**2/gpl).sum())/n)
print("empirical weight profile w(u), u=log p/log t0 (bin width 0.02):")
for b in range(nb):
    if gb[b]>0: print(f"  u in [{0.02*b:.2f},{0.02*(b+1):.2f}): w={cb[b]/gb[b]:+.3f}")
