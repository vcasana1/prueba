import numpy as np, sys
from rho import model
T=float(sys.argv[1]); h=0.002
om,c=model(T)
tmax=2*T+60
M=int(tmax/h)
F=np.empty(M)
ch=2_000_000
for s in range(0,M,ch):
    t=h*(np.arange(s,min(s+ch,M))+0.5)
    Phi=np.zeros_like(t)
    for o,cc in zip(om,c): Phi-=cc*np.sin(o*t)
    F[s:s+len(t)]=t/(2*np.pi)*np.log(t/(2*np.pi*np.e))+7/8+Phi
Fs=np.sort(F)
N0=lambda t: t/(2*np.pi)*np.log(t/(2*np.pi*np.e))+7/8
# RvM tau_k and tauC_k for k with tau in (T,2T]
ks=np.arange(int(np.ceil(N0(T)+.5)),int(np.floor(N0(2*T)+.5))+1)
cs=ks-0.5
tauC=h*np.searchsorted(Fs,cs)
# tau by bisection
lo=np.full(len(cs),1.0);hi=np.full(len(cs),3*T)
for _ in range(80):
    mid=(lo+hi)/2; m=N0(mid)<cs; lo=np.where(m,mid,lo); hi=np.where(m,hi,mid)
tau=(lo+hi)/2
sel=(tau>T)&(tau<=2*T); n=sel.sum()
print("T=%g n=%d  mean(tauC-tau)=%.4f rms=%.4f max=%.3f"%(T,n,np.mean((tauC-tau)[sel]),np.sqrt(np.mean((tauC-tau)[sel]**2)),np.max(abs(tauC-tau)[sel])))
def Lam(x):
    for p in [2,3,5,7,11,13]:
        k=0;y=x
        while abs(y-round(y))<1e-9 and round(y)%p==0: y=round(y)/p;k+=1
        if k>0 and abs(y-1)<1e-9: return np.log(p)
    return 0.0
for name,x in [("log2",2),("log3",3),("log4",4),("log5",5),("log6",6),("log7",7),("log8",8),("log9",9),("log10",10)]:
    u=np.log(x)
    sT=np.mean(np.exp(1j*u*tau[sel])); sC=np.mean(np.exp(1j*u*tauC[sel]))
    pred=-(T/(2*np.pi*n))*Lam(x)*x**-0.5
    print("%-6s Lam=%.3f pred=%+.5f | tau: %+.5f%+.5fj | tauC: %+.5f%+.5fj"%(name,Lam(x),pred,sT.real,sT.imag,sC.real,sC.imag))
for u in [0.5,1.0,2.0]:
    sT=np.mean(np.exp(1j*u*tau[sel])); sC=np.mean(np.exp(1j*u*tauC[sel]))
    print("u=%.1f  tau: %+.5f%+.5fj | tauC: %+.5f%+.5fj"%(u,sT.real,sT.imag,sC.real,sC.imag))
