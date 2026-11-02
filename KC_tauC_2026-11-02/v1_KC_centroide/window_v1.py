# windowed difference D_phi(N) = sum_k phi(k/N) (e^{i gamma_k} - e^{i x_k}),  phi(x)=exp(4-1/(x(1-x))) in C_c^inf(0,1)
import numpy as np, sys
g=np.fromfile('zeross.bin'); K=len(g)
sets={'tauC':np.fromfile('run_h32_tauC.bin'),'tau':np.fromfile('run_h32_tau.bin'),
      'tauC_M25':np.fromfile('cap25_tauC.bin'),'tauC_M168':np.fromfile('cap168_tauC.bin')}
eg=np.exp(1j*g.astype(np.longdouble))
def phi(x):
    x=np.asarray(x,dtype=np.longdouble); out=np.zeros_like(x); m=(x>0)&(x<1)
    out[m]=np.exp(4-1/(x[m]*(1-x[m]))); return out
print("N, then |D_phi| for: "+", ".join(sets))
res={}
for N in [10**3,3*10**3,10**4,3*10**4,10**5,3*10**5,10**6,K]:
    w=phi(np.arange(1,N+1)/N); row=[]
    for nm,x in sets.items():
        d=(w*(eg[:N]-np.exp(1j*x[:N].astype(np.longdouble)))).sum(); row.append(complex(d))
    res[N]=row
    print(f"{N:8d} "+"  ".join(f"{r.real:+8.4f}{r.imag:+8.4f}i |{abs(r):7.4f}|" for r in row))
