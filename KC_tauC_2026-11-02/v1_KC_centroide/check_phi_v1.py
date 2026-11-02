# independent check of F, F' from tauC_v1 'dump' against direct numpy evaluation of the paper's definition
import numpy as np, subprocess, math
def primes(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
P=primes(20000).astype(np.longdouble)
def phi_m(t,m):
    p=P[:m]; z=np.exp(-1j*(np.longdouble(t)*np.log(p)))/np.sqrt(p)
    return float(-np.sum(np.angle(1-z.astype(np.clongdouble)))/np.pi)
def F(t):
    j=int(math.floor(math.log2(t))); mj=lambda j: 0 if j<0 else int(math.floor(2**(j/2)))
    u=t/2**j
    po=(2-u)*phi_m(t,mj(j))+(u-1)*phi_m(t,mj(j+1))
    hc=1.0 if t<=math.exp(math.e) else math.log(math.log(t))
    ph=max(-hc,min(hc,po))
    x=t/(2*math.pi)
    return x*(math.log(x)-1)+0.875+ph, ph
for g0 in [3200, 6400*5+13, 29000000, 33000000]:
    out=subprocess.run(['./tauC_v1','dump',str(g0),'3'],capture_output=True,text=True).stdout.split('\n')[:4]
    for ln in out:
        t,Fc,Dc=map(float,ln.split()); Fr,ph=F(t)
        h=1e-5; Dr=(F(t+h)[0]-F(t-h)[0])/(2*h)
        print(f"t={t:.5f}  F_C={Fc:.9f} F_py={Fr:.9f} dF={Fc-Fr:+.2e} | F'_C={Dc:.8f} F'_py(fd)={Dr:.8f} d={Dc-Dr:+.1e}")
