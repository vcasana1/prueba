import numpy as np, sys
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
def model(T):
    X=np.sqrt(1.5*T)
    ps=primes_upto(int(X)+1)
    terms=[]
    for p in ps:
        k=1;n=p
        while n<=X:
            if n<=X**0.5: w=1.0
            else: w=np.log(X/n)/np.log(X**0.5)
            terms.append((k*np.log(p), w/(k*p**(k/2))/np.pi))
            k+=1;n*=p
    a=np.array(terms)
    return a[:,0],a[:,1]
def run(T,us,h=0.002,chunk=2_000_000):
    om,c=model(T)
    acc={u:0j for u in us}; n=0
    t0=T; M=int(T/h)
    for s in range(0,M,chunk):
        t=T+h*(np.arange(s,min(s+chunk,M))+0.5)
        Phi=np.zeros_like(t)
        for o,cc in zip(om,c): Phi-=cc*np.sin(o*t)
        N0=t/(2*np.pi)*np.log(t/(2*np.pi*np.e))+7/8
        F=N0+Phi
        sw=F-np.floor(F+0.5)
        for u in us: acc[u]+=np.sum(np.exp(1j*u*t)*sw)
        n+=len(t)
    return {u:-2*np.pi*acc[u]/n for u in us}, Phi.std()
if __name__=="__main__":
    T=float(sys.argv[1]); h=float(sys.argv[2])
    us=[0,0.5,np.log(2),1,np.log(3),2,3]
    r,sd=run(T,us,h)
    print(T,h,"sdPhi(last chunk)=%.3f"%sd)
    for u in us: print("  u=%.4f rho=%.5f%+.5fj |rho|=%.5f"%(u,r[u].real,r[u].imag,abs(r[u])))
