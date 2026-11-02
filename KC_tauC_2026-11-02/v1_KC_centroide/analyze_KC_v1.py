# analyze_KC_v1.py : K^C = lim_N [C_N(gamma) - C_N(tau^C)],  C_N = (1/N) sum_{k<=N} (N-k+1) e^{i x_k} = (1/N) sum_{n<=N} S_n
import numpy as np, sys
pref = sys.argv[1]
g  = np.fromfile('zeross.bin').astype(np.longdouble)
tc = np.fromfile(pref+'_tauC.bin').astype(np.longdouble)
t0 = np.fromfile(pref+'_tau.bin').astype(np.longdouble)
K = len(g)
def cesaro(x):
    e = np.exp(1j*x.astype(np.clongdouble) if False else (1j*x).astype(np.clongdouble))
    S = np.cumsum(e); T = np.cumsum(S)
    return S, T/np.arange(1,len(x)+1,dtype=np.longdouble)
Sg, Cg = cesaro(g); Sc, Cc = cesaro(tc); S0, C0 = cesaro(t0)
DC, D0 = Cg-Cc, Cg-C0          # Cesaro differences (gamma - tau^C), (gamma - tau)
BC, B0 = Sg-Sc, Sg-S0          # partial-sum differences  B - B^C
Ns = sorted(set([10**e for e in range(2,7)] + [int(a*10**e) for e in range(3,7) for a in (2,3,5,7)] + [K]))
Ns = [n for n in Ns if n<=K]
f = lambda z: f"{float(z.real):+.6f}{float(z.imag):+.6f}i"
rows=[]
hdr="N, Re C_N(gamma), Im C_N(gamma), Re C_N(tauC), Im C_N(tauC), Re D^C_N, Im D^C_N, |D^C_N|, Re D^0_N, Im D^0_N, |D^0_N|"
print(hdr)
for N in Ns:
    i=N-1
    r=[N,Cg[i].real,Cg[i].imag,Cc[i].real,Cc[i].imag,DC[i].real,DC[i].imag,abs(DC[i]),D0[i].real,D0[i].imag,abs(D0[i])]
    rows.append(r)
    print(f"{N:8d}  C(g)={f(Cg[i])}  C(tC)={f(Cc[i])}  D^C={f(DC[i])} |.|={float(abs(DC[i])):.5f}   D^0={f(D0[i])} |.|={float(abs(D0[i])):.5f}")
np.savetxt(pref+'_KC_table.csv', np.array(rows,dtype=float), delimiter=',', header=hdr, fmt='%.10g')
# stability of the Cesaro limit: statistics of D_N over windows [N/2,N] and last-decile
print("\nwindow statistics of D^C_N (mean +- std over n in [N/2,N]) and of the partial-sum difference B-B^C")
for N in [10**4,10**5,5*10**5,10**6,K]:
    w=slice(N//2-1,N)
    m=DC[w].mean(); s=np.sqrt(((DC[w]-m).real**2+(DC[w]-m).imag**2).mean())
    mb=BC[w].mean(); sb=np.sqrt(((BC[w]-mb).real**2+(BC[w]-mb).imag**2).mean())
    m0=D0[w].mean()
    print(f"N={N:8d}: <D^C>={f(m)} std={float(s):.5f} | <D^0>={f(m0)} | <B-B^C>={f(mb)} std={float(sb):.4f}")
# Cesaro of order 2 : mean of D_N over dyadic ranges
print("\nlast values: B-B^C (N=K) =",f(BC[-1])," C_N(g)(K)=",f(Cg[-1])," C_N(tauC)(K)=",f(Cc[-1]))
np.save(pref+'_DC.npy', np.stack([DC.real,DC.imag]).astype(np.float64)[:, ::100])
