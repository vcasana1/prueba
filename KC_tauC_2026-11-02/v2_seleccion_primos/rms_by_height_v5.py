# RMS of (tau^C - gamma) in mean spacings, by height ranges, for full-range runs
import numpy as np, sys
g=np.fromfile('zeross.bin'); sp=2*np.pi/np.log(np.maximum(g,20)/(2*np.pi))
runs=[('tau (sin primos)','/home/claude/KC_v1/run_h32_tau.bin'),('articulo m=2^(j/2)','/home/claude/KC_v1/run_h32_tauC.bin')]+[(f'p<=2^(L*{th})',f'full_size{th}_tauC.bin') for th in sys.argv[1:]]
edges=[1,10**3,10**4,10**5,10**6,len(g)+1]
print(f"{'regla':22s}"+"".join(f"  k in [{a:.0e},{b:.0e})" for a,b in zip(edges[:-1],edges[1:]))+"   max|.|/sp   global")
for nm,f in runs:
    x=np.fromfile(f); d=(x-g)/sp; row=[]
    for a,b in zip(edges[:-1],edges[1:]): row.append(np.sqrt((d[a-1:b-1]**2).mean()))
    print(f"{nm:22s}"+"".join(f"  {r:17.5f}" for r in row)+f"   {np.abs(d[1:]).max():9.4f}   {np.sqrt((d[1:]**2).mean()):.5f}")
