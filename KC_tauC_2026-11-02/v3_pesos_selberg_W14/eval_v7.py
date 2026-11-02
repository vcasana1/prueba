# eval_v7.py : comparacion de taucs por reglas (RMS por alturas en espaciados, parejas cercanas, centroide, Cauchy)
import numpy as np, re, os
g=np.fromfile('zeross.bin'); sp=2*np.pi/np.log(np.maximum(g,20)/(2*np.pi))
R=[('gamma^RvM (zLamb)','/home/claude/KC_v1/run_h32_tau.bin',None),
   ('articulo m=2^(L/2)','/home/claude/KC_v1/run_h32_tauC.bin','/home/claude/KC_v1/analysis_h32_v1.txt'),
   ('agudo p<=2^(1.4L)','full_size1.4_tauC.bin','analysis_size1.4_v5.txt')]
for tag,nm in [('selb_X200','Selberg X=200'),('w14_X200','W1/4 X=200'),('selb_th0.5','Selberg X=2^(0.5L)'),('w14_th0.5','W1/4 X=2^(0.5L)'),('selb_th0.7','Selberg X=2^(0.7L)'),('w14_th0.7','W1/4 X=2^(0.7L)')]:
    if os.path.exists(f'full7_{tag}_tauC.bin'): R.append((nm,f'full7_{tag}_tauC.bin',f'analysis7_{tag}.txt'))
edges=[2,10**3,10**4,10**5,10**6,len(g)+1]
print(f"{'regla':22s}"+"".join(f" [{a:.0e},{b:.0e})" for a,b in zip(edges[:-1],edges[1:]))+"  global   max  gap4765 gap6709 |D_2e6|")
X={}
for nm,f,an in R:
    x=np.fromfile(f); X[nm]=x; d=(x-g)/sp
    rows=[np.sqrt((d[a-1:b-1]**2).mean()) for a,b in zip(edges[:-1],edges[1:])]
    gp=[(x[k]-x[k-1])/sp[k-1] for k in (4765,6709)]
    D='-'
    if an and os.path.exists(an):
        m=[l for l in open(an) if l.strip().startswith('2001052')]
        if m: D=re.search(r'\|\.\|=([0-9.]+)',m[0]).group(1)
    print(f"{nm:22s}"+"".join(f" {r:15.5f}" for r in rows)+f"  {np.sqrt((d[1:]**2).mean()):.5f} {np.abs(d[1:]).max():.4f} {gp[0]:.4f}  {gp[1]:.4f}  {D}")
print("gap verdadero: 4765-4766 = %.4f, 6709-6710 = %.4f espaciados"%((g[4765]-g[4764])/sp[4764],(g[6709]-g[6708])/sp[6708]))
def c(a,b): return np.sqrt((((X[a]-X[b])/sp)[1:]**2).mean())
for a,b in [('Selberg X=2^(0.5L)','Selberg X=2^(0.7L)'),('W1/4 X=2^(0.5L)','W1/4 X=2^(0.7L)')]:
    if a in X and b in X: print(f"Cauchy (solo primos) ||{a} - {b}|| = {c(a,b):.5f} espaciados")
