# celda de gamma: en que intervalo entre ceros [g_j, g_{j+1}) cae tauc_k ; j = N(tauc_k)
import numpy as np, sys
g=np.fromfile(sys.argv[1]); x=np.fromfile(sys.argv[2]); tag=sys.argv[3]; K=g.size
j=np.searchsorted(g,x,side='right')      # N(tauc_k)
k=np.arange(1,K+1); off=j-k
u,c=np.unique(off,return_counts=True)
print(f"[{tag}] N(tauc_k) - k :", dict(zip(u.tolist(),c.tolist())))
print("   -1: tauc_k en [g_(k-1), g_k)   0: tauc_k en [g_k, g_(k+1))")
bad=np.nonzero((off!=0)&(off!=-1))[0]
print("   fuera de las dos celdas adyacentes a gamma_k:", bad.size, (bad[:20]+1).tolist())
