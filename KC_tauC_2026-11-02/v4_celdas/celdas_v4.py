# celdas_v4.py : comprueba si cada tauc_k cae en la celda de gamma_k
#   celda de media separacion (Voronoi): ((g_{k-1}+g_k)/2, (g_k+g_{k+1})/2)   [g_0 := 0 ; k=K: lado derecho abierto]
#   celda ancha: (g_{k-1}, g_{k+1})
# usage: python3 celdas_v4.py zeross.bin taucs.bin etiqueta
import numpy as np, sys
g=np.fromfile(sys.argv[1]); x=np.fromfile(sys.argv[2]); tag=sys.argv[3]; K=g.size
gl=np.concatenate([[0.0],g[:-1]]); gr=np.concatenate([g[1:],[np.inf]])
lo=(gl+g)/2; hi=np.where(np.isfinite(gr),(g+gr)/2,np.inf)
sp=2*np.pi/np.log(np.maximum(g,20)/(2*np.pi))
outV=(x<=lo)|(x>=hi); outW=(x<=gl)|(x>=gr)
k=np.arange(1,K+1)
print(f"[{tag}] K={K}: fuera de la celda de media separacion: {outV.sum()}   fuera de (g_(k-1), g_(k+1)): {outW.sum()}")
idx=np.nonzero(outV)[0]
with open(f'fuera_celda_{tag}.csv','w') as f:
    f.write("k,gamma_k,tauc_k,err_espaciados,gap_izq_esp,gap_der_esp,lado,fuera_celda_ancha\n")
    for i in idx:
        lado='izq' if x[i]<=lo[i] else 'der'
        gi=(g[i]-gl[i])/sp[i] if i>0 else np.nan; gd=(gr[i]-g[i])/sp[i] if i<K-1 else np.nan
        f.write(f"{i+1},{g[i]:.9f},{x[i]:.9f},{(x[i]-g[i])/sp[i]:+.5f},{gi:.5f},{gd:.5f},{lado},{int(outW[i])}\n")
if idx.size:
    print("  primeros/ultimos k fuera:", (idx[:12]+1).tolist(), "...", (idx[-5:]+1).tolist())
    # las parejas: k fuera suele ir con k+1 fuera
    pares=np.intersect1d(idx, idx-1).size
    print(f"  casos con k y k+1 ambos fuera: {pares}")
    gmin=np.minimum(np.concatenate([[np.inf],np.diff(g)]), np.concatenate([np.diff(g),[np.inf]]))/sp
    print(f"  gap minimo vecino (espaciados) de los fuera: mediana {np.median(gmin[idx]):.4f}, max {gmin[idx].max():.4f}")
