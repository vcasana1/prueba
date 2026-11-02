# condicion simetrica de media separacion: |tauc_k - g_k| < min(g_k - g_{k-1}, g_{k+1} - g_k)/2
import numpy as np, sys
g=np.fromfile(sys.argv[1]); x=np.fromfile(sys.argv[2]); tag=sys.argv[3]; K=g.size
sp=2*np.pi/np.log(np.maximum(g,20)/(2*np.pi))
gl=np.concatenate([[np.inf],np.diff(g)]); gr=np.concatenate([np.diff(g),[np.inf]]); m=np.minimum(gl,gr)
bad=np.abs(x-g)>=m/2; idx=np.nonzero(bad)[0]
print(f"[{tag}] incumplen |tauc-gamma| < min(gap)/2: {idx.size} de {K}")
with open(f'media_separacion_{tag}.csv','w') as f:
    f.write("k,gamma_k,tauc_k,err_esp,gap_min_esp,ratio_err_sobre_mitad_gap,lado_del_gap_corto\n")
    for i in idx:
        lado='izq' if gl[i]<gr[i] else 'der'
        f.write(f"{i+1},{g[i]:.9f},{x[i]:.9f},{(x[i]-g[i])/sp[i]:+.5f},{m[i]/sp[i]:.5f},{abs(x[i]-g[i])/(m[i]/2):.3f},{lado}\n")
for i in idx[:20]: print(f"  k={i+1:8d} gamma={g[i]:.6f} err={(x[i]-g[i])/sp[i]:+.4f} esp  gap_min={m[i]/sp[i]:.4f} esp")
if idx.size>20: print("  ...", idx.size-20, "mas; lista completa en el csv")
# hacia donde se mueve: hacia el vecino cercano o alejandose
if idx.size:
    toward=np.sign(x[idx]-g[idx])==np.where(gl[idx]<gr[idx],-1,1)
    print(f"  de ellos, desplazados HACIA el vecino cercano: {toward.sum()}, ALEJANDOSE: {(~toward).sum()}")
