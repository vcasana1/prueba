# margenes_v4.py : margen de cada tauc dentro de (g_{k-1}, g_{k+1}) en espaciados, y ternas de ceros mas apretadas
import numpy as np, sys
g=np.fromfile(sys.argv[1]); x=np.fromfile(sys.argv[2]); K=g.size; sp=2*np.pi/np.log(np.maximum(g,20)/(2*np.pi))
gl=np.r_[0.0,g[:-1]]; gr=np.r_[g[1:],np.inf]; m=np.minimum(x-gl,gr-x)/sp
o=np.argsort(m)[:20]; print("20 menores margenes (k, margen):",[(int(i+1),round(float(m[i]),4)) for i in o])
w=(g[2:]-g[:-2])/sp[:-2]; o=np.argsort(w)[:20]
print("20 ternas mas apretadas (k, anchura, margen del tauc central):",[(int(i+1),round(float(w[i]),4),round(float(m[i+1]),4)) for i in o])
