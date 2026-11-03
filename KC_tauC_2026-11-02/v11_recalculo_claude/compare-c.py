import numpy as np, sys
A='/mnt/user-data/attachments/17d8d91c-7380-5c16-ac0b-be55137c94f8/'
z=np.fromfile(A+'a4b326cf-zeross.bin','<f8'); zl=np.fromfile(A+'2fff6892-zLamb.bin','<f8')
ref=np.fromfile(A+'2a348cb2-zTaucs2e6.bin','<f8'); mine=np.fromfile('full_tauC-c.bin','<f8'); nw=np.fromfile('full_newton-c.bin','<f8')
def cen(x):
    P=np.cumsum(np.exp(1j*x)); return P, np.cumsum(P)/np.arange(1,len(x)+1)
out=[]
def p(*a):
    s=' '.join(str(x) for x in a); print(s); out.append(s)
d=mine-ref
p('N',len(mine),'monotone mine',bool(np.all(np.diff(mine)>0)))
p('max|mine-file|',np.abs(d).max(),'at k=',int(np.abs(d).argmax())+1,' rms',np.sqrt(np.mean(d**2)),' count>1e-6:',int((np.abs(d)>1e-6).sum()),' count>1e-3:',int((np.abs(d)>1e-3).sum()))
p('max|newton-zLamb|',np.abs(nw-zl).max())
p('rms(tauC-gamma) mine',np.sqrt(np.mean((mine-z)**2)),' file',np.sqrt(np.mean((ref-z)**2)))
Pm,Cm=cen(mine); Pr,Cr=cen(ref); Pg,Cg=cen(z); Pl,Cl=cen(nw)
for n in (220,10000,100000,1000000,2001052):
    p(n,'Cen mine',Cm[n-1],'file',Cr[n-1],'|diff|',abs(Cm[n-1]-Cr[n-1]),'| Cen(gamma)-Cen(mine)',abs(Cg[n-1]-Cm[n-1]),'| Cen(gamma)-Cen(file)',abs(Cg[n-1]-Cr[n-1]),'| Cen newton',Cl[n-1],'Cen gamma',Cg[n-1])
open('comparacion-c.txt','w').write('\n'.join(out)+'\n')
