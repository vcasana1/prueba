import numpy as np
g=np.fromfile('zeross.bin'); tc=np.fromfile('run_h32_tauC.bin'); t0=np.fromfile('run_h32_tau.bin'); K=len(g)
d=tc-t0                      # correction delta_k = tau^C - tau
e0=np.exp(1j*t0.astype(np.longdouble)).astype(np.clongdouble)
eC=np.exp(1j*tc.astype(np.longdouble)).astype(np.clongdouble)
eG=np.exp(1j*g.astype(np.longdouble)).astype(np.clongdouble)
Sg=np.cumsum(eG); S0=np.cumsum(e0); SC=np.cumsum(eC)
first=np.cumsum(1j*d*e0)     # first-order proxy for S^C - S^0
act=SC-S0
print("S_n(gamma)  max|.| =",float(np.abs(Sg).max()), "  S_n(tau) max|.| =",float(np.abs(S0).max()), "  S_n(tauC) max|.| =",float(np.abs(SC).max()))
for n in [10**3,10**4,10**5,5*10**5,10**6,K]:
    i=n-1; print(f"n={n:8d}  S(tauC)-S(tau)={act[i].real:+9.3f}{act[i].imag:+9.3f}i   first-order={first[i].real:+9.3f}{first[i].imag:+9.3f}i   |S(g)|={float(abs(Sg[i])):.3f} |S(tau)|={float(abs(S0[i])):.3f}")
# local statistics of delta
for a,b in [(0,10**4),(10**5,10**5+10**4),(10**6,10**6+10**4),(K-10**4,K)]:
    print(f"k in [{a},{b}): mean delta={d[a:b].mean():+.5f}  rms delta={np.sqrt((d[a:b]**2).mean()):.5f}  rms(tauC-g)={np.sqrt(((tc-g)[a:b]**2).mean()):.5f}  rms(tau-g)={np.sqrt(((t0-g)[a:b]**2).mean()):.5f}")
# spectral content at frequency 1 of the correction: sum e^{i tau_k} delta_k over dyadic ranges
print("\nhalf-range sums  sum_k delta_k e^{i tau_k}:")
for a,b in [(0,10**4),(10**4,10**5),(10**5,10**6),(10**6,K)]:
    z=(d[a:b]*e0[a:b]).sum(); print(f"  k in [{a},{b}): {complex(z):.3f}   (rms term*sqrt(n) = {np.sqrt((d[a:b]**2).mean())*np.sqrt(b-a):.3f})")
# slow beating: FFT of act over n
x=act[::50].astype(np.complex128); X=np.abs(np.fft.fft(x-x.mean())[:len(x)//2]); fr=np.fft.fftfreq(len(x),d=50)[:len(x)//2]
idx=np.argsort(X)[-5:][::-1]; print("\ndominant slow frequencies of S^C-S^0 (cycles per zero):",[(float(fr[i]),float(X[i])) for i in idx])
