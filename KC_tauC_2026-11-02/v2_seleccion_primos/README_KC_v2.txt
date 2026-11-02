KC_v2 -- mejora de la seleccion de primos en tau_k^C (IJNT paper, eq. xk_def) -- 2026-11-02
Datos: zeross.bin propio (2 001 052 ceros). Los .bin de datos (zeross.bin, phi/dphi, tauC generados) NO van en el zip;
el tau^C mejorado se entrega aparte como zTaucs2e6_v2.bin.

Reglas comparadas (todas con la interpolacion diadica del paper sobre niveles L, y el recorte |Phi|<=loglog t):
  paper      : P_L = primeros floor(2^{L/2}) primos
  count c    : P_L = primeros floor(c 2^{L/2}) primos
  size theta : P_L = { p <= 2^{L theta} }
  (criba: tambien pesos tipo Selberg y corte p<=t^theta directo, solo en screen_primes_v2.py)

Ficheros
  screen_primes_v2.py         criba: E[(Phi-S)^2], E[(Phi-S)^4] en ventanas, S(t)=N(t)-N0(t) contado sobre zeross.bin,
                              y perfil empirico de pesos por minimos cuadrados.  Salidas screen_1e4/1e5/1e6_v2.txt
  tauC_v3.c / tauC_v4.c       motor directo con pesos por nivel (v3) y modo ventana exacto WIN=TA,TB (v4).
                              Salidas w*_win.txt evaluadas con eval_win_v4.py -> eval_win_v4.txt ; winsum_v4.py
  nufft_phi_v5.py             Phi^o y Phi^o' en toda la malla t=g/32 por NUFFT tipo 1 (FINUFFT 2.5.1, eps 1e-13)
  tauC_fromphi_v5.c           mismo motor de subnivel que tauC_v1.c leyendo Phi de nufft_phi_v5.py
  full_runs_v5.sh             rango completo para theta = 1.0, 1.2, 1.4 ; full_size*.out ; analysis_size*_v5.txt
  rms_by_height_v5.py         RMS(tau^C-gamma)/espaciado por rangos de k -> rms_by_height_v5.txt

Validaciones
  modo ventana vs rango completo (regla paper): max diferencia 1.5e-11
  NUFFT vs evaluacion directa de Phi: 4.6e-10 (redondeo de N0 en double); tau^C NUFFT vs tauC_v1: 9.3e-10

Orden: gcc -O3 -march=native -fopenmp -o tauC_fromphi_v5 tauC_fromphi_v5.c -lm ;
       python3 nufft_phi_v5.py zeross.bin nf_size1.2 size 1.2 32 ; ./tauC_fromphi_v5 zeross.bin nf_size1.2 full_size1.2
  falpha_v6.py                V(theta)=E[(S-Phi_X)^2] con corte agudo X=Tc^theta en ventana; F_est=-2pi^2 theta^2 V'(theta)
                              -> falpha_1e4_v6.txt, falpha_1e6_v6.txt
  cauchy_theta_v6.txt         ||tauC(theta)-tauC(theta')|| (solo primos) y RMS(tauC-gamma)*theta^2
  zTaucs2e6_v2.bin (aparte)   tau_k^C, k=1..2 001 052, regla size theta=1.4 (P_L = {p <= 2^{1.4 L}}), float64 LE, h=1/32
