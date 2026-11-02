KC_v1 -- calculo numerico de K^C = lim [C_N(gamma) - C_N(tau^C)]   (C_N = (1/N) sum_k (N-k+1) e^{i x_k})
Datos: zeross.bin propio (2 001 052 ceros, float64).  Los .bin (zeross.bin y los tauC/tau generados, ~16 MB c/u) NO van en el zip.
Orden de ejecucion:
  gcc -O3 -march=native -ffast-math -fopenmp -o tauC_v2 tauC_v2.c -lm
  ./tauC_v2 zeross.bin run_h32 32            # tau^C (definicion del paper, h=1/32); run_h64 con 64
  ./tauC_v2 zeross.bin test_nophi 32 nophi   # validacion: Phi=0 -> tau^C = tau (error 4.7e-10)
  MCAP=25 ./tauC_v2 zeross.bin cap25 32      # presupuesto de primos truncado a M=25 (tambien 168)
  python3 check_phi_v1.py                    # Phi, F, F' contra evaluacion independiente en numpy
  python3 analyze_KC_v1.py run_h32           # tabla C_N, D_N = C_N(gamma)-C_N(x)
  python3 diag_v1.py ; python3 window_v1.py  # diagnosticos y diferencia con ventana suave
tauC_v1.c es la version usada en los runs h32/h64; tauC_v2.c es identica salvo la variable MCAP.
