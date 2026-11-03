Recalculo independiente de tau_k^C (Claude, 2026-11-03), regla del articulo (m_j=floor(2^{j/2}), clip loglog), malla h=1/32.
Entrada: zeross.bin (2 001 052 ceros; incluido en datos_tauC), zLamb.bin (incluido en datos_tauC; solo para Newton de las RvM).
Compilar: gcc -O3 -march=native -fopenmp -o tauC-c tauC-c.c -lm
Ejecutar: ./tauC-c zeross.bin full 32 phi zLamb.bin   # -> full_tauC-c.bin (= datos_tauC/zTaucs2e6c.bin), full_newton-c.bin (= datos_tauC/zLamb_newton-c.bin)
Comparar con zTaucs2e6.bin: python3 compare-c.py  (resultado en comparacion-c.txt)
