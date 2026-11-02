KC_v9 -- identidad exacta del defecto (IJNT v55, Thm identity) evaluada en las ordenadas de zeross.bin (k <= 2 001 052)
  u_k = k - 1/2 - F(gamma_k),  F = N0 + clip(Phi^o, +-loglog t)  (misma F que los taucs de cada regla)
  reglas: draft  = borrador IJNT (floor(2^{j/2}) primos, sin pesos)
          w14_th0.5, w14_th0.7 = peso W_{1/4}, X_L = 2^{theta L}, n <= X_L^2
uk_at_zeros_v9.py   calcula u_k (NUFFT tipo 3, FINUFFT 2.5.1)          -> *_uk.bin (float64, k=1..K)
analyze_uk_v9.py    j=1 (eta = media u^2, cota 1-4 eta), certificacion por ventanas (j*), peores puntos -> *_analysis.txt, *_windows_W1000.txt
gap_margin_v9.py    margen delta = 1/2-|u_k| frente al hueco al vecino mas cercano               -> gap_margin_v9.txt
decomp_j1_v9.py     descompone la identidad j=1: A = int (N-F)^2 N0', B = int (N-F)^2 (clip Phi)'; comprueba A+B-1/12 = eta -> decomp_w14_th0.5.txt
budget_j1_v9.py     presupuesto explicito: mean D^2, mean D^4, RMS Phi'/N0', cota de Cauchy-Schwarz de B         -> budget_w14_th0.5_j19.txt
mv_dphi_v9.py       valor Montgomery-Vaughan (incondicional) de mean Phi'^2 con constante de Preissmann          -> mv_dphi_th0.5_j19.txt
(g_w14_th0.5_phi.bin / _dphi.bin: rejilla de Phi^o con h=1/32, 290 MB cada uno, no se incluyen; se regeneran con KC_v2/nufft_phi_v7.py zeross.bin g_w14_th0.5 wsize w14 0.5 32)
