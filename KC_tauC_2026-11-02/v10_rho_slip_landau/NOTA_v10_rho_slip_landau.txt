# v10 — coeficiente de Fourier del diente de sierra, cancelación de Landau y líneas de tau^C (2 nov 2026)

Estado: derivaciones a mano (DED); fórmula de Landau para x fijo citada, sin releer (REC); numérico = modelo F, no el F exacto de IJNT.

## 1. Reducción sin ceros
N entera => <D> = <N-F> = -<F> c.p.t. Por tanto, con F = N_0 + Phi°,
  rho_u(T) = (2pi/T) ∫_T^{2T} e^{iut} <D(t)> dt = -(2pi/T) ∫_T^{2T} e^{iut} <F(t)> dt.
Media del slip: <eps> = rho_0(T)/log T + o(1/log T).

## 2. Cancelación de la línea de Landau
Con G = N - N_0 = S + varpi/pi y D = G - Phi°:
  ∫_T^{2T} e^{iut} G dt = -(1/(iu)) Σ_{T<γ<=2T} e^{iuγ} + B,  |B| <= [2(0.08382 log 2T + 5.326) + O(L)]/u   (cota MCOM v7 solo en el borde).
Landau (x=e^u fijo): Σ e^{iuγ} = -(T/2pi) Λ(x) x^{-1/2} + O_x(log T)  [REC].
=> (2pi/T)∫ e^{iut} G dt = -i Λ(x)/(sqrt(x) log x) + o(1) = -i/(k p^{k/2}) para x=p^k.
Los coeficientes de Fourier de Phi° = -(1/pi) Σ sin(kt log p)/(k p^{k/2}) dan lo mismo => (2pi/T)∫ e^{iut} D dt = o(1) y
  rho_u(T) = -(2pi/T) ∫ e^{iut} Z dt + o_u(1),  Z = floor(D+1/2) (slip entero).
Cota trivial |rho_u| <= 2pi mean|Z|; heurística: Var D ≈ log2/(2pi^2)=0.035, P(|D|>=1/2)≈0.008 => |rho_u| ≲ 0.05.

## 3. tau^C hereda las líneas de Landau
Con b_k ≈ -Phi(tau_k)/L(tau_k) y Σ_k g(tau_k) ≈ ∫ g L dt: Σ e^{iu tau^C_k} ≈ -(T/2pi) Λ(x) x^{-1/2}. tau (RvM) no las tiene.
Delta_C(u) = (iu/log T) rho_u(T) + o(1/log T).

## 4. Numérico (modelo): ver output_rho_u.txt y output_lines.txt
- |rho_u| en u=log x (x potencia de primo) decae ~T^{-0.27}; fuera de ellas, nivel de ruido (~1e-3); rho_0 ~1e-4.
- tau^C modelo reproduce -(T/2 pi n)Λ(x)x^{-1/2} con desviación relativa +1% (log2) a +6% (log7) en T=1e4.
Modelo: primos p^k <= sqrt(1.5T), peso 1 hasta X^{1/2}, taper log-lineal hasta X, fijo por bloque. NO es la interpolación diádica m_j=floor(2^{j/2}) de IJNT.

## 5. Contraste con el texto de ChatGPT (centroides_tau_tauc_taucm)
- Identidades W_{T,x}(u,v), coherencia espectral Psi_C(u_T(v))->1: correctas; por Cauchy-Schwarz Psi_C(u)->1 para todo u=o(log T).
- Discrepancia de definición: la nota v5 (eq. old) dice tau^C = u_k(F_h) con clip_h, h=loglog t; IJNT v50 (eqs. F_def, xk_def) define tau^C con N_0+Phi sin truncar (0 apariciones de "clip"). Con la sin truncar, la cota Δ_Cm-Δ_C ≪ u (loglog T)^2/log^4 T se deduce del cuarto momento; la cota T/loglog T solo vale para el par (F_h,F_H). Revisar v55.
- En u=1, Λ(e)=0: por Landau (x fijo, error O(log T)) Σ_{k<=n} e^{iγ_k}=O(log n) y lo mismo para tau; C_N(1)=O(log N) en ambos. El término lineal b(u)N/(2 log N) no existe en u=1; el coeficiente del diente de sierra no decide el límite del centroide en u=1, solo en u=r log p.
- Abierto: demostrar rho_u(T)->0; restos O(log N) en u=1.

Scripts: rho.py (rho_u), rho2.py (tabla de escalado), lines.py (líneas de tau^C). Resultados: output_rho_u.txt, output_lines.txt.
