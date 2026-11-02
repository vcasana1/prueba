# Notas de la investigación (2026-11-01/02): taucs, centroide, selección de primos, DBSRH

Nomenclatura: **tauc** (singular), **taucs** (plural) = τ_k^C = meas{t : F(t) < k − ½};
**γ_k^RvM** = ordenadas del reloj sin corrección prima (`zLamb.bin`, coincide con Newton a 2·10⁻¹⁰).
Etiquetas: [MED] medido, [DED] deducido aquí, [LIT] localizado, [REC] de memoria sin comprobar.

## 1. Centroide (C,1)

C_N(x) = Σ_{k≤N} w_{k,N} e^{i x_k}, w_{k,N} = (N−k+1)/N (taper móvil: ningún cero se suprime).

* [MED] D_N = C_N(γ) − C_N(γ^RvM) = 0.735337 − 0.989697 i en N = 2 001 052 (es el D_N del corpus).
* [MED] C_N(γ) − C_N(taucs) no converge para ninguna regla probada hasta 2·10⁶.
* [DED] C_N(γ) converge ⇔ D_N converge (la convergencia RVM está probada). La parte no trivial de D_N es
  la transformada de Cesàro en frecuencia 1 de S, es decir, el paseo B(T).
* [LIT] Guinand: e^{−t/2}V(t) + e^{t/2}V(−t) = e^{t/2} + e^{−t/2} − e^{−t/2}/(1 − e^{−2t}), V(t) = Σ_{Im ρ>0} e^{ρt}.
  [DED, verificado algebraicamente] la mitad del lado derecho en t = 1 es Re C_∞ = (1/(2√e))(e − 1/(e²−1)).
* [DED] Teorema abeliano: si el centroide de Cesàro en k existe, vale el de Guinand. La existencia es un
  problema tauberiano para (e^{iγ_k}), abierto incluso bajo RH (Landau da S_n = O(log n)); medido max|S_n| = 4.74.
* Una sola frecuencia no es criterio de RH (el defecto 2(cosh d − 1)e^{iγ} es complejo y puede cancelarse).
  [DED] Familia en u: Δ(z) = V^ord(z) − V^hol(z) = Σ_pares 2(1 − cosh(d_k z)) e^{iγ_k z}, analítica en Im z > 0;
  Δ(iy) = Σ 2(1 − cos d_k y) e^{−γ_k y} ⇒ Δ ≡ 0 ⇔ RH; anulación del valor frontera en un intervalo sin
  {0, ±log p^m} basta por Luzin–Privalov [REC]. Pendiente: paso Cesàro-en-k → Abel-en-γ (Hardy–Riesz [REC]).

## 2. Selección de primos para los taucs

| regla | RMS(tauc − γ) global, espaciados |
|---|---|
| γ^RvM | 0.332 |
| artículo, m = ⌊2^{L/2}⌋ | 0.0557 |
| agudo p ≤ 2^{θL}, θ = 1.0 / 1.2 / 1.4 | 0.0281 / 0.0197 / 0.0146 |

* [MED] RMS·θ² = 0.0281, 0.0283, 0.0285 (θ ∈ [1, 1.4]).
* [DED, heurístico] S − Φ_X ≈ Σ_γ g((γ−t) log X), g(x) = Si(x)/π − ½ sgn x (corte agudo);
  tauc_k − γ_k ≈ π R_k / log X, R_k = suma de núcleos de los demás ceros ⇒ O(θ⁻²).
* Parejas cercanas (las dos menores bajo k = 10⁴): 4765–4766 (0.046 esp.) y 6709–6710 (0.042 esp.).
  Los taucs se repelen; separación 0.31 (artículo) → 0.15 (θ = 1.4). [DED] error de parejas ~ 1/log X.

## 3. Línea GUE (perfil casi GUE) — no abandonar

* [DED, heurístico; ADAPT de Montgomery 1973 y Goldston 1987, ambos REC]
  V(θ) = E[(S − Φ_{T^θ})²] ≈ (1/2π²) ∫_θ^∞ F(α) α⁻² dα, de donde F(θ) = −2π²θ² V′(θ).
* [MED] T ≈ 1.05·10⁶: F_est(α) = α con ±0.02 para α ≤ 1; meseta ≈ 0.87 para α ∈ (1, 1.4).
  T ≈ 10⁴: meseta ~0.8 (ruidosa). Ajuste propio del déficit: ≈ 1.5/log(T/2π) (dos alturas, sin identificar).
* [DED, heurístico] Con perfil W: V_W(θ) ≈ (1/2π²) ∫ (1 − W(α/θ))² F(α) α⁻² dα.
* Importante: la varianza en t uniforme está dominada por el salto propio del núcleo; para los taucs
  importa E[R_k²] (promedio de Palm, condicionado a un cero). Por eso Selberg y W_{1/4} ganan al corte agudo
  en puntos intermedios. Objetivo propuesto: minimizar E[R_k²] con correlación tipo GUE → perfil óptimo analítico.

## 4. Pesos de Doktor

W_b(u) = 1 (u ≤ 1); 1 − v + b v(1 − v), v = u − 1 (1 < u < 2); 0 (u ≥ 2). Selberg: b = 0; candidato b = 1/4.
[DED, verificado] W(u) = ∫_1^2 q(v) 1{u ≤ v} dv con q = −W′; W_{1/4}: q(v) = v/2 + 1/4 > 0, ∫q = 1.
[DED] Para W continuo el término del polo es O(X/(t² log² X)), no √X/(t log X) como en el corte agudo.

## 5. DBSRH

DBSRH(W): lim_{X→∞} tauc_k(X) = γ_k para todo k, con recorte H(t) ≥ |S(t)| (H = 0.08382 log t + 4.444).
* [DED, riguroso] Φ^W_X → S c.t.p. ⇒ DBSRH(W) (convergencia dominada; N = N_0 + S + ϖ/π, |ϖ|/π < ½).
* [DED sobre fórmula de Selberg REC] con pesos suaves basta RH; con corte agudo hace falta fuerza DRH
  (Conrad, Canad. J. Math. [LIT]; DRH: Kaneko–Koyama–Kurokawa arXiv:2206.02612 [LIT]; la DRH no usa RvM).
* Falta: DBSRH(W) ⇒ RH (que los términos x^{β−1/2} impidan la convergencia de las medidas de subnivel).
