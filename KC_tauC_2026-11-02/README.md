# τ_k^C: centroide K^C y mejora de la selección de primos (2026-11-02)

Material numérico sobre los puntos corregidos de Riemann–von Mangoldt
τ_k^C = meas{t > 0 : F(t) < k − ½}, con F = N₀ + Φ (artículo IJNT, ec. `xk_def`),
calculados para los 2 001 052 ceros de `zeross.bin` (γ ≤ 1 132 490.66).
`zeross.bin` no se incluye. Todos los `.bin` son float64 little-endian, k = 1..2 001 052.

## `v1_KC_centroide/` — límite de Cesàro K^C

C_N(x) = (1/N) Σ_{k≤N} (N−k+1) e^{i x_k}.

* Sin primos (τ_k): C_N(γ) − C_N(τ) → **K⁰ = 0.735337 − 0.989697 i** (estable a 6 cifras desde N ≈ 10⁵).
  Con ventana suave φ(k/N) la diferencia es ≈ 10⁻⁴: K⁰ procede íntegramente de la base.
* Con la regla del artículo (m_j = ⌊2^{j/2}⌋): C_N(γ) − C_N(τ^C) **no converge** hasta N = 2·10⁶
  (|D_N| = 1.39, 1.97, 10.8, 14.7 en N = 10⁴, 10⁵, 10⁶, 2·10⁶), y la versión con ventana no tiende a 0.
* Validaciones: Φ = 0 → τ^C = τ con error 4.7·10⁻¹⁰; malla h = 1/32 frente a 1/64 → 3·10⁻⁶.

## `v2_seleccion_primos/` — selección por tamaño p ≤ 2^{θL}

RMS(τ^C − γ) en espaciados medios, en todo el rango:

| regla | global | máx |
|---|---|---|
| τ sin primos | 0.332 | 1.33 |
| artículo, m = ⌊2^{L/2}⌋ | 0.0557 | 0.25 |
| θ = 1.0 | 0.0281 | 0.13 |
| θ = 1.2 | 0.0197 | 0.11 |
| θ = 1.4 | 0.0146 | 0.089 |

* Ley medida: RMS · θ² = 0.0281, 0.0283, 0.0285 (θ = 1.0, 1.2, 1.4).
* V(θ) = E[(S − Φ_{T^θ})²] con corte agudo: F_est(θ) = −2π²θ²V′(θ) reproduce F(α) = α de
  Montgomery a ±0.02 para α ≤ 1 (T ≈ 10⁶); para α > 1 da una meseta ≈ 0.87 (`falpha_*_v6.txt`).
  La relación es heurística, no está demostrada.
* El centroide sigue sin converger: |D_N| en N = 2·10⁶ = 5.96, 7.02, 3.95 para θ = 1.0, 1.2, 1.4.
* Φ se calcula en la malla t = g/32 por NUFFT de tipo 1 (FINUFFT 2.5.1, eps 10⁻¹³); con la regla
  del artículo reproduce el motor directo a 9·10⁻¹⁰.

## `datos_tauC/`

| fichero | contenido |
|---|---|
| `zTaucs2e6.bin` | τ^C, regla del artículo, malla h = 1/64 |
| `zTaucs2e6_v2.bin` | τ^C, regla p ≤ 2^{1.4 L} (la mejor calculada) |
| `tauC_paper_h32.bin` | τ^C, regla del artículo, h = 1/32 |
| `tauC_paper_MCAP25.bin`, `tauC_paper_MCAP168.bin` | regla del artículo con tope de 25 / 168 primos |
| `tauC_size1.0.bin`, `tauC_size1.2.bin` | p ≤ 2^{θL}, θ = 1.0 y 1.2 |
| `tau_RvM.bin` | τ_k con N₀(τ_k) = k − ½ (Newton) |

## Reproducir

```
gcc -O3 -march=native -fopenmp -o tauC_fromphi_v5 v2_seleccion_primos/tauC_fromphi_v5.c -lm
python3 v2_seleccion_primos/nufft_phi_v6.py zeross.bin nf size 1.4 32
./tauC_fromphi_v5 zeross.bin nf out          # -> out_tauC.bin
python3 v1_KC_centroide/analyze_KC_v1.py out  # tabla de C_N y D_N
```
