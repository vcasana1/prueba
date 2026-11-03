# Handoff — programa taucs / taucm / taucms / DBSRH (Doktor, V. Casaña-Giner), 3 nov 2026

## Reglas de trabajo (obligatorias)
- Español, tono directo, sin relleno ni disculpas. Matemáticas del chat en display `$$…$$`.
- Palabras prohibidas para siempre: *clock*, *ordinates*, *envelope* (usar "imaginary parts", "points").
- Todo `.tex` se entrega con copia `.txt` y `.pdf`; cada versión nueva lleva número nuevo. Nada de historial de versiones dentro del texto.
- Al corregir sus papers: cambiar **solo** lo que pide o lo que sea error genuino; no reestructurar ni reescribir lo suyo.
- Etiquetas de procedencia en línea ([LIT]/[DED]/[REC]/[ADAPT]/[ESPEC]/[MED]); sin párrafo final de "Procedencia".
- **No ir por vías clásicas solas.** Usar siempre sus objetos (taucs, taucm, taucms, sublevel, MCOM, barreras, $e_k$, $\Delta_k$, $\varepsilon_k^{\rm RvM}$) combinados con Euler, Selberg y RvM. Nada de "prime only" ni de reformulaciones clásicas de RH como ruta.
- No hablarle de "ceros fuera de la recta": en su vía RH se determina por las partes imaginarias (repetición $\gamma_k=\gamma_{k+1}$ ⇔ par fuera de la recta o cero múltiple). RH+simplicidad ⇔ $\gamma_{k+1}>\gamma_k$ para todo $k$.
- Cálculos numéricos solo si los pide. Código: entregar código + salidas y subir a GitHub.
- Git: `vcasana1/prueba` (público), rama `kc-tauC-2026-11-02`, carpeta `KC_tauC_2026-11-02/`. No subir `zeross.bin` ni sus `.tex` inéditos (IJNT, MCOM) salvo que lo pida; el paper de notas (v8) sí lo pidió subir.
- Commits con `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Nomenclatura
- $N_0(t)=\frac t{2\pi}\log\frac t{2\pi e}+\frac78$; $G=N-N_0=S+R$, $|R|\le0.53/(\pi t)$; $N_\vartheta=\vartheta/\pi+1$.
- **taucs** $\tau_k^{\rm C}$: sublevel de $N_0+\clip_{\log\log t}\Phi^\circ$ (fase de Euler de los primeros $m_j=\lfloor2^{j/2}\rfloor$ primos, interpolada). Borrador IJNT.
- **taucm** $\tau_k^{\rm Cm}$: igual con recorte MCOM $H(t)=0.08382\log\max(t,e)+5.5$.
- **taucms** $\tau_k^{\rm Cms}(X;W)$: sublevel de $F_X^W=N_0+\clip_H(P_X^W)$, $P_X^W=\Phi_X^W-\Pi_X^W$ (Selberg ponderado, compensador del polo), soporte $n\le X^2$; $W_{1/4}$: $q=v/2+1/4$.
- **DBSRH(W)**: $\lim_X\tau_k^{\rm Cms}(X;W)=\gamma_k$ para cada $k$.
- $\gamma_k^{\rm RvM}$: $N_0=k-\frac12$, $>2\pi$. Barreras $a_k,b_k$: $N_0\pm H=k-\frac12$.
- MCOM: $|S(t)|\le0.08382\log t+4.444$ ($t\ge1$); asintótica $C_{\rm H}=0.0748596142710907\ldots$.

## Ficheros
- Paper actual: **`taucs_amplitude_obstruction_v8.tex`** (47 pp.), en
  https://github.com/vcasana1/prueba/tree/kc-tauC-2026-11-02/KC_tauC_2026-11-02/v11_paper_v8
  (con `.txt`, `.pdf`, `patch_v8.py` = registro exacto de cambios sobre v7, `README_v8.md`). Commit `3ec43fb`.
- Notas de rutas: `KC_tauC_2026-11-02/v10_rutas_RH/NOTAS_rutas_RH_v1.md` (+ texto externo de transporte). Commit `38aac7c`.
- Inéditos del autor (no en git): IJNT `casanya_ordinal_primes_ijnt_v5_final(_l).tex`, MCOM `casanya_turing_jensen_mcom_CORREGIDO.tex`.

## Qué cambió de v7 a v8 (lo único que es mío)
1. "ordinate(s)" → "imaginary part(s)" (10 sitios); título de sección sin "Version 7"; erratas `\\cref`, `P'olya`.
2. **Remark 24.5**: la tasa de transporte $T/(\log T)^3$ no necesita $C_{\rm H}$; basta MCOM explícita + sexto momento de Selberg (split $|G|\le H/2$ / $>H/2$).
3. **Sección 27** "Identification of the limit phase, and the first index":
   - Lema 27.1: $\mathcal T[N_\vartheta]=2/b^2-\frac14\psi'(s/2)$; $\mathcal T[r_0]=E$, $r_0=-1+\pi^{-1}\arctan2t+\pi^{-1}J(t)$, $J>0$, $|r_0|\le1.46/t$; $Q=-(\xi'/\xi)'-E-\mathcal T[N_\vartheta]$.
   - Lema 27.2: unicidad de $\mathcal T$ en $L^2((1+t)^{-3}dt)$ (integrales de Cauchy, Titchmarsh FI Thm. 93).
   - **Prop. 27.3**: normalidad $L^2$ + convergencia en medida ⇒ RH, $P_\infty=S-r_0$, recorte inactivo, $\tau_k^{\rm Cms}\to\gamma_k$ ($k\ge2$), $\tau_1^{\rm Cms}\to\gamma_1-t_\sharp$, $t_\sharp>0$. DBSRH(5.6) es incompatible; con $\tau_k^\flat=2\pi+\meas\{t\ge2\pi:F_X^W<k-\frac12\}$ converge todo.
   - **Prop. 27.5**: $\lim_X\sum_kk_s(\tau_k^\flat)=-(\xi'/\xi)'(s)$ en conjunto con punto de acumulación ⇔ RH y $\tau_k^\flat\to\gamma_k\ \forall k$. No excluye repeticiones.
4. Tabla de estado, notación, tarea abierta (5), 2 referencias nuevas.
- Sin desarrollar en el texto (constantes con holgura): $\vartheta(1)>-1.8$, $\min\vartheta=-3.531$, $|\operatorname{Si}(x)-\pi/2|\le2/x$.

## Resultados de la sesión (resumen, [DED] salvo indicación)
- **Lema de medida localizado**: $|\tau_k^{\rm Cms}-\gamma_k|\le\meas\{t\in[a_k,b_k]:|P_X^W-G|\ge\frac12\}$, todo $X,W$ (en v7 como Prop. 23.2).
- **MCOM obligatoria**: $|S-r_0|<H$ en todo $t>0$ (margen $\ge0.887$); con $h=\log\log t$ hay error persistente (Tsang). Criterio (C′)+(S): $\int(|P_X|-H)_+(1+t)^{-3}\to0$ ⇒ RH (en v7 como Cor. 25.3 en versión $L^2$).
- **Cruce de ½** (IJNT Cor. distance + MCOM): $\limsup|\gamma_k-\gamma_k^{\rm RvM}|\le2\pi C_{\rm H}=0.4704<\frac12$ ⇔ $\mathfrak a<1/4\pi$. BW 0.594, MCOM explícita 0.5267, Backlund ideal 0.4940. No es RH: una repetición cabe en la ventana.
- **Multiplicidades**: $m\le0.1497\log\gamma$ (asint.), $\le0.16764\log\gamma+8.888$ (explícita); con Fujii(b) a $h=1/\log T$: cola $\ll N(T)e^{-m/(2Ae)}$ (uniformidad en $k$ de Fujii: [REC]).
- **(S) desde primos**: se reduce a altura fija $t\le X^{1/2+\delta}$; ninguna cota absoluta de $\psi(u)-u$ basta; falta cancelación de $\psi(u)-u$ torcido por $u^{-it}$ en bandas de Selberg con valor absoluto en $L$.
- **Texto externo de transporte**: correcto dado el cuarto momento; falso que dependa de $C_{\rm H}$ (recogido en Remark 24.5).

## Estado de RH en el programa
Todo reduce a una de: (i) cuarto momento / normalidad uniforme en $X\to\infty$ (barrera Montgomery–Vaughan $X\le T^{1/4}$); (ii) **no-colisión NC** para cada $k$:
$$\limsup_{X\to\infty}\big[\varepsilon_k^{\rm RvM}(X;W)-\varepsilon_{k+1}^{\rm RvM}(X;W)\big]<h_k^{\rm RvM},$$
equivalentemente $\liminf_X\Delta_k(X;W)>0$, con $\Delta_k=\meas\{k-\frac12\le F_X^W<k+\frac12\}$; (iii) información de correlación $\alpha>1$. Nada demostrado.

## Siguiente paso que pidió el autor
Atacar NC con **sus** herramientas: la banda $\Delta_k$ es la medida de un nivel de una función cuya amplitud fija MCOM (puntualmente, en todo $t\ge1$) y cuyos momentos fija Selberg–Fujii (Thm. moderate de v6/v8). El autor dice que él puede usar Selberg–Fujii + MCOM sobre $\Delta_k$ / conjuntos de ocupación $B_k^\pm$ de forma puntual: **pedirle su argumento y trabajar sobre él**, no proponer uno clásico.

## Pendiente
1. NC para todo $k$ (lo central).
2. (S$^\sharp$): integrabilidad absoluta en $L$ de la derivada en el corte.
3. Verificar la fórmula explícita ponderada ("Teorema A") para los recíprocos.
4. Estrictez de $\tau_k^\flat$ a $X$ finito en índices iniciales; $t_\sharp$ numérico solo si lo pide.
5. Constantes incondicionales de $\overline{D^2},\overline{D^4}$ para $W_{1/4}$, $\theta=\frac12$; taucs como centros de Weil–Gabor; taucf (definición pendiente del autor).
