# Rutas hacia RH desde taucm / taucms — notas de sesión (2 nov 2026), v1

Programa de V. Casaña-Giner. Objetos: $F_X^W=N_0+\operatorname{clip}_H(\Phi_X^W-\Pi_X^W)$,
$H(t)=0.08382\log\max(t,e)+5.5$ (MCOM), $\tau_k^{\rm Cms}(X;W)=\operatorname{meas}\{F_X^W<k-\tfrac12\}$,
$G=N-N_0=S+\varpi/\pi$, $\gamma_k^{\rm RvM}$ con $N_0(\gamma_k^{\rm RvM})=k-\tfrac12$.

Etiquetas: [DED] deducido en la sesión (sin comprobar numéricamente salvo que se diga);
[LIT] literatura; [REC] de memoria; [ESPEC] especulativo. Los resultados citados de IJNT, MCOM y v6 se toman tal cual.

---

## 1. Lema de medida (elimina recorte, reordenamiento y colisiones) [DED]

Si en $t$ los valores $F$ y $N$ quedan a lados opuestos de $\lambda=k-\tfrac12$, entonces $|F-N|\ge\tfrac12$ ($N$ es entero).
Con la contractividad de v6 y $|G|\le H$ (MCOM), para todo $k$, $X$, $W$:

$$|\tau_k^{\rm Cms}(X;W)-\gamma_k|\le\operatorname{meas}\{t\in(0,B_k]:\ |P_X^W(t)-G(t)|\ge\tfrac12\}.$$

**Versión localizada** ($k\ge k_0$): el cruce exige $N_0-H<\lambda<N_0+H$, luego basta medir en $(a_k,b_k)$ (barreras de v6).
DBSRH se reduce a $P_X^W\to G$ en medida en compactos.

## 2. Familia Laguerre–Pólya de las taucs [DED]

$E_X(u)=\prod_{k<k_0}(\tau_k^2-u^2)\prod_{k\ge k_0}(1-u^2/\tau_k^2)$. Por las barreras MCOM la familia es normal, $E_X(i)\ge1$,
y todo límite es el producto con $\lambda_k=\lim\tau_k(X_n)$ reales. Las colisiones solo dan ceros reales dobles.

**Criterio localizado (núcleos $j=2$ de MCOM, Lemma *kernel sums*).** Para $d>\tfrac12$ e intervalo $I$:
$$\lim_X\sum_k\big[F_\iota(d,T-\tau_k)+F_\iota(d,T+\tau_k)\big]=\mathfrak M_\iota(d,T),\ \iota\in\{(\mathrm R,2),(\mathrm I,2)\},\ T\in I
\iff \mathrm{RH}\wedge\mathrm{DBSRH}(W),$$
con $\mathfrak M_\iota(d,T)=\Phi_{\rm X}(D_j(s_1))-\frac1{(j-1)!}\sum\Lambda(n)(\log n)^{j-1}\Phi_{\rm X}(n^{-s_1})$.
**Corrección (§4):** con la definición de v6 la hipótesis nunca se cumple porque $\lambda_1\ne\gamma_1$; usar $\tau^\flat$.

## 3. Transformada $\mathcal T$ y teorema de convergencia en medida [DED]

$k_s(t)=(s-\tfrac12-it)^{-2}+(s-\tfrac12+it)^{-2}$, $\ \mathcal T[g](s)=-\int_0^\infty g\,\partial_tk_s\,dt$ ($\operatorname{Re}s>\tfrac12$).

- $\mathcal T[N]=\sum_{\gamma>0}[(s-\tfrac12-i\gamma)^{-2}+(s-\tfrac12+i\gamma)^{-2}]$.
- $\mathcal T[1]=2/(s-\tfrac12)^2$, $\ \mathcal T[\vartheta/\pi]=-\tfrac14\psi_\Gamma'(s/2)$ (residuos; coincide con Malmstén).
- $\mathcal T[P_X^W]=-\sum_nW\Lambda(n)\log n\,n^{-s}+\int_2^\infty W\log u\,u^{-s}du$.
- Para $\operatorname{Re}s>1$: $\lim_X\mathcal T[N_\vartheta+P_X^W]=-(\xi'/\xi)'(s)-E(s)$,
  $E(s)=s^{-2}-2(s-\tfrac12)^{-2}+\int_1^2\log u\,u^{-s}du$ (holomorfa en $\operatorname{Re}s>\tfrac12$).
- $\mathcal T[r_0]=E$ con $r_0(t)=-1+\frac1\pi\arctan2t+\frac1\pi\int_0^{\log2}e^{v/2}\frac{\sin tv}v\,dv$, $r_0(0)=-1$, $r_0=O(1/t)$.

**Teorema.** (C) $P_X^W$ Cauchy en medida en compactos y (UI$_2$) $\sup_X\int|P_X^W|^2(1+t)^{-3}<\infty$ $\Rightarrow$ RH, y $P_\infty=S-r_0$.
Comprobación: $S(0^+)=-1=r_0(0)$, luego $P_\infty(0^+)=0=P_X(0)$. Mecanismo tipo Landau [ADAPT].
Recíproco vía la fórmula explícita ponderada ("Teorema A" de la sesión), **sin comprobar**.

## 4. Corrección a DBSRH de v6 [DED]

Bajo (C)+(UI$_2$): $F_\infty=N-\varpi/\pi-r_0$. Para $k\ge2$, $\tau_k\to\gamma_k$. Para $k=1$, $F_\infty(0^+)=7/8>\tfrac12$ y
$$\tau_1\to\gamma_1-t_\sharp,\qquad t_\sharp=\operatorname{meas}\{t<\gamma_1:-\varpi/\pi-r_0\ge\tfrac12\}>0\quad(t_\sharp\approx0.7,\ \text{a mano}).$$
DBSRH "para todo $k\ge1$" es falsa en el escenario en que debería valer. Corrección mínima:
$$\tau_k^\flat=2\pi+\operatorname{meas}\{t\ge2\pi:F_X^W<k-\tfrac12\}$$
(misma convención que $\gamma_k^{\rm RvM}>2\pi$). Crecimiento estricto de $\tau_k^\flat$ a $X$ finito para $k$ pequeño: pendiente.

## 5. Ruta con MCOM obligatoria [DED]

Margen: $|S-r_0|\le0.08382\log t+4.444+0.39<H$, $|G|\le H-0.887$ en **todo** $t>0$.
$\mathcal E_h=\{|G|\ge h+\tfrac12\}$ es error de taucs independiente de $X$; $\mathcal E_H=\varnothing$, $\mathcal E_{\log\log}$ no (Tsang).

- Desplazamiento de índice: $\gamma_{k-m_k}\le\tau_k^{\rm Cms}\le\gamma_{k+m_k}$, $m_k=\lceil0.16764\log b_k+9.46\rceil$, todo $k$, $X$.
- $\limsup_k\sup_X|\tau_k^{\rm Cms}-\gamma_k|\le2\pi(0.08382+0.0748596)=0.99702<1$.
- Pendiente mínima admisible del recorte $=\mathfrak a=\limsup|\tau_k-\gamma_k|/2\pi$ (IJNT).

**Teorema.** (C′) $\operatorname{clip}_HP_X^W$ Cauchy en medida en compactos y
(S) $\int_0^\infty(|P_X^W|-H)_+(1+t)^{-3}dt\to0$ $\Rightarrow$ RH. Con cualquier recorte que no domine $|S-r_0|$ en todo $t$, (S) es falsa incluso bajo RH.
Bajo (C′), un cero con $\operatorname{Re}\rho>\tfrac12$ aparece como polo doble de $\lim\mathcal T[\operatorname{clip}_HP_X-P_X]$.

## 6. El cruce de ½ [DED sobre IJNT + MCOM]

$\limsup_k|\gamma_k-\gamma_k^{\rm RvM}|=2\pi\mathfrak a<2\pi\cdot0.0748596142711=0.470357<\tfrac12$ (IJNT, Cor. *distance* (a)).
Equivale a $\mathfrak a<1/4\pi=0.0795775$. BW $0.0946\to0.594$; MCOM explícita $\to0.5267$; suelo de Backlund $\to0.4940$.

- Multiplicidad $\le2C_{\rm H}\log\gamma+o=0.1497\log\gamma$ (explícita: $\le0.16764\log\gamma+8.888$).
- Con recorte asintótico de pendiente $c\in(C_{\rm H},1/4\pi)$, la ventana $(a_k,b_k)$ mide $4\pi c<1$ y contiene a $\gamma_k$ y a todos los $\tau_k(X)$.
- **No es RH**: un índice repetido $\gamma_k=\gamma_{k+1}$ cabe en la ventana; hace falta NC de v6:
  $\limsup_X[\varepsilon_k^{\rm RvM}-\varepsilon_{k+1}^{\rm RvM}]<h_k^{\rm RvM}$ para todo $k$.

## 7. (S) desde los primos [DED]

- Reducción incondicional: $|P_X^W|\ll X/\log X+X/t$ ⇒ (S) ⟺ misma integral en $t\le X^{1/2+\delta}$ (altura fija, $X\to\infty$).
- Flujo del exceso: $\frac d{dL}\int_K(|P_L|-H)_+w=\int_{K\cap A_L}\operatorname{sgn}(P_L)\dot P_L\,w$,
  $L^2\dot P_L=\frac1\pi\operatorname{Im}\int_{e^L}^{e^{2L}}q(\log u/L)u^{-1/2-it}d(\psi(u)-u)$ (v6, *cutoff derivative*).
- Negativo: ninguna cota absoluta de $\psi(u)-u$ basta (con fuerza RH, $|\dot P_L|\ll1+|t|$, no integrable).
- Falta (S$^\sharp$): cancelación del error de Chebyshev torcido por $u^{-it}$ en bandas de Selberg, con valor absoluto en $L$.

## 8. Selberg–Fujii + MCOM sobre repeticiones [DED; uniformidad en $k$ de Fujii: REC]

$$\#\{\text{grupos en }(T,2T]\text{ con }m\ge m_0\}\ll N(T)e^{-m_0/(2Ae)},\qquad m\le0.1497\log\gamma.$$
Cola exponencial con corte duro. No da $R(T)=\#\{\gamma_{k+1}=\gamma_k\}=0$: eso requiere la firma aritmética de la repetición en $\Delta_k$.

## 9. Texto externo "transporte ordinal" (`transporte_ordinal_texto_externo.txt`)

Veredicto: (A)–(7) correctos como deducción dados los cuartos momentos de IJNT/v6. Falso que dependa del margen asintótico de MCOM:
basta cualquier $H\asymp\log t$, separando $|G|\le H/2$ (Markov) y $|G|>H/2$ (momentos de Selberg + Hölder).
Compara $\tau^{\rm Cm}$ con $u^\circ$ (fase sin recortar), no con $\gamma_k$; no toca $\Delta_k$ ni $X\to\infty$ con $k$ fijo.
Con momentos $2k$ (Fujii (a)): $\int D\ll_kT/(\log T)^{2k-1}$ para la fase ponderada con $X\le T^{1/(2k)}$; no para $\Phi^\circ$ con presupuesto IJNT.
MCOM sí cuenta frente a $\gamma_k$ (§1, §5): un solo punto con $|G|>H+\tfrac12$ deja error que no desaparece con $X$.

## Pendiente

1. NC de v6 para todo $k$ (resolución puntual de pares).
2. (S$^\sharp$) / (A): integrabilidad absoluta en $L$ de la banda de Selberg del error de Chebyshev torcido.
3. Verificar "Teorema A" (fórmula explícita ponderada) para los recíprocos.
4. Calcular $t_\sharp$ con precisión (solo si se pide) y aplicar $\tau^\flat$ en v7.
