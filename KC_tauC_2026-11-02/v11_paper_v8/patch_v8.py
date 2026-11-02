#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Apply the v8 edits to a copy of taucs_amplitude_obstruction_v7.tex.
Only: (i) the author's standing word rules, (ii) two LaTeX typos,
(iii) new material derived in the 2 Nov 2026 session, (iv) header/version."""
import sys, re
src = sys.argv[1]; dst = sys.argv[2]
s = open(src, encoding='utf-8').read()

def rep(old, new, count=1):
    global s
    assert s.count(old) == count, (s.count(old), old[:70])
    s = s.replace(old, new)

# ---------- header / version ----------
rep("pdftitle={Amplitude and fourth-moment obstructions for prime-corrected sublevel approximations - Research Note v7}",
    "pdftitle={Amplitude and fourth-moment obstructions for prime-corrected sublevel approximations - Research Note v8}")
rep(r"\fancyhead[R]{\small Research note v7 \textbullet\ 2 Nov 2026}",
    r"\fancyhead[R]{\small Research note v8 \textbullet\ 2 Nov 2026}")
rep(r"{\small Research note v7 --- 2 November 2026}\par",
    r"{\small Research note v8 --- 2 November 2026}\par")

# ---------- author's word rule (no "ordinate(s)") ----------
rep("We study the author's Euler-prime sublevel ordinates $\\tau_k^{\\mathrm C}$,",
    "We study the author's Euler-prime sublevel points $\\tau_k^{\\mathrm C}$,")
rep("corresponding zero ordinates (DBSRH).", "corresponding imaginary parts of zeros (DBSRH).")
rep("$\\mathbf1_{\\{t<\\gamma_k\\}}$ a.e., even for repeated ordinates,",
    "$\\mathbf1_{\\{t<\\gamma_k\\}}$ a.e., even for repeated imaginary parts,")
rep("statement in \\emph{physical ordinate units}; it is not a bound by",
    "statement in \\emph{physical units of the imaginary part}; it is not a bound by")
rep("$u_k^\\circ$ are \\emph{not} the true ordinates, and these results",
    "$u_k^\\circ$ are \\emph{not} the true imaginary parts, and these results")
rep("Let $0<a<b$ not be ordinates of zeros. Suppose $F$ is",
    "Let $0<a<b$ not be imaginary parts of zeros. Suppose $F$ is")
rep("The sum counts zero ordinates with multiplicity.",
    "The sum counts imaginary parts of zeros with multiplicity.")
rep("At an ordinate with multiplicity $m$, the missing telescopic",
    "At an imaginary part with multiplicity $m$, the missing telescopic")
rep("zeros with the same ordinate. It is important not to import",
    "zeros with the same imaginary part. It is important not to import")
rep("real pair has $h^2\\ge0$, a double ordinate has $h^2=0$,",
    "real pair has $h^2\\ge0$, a repeated imaginary part has $h^2=0$,")

# ---------- author's rule: no version history in the text ----------
rep("\\section{Version 7: arithmetic scope and source attribution}\\label{sec:v7-scope}",
    "\\section{Arithmetic scope and source attribution}\\label{sec:v7-scope}")

# ---------- two LaTeX typos ----------
rep("The inverse estimator from \\\\cref{sec:identify-v4} measures a",
    "The inverse estimator from \\cref{sec:identify-v4} measures a")
rep("A separate Laguerre--P'olya proposal in the public notes",
    "A separate Laguerre--P\\'olya proposal in the public notes")

# ---------- notation table: new rows ----------
rep("$\\Psi_g$ & finite-interval box convolution of a test function & \\cref{thm:v7-box}\\\\",
    "$\\Psi_g$ & finite-interval box convolution of a test function & \\cref{thm:v7-box}\\\\\n"
    "$J(t)$, $r_0(t)$ & Dirichlet-type integral and the limit correction $r_0=-1+\\pi^{-1}\\arctan(2t)+\\pi^{-1}J(t)$ & \\cref{lem:v8-smooth}\\\\\n"
    "$E(s)$ & $s^{-2}-2(s-\\frac12)^{-2}+\\int_1^2(\\log u)u^{-s}du$, the transform of $r_0$ & \\cref{lem:v8-smooth}\\\\\n"
    "$P_\\infty$, $F_\\infty$ & limit in measure of $P_{X_j}^W$; $N_0+P_\\infty$ & \\cref{prop:v8-first}\\\\\n"
    "$t_\\sharp$ & $\\meas\\{0<t<\\gamma_1:R(t)+r_0(t)\\le-\\frac12\\}$, the first-index defect & \\cref{prop:v8-first}\\\\\n"
    "$\\tau_k^\\flat(X;W)$ & $2\\pi+\\meas\\{t\\ge2\\pi:F_X^W<k-\\frac12\\}$, sublevel points measured from $2\\pi$ & \\cref{prop:v8-first}\\\\\n"
    "$\\mathcal K_X(s)$ & $\\sum_k k_s(\\tau_k^\\flat(X;W))$ & \\cref{prop:v8-criterion}\\\\")

# ---------- transport: the asymptotic MCOM estimate is not needed ----------
old_rem_end = ("$u_k^\\circ$ are \\emph{not} the true imaginary parts, and these results\n"
               "do not imply noncollapse for a fixed $k$ as $X\\to\\infty$.\n"
               "\\end{remark}\n")
new_rem = old_rem_end + r"""
\begin{remark}[The asymptotic estimate \eqref{eq:asymptMCOM} is not needed for the rate]
\label{rem:v8-transport-noasympt}
The conclusions of \cref{thm:v7-transport} and
\cref{cor:v7-weighted-transport} hold with the explicit bound
\eqref{eq:MCOM} alone, provided Selberg's unconditional moment
estimate $\int_T^{3T}|S(t)|^{6}\,dt\ll T(\log\log T)^3$
\cite{Selberg1946} is added. Split $[T/2,3T]$ into
$B_1=\{|G|\le H/2\}$ and $B_2=\{|G|>H/2\}$. On $E_T\cap B_1$ one has
$|P-G|\ge D+H/2\ge D+\tfrac12A_{\mathrm M}\log(T/2)$, and the Markov
step of the proof applies verbatim with $c=\tfrac12A_{\mathrm M}$.
On $B_2$, $|S|\ge H/2-|R|\ge\tfrac12A_{\mathrm M}\log T+2$ for large
$T$, so $\meas(B_2)\ll T(\log\log T)^3(\log T)^{-6}$. Since $|G|\le H$,
$D\le|P-G|$ pointwise, and H\"older's inequality with the fourth
moment gives
$\int_{B_2}D\,dt\le\bigl(\int|P-G|^4\bigr)^{1/4}\meas(B_2)^{3/4}
\ll T(\log\log T)^{9/4}(\log T)^{-9/2}\ll T(\log T)^{-3}$,
while $\meas(E_T\cap B_2)\le\meas(B_2)\ll T(\log T)^{-4}$.
Thus the only property of the truncation used by the transport rate
is $H\asymp\log t$ together with the covering $|G|\le H$; the
asymptotic coefficient $C_{\mathrm H}$ removes the factor
$(\log\log T)^{O(1)}$ from the exceptional set but does not change
the order $T(\log T)^{-3}$.
\end{remark}
"""
rep(old_rem_end, new_rem)

# ---------- new section after the transform section ----------
anchor = "\\section{A box-transform and a parity-sensitive multiplicity test}\n\\label{sec:v7-parity}"
new_sec = r"""\section{Identification of the limit phase, and the first index}
\label{sec:v8-limit}
The criterion of \cref{thm:v7-normality} uses no convergence of the
prime phases. If convergence in measure is added, the limit phase is
identified exactly, and the first sublevel point of \eqref{eq:Cms} is
seen \emph{not} to converge to $\gamma_1$. The deductions of this
section were made in the research conversation of 2 November 2026;
their external inputs are Malmst\'en's integral for $\log\Gamma$
\cite[\S12.31]{WhittakerWatson}, the boundary formula for Cauchy
integrals of $L^2$ functions \cite[Thm.~93]{TitchmarshFI}, the Hadamard
product of $\xi$, and the explicit bound $|R(t)|\le0.53/(\pi t)$
already used in \eqref{eq:MCOMH}. Throughout, $N_\vartheta=\vartheta/\pi+1$,
so that $N=N_\vartheta+S$ for every $t>0$ by \eqref{eq:RVM}, $b:=s-\frac12$,
and $\psi=\Gamma'/\Gamma$.

\begin{lemma}[Transforms of the smooth terms]\label{lem:v8-smooth}
Let $\Re s>\frac12$.
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item For $a>0$, $\mathcal T[\sin(a\,\cdot)](s)=\pi a^2e^{-ab}$;
 moreover $\mathcal T[1](s)=2/b^2$ and $\mathcal T[t](s)=0$.
\item $\mathcal T[N_\vartheta](s)=2/b^2-\tfrac14\psi'(s/2)$.
\item Put
\begin{equation}\label{eq:v8-r0}
 J(t)=\int_0^{\log2}e^{v/2}\frac{\sin(tv)}{v}\,dv,\qquad
 r_0(t)=-1+\frac1\pi\arctan(2t)+\frac1\pi J(t).
\end{equation}
Then $r_0$ is real-analytic, $r_0(0)=-1$, $J(t)>0$ and
$|r_0(t)|\le1.46/t$ for $t>0$, and
\begin{equation}\label{eq:v8-E}
 \mathcal T[r_0](s)=E(s):=\frac1{s^2}-\frac2{b^2}
 +\int_1^2(\log u)\,u^{-s}\,du .
\end{equation}
\item For $\Re s>1$ the limit $Q$ of \eqref{eq:v7-meromorphic-limit}
 satisfies $Q(s)=-(\xi'/\xi)'(s)-E(s)-\mathcal T[N_\vartheta](s)$.
\end{enumerate}
\end{lemma}
\begin{proof}
(i) Integrating by parts,
$\mathcal T[\sin(a\,\cdot)](s)=a\int_0^\infty\cos(at)k_s(t)\,dt
=\frac a2\int_{\mathbb R}e^{iat}k_s(t)\,dt$, since $k_s$ is even.
For $a>0$ and $\Re b>0$, residues give
$\int_{\mathbb R}e^{iat}(b+it)^{-2}dt=2\pi ae^{-ab}$ and
$\int_{\mathbb R}e^{iat}(b-it)^{-2}dt=0$. Next,
$\mathcal T[1]=k_s(0)-k_s(\infty)=2/b^2$, and
$\mathcal T[t]=\int_0^\infty k_s\,dt=\frac12\int_{\mathbb R}k_s\,dt=0$
because each of $(b\pm it)^{-2}$ integrates to zero over $\mathbb R$.

(ii) Malmst\'en's formula gives, for $\Re z>0$,
$\log\Gamma(z)=\int_0^\infty\bigl[(z-1)e^{-u}
-\frac{e^{-u}-e^{-zu}}{1-e^{-u}}\bigr]\frac{du}u$. With
$z=\frac14+\frac{it}2$,
$\Im\log\Gamma(z)=\int_0^\infty\bigl[\frac t2e^{-u}
-\frac{e^{-u/4}\sin(tu/2)}{1-e^{-u}}\bigr]\frac{du}u$,
whose integrand is $O(t+ut^3)$ as $u\to0^+$; hence Fubini's theorem
allows $\mathcal T$ to be applied under the integral. By (i),
the term linear in $t$ contributes nothing and
$\mathcal T[\sin(\tfrac u2\,\cdot)]=\pi\tfrac{u^2}4e^{-bu/2}$, so
$\mathcal T[\Im\log\Gamma(\tfrac14+\tfrac{i\cdot}2)]
=-\frac\pi4\int_0^\infty\frac{u\,e^{-su/2}}{1-e^{-u}}\,du
=-\frac\pi4\psi'(s/2)$. Since
$\vartheta(t)=\Im\log\Gamma(\frac14+\frac{it}2)-\frac t2\log\pi$ and
$\mathcal T[t]=0$, (ii) follows.

(iii) From $\arctan(2t)=\int_0^\infty e^{-v/2}\sin(tv)\,dv/v$ and (i),
$\mathcal T[\pi^{-1}\arctan(2\,\cdot)]=\int_0^\infty ve^{-sv}dv=s^{-2}$
and $\mathcal T[\pi^{-1}J]=\int_0^{\log2}ve^{(1-s)v}dv
=\int_1^2(\log u)u^{-s}du$; with $\mathcal T[-1]=-2/b^2$ this is
\eqref{eq:v8-E}. Analyticity and $r_0(0)=-1$ are clear. For the
bounds write $J(t)=\operatorname{Si}(t\log2)+\int_0^{\log2}\varphi(v)\sin(tv)\,dv$
with $\varphi(v)=(e^{v/2}-1)/v$ increasing from $\frac12$ to $0.598$.
Integration by parts gives $|\int_0^{\log2}\varphi\sin(tv)\,dv|\le2\varphi(\log2)/t\le1.2/t$,
and also $\le\varphi(\log2)\,t(\log2)^2/2\le0.144\,t$. For $t\le2.9$,
$\operatorname{Si}(t\log2)\ge0.79\,t\log2$ by concavity of $\operatorname{Si}$ on
$[0,\pi]$, so $J(t)\ge0.40\,t>0$; for $t\ge2.9$,
$\operatorname{Si}(t\log2)\ge\operatorname{Si}(2\pi)>1.41$ and $1.2/t\le0.42$, so $J(t)>0$.
Finally $r_0=\pi^{-1}\bigl[J-\frac\pi2-\arctan\frac1{2t}\bigr]$ and
$|\operatorname{Si}(x)-\frac\pi2|\le2/x$ give
$|r_0(t)|\le\pi^{-1}(2/\log2+1.2+\tfrac12)/t\le1.46/t$.

(iv) From $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$,
$(\xi'/\xi)'=(\zeta'/\zeta)'-s^{-2}-(s-1)^{-2}+\frac14\psi'(s/2)$.
Using (ii), (iii) and $\int_1^\infty(\log u)u^{-s}du=(s-1)^{-2}$,
\[
-(\xi'/\xi)'-E-\mathcal T[N_\vartheta]
=-(\zeta'/\zeta)'+\int_1^\infty(\log u)u^{-s}du-\int_1^2(\log u)u^{-s}du
=Q(s).\qedhere
\]
\end{proof}

\begin{lemma}[Uniqueness for the transform]\label{lem:v8-unique}
Let $g$ be measurable on $(0,\infty)$ with
$\int_0^\infty|g(t)|^2(1+t)^{-3}dt<\infty$ and $\mathcal T[g](s)=0$
for all $\Re s>\frac12$. Then $g=0$ almost everywhere.
\end{lemma}
\begin{proof}
Let $g_o$ be the odd extension of $g$ to $\mathbb R$. Since
$\partial_tk_s$ is odd, $\mathcal T[g](s)=-\frac12\int_{\mathbb R}g_o\,\partial_tk_s\,dt$,
and $\partial_tk_s(t)=2(t+ib)^{-3}+2(t-ib)^{-3}$. Put
$C(z)=\int_{\mathbb R}g_o(t)(t-z)^{-3}dt$ for $z\notin\mathbb R$
(absolutely convergent and holomorphic). Oddness gives $C(-z)=C(z)$,
so $\mathcal T[g](s)=-2C(ib)$, and the hypothesis says $C\equiv0$ on
$\{\Im z>0\}$, hence on $\{\Im z<0\}$ as well. Let
$C_1(z)=\int_{\mathbb R}g_o(t)\bigl[(t-z)^{-1}-(t-i)^{-1}-(z-i)(t-i)^{-2}\bigr]dt$;
the bracket is $O(|t|^{-3})$, $C_1$ is holomorphic off $\mathbb R$ and
$C_1''=2C\equiv0$, so $C_1$ is a quadratic polynomial $p_\pm$ on each
half-plane. The bracket equals $(z-i)^2(t-z)^{-1}(t-i)^{-2}$, so
$C_1(z)=(z-i)^2\int_{\mathbb R}w(t)(t-z)^{-1}dt$ with
$w(t)=g_o(t)(t-i)^{-2}\in L^2(\mathbb R)$. By the boundary formula for
Cauchy integrals of $L^2$ functions \cite[Thm.~93]{TitchmarshFI},
$2\pi i\,w(x)=\lim_{y\to0^+}[C_1(x+iy)-C_1(x-iy)](x-i)^{-2}
=(p_+(x)-p_-(x))(x-i)^{-2}$ for almost every $x$. Hence $g_o$ agrees
a.e. with a polynomial of degree at most two; being odd, $g_o(t)=ct$,
and $\int_0^\infty c^2t^2(1+t)^{-3}dt<\infty$ forces $c=0$.
\end{proof}

\begin{proposition}[Limit phase and the first index]\label{prop:v8-first}
Fix $W\in\mathcal W$. Suppose that a sequence $X_j\to\infty$
satisfies \eqref{eq:v7-normality-assumption} and that
$P_{X_j}^W$ converges in measure on every compact subset of
$(0,\infty)$ to a function $P_\infty$. Then:
\begin{enumerate}[label=(\alph*),leftmargin=*]
\item RH holds and $P_\infty=S-r_0$ almost everywhere on $(0,\infty)$,
 with $r_0$ from \eqref{eq:v8-r0};
\item $|P_\infty(t)|<H(t)$ for every $t>0$: the MCOM truncation is
 inactive in the limit, and $F_\infty:=N_0+\clip_H(P_\infty)=N-R-r_0$;
\item $\lim_{j\to\infty}\tau_k^{\Cms}(X_j;W)=\gamma_k$ for every $k\ge2$;
\item $\lim_{j\to\infty}\tau_1^{\Cms}(X_j;W)=\gamma_1-t_\sharp$, where
\begin{equation}\label{eq:v8-tsharp}
 t_\sharp:=\meas\{0<t<\gamma_1:\ R(t)+r_0(t)\le-\tfrac12\}>0 .
\end{equation}
\end{enumerate}
Consequently $\mathrm{DBSRH}(W)$ in the form \eqref{eq:DBSRH} is
incompatible with the two hypotheses. With the sublevel points
measured from $2\pi$, the minimum of $N_0$,
\begin{equation}\label{eq:v8-flat}
 \tau_k^\flat(X;W):=2\pi+\meas\{t\ge2\pi:\ F_X^W(t)<k-\tfrac12\},
\end{equation}
one has $\lim_j\tau_k^\flat(X_j;W)=\gamma_k$ for every $k\ge1$.
\end{proposition}
\begin{proof}
(a) RH is \cref{thm:v7-normality}. By Fatou's lemma $P_\infty$
satisfies the weighted $L^2$ bound; the bound gives uniform
integrability on compact sets and, by Cauchy--Schwarz, uniformly small
contributions of $(0,\delta)\cup(B,\infty)$ to
$\int|P_{X_j}^W|(1+t)^{-3}dt$, so $\mathcal T[P_{X_j}^W](s)\to\mathcal T[P_\infty](s)$
for every $\Re s>\frac12$. Hence $\mathcal T[P_\infty]$ is the
holomorphic continuation of $Q$ obtained in the proof of
\cref{thm:v7-normality}. Integrating by parts,
$\mathcal T[N](s)=\sum_kk_s(\gamma_k)=\sum_\rho(s-\rho)^{-2}$ under RH,
which equals $-(\xi'/\xi)'(s)$ by the Hadamard product. By
\cref{lem:v8-smooth}(iv), on $\Re s>\frac12$ (where both sides are
holomorphic under RH),
$\mathcal T[S-r_0]=\mathcal T[N]-\mathcal T[N_\vartheta]-\mathcal T[r_0]
=-(\xi'/\xi)'-\mathcal T[N_\vartheta]-E=Q=\mathcal T[P_\infty]$.
\Cref{lem:v8-unique} applied to $g=P_\infty-S+r_0$, which satisfies
the weighted $L^2$ bound since $|S|\ll\log(2+t)$ and $r_0$ is bounded,
proves (a).

(b) For $0<t<\gamma_1$ we have $N=0$, so $S=-1-\vartheta/\pi\in[-1,0.13]$
because $\vartheta\ge-3.531$ on $(0,\infty)$; also
$|J(t)|\le t\int_0^{\log2}e^{v/2}dv\le0.83\,t$, so for $t<2$,
$|r_0|\le1.53$ and $|S-r_0|\le2.6<5.5\le H(t)$. For $t\ge2$,
\eqref{eq:MCOM} (extended to zero heights by right limits) and
$|r_0|\le0.73$ give $|S-r_0|\le0.08382\log t+5.18<H(t)$.
Thus $\clip_H(P_\infty)=P_\infty$ and $F_\infty=N_0+S-r_0=N-R-r_0$.

(c), (d) Since $F_X^W\ge N_0-H$, every sublevel set at level
$\lambda_k=k-\frac12$ lies in $(0,B_k]$ with
$B_k:=\sup\{N_0-H<\lambda_k\}$ independent of $X$. The truncation is
$1$-Lipschitz, so $F_{X_j}^W\to F_\infty$ in measure on compact subsets
of $(0,B_k]$. On each interval between consecutive distinct zero heights,
$F_\infty-\lambda_k$ is real-analytic and not constant (as
$R+r_0\to0$ at infinity and $R+r_0\to\frac18-1$ as $t\to0^+$), so
$\meas\{F_\infty=\lambda_k\}=0$. Therefore, for each $\delta>0$,
$\meas(\{F_{X_j}^W<\lambda_k\}\cap[\delta,B_k])\to\meas(\{F_\infty<\lambda_k\}\cap[\delta,B_k])$,
while the part in $(0,\delta)$ has measure at most $\delta$; letting
$\delta\to0$, $\lim_j\tau_k^{\Cms}(X_j;W)=\meas\{F_\infty<\lambda_k\}$.
We compute this set. For $t\ge\gamma_1$, $|R|\le0.53/(\pi\gamma_1)<0.012$
and $|r_0|\le1.46/\gamma_1<0.104$, so $|R+r_0|<\frac12$ and, $N$ being
an integer, $N-R-r_0<k-\frac12$ holds exactly when $N\le k-1$, that is
for $t<\gamma_k$. For $0<t<\gamma_1$, $F_\infty=-R-r_0$; since $J>0$,
$-r_0\le1-\pi^{-1}\arctan(2t)<1$, and $R>-\frac12$ there (for
$t\ge1$ by $|R|\le0.53/\pi$; for $t<1$ by $\vartheta(t)\ge\vartheta(1)>-1.8$
and $N_0(t)\le\frac78$, which give $R(t)>-0.45$), so $F_\infty<\frac32$.
Hence $\{F_\infty<\lambda_k\}=(0,\gamma_k)$ up to a null set for every
$k\ge2$, proving (c), while
$\{F_\infty<\frac12\}=\{t<\gamma_1:\ R+r_0>-\frac12\}$, whose measure is
$\gamma_1-t_\sharp$. Since $R(t)+r_0(t)\to-\frac78<-\frac12$ as
$t\to0^+$ and both functions are continuous on $(0,\gamma_1)$, the set in
\eqref{eq:v8-tsharp} contains an interval $(0,\epsilon)$, so $t_\sharp>0$.
For \eqref{eq:v8-flat}: on $[2\pi,\gamma_1)$, $|R|\le0.53/(2\pi^2)<0.03$
and $-r_0\le\pi^{-1}\arctan(1/(2t))\le1/(4\pi^2)$, so $F_\infty<\frac12$
there, giving $\lim_j\tau_1^\flat=2\pi+(\gamma_1-2\pi)=\gamma_1$; for
$k\ge2$ the interval $(0,2\pi)$ lies inside $\{F_\infty<\lambda_k\}$, so
$\tau_k^\flat$ and $\tau_k^{\Cms}$ have the same limit.
\end{proof}

\begin{remark}[Scope]
The two hypotheses of \cref{prop:v8-first} concern only the compensated
weighted prime phase. The proposition does not prove RH; it shows that,
in the only scenario in which the transform criterion applies, the
limit of the prime phase is $S-r_0$ and not $S$, that the correction
$r_0$ is invisible at every zero height above $\gamma_1$ (where
$|R+r_0|<\frac12$), and that it displaces the first sublevel point by a
positive amount $t_\sharp$. The numerical observation $P_X^W(0)=0$ versus
$N_0(0^+)=\frac78$ in \cref{sec:v7-complex} is the finite-cutoff trace of
the same fact. The conditional statements of
\cref{cor:rh-rvm-v5} and \cref{sec:v7-open} remain as formulated;
whenever DBSRH is assumed there, the first index should be read with
\eqref{eq:v8-flat}.
\end{remark}

\begin{proposition}[Transform criterion for the normalized sublevel points]
\label{prop:v8-criterion}
Fix $W\in\mathcal W$ and put, for $\Re s>\frac12$,
\begin{equation}\label{eq:v8-KX}
 \mathcal K_X(s):=\sum_{k\ge1}k_s\bigl(\tau_k^\flat(X;W)\bigr),
\end{equation}
which converges absolutely and locally uniformly, uniformly in $X\ge2$,
because $\tau_k^\flat\ge a_k$ for $k\ge k_0$ (\cref{thm:block-sep-v3})
and $\sum_ka_k^{-2}<\infty$, while $2\pi\le\tau_k^\flat\le B_k$ for
$k<k_0$. Let $\Sigma\subset\{\Re s>\frac12\}$ have an accumulation point
in that half-plane. Then
\begin{equation}\label{eq:v8-criterion}
 \lim_{X\to\infty}\mathcal K_X(s)=-\Bigl(\frac{\xi'}{\xi}\Bigr)'(s)
 \ \ \text{for every }s\in\Sigma
 \quad\Longleftrightarrow\quad
 \mathrm{RH}\ \text{and}\ \lim_{X\to\infty}\tau_k^\flat(X;W)=\gamma_k
 \ \text{for every }k\ge1 .
\end{equation}
For $\Sigma\subset\{\Re s>1\}$ the right side of the limit is the
convergent prime expression
$-(\zeta'/\zeta)'(s)-s^{-2}-(s-1)^{-2}+\frac14\psi'(s/2)$, so the
criterion compares the sublevel points with Euler, $\Gamma$ and the
pole only; for $s=\frac12+d+iT$ the real and imaginary parts of
$k_s(\tau)$ are the kernels of order $j=2$ used in \cite{CasanaMCOM}.
\end{proposition}
\begin{proof}
($\Leftarrow$) Under RH, $-(\xi'/\xi)'(s)=\sum_\rho(s-\rho)^{-2}=\sum_kk_s(\gamma_k)$,
and dominated convergence with the bound
$|k_s(\tau)|\le C_s\min(1,\tau^{-2})$ gives the limit.
($\Rightarrow$) Let $X_n\to\infty$. As $\tau_k^\flat(X_n)\in[2\pi,B_k]$
is nondecreasing in $k$, a diagonal subsequence satisfies
$\tau_k^\flat(X_{n_i})\to\lambda_k$ for every $k$, with
$2\pi\le\lambda_k\le\lambda_{k+1}$ and $\lambda_k\ge a_k$ for $k\ge k_0$.
By dominated convergence $\mathcal K_{X_{n_i}}\to\mathcal K_\infty:=\sum_kk_s(\lambda_k)$
locally uniformly on $\Re s>\frac12$; the same series converges locally
uniformly on $\mathbb C$ away from the points $\frac12\pm i\lambda_k$, so
$\mathcal K_\infty$ is meromorphic on $\mathbb C$ with double poles
exactly at $\frac12\pm i\lambda_k$, the leading coefficient at
$\frac12+i\mu$ being the number of $k$ with $\lambda_k=\mu$. The
hypothesis and the identity theorem for meromorphic functions give
$\mathcal K_\infty=-(\xi'/\xi)'$ on $\Re s>\frac12$; as
$\mathcal K_\infty$ is holomorphic there, $\xi$ has no zero in that
half-plane, which is RH. Then $-(\xi'/\xi)'=\sum_kk_s(\gamma_k)$, and
comparing the double poles and their coefficients of the two meromorphic
functions shows that the multisets $\{\lambda_k\}$ and $\{\gamma_k\}$
coincide; both being enumerated in nondecreasing order,
$\lambda_k=\gamma_k$ for every $k$. Since this holds along a subsequence
of every sequence $X_n\to\infty$, the full limits exist.
\end{proof}

\begin{remark}
\Cref{prop:v8-criterion} is the precise content of the
Laguerre--P\'olya proposal in \cite{CasanaV10Notes}: the finite-cutoff
points are real by construction, so the only way the limit identity can
fail is through a defect in the \emph{positions} of the $\lambda_k$,
never through a non-real limit. A repeated imaginary part is
\emph{not} excluded by the right side of \eqref{eq:v8-criterion}: if
$\gamma_k=\gamma_{k+1}$, RH still forces $\lambda_k=\lambda_{k+1}$ and the
identity holds. Simplicity therefore has to come from the arithmetic
noncollapse of \cref{prop:correlated-v5}, exactly as in
\cref{cor:rh-rvm-v5}.
\end{remark}

"""
rep(anchor, new_sec + anchor)

# ---------- paragraph on the first-index normalization (sec:v7-complex) ----------
old_par = r"""The public notes also flag a possible lower-endpoint defect
for the naive ``DBSRH for every $k\ge1$'' statement: the
Euler phase satisfies $P_X^W(0)=0$, whereas the fixed main
term extends to $N_0(0^+)=7/8>1/2$. The formal proposal
$\tau_k^\flat=2\pi+\meas\{t\ge2\pi:
F_X^W(t)<k-\frac12\}$ is an \emph{alternative proposed
normalization}; its convergence and finite-cutoff strictness
in all indices are not established here. We therefore retain
the original definitions for existing theorems and do not
silently replace $\mathrm{DBSRH}(W)$ by its corrected variant."""
new_par = r"""The public notes also flag a lower-endpoint defect
for the naive ``DBSRH for every $k\ge1$'' statement: the
Euler phase satisfies $P_X^W(0)=0$, whereas the fixed main
term extends to $N_0(0^+)=7/8>1/2$. \Cref{prop:v8-first} makes
this precise: under the hypotheses of the transform criterion
together with convergence in measure, $\tau_1^{\Cms}$ converges to
$\gamma_1-t_\sharp$ with $t_\sharp>0$, while the points
$\tau_k^\flat=2\pi+\meas\{t\ge2\pi:F_X^W(t)<k-\frac12\}$ of
\eqref{eq:v8-flat} converge to $\gamma_k$ for every $k\ge1$.
Finite-cutoff strictness of $\tau_k^\flat$ in the initial indices is
not established here. We therefore retain the original definitions
for existing theorems and do not silently replace
$\mathrm{DBSRH}(W)$ by its corrected variant."""
rep(old_par, new_par)

# ---------- open task (5): first-index normalization now settled conditionally ----------
rep(r""" Verify the weighted zero-explicit formula including
 branch conventions, the exact IJNT/Fujii ranges, and the
 proposed first-index normalization; no claims about""",
    r""" Verify the weighted zero-explicit formula including
 branch conventions and the exact IJNT/Fujii ranges (the
 first-index normalization is settled, under the hypotheses of
 \cref{prop:v8-first}, by \eqref{eq:v8-flat}); no claims about""")

# ---------- status table ----------
rep("Fourier-box identity and parity test & Exact & Finite-interval integration by parts and endpoint terms; no simplicity claim.\\\\",
    "Fourier-box identity and parity test & Exact & Finite-interval integration by parts and endpoint terms; no simplicity claim.\\\\\n"
    "Limit phase $P_\\infty=S-r_0$ and first-index defect $t_\\sharp>0$ & Conditional on weighted-$L^2$ bound and convergence in measure & Transforms of $N_\\vartheta$ and $r_0$, Hadamard product, uniqueness via Cauchy integrals; MCOM bound for inactivity of the truncation.\\\\\n"
    "Transform criterion for $\\tau_k^\\flat$ & Proved equivalence; hypothesis open & Normality from the MCOM barriers, identity theorem, pole comparison; RH together with convergence of all $\\tau_k^\\flat$.\\\\")

# ---------- bibliography ----------
rep(r"""\bibitem{Titchmarsh}
E.~C.~Titchmarsh, \emph{The Theory of the Riemann Zeta-Function}, 2nd edition, revised by D.~R.~Heath-Brown, Oxford University Press (1986).""",
    r"""\bibitem{Titchmarsh}
E.~C.~Titchmarsh, \emph{The Theory of the Riemann Zeta-Function}, 2nd edition, revised by D.~R.~Heath-Brown, Oxford University Press (1986).

\bibitem{TitchmarshFI}
E.~C.~Titchmarsh, \emph{Introduction to the Theory of Fourier Integrals}, 2nd edition, Oxford University Press (1948); Theorem~93 (boundary values of the Cauchy integral of an $L^2$ function).

\bibitem{WhittakerWatson}
E.~T.~Whittaker and G.~N.~Watson, \emph{A Course of Modern Analysis}, 4th edition, Cambridge University Press (1927), \S12.31 (Binet's and Malmst\'en's integrals for $\log\Gamma$), \href{https://doi.org/10.1017/CBO9780511608759}{doi:10.1017/CBO9780511608759}.""")

open(dst, 'w', encoding='utf-8').write(s)
print("ok", len(s))
