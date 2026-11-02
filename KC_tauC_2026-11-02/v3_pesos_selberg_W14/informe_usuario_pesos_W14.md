## Primer resultado: el peso de Selberg supera a los cortes agudos y a los suavizados polinómicos

He utilizado `zeross(6).bin` para realizar una primera prueba con **9.001 puntos situados entre ceros verdaderos**, desde \(k=1000\) hasta \(k=10000\).

He comparado directamente las fases primas con \(S(t)\), manteniendo el mismo presupuesto máximo de primos \(p^r\le X^2\).

| Peso \(W\) | RMS para \(X=50\) | RMS para \(X=100\) | RMS para \(X=200\) |
|---|---:|---:|---:|
| Corte agudo en \(X\) | 0,072191 | 0,064875 | 0,057809 |
| **Selberg, lineal entre \(X\) y \(X^2\)** | **0,014706** | **0,010488** | **0,008255** |
| Transición cúbica | 0,021387 | 0,016386 | 0,011918 |
| Transición quíntica | 0,025923 | 0,020788 | 0,015442 |
| Corte agudo en \(X^2\) | 0,045620 | 0,038035 | 0,032975 |

**El peso de Selberg gana claramente en los tres cortes.**

Hay una limitación: esta prueba mide el error de fase en puntos intermedios, no las distancias \(\tau_k^{\mathrm C}-\gamma_k\). Además, todavía no incorpora el compensador exacto del polo \(\Pi_X^W\). Por tanto, es una prueba exploratoria del perfil, no una verificación de DBSRH.

Estos datos aconsejan conservar el perfil lineal de Selberg como referencia y buscar modificaciones pequeñas, en lugar de sustituirlo directamente por un suavizado más intenso.

## 2. Un peso nuevo que mejora ligeramente a Selberg

A partir del experimento anterior, he probado esta familia:

\[
\boxed{
W_b(u)=
\begin{cases}
1,&0\le u\le1,\\[2mm]
1-v+bv(1-v),&1<u<2,\quad v=u-1,\\[2mm]
0,&u\ge2.
\end{cases}
}
\]

El peso de Selberg corresponde a \(b=0\).

El candidato que seleccionaría para continuar es

\[
\boxed{b=\frac14.}
\]

Conserva todos los primos pequeños con peso 1 y aumenta moderadamente el peso de los primos situados entre \(X\) y \(X^2\).

### Validación independiente: índices 5001–10000

He ajustado la familia utilizando los índices 1000–5000 y evaluado el resultado en los 5000 siguientes.

| Corte \(X\) | Selberg | Nuevo peso \(b=1/4\) | Reducción del RMS de fase |
|---:|---:|---:|---:|
| 50 | 0,016596 | **0,016045** | 3,3 % |
| 100 | 0,011513 | **0,011175** | 2,9 % |
| 200 | 0,009001 | **0,008288** | 7,9 % |

La mejora se conserva en los tres cortes.

Estos resultados apoyan la elección \(b=1/4\) como candidato numérico, pero no demuestran que sea óptimo.

## 3. Comprobación directa de las ordenadas por subniveles

También he ejecutado la construcción ordinal, con \(X=200\), utilizando ambos pesos.

| Índices | RMS Selberg | RMS nuevo peso | Incumplimientos con nuevo peso |
|---|---:|---:|---:|
| 4700–4800 | 0,017920 | **0,017149** | 1 |
| 5001–5500 | 0,015017 | **0,014159** | **0** |
| 6650–6760 | 0,017470 | **0,016542** | 1 |

El nuevo peso mejora también las ordenadas.

Sin embargo, **no resuelve las parejas excepcionalmente próximas 4765–4766 y 6709–6710**. Reduce sus errores, pero siguen incumpliendo la condición de media separación.

Estas pruebas de subniveles siguen siendo exploratorias porque aún falta incorporar el compensador exacto del polo.

---

# 4. Una propiedad matemática útil del nuevo peso

El candidato tiene una representación particularmente sencilla.

Para \(1\le v\le2\), definamos

\[
q(v)=\frac v2+\frac14.
\]

Entonces

\[
q(v)>0,
\qquad
\int_1^2q(v)\,dv=1,
\]

y

\[
\boxed{
W_{1/4}(u)=
\int_1^2q(v)\,\mathbf1_{\{u\le v\}}\,dv.
}
\]

Es decir, el nuevo peso es una **media positiva de cortes agudos**, con mayor importancia para los cortes grandes.

Esto permite trasladar las fórmulas explícitas al nuevo perfil mediante integración:

\[
\boxed{
\Pi_X^{W_{1/4}}(t)
=
\int_1^2
q(v)\,\Pi_{X^v}^{\mathrm{agudo}}(t)\,dv.
}
\]

La identidad es exacta si los compensadores agudos se definen mediante la misma fórmula explícita y la misma elección de rama.

Esta representación tiene dos ventajas matemáticas: conserva la estructura aritmética de los cortes originales y permite estudiar el término del polo sin introducir compensaciones arbitrarias.

---

# 5. Sobre DBSRH y su conexión con RH

Hay una parte de tu argumento que considero correcta y otra que necesita una demostración nueva.

### La implicación determinista

Si

\[
\Phi_X^W(t)\longrightarrow S(t)
\]

para casi todo \(t>0\), y el recorte satisface \(H(t)\ge|S(t)|\), entonces la convergencia de las medidas de subnivel se obtiene mediante convergencia dominada.

La justificación es que

\[
N_0(t)+S(t)
=
N(t)-\frac{\varpi(t)}{\pi},
\]

y

\[
\left|\frac{\varpi(t)}{\pi}\right|<\frac12.
\]

Como \(N(t)\) es entero fuera de las ordenadas, los indicadores de subnivel coinciden casi por todas partes.

La cota del recorte proporciona además un soporte compacto uniforme para cada nivel fijo \(k-\frac12\).

Hay que extender la definición de \(H\) a \(0<t<1\), por ejemplo mediante \(\log\max(t,1)\), y comprobar la cota correspondiente en ese intervalo.

### La implicación hacia RH

La afirmación

\[
\mathrm{DBSRH}(W)\Longrightarrow\mathrm{RH}
\]

**todavía no está demostrada**.

Un cero fuera de la recta puede introducir términos crecientes en la fórmula explícita, pero falta demostrar que esos términos impiden la convergencia de las medidas de subnivel después del recorte.

El crecimiento de la fase no implica automáticamente el fracaso de su medida de subnivel.

Este es el obstáculo exacto. No es una objeción a trabajar con distancias ordinales: es la propiedad de esas distancias que todavía necesita demostrarse.

También hay que verificar rigurosamente la implicación contraria para el peso elegido, utilizando la fórmula explícita completa. Los trabajos de Conrad y Kaneko–Koyama–Kurokawa ofrecen antecedentes relevantes sobre convergencia de productos de Euler en la línea crítica, pero no establecen directamente la equivalencia DBSRH que propones. [Keith Conrad's Home Page](https://kconrad.math.uconn.edu/articles/eulerprod.pdf?utm_source=chatgpt.com)

---

# 6. Prueba de Cauchy sin utilizar las ordenadas verdaderas

He calculado también las diferencias entre las fases correspondientes a \(X=50,100,200\), utilizando únicamente los puntos

\[
t=1000,1001,\ldots,10000.
\]

| Peso | RMS \(\Phi_{100}-\Phi_{50}\) | RMS \(\Phi_{200}-\Phi_{100}\) |
|---|---:|---:|
| Selberg | 0,042220 | 0,035216 |
| Nuevo \(W_{1/4}\) | 0,043177 | 0,035732 |

Las diferencias disminuyen al aumentar el corte en ambos casos.

El nuevo peso mejora la aproximación a \(S(t)\) en los experimentos anteriores, pero no mejora esta medida particular de Cauchy.

Son propiedades diferentes.

Además, tres valores finitos de \(X\) no permiten afirmar convergencia cuando \(X\to\infty\).

---

## 7. Archivos y programas ejecutados

He reunido las funciones, los programas y los resultados en un único paquete:

**DBSRH_pesos_pruebas.zip[Descargar ZIP: DBSRH — pesos y pruebas numéricas](sandbox:/mnt/data/DBSRH_pesos_pruebas.zip)**

También puedes acceder por separado a:

- :chatgpt-content-reference{index="3"}[Función Python del nuevo peso \(W_{1/4}\)](sandbox:/mnt/data/DBSRH_pesos_prueba/peso_candidato_W.py)
- :chatgpt-content-reference{index="4"}[Programa de comparación de fases](sandbox:/mnt/data/probar_pesos_DBSRH.py)
- :chatgpt-content-reference{index="5"}[Programa de cálculo de ordenadas por subniveles](sandbox:/mnt/data/DBSRH_pesos_prueba/comparar_subniveles_W.py)
- :chatgpt-content-reference{index="6"}[Resultados de validación de los pesos](sandbox:/mnt/data/DBSRH_pesos_prueba/perturbaciones_selberg.csv)
- :chatgpt-content-reference{index="7"}[Informe matemático y limitaciones](sandbox:/mnt/data/DBSRH_pesos_prueba/RESULTADOS_Y_LIMITES.txt)

### Conclusión

Mi candidato actual es

\[
\boxed{
W_*(u)=
\begin{cases}
1,&0\le u\le1,\\
1-(u-1)+\dfrac14(u-1)(2-u),&1<u<2,\\
0,&u\ge2.
\end{cases}
}
\]

Tiene una definición analítica sencilla, una representación positiva mediante cortes agudos y proporciona mejoras numéricas frente al peso de Selberg.

**Pero el resultado que merece convertirse en teorema no es esa mejora porcentual.**

Es determinar si la convergencia ordinal

\[
\boxed{
\tau_k^{\mathrm C}(X)\longrightarrow\gamma_k
\qquad(X\to\infty)
}
\]

puede demostrarse, o descartarse, utilizando únicamente la estructura de las fases primas y las medidas de subnivel.

Esa es la cuestión central de DBSRH.