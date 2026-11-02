# taucs_amplitude_obstruction_v8 (2 nov 2026)

Cambios respecto a v7 (todos aplicados por `patch_v8.py` sobre el v7 del autor; nada más se ha tocado):

1. Palabras prohibidas: "ordinate(s)" -> "imaginary part(s)" / "points" (10 apariciones). Título de sección sin número de versión.
2. Dos erratas LaTeX: `\\cref{sec:identify-v4}` y `P'olya`.
3. Nueva Remark 24.5: la tasa de transporte T/(log T)^3 no necesita la cota asintótica MCOM; basta la cota explícita y el sexto momento de Selberg.
4. Nueva sección 27 "Identification of the limit phase, and the first index":
   - Lema 27.1: transformadas de N_theta y de r_0 (T[r_0] = E), identidad Q = -(xi'/xi)' - E - T[N_theta].
   - Lema 27.2: unicidad de la transformada en L^2 ponderado (Plemelj).
   - Prop. 27.3: bajo normalidad L^2 + convergencia en medida: RH, P_inf = S - r_0, tau_k -> gamma_k (k>=2) y tau_1 -> gamma_1 - t_sharp con t_sharp > 0; corrección tau^flat medido desde 2pi.
   - Prop. 27.5: criterio por transformada para tau^flat: identidad límite en un conjunto con punto de acumulación <=> RH y tau_k^flat -> gamma_k para todo k.
5. Párrafo de la sección "Complex signed-pair kernel" sobre tau^flat remitido a la Prop. 27.3; tarea abierta (5) actualizada; dos filas nuevas en la tabla de estado; filas nuevas en la tabla de notación; dos referencias nuevas (Titchmarsh, Fourier Integrals; Whittaker-Watson).
