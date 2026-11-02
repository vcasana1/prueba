/* subnivel_rel_v8.c -- taucs en una ventana [T0, T0+W] a cualquier altura, con N0 relativo:
 *   F_rel(t) = F(t) - K0 = fr0 + dN0(u) + clip(Phi^o(t), +-loglog t),  u = t - T0,
 *   dN0(u) = [u log(T0/(2 pi e)) + (T0+u) log1p(u/T0)] / (2 pi),   N0' = log(t/2pi)/(2pi).
 *   tauc_k = T0 + meas{t in [T0,T0+W] : F_rel < (k-K0) - 1/2}, valido si fr0 + hc(T0) < (k-K0) - 1/2 < dN0(W)+fr0 - hc(T0+W).
 *   Mismo motor de celdas (Hermite cubico, piezas monotonas, biseccion) que tauC_v1.c.
 * usage: ./subnivel_rel_v8 phipref T0 W h_inv K0 fr0 out.txt      (out: k  tauc_k-T0)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979323846
static double H;
static void piece(double F0, double c1, double c2, double c3, double sa, double sb, double Fa, double Fb, int K, double *hist, double *add) {
    if (sb - sa <= 1e-15) return;
    double len = (sb - sa) * H, lo = Fa < Fb ? Fa : Fb, hi = Fa < Fb ? Fb : Fa;
    int kf = (int)floor(hi + 0.5) + 1; if (kf < 0) kf = 0; if (kf <= K) hist[kf] += len;
    int k1 = (int)floor(lo + 0.5) + 1, k2 = (int)floor(hi + 0.5); if (k1 < 0) k1 = 0; if (k2 > K) k2 = K;
    int inc = Fa < Fb;
    for (int k = k1; k <= k2; k++) {
        double L = (k - 0.5) - F0, a = sa, b = sb;
        for (int it = 0; it < 60; it++) { double m = 0.5 * (a + b), pm = ((c3 * m + c2) * m + c1) * m; if ((pm < L) == inc) a = m; else b = m; }
        double s = 0.5 * (a + b); add[k] += (inc ? (s - sa) : (sb - s)) * H;
    }
}
int main(int argc, char **argv) {
    const char *pp = argv[1]; double T0 = atof(argv[2]), W = atof(argv[3]); long HI = atol(argv[4]); long long K0 = atoll(argv[5]); double fr0 = atof(argv[6]);
    H = 1.0 / HI; long M = (long)ceil(W * HI) + 1;
    double *phi = malloc(8 * M), *dphi = malloc(8 * M); char fn[512];
    snprintf(fn, 512, "%s_phi.bin", pp); FILE *f = fopen(fn, "rb"); if (fread(phi, 8, M, f) != (size_t)M) return 2; fclose(f);
    snprintf(fn, 512, "%s_dphi.bin", pp); f = fopen(fn, "rb"); if (fread(dphi, 8, M, f) != (size_t)M) return 2; fclose(f);
    double l0 = log(T0 / (2 * PI * exp(1.0)));
    double *Fv = malloc(8 * M), *Dv = malloc(8 * M); long clip = 0;
    for (long j = 0; j < M; j++) {
        double u = j * H, t = T0 + u, hc = log(log(t)), ph = phi[j], dp = dphi[j];
        if (ph > hc) { ph = hc; dp = 1.0 / (t * log(t)); clip++; } else if (ph < -hc) { ph = -hc; dp = -1.0 / (t * log(t)); clip++; }
        Fv[j] = fr0 + (u * l0 + t * log1p(u / T0)) / (2 * PI) + ph;
        Dv[j] = log(t / (2 * PI)) / (2 * PI) + dp;
    }
    int K = (int)ceil(Fv[M - 1]) + 5; double *hist = calloc(K + 3, 8), *add = calloc(K + 3, 8);
    for (long c = 0; c + 1 < M; c++) {
        double F0 = Fv[c], F1 = Fv[c + 1], A = H * Dv[c], B = H * Dv[c + 1], Dl = F1 - F0;
        double c1 = A, c2 = 3 * Dl - 2 * A - B, c3 = A + B - 2 * Dl, r[4]; int nr = 0; r[nr++] = 0;
        double cr[2]; int ncr = 0;
        if (fabs(c3) < 1e-13) { if (fabs(c2) > 1e-14) cr[ncr++] = -c1 / (2 * c2); }
        else { double d = 4 * c2 * c2 - 12 * c3 * c1; if (d > 0) { double q = sqrt(d), r1 = (-2 * c2 - q) / (6 * c3), r2 = (-2 * c2 + q) / (6 * c3); if (r1 > r2) { double tt = r1; r1 = r2; r2 = tt; } cr[ncr++] = r1; cr[ncr++] = r2; } }
        for (int i = 0; i < ncr; i++) if (cr[i] > 0 && cr[i] < 1) r[nr++] = cr[i];
        r[nr++] = 1;
        double sp = 0, Fp = F0;
        for (int i = 1; i < nr; i++) { double s = r[i], Fs = (i == nr - 1) ? F1 : F0 + ((c3 * s + c2) * s + c1) * s; piece(F0, c1, c2, c3, sp, s, Fp, Fs, K, hist, add); sp = s; Fp = Fs; }
    }
    double hcA = log(log(T0)), hcB = log(log(T0 + W));
    int k1 = (int)ceil(fr0 + hcA + 0.5 + 1e-9), k2 = (int)floor(Fv[M - 1] - (Fv[M - 1] - (fr0 + ((W) * l0 + (T0 + W) * log1p(W / T0)) / (2 * PI))) - hcB + 0.5 - 1e-9);
    f = fopen(argv[7], "w"); long double cum = 0;
    for (int k = 0; k <= K && k <= k2; k++) { cum += hist[k]; if (k >= k1) fprintf(f, "%lld %.9Lf\n", K0 + k, cum + add[k]); }
    fclose(f);
    fprintf(stderr, "M=%ld clip=%ld  niveles validos k-K0 en [%d,%d]\n", M, clip, k1, k2);
    return 0;
}
