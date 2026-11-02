/* tauC_v2.c  --  (v2: adds env MCAP = cap on the number of primes m_j; default = paper definition)    prime-corrected Riemann-von Mangoldt points tau_k^C
 *
 *   F(t) = N0(t) + Phi(t),   N0(t) = (t/2pi) log(t/(2 pi e)) + 7/8
 *   Phi_m(t) = -(1/pi) sum_{nu<=m} arg(1 - p_nu^{-1/2-it})
 *   Phi^o(t) = (2-t/2^j) Phi_{m_j}(t) + (t/2^j-1) Phi_{m_{j+1}}(t),  2^j<=t<=2^{j+1},  m_j=floor(2^{j/2})
 *   Phi = clip(Phi^o, -hc, hc),  hc = loglog max(t,e^e)
 *   tau_k^C = meas{ t>0 : F(t) < k-1/2 }          (IJNT paper, eq. xk_def)
 *
 * Method: F and F' are evaluated on the dyadic-aligned grid t = g*h (h = 1/HINV, power of 2),
 * cell-wise cubic Hermite interpolation, each cell split at the critical points of the cubic into
 * monotone pieces.  meas{F<lam_k} = sum over pieces:  (full length if lam_k > max)  +  (partial via bisection).
 * Full lengths are accumulated in a histogram by k_full and cumulated afterwards.
 *
 * usage:  ./tauC_v1 zeross.bin outprefix [HINV=32] [nophi]
 *         ./tauC_v1 dump g0 n [HINV]       (debug: F, F' at grid points g0..g0+n)
 * outputs: outprefix_tauC.bin, outprefix_tau.bin (float64, k=1..K)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265358979323846264338327950288
#define PIL 3.14159265358979323846264338327950288L
#define NEX 25      /* primes < 100: exact atan2 */
#define N1 168      /* primes < 1000 */
#define N2 1229     /* primes < 10000 */

static int NP;
static double *pa, *plg, *prc, *prs;  /* a=p^-1/2, log p, cos(h lg), sin(h lg) */
static long double *plgl;
static double H;
static long HINV;
static int NOPHI = 0;

static void sieve(int lim, int need) {
    char *c = calloc(lim + 1, 1);
    int *pr = malloc(sizeof(int) * (lim + 1));
    int n = 0;
    for (int i = 2; i <= lim; i++) {
        if (!c[i]) { pr[n++] = i; for (long q = (long)i * i; q <= lim; q += i) c[q] = 1; }
    }
    if (n < need) { fprintf(stderr, "sieve too small\n"); exit(1); }
    NP = need;
    pa = malloc(sizeof(double) * NP); plg = malloc(sizeof(double) * NP);
    prc = malloc(sizeof(double) * NP); prs = malloc(sizeof(double) * NP);
    plgl = malloc(sizeof(long double) * NP);
    for (int i = 0; i < NP; i++) {
        long double l = logl((long double)pr[i]);
        plgl[i] = l; plg[i] = (double)l; pa[i] = 1.0 / sqrt((double)pr[i]);
        long double ang = (long double)H * l;
        prc[i] = (double)cosl(ang); prs[i] = (double)sinl(ang);
    }
    free(c); free(pr);
}

/* weighted sums over prime indices [0,i1):  PA = sum wA_i pi*Phi_p, PB = sum wB_i pi*Phi_p ; QA,QB likewise with -pi*Phi_p' */
#define TIERW(IA, IB, R)                                                        \
    { int a_ = (IA); int b_ = (IB) < i1 ? (IB) : i1;                            \
      for (int i = a_; i < b_; i++) {                                           \
        double zr = pa[i] * wr[i], zi = pa[i] * wi[i];                          \
        double qr = zr, qi = zi, sp = zi, sd = zr;                              \
        _Pragma("GCC unroll 10")                                                \
        for (int r = 2; r <= (R); r++) {                                        \
            double t1 = qr * zr - qi * zi; qi = qr * zi + qi * zr; qr = t1;     \
            sp += qi / r; sd += qr; }                                           \
        double q_ = plg[i] * sd;                                                \
        PA += wA[i] * sp; PB += wB[i] * sp; QA += wA[i] * q_; QB += wB[i] * q_; } }

static inline void accw(int i1, const double *wr, const double *wi, const double *wA, const double *wB,
                        double *PAp, double *PBp, double *QAp, double *QBp) {
    double PA = 0, PB = 0, QA = 0, QB = 0;
    { int b_ = NEX < i1 ? NEX : i1;
      for (int i = 0; i < b_; i++) {
        double a = pa[i], zr = a * wr[i], zi = a * wi[i];
        double x = 1.0 - zr, den = x * x + zi * zi;
        double v = atan2(zi, x), q = plg[i] * (zr - a * a) / den;
        PA += wA[i] * v; PB += wB[i] * v; QA += wA[i] * q; QB += wB[i] * q; } }
    TIERW(NEX, N1, 9)
    TIERW(N1, N2, 6)
    TIERW(N2, 1 << 30, 5)
    *PAp = PA; *PBp = PB; *QAp = QA; *QBp = QB;
}

/* scheme */
static int SCH = 0; static double PARAM = 1.0;   /* 0 paper, 1 count, 2 size, 3 taper */
static double wlevel(int L, int i) {             /* W_L(i) */
    if (L < 0) return 0.0;
    if (SCH == 0) return i < (int)floor(pow(2.0, L / 2.0) + 1e-12) ? 1.0 : 0.0;
    if (SCH == 1) return i < (int)floor(PARAM * pow(2.0, L / 2.0) + 1e-12) ? 1.0 : 0.0;
    double lx = PARAM * L * log(2.0);
    if (lx <= 0) return 0.0;
    if (SCH == 2) return plg[i] <= lx ? 1.0 : 0.0;
    double w = 2.0 - plg[i] / lx; return w < 0 ? 0 : (w > 1 ? 1 : w);
}
static int support(int L) {                      /* number of indices with W_L(i) != 0 */
    int m = 0; while (m < NP && wlevel(L, m) != 0.0) m++; return m;
}
static int MCAP = 1<<30;
static inline int mj_of(int j) { int m = j < 0 ? 0 : (int)floor(pow(2.0, j / 2.0) + 1e-12); return m < MCAP ? m : MCAP; }

static long nclip = 0;

/* evaluate F and F' at t = (g0 + k) h, k = 0..nc, all in block j */
static void eval_chunk(int j, long g0, int nc, double *Fv, double *Dv, double *wr, double *wi, long *clipcount,
                       double *wA, double *wB) {
    int m2 = 0;
    if (!NOPHI) {
        int s1 = support(j), s2 = support(j + 1); m2 = s1 > s2 ? s1 : s2;
        for (int i = 0; i < m2; i++) { double a = wlevel(j, i), b = wlevel(j + 1, i); wA[i] = a; wB[i] = b - a; }
    }
    double t0 = (double)g0 * H;
    for (int i = 0; i < m2; i++) {
        long double th = (long double)t0 * plgl[i];
        th = fmodl(th, 2.0L * PIL);
        wr[i] = (double)cosl(th); wi[i] = -(double)sinl(th);
    }
    double pj = ldexp(1.0, j);
    const double ee = exp(exp(1.0));
    for (int k = 0; k <= nc; k++) {
        double t = (double)(g0 + k) * H;
        double P1 = 0, Q1 = 0, P2 = 0, Q2 = 0;
        if (m2 > 0) accw(m2, wr, wi, wA, wB, &P1, &P2, &Q1, &Q2);
        double u = t / pj;
        double Phio = (P1 + (u - 1.0) * P2) / PI;
        double dPhio = (P2 / pj - Q1 - (u - 1.0) * Q2) / PI;
        double hc, dhc;
        if (t <= ee) { hc = 1.0; dhc = 0; } else { double lt = log(t); hc = log(lt); dhc = 1.0 / (t * lt); }
        double Ph = Phio, dPh = dPhio;
        if (Phio > hc) { Ph = hc; dPh = dhc; (*clipcount)++; }
        else if (Phio < -hc) { Ph = -hc; dPh = -dhc; (*clipcount)++; }
        double x = t / (2.0 * PI);
        double lx = log(x);
        double N0 = x * (lx - 1.0) + 0.875;
        double dN0 = lx / (2.0 * PI);
        Fv[k] = N0 + Ph; Dv[k] = dN0 + dPh;
        if (k < nc) {
            for (int i = 0; i < m2; i++) {
                double a = wr[i] * prc[i] + wi[i] * prs[i];
                double b = wi[i] * prc[i] - wr[i] * prs[i];
                wr[i] = a; wi[i] = b;
            }
        }
    }
}

static inline int kfloor(double v) { return (int)floor(v + 0.5); }

/* process one monotone piece [sa,sb] of the cell cubic (offset F0, scale h) */
static void piece(double F0, double c1, double c2, double c3, double sa, double sb,
                  double Fa, double Fb, int K, double *hist, double *add) {
    if (sb - sa <= 1e-15) return;
    double len = (sb - sa) * H;
    double lo = Fa < Fb ? Fa : Fb, hi = Fa < Fb ? Fb : Fa;
    int kf = kfloor(hi) + 1;
    if (kf < 1) kf = 1;
    if (kf <= K) hist[kf] += len;
    int k1 = kfloor(lo) + 1, k2 = kfloor(hi);
    if (k1 < 1) k1 = 1;
    if (k2 > K) k2 = K;
    int inc = Fa < Fb;
    for (int k = k1; k <= k2; k++) {
        double L = (k - 0.5) - F0;
        double a = sa, b = sb;
        for (int it = 0; it < 60; it++) {
            double m = 0.5 * (a + b);
            double pm = ((c3 * m + c2) * m + c1) * m;
            if ((pm < L) == (inc ? 1 : 0)) a = m; else b = m;
            /* inc: p<L  => root to the right => a=m ; dec: p<L => root to the left => b=m */
        }
        double s = 0.5 * (a + b);
        double part = inc ? (s - sa) : (sb - s);
        add[k] += part * H;
    }
}

typedef struct { int j; long g0; int nc; } Task;

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "dump")) {
        long g0 = atol(argv[2]); int n = atoi(argv[3]);
        HINV = argc > 4 ? atol(argv[4]) : 32; H = 1.0 / HINV;
        sieve(60000, 4000);
        double t = g0 * H; int j = (int)floor(log2(t)); if (j < -2) j = -2;
        double *Fv = malloc(8 * (n + 1)), *Dv = malloc(8 * (n + 1)), *wr = malloc(8 * NP), *wi = malloc(8 * NP);
        long cc = 0;
        double *wA = malloc(8 * NP), *wB = malloc(8 * NP);
        eval_chunk(j, g0, n, Fv, Dv, wr, wi, &cc, wA, wB);
        for (int k = 0; k <= n; k++) printf("%.17g %.17g %.17g\n", (g0 + k) * H, Fv[k], Dv[k]);
        return 0;
    }
    if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
    const char *zf = argv[1], *pref = argv[2];
    HINV = argc > 3 ? atol(argv[3]) : 32; H = 1.0 / HINV;
    NOPHI = (argc > 4 && !strcmp(argv[4], "nophi"));
    if (getenv("MCAP")) MCAP = atoi(getenv("MCAP"));
    { const char *sc = getenv("SCHEME"); if (sc) { if (!strcmp(sc,"count")) SCH=1; else if (!strcmp(sc,"size")) SCH=2; else if (!strcmp(sc,"taper")) SCH=3; }
      if (getenv("PARAM")) PARAM = atof(getenv("PARAM")); }
    FILE *f = fopen(zf, "rb"); fseek(f, 0, SEEK_END); long nb = ftell(f); fseek(f, 0, SEEK_SET);
    int K = (int)(nb / 8); double *gam = malloc(nb); if (fread(gam, 8, K, f) != (size_t)K) return 2; fclose(f);

    /* Tmax: N0(Tmax) = K + 12 */
    long double lam = K + 12.0L, T = 2.0L * PIL * (K + 12.0L) / logl((K + 12.0L)) + 50;
    for (int it = 0; it < 100; it++) {
        long double x = T / (2 * PIL), lx = logl(x);
        long double F = x * (lx - 1) + 0.875L - lam, dF = lx / (2 * PIL);
        T -= F / dF;
    }
    double Tmax = (double)T;
    long gmax = (long)ceil(Tmax / H);
    int jmax = (int)floor(log2(gmax * H));
    sieve(8000000, 500000);            /* provisional: all primes < 8e6 */
    int mneed = support(jmax + 1) + 2;
    if (mneed > NP) { fprintf(stderr, "need more primes\n"); return 3; }
    { int keep = mneed; NP = keep; }
    fprintf(stderr, "K=%d  Tmax=%.6f  gmax=%ld  h=1/%ld  max primes=%d  threads=%d  nophi=%d  scheme=%d param=%g\n", K, Tmax, gmax, HINV, mneed, omp_get_max_threads(), NOPHI, SCH, PARAM);

    /* tasks */
    int CH = 2048;
    long ntmax = (gmax / CH) + 200; Task *tasks = malloc(sizeof(Task) * ntmax); long nt = 0;
    { /* block j=-2 : g from 1 to 0.5/h */
        long ga = 1, gb = (long)(0.5 * HINV);
        for (long g = ga; g < gb; g += CH) { long e = g + CH < gb ? g + CH : gb; tasks[nt].j = -2; tasks[nt].g0 = g; tasks[nt].nc = (int)(e - g); nt++; }
        for (int j = -1; j <= jmax; j++) {
            long ga = (long)ldexp((double)HINV, j), gb = 2 * ga; if (gb > gmax) gb = gmax;
            if (ga >= gb) break;
            for (long g = ga; g < gb; g += CH) { long e = g + CH < gb ? g + CH : gb; tasks[nt].j = j; tasks[nt].g0 = g; tasks[nt].nc = (int)(e - g); nt++; }
        }
    }
    int nth = omp_get_max_threads();
    double **hist = malloc(sizeof(double *) * nth), **add = malloc(sizeof(double *) * nth);
    for (int i = 0; i < nth; i++) { hist[i] = calloc(K + 3, 8); add[i] = calloc(K + 3, 8); }
    long totclip = 0; double minmargin = 1e9;
    long done = 0;

    #pragma omp parallel
    {
        int th = omp_get_thread_num();
        double *Fv = malloc(8 * (CH + 2)), *Dv = malloc(8 * (CH + 2)), *wr = malloc(8 * NP), *wi = malloc(8 * NP), *wA = malloc(8 * NP), *wB = malloc(8 * NP);
        double *hs = hist[th], *ad = add[th];
        long cc = 0; double mm = 1e9;
        #pragma omp for schedule(dynamic, 1)
        for (long q = 0; q < nt; q++) {
            Task *tk = &tasks[q];
            eval_chunk(tk->j, tk->g0, tk->nc, Fv, Dv, wr, wi, &cc, wA, wB);
            for (int c = 0; c < tk->nc; c++) {
                double F0 = Fv[c], F1 = Fv[c + 1];
                double A = H * Dv[c], B = H * Dv[c + 1], Dl = F1 - F0;
                double c1 = A, c2 = 3 * Dl - 2 * A - B, c3 = A + B - 2 * Dl;
                double rts[4]; int nr = 0; rts[nr++] = 0.0;
                double cr[2]; int ncr = 0;
                if (fabs(c3) < 1e-13) { if (fabs(c2) > 1e-14) cr[ncr++] = -c1 / (2 * c2); }
                else { double disc = 4 * c2 * c2 - 12 * c3 * c1;
                       if (disc > 0) { double sq = sqrt(disc); double r1 = (-2 * c2 - sq) / (6 * c3), r2 = (-2 * c2 + sq) / (6 * c3);
                                       if (r1 > r2) { double tt = r1; r1 = r2; r2 = tt; } cr[ncr++] = r1; cr[ncr++] = r2; } }
                for (int i = 0; i < ncr; i++) if (cr[i] > 0 && cr[i] < 1) rts[nr++] = cr[i];
                rts[nr++] = 1.0;
                double sprev = rts[0], Fprev = F0;
                for (int i = 1; i < nr; i++) {
                    double s = rts[i];
                    double Fs = (i == nr - 1) ? F1 : F0 + ((c3 * s + c2) * s + c1) * s;
                    piece(F0, c1, c2, c3, sprev, s, Fprev, Fs, K, hs, ad);
                    sprev = s; Fprev = Fs;
                }
            }
            if (tk->j == jmax && tk->g0 + tk->nc >= gmax) { /* margin check at the end */
                double m = Fv[tk->nc] - (K - 0.5); if (m < mm) mm = m;
            }
            #pragma omp atomic
            done++;
            if (th == 0 && (done % 2000) == 0) fprintf(stderr, "  tasks %ld / %ld\n", done, nt);
        }
        #pragma omp critical
        { totclip += cc; if (mm < minmargin) minmargin = mm; }
        free(Fv); free(Dv); free(wr); free(wi);
    }
    /* region [0, g=1]: F < lam for all k>=2 */
    double *tauC = malloc(8 * (K + 1)), *tau = malloc(8 * (K + 1));
    long double cum = 0, base = 1.0L * H;  /* measure of [0,h] */
    long double *Hs = calloc(K + 3, sizeof(long double)), *As = calloc(K + 3, sizeof(long double));
    for (int i = 0; i < nth; i++) for (int k = 0; k <= K + 1; k++) { Hs[k] += hist[i][k]; As[k] += add[i][k]; }
    Hs[2] += base;
    for (int k = 1; k <= K; k++) { cum += Hs[k]; tauC[k - 1] = (double)(cum + As[k]); }
    /* Newton for tau_k */
    for (int k = 1; k <= K; k++) {
        long double l = k - 0.5L, t = 2 * PIL * (k + 1) / logl(k + 3.0L) + 20;
        for (int it = 0; it < 100; it++) { long double x = t / (2 * PIL), lx = logl(x); long double Fv = x * (lx - 1) + 0.875L - l; t -= Fv / (lx / (2 * PIL)); }
        tau[k - 1] = (double)t;
    }
    char fn[512];
    snprintf(fn, 512, "%s_tauC.bin", pref); f = fopen(fn, "wb"); fwrite(tauC, 8, K, f); fclose(f);
    snprintf(fn, 512, "%s_tau.bin", pref); f = fopen(fn, "wb"); fwrite(tau, 8, K, f); fclose(f);
    /* stats */
    long double s1 = 0, s2 = 0, u1 = 0, u2 = 0; double maxd = 0;
    for (int k = 0; k < K; k++) { double d = tauC[k] - gam[k], e = tau[k] - gam[k]; s1 += d; s2 += d * d; u1 += e; u2 += e * e; if (fabs(d) > maxd) maxd = fabs(d); }
    printf("scheme=%d param=%g  K=%d  h=1/%ld  nophi=%d  clipped_points=%ld  min margin F(Tmax)-(K-1/2)=%.4f\n", SCH, PARAM, K, HINV, NOPHI, totclip, minmargin);
    printf("mean(tauC-gamma)=%.6Lf rms=%.6Lf max|.|=%.4f   mean(tau-gamma)=%.6Lf rms=%.6Lf\n", s1 / K, sqrtl(s2 / K), maxd, u1 / K, sqrtl(u2 / K));
    printf("tauC_1=%.8f tauC_2=%.8f tauC_K=%.8f  tau_K=%.8f gamma_K=%.8f\n", tauC[0], tauC[1], tauC[K - 1], tau[K - 1], gam[K - 1]);
    return 0;
}
