/* tauC-c.c  (Claude) -- tau_k^C = meas{t>0 : F(t) < k-1/2},  F = N0 + clip_h(Phi^o)
 * N0(t) = (t/2pi) log(t/(2 pi e)) + 7/8
 * Phi_m(t) = -(1/pi) sum_{nu<=m} arg(1 - p_nu^{-1/2-it})
 * Phi^o(t) = (2-t/2^j) Phi_{m_j} + (t/2^j-1) Phi_{m_{j+1}},  2^j<=t<2^{j+1},  m_j = floor(2^{j/2})
 * clip level hc = loglog max(t,e^e).
 * Method: F and F' on the grid t=g/HINV (t>=1), cubic Hermite per cell, split at critical points,
 * exact measure of {F<lambda_k} by root finding.  t<1 handled by direct sampling.
 * usage: ./tauC-c zeross.bin out_prefix [HINV=32] [nophi]
 * Also writes out_prefix_newton-c.bin: Newton solution of N0(t)=k-1/2 started from the file given as argv[5] (zLamb.bin), if present.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define PI 3.14159265358979323846
#define PIL 3.14159265358979323846264338327950288L
static int NP; static double *pa,*plg,*prc,*prs; static long double *plgl;
static long HINV; static double H; static int NOPHI=0;
static long isqrt_l(long x){ long r=(long)sqrt((double)x); while(r*r>x) r--; while((r+1)*(r+1)<=x) r++; return r; }
static long mj(int j){ if(j<0) return 0; if(j%2==0) return 1L<<(j/2); return isqrt_l(1L<<j); }
static void sieve(int lim,int need){
  char*c=calloc(lim+1,1); int*pr=malloc(sizeof(int)*(lim+1)); int n=0;
  for(int i=2;i<=lim;i++){ if(!c[i]){pr[n++]=i; for(long q=(long)i*i;q<=lim;q+=i)c[q]=1;} }
  if(n<need){fprintf(stderr,"sieve too small %d<%d\n",n,need);exit(1);}
  NP=need; pa=malloc(8*NP);plg=malloc(8*NP);prc=malloc(8*NP);prs=malloc(8*NP);plgl=malloc(sizeof(long double)*NP);
  for(int i=0;i<NP;i++){ long double l=logl((long double)pr[i]); plgl[i]=l; plg[i]=(double)l; pa[i]=1.0/sqrt((double)pr[i]);
    long double ang=(long double)H*l; prc[i]=(double)cosl(ang); prs[i]=(double)sinl(ang);} free(c);free(pr);
}
/* arg part and derivative part of one prime at state w=e^{-i t lg}: returns arg(1-z) and lg*Re(z/(1-z)) */
static inline void onep(int i,double wr,double wi,double*Ap,double*Qp){
  double a=pa[i], zr=a*wr, zi=a*wi;   /* z = a w, w=e^{-i theta}: Im z = -a sin theta */
  if(a>0.1){ /* exact */
    double x=1.0-zr; double den=x*x+zi*zi;
    *Ap=atan2(zi,x)*0+atan2(-zi,x)*0; /* placeholder, replaced below */
    *Ap=atan2(zi,x); /* arg(1-z)=atan2(-Im z, 1-Re z); Im z=zi => -zi */
    *Ap=atan2(-zi,x);
    *Qp=plg[i]*(zr-a*a)/den;
  } else {
    int R=(a>0.0316)?8:5;
    double qr=zr,qi=zi,s=-zi,d=zr;  /* arg(1-z) = -Im log... = sum Im(z^r)/r with sign: log(1-z)=-sum z^r/r, arg = -sum Im(z^r)/r */
    for(int r=2;r<=R;r++){ double t1=qr*zr-qi*zi; qi=qr*zi+qi*zr; qr=t1; s+=-qi/r; d+=qr; }
    *Ap=s; *Qp=plg[i]*d;
  }
}
static long double N0l(long double t){ return t/(2*PIL)*logl(t/(2*PIL*expl(1.0L)))+7.0L/8; }
static long double dN0l(long double t){ return logl(t/(2*PIL))/(2*PIL); }
/* direct point evaluation of F for t<1 (and for tests) */
static double F_direct(double t,double*dF){
  int j=(int)floor(log2(t)); long m0=mj(j),m1=mj(j+1);
  double A=0,B=0,dA=0,dB=0;
  for(long i=0;i<m1&&i<NP;i++){ double th=(double)((long double)t*plgl[i]); double wr=cos(th),wi=-sin(th),ap,q; onep((int)i,wr,wi,&ap,&q);
    if(i<m0){A+=ap;dA+=q;} B+=ap;dB+=q; }
  A*=-1/PI;B*=-1/PI;dA*=-1/PI;dB*=-1/PI;
  double p2=ldexp(1.0,j); double w1=2-t/p2,w2=t/p2-1;
  double ph=w1*A+w2*B, dph=w1*dA+w2*dB+(B-A)/p2;
  double tt=t>exp(exp(1.0))?t:exp(exp(1.0)); double hc=log(log(tt)); double dhc=(t>exp(exp(1.0)))?1.0/(t*log(t)):0;
  if(NOPHI){ph=0;dph=0;}
  if(ph>hc){dph=dhc;ph=hc;} else if(ph<-hc){dph=-dhc;ph=-hc;}
  if(dF)*dF=(double)dN0l(t)+dph; return (double)N0l(t)+ph;
}
static double *full_g,*part_g; static long K;
static inline double pcub(const double*c,double s){ return ((c[0]*s+c[1])*s+c[2])*s+c[3]; }
static void piece(const double*c,double sa,double sb,double h,double*full,double*part){
  double pa_=pcub(c,sa),pb_=pcub(c,sb); int inc=pb_>=pa_; double lo=inc?pa_:pb_,hi=inc?pb_:pa_;
  double L=(sb-sa)*h; long kf=(long)floor(hi+0.5)+1; if(kf<1)kf=1; if(kf<=K) full[kf]+=L;
  long k0=(long)floor(lo+0.5)+1,k1=(long)floor(hi+0.5); if(k0<1)k0=1;
  for(long k=k0;k<=k1&&k<=K;k++){ double lam=k-0.5; double a=sa,b=sb;
    for(int it=0;it<60;it++){ double m=0.5*(a+b); double v=pcub(c,m); if((v<lam)==inc) a=m; else b=m; }
    double s=0.5*(a+b); part[k]+= inc?(s-sa)*h:(sb-s)*h; }
}
static void cell(double F0,double F1,double d0,double d1,double h,double*full,double*part){
  double m0=h*d0,m1=h*d1; double c[4]; c[0]=2*(F0-F1)+m0+m1; c[1]=3*(F1-F0)-2*m0-m1; c[2]=m0; c[3]=F0;
  double sp[4]; int n=0; sp[n++]=0;
  double A=3*c[0],B=2*c[1],C=c[2];
  if(fabs(A)>1e-300){ double D=B*B-4*A*C; if(D>0){ double sq=sqrt(D); double r1=(-B-sq)/(2*A),r2=(-B+sq)/(2*A); if(r1>r2){double t=r1;r1=r2;r2=t;}
      if(r1>0&&r1<1)sp[n++]=r1; if(r2>0&&r2<1)sp[n++]=r2; } }
  else if(fabs(B)>1e-300){ double r=-C/B; if(r>0&&r<1)sp[n++]=r; }
  sp[n++]=1;
  for(int i=0;i+1<n;i++) piece(c,sp[i],sp[i+1],h,full,part);
}
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage\n");return 1;}
  HINV=argc>3?atol(argv[3]):32; H=1.0/HINV; if(argc>4&&!strcmp(argv[4],"nophi"))NOPHI=1;
  FILE*f=fopen(argv[1],"rb"); fseek(f,0,SEEK_END); long nb=ftell(f); fseek(f,0,SEEK_SET); K=nb/8; double*z=malloc(nb); fread(z,8,K,f); fclose(f);
  double Tmax=z[K-1]+60.0; int jmax=(int)floor(log2(Tmax)); long mneed=mj(jmax+1)+5;
  sieve(40000,(int)mneed); fprintf(stderr,"K=%ld Tmax=%.3f primes=%d\n",K,Tmax,NP);
  long G0=HINV*1,G1=(long)ceil(Tmax*HINV); const long B=2048; long nch=(G1-G0+B-1)/B;
  int nth=omp_get_max_threads(); double**fu=malloc(sizeof(double*)*nth),**pt=malloc(sizeof(double*)*nth);
  for(int i=0;i<nth;i++){fu[i]=calloc(K+3,8);pt[i]=calloc(K+3,8);}
  long done=0;
  #pragma omp parallel
  { int th=omp_get_thread_num(); double*full=fu[th],*part=pt[th];
    double*wr=malloc(8*NP),*wi=malloc(8*NP); double*Fv=malloc(8*(B+1)),*dFv=malloc(8*(B+1));
    #pragma omp for schedule(dynamic,1)
    for(long ch=0;ch<nch;ch++){
      long g0=G0+ch*B; long np=B+1;
      long jtop=ilogb((double)(g0+B)/HINV); long mtop=mj((int)jtop+1); if(mtop>NP)mtop=NP;
      for(long i=0;i<mtop;i++){ long double th0=(long double)(g0)/HINV*plgl[i]; th0=fmodl(th0,2*PIL); wr[i]=(double)cosl(th0); wi[i]=-(double)sinl(th0);}
      for(long q=0;q<np;q++){
        long g=g0+q; double t=(double)g/HINV; int j=ilogb(t); long m0=mj(j),m1=mj(j+1); if(m1>NP)m1=NP; if(m0>NP)m0=NP;
        double A=0,Bs=0,dA=0,dB=0;
        for(long i=0;i<m0;i++){ double ap,qq; onep((int)i,wr[i],wi[i],&ap,&qq); A+=ap; dA+=qq; }
        Bs=A;dB=dA;
        for(long i=m0;i<m1;i++){ double ap,qq; onep((int)i,wr[i],wi[i],&ap,&qq); Bs+=ap; dB+=qq; }
        A*=-1/PI;Bs*=-1/PI;dA*=-1/PI;dB*=-1/PI;
        double p2=ldexp(1.0,j); double w1=2-t/p2,w2=t/p2-1;
        double ph=w1*A+w2*Bs, dph=w1*dA+w2*dB+(Bs-A)/p2;
        double tt=t>15.154262241479264?t:15.154262241479264; double hc=log(log(tt)); double dhc=(t>15.154262241479264)?1.0/(t*log(t)):0;
        if(NOPHI){ph=0;dph=0;}
        if(ph>hc){dph=dhc;ph=hc;} else if(ph<-hc){dph=-dhc;ph=-hc;}
        Fv[q]=(double)N0l(t)+ph; dFv[q]=(double)dN0l(t)+dph;
        for(long i=0;i<mtop;i++){ double a=wr[i]*prc[i]+wi[i]*prs[i]; double b=wi[i]*prc[i]-wr[i]*prs[i]; wr[i]=a; wi[i]=b; }
      }
      for(long q=0;q<B;q++) cell(Fv[q],Fv[q+1],dFv[q],dFv[q+1],H,full,part);
      #pragma omp atomic
      done++;
      if(done%2000==0) fprintf(stderr,"chunks %ld/%ld\n",done,nch);
    }
  }
  /* region (0,1): F<lam_k */
  double *full=calloc(K+3,8),*part=calloc(K+3,8);
  for(int i=0;i<nth;i++)for(long k=0;k<=K+2;k++){full[k]+=fu[i][k];part[k]+=pt[i][k];}
  { long n=4000000; double m01=0; for(long i=0;i<n;i++){ double t=(i+0.5)/n; if(F_direct(t,NULL)<0.5) m01+=1.0/n; }
    part[1]+=m01; full[2]+=1.0; fprintf(stderr,"measure of F<1/2 in (0,1): %.9f\n",m01); }
  double*tau=malloc(8*K); double cum=0; for(long k=1;k<=K;k++){ cum+=full[k]; tau[k-1]=cum+part[k]; }
  char fn[512]; snprintf(fn,512,"%s_tauC-c.bin",argv[2]); FILE*o=fopen(fn,"wb"); fwrite(tau,8,K,o); fclose(o);
  /* Newton for N0(t)=k-1/2 from zLamb.bin if given */
  if(argc>5){ FILE*g=fopen(argv[5],"rb"); double*zl=malloc(8*K); fread(zl,8,K,g); fclose(g); double*nw=malloc(8*K); double maxdiff=0;
    for(long k=1;k<=K;k++){ long double t=zl[k-1],lam=k-0.5L; for(int it=0;it<20;it++){ long double d=(N0l(t)-lam)/dN0l(t); t-=d; if(fabsl(d)<1e-15L)break; }
      nw[k-1]=(double)t; double df=fabs(nw[k-1]-zl[k-1]); if(df>maxdiff)maxdiff=df; }
    snprintf(fn,512,"%s_newton-c.bin",argv[2]); o=fopen(fn,"wb"); fwrite(nw,8,K,o); fclose(o); fprintf(stderr,"max |Newton - zLamb| = %.3e\n",maxdiff); }
  return 0;
}
