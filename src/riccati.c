/******************
  Ricatti Eq
  2016/09/06
*******************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "libLBFGS\lbfgs.h"

#include "mxlib.h"
#include "usrlib.h"

#define EPS 0.00001

#define LBFG

/* protoType Declare */
double riccati(int,int,double **,double **,double **,double **,double **,double **);


/*************
data format 
A    0
Q    0
B'   R
--------------
N*N 0*M
N*N 0*M
M*N M*M
**************/

#ifdef LBFG
struct riccatiPrm {
  int iter;

  int n;
  int m;
  double **A;
  double **B;
  double **Q;
  double **R;

  double **P;
  double **Z;

  double **WP;
  double **WZ;
};
static double xGrant();
static double xFunc();

static double **g_P;
static double **g_Z;
/****************
  libLBFGS evaluate 
*****************/
static lbfgsfloatval_t evaluate (
  void *instance,
  const lbfgsfloatval_t *x,
  lbfgsfloatval_t       *g,
  const int n,
  const lbfgsfloatval_t step
)
{
   int i;
   struct riccatiPrm *prm;
   double delta;
   double fx;

   prm=(struct riccatiPrm *)instance;

   for(i=0;i<n;i++) {
     //delta=fabs(x[i])/100.0;
     delta=0.0001;
     g[i]=xGrant(i,prm,n,x,delta);
   }
   fx=xFunc(prm,n,x);
   return(fx);
}
/*****************
  libLBFGS process
******************/
static int progress(
  void *instacne,
  const lbfgsfloatval_t *x,
  const lbfgsfloatval_t *g,
  const lbfgsfloatval_t fx,
  const lbfgsfloatval_t xnorm,
  const lbfgsfloatval_t gnorm,
  const lbfgsfloatval_t step,
  int n,
  int k,
  int ls
) 
{

   int i;

   printf("Iteration %d:\n", k);
   for(i=0;i<n;i++) {
     printf("x[%2d]=%lf\n",i,x[i]);
   }
   printf("  fx = %lf\n",fx);

   printf("xnorm=%lf gnorm=%lf step=%lf\n",xnorm,gnorm,step);

   return(0);
}
/************************
  目的関数
*************************/
static double xFunc(prm,n,x)
struct riccatiPrm *prm;
lbfgsfloatval_t *x;
int n;
{
     double fx;
     int i,j;
     //double **P,**Z;


     for(i=0;i<prm->n;i++) {
       for(j=0;j<prm->n;j++) {
         g_P[i][j] = x[i*prm->n+j];
       }
     }
     fx=riccati(prm->n,prm->m,prm->A,prm->B,prm->Q,prm->R,g_P,g_Z);


	 return(fx);   

}
/*************************
  偏微分（擬似）
**************************/
static double xGrant(id,prm,n,x,step)
int id;  /* 母数のid */
struct riccatiPrm *prm;
int n;
lbfgsfloatval_t *x;
double step;
{

   double fx1,fx2,grant;
   lbfgsfloatval_t *xx;
   int i;
   
   /* 領域 */
   xx=(lbfgsfloatval_t *)malloc(sizeof(lbfgsfloatval_t)*n);
   for(i=0;i<n;i++) {
	 xx[i]=x[i];
	 if(i == id) xx[i] += step;
   }
   
   fx1=xFunc(prm,n,x);
   fx2=xFunc(prm,n,xx);
   
   if(fabs(fx1 - fx2) > EPS) grant=(fx2 - fx1)/step;
   else                      grant=0;
   
   /* 領域の開放 */
   free(xx);
   
   return(grant);
}   
#endif
/*********************
    Main
**********************/
#define MAX_REC 2045

static double ***g_wk;


int main(argc,argv)
int argc;
char *argv[];
{
    char record[MAX_REC],*pc;
    int i,j,k;
    int n,m;
    double **data;
    double **A,**B,**Q,**R,**P;
    //double ***wk;
    //int N;
    double **PP,**df,**dx,**Idf;
    double r0;
    double **Z,**ZZ;

    FILE *fp,*fw;

    lbfgs_parameter_t param;
    //lbfgsfloatval_t fx;

    lbfgsfloatval_t *x;   /* n次元 独立変数 */
    double *xold;

    struct riccatiPrm *parm;
    double sum;
    int maxLoop,loop,ret;

    int N;



    if(argc < 5) {
      fprintf(stderr,"USAGE Xsiae Usize readFile outFile!\n");
      exit(-9);
    }
    n=atoi(argv[1]);
    m=atoi(argv[2]);
    if(n * m <= 0) {
      fprintf(stderr,"Wrong Matrix size Xsize=%d Usize=%d\n",n,m);
      exit(-9);
    }
    if(!(fp=fopen(argv[3],"r"))) {
      fprintf(stderr,"Cannot Read file=[%s]\n",argv[3]);
      exit(-1);
    }
    if(!(fw=fopen(argv[4],"w"))) {
      fprintf(stderr,"Cannot Write file=[%s]\n",argv[4]);
      exit(-1);
    }

    i=0;
    k=0;
    while(fgets(record,MAX_REC,fp)) {
      if(record[0] == '#') continue;
      if(record[0] == '$') continue;
      record[strlen(record)-1]='\0';

      pc=strtok(record,", ");
      j=0;
      while(pc) {
        j++;
        pc=strtok(NULL,", ");
      }
      if(k < j) k=j;
      i++;
    }
    /* record size check */
    if(k != n+m) {
      fprintf(stderr,"Unmatch col size read=%d correct=%d\n",k,n+m);
      exit(-9);
    }
    if(i != 2*n+m) {
      fprintf(stderr,"Unmatch row size read=%d correct=%d\n",k,2*n+m);
      exit(-9);
    }

    if(n < m) N=m;
    else      N=n;
    g_wk=(double ***)comMalloc(sizeof(double **)*15);
    for(i=0;i<15;i++) {
      g_wk[i] = (double **)comMalloc(sizeof(double *)*N);
      for(j=0;j<N;j++) {
        g_wk[i][j]=(double *)comMalloc(sizeof(double)*N);
      }
    }
    g_P=(double **)comMalloc(sizeof(double *)*n);
    g_Z=(double **)comMalloc(sizeof(double *)*n);
    for(i=0;i<n;i++) {
      g_P[i]=(double *)comMalloc(sizeof(double)*n);
      g_Z[i]=(double *)comMalloc(sizeof(double)*n);
    }
    /* 領域の獲得 */
    A = (double **)comMalloc(sizeof(double *)*n);
    Q = (double **)comMalloc(sizeof(double *)*n);
    B = (double **)comMalloc(sizeof(double *)*n);
    P = (double **)comMalloc(sizeof(double *)*n);
    PP= (double **)comMalloc(sizeof(double *)*n);
    df= (double **)comMalloc(sizeof(double *)*n);
    dx= (double **)comMalloc(sizeof(double *)*n);
    Idf=(double **)comMalloc(sizeof(double *)*n);
    Z =(double **)comMalloc(sizeof(double *)*n);
    ZZ=(double **)comMalloc(sizeof(double *)*n);

    for(i=0;i<n;i++) {
      A[i]=(double *)comMalloc(sizeof(double)*n);
      Q[i]=(double *)comMalloc(sizeof(double)*n);
      P[i]=(double *)comMalloc(sizeof(double)*n);
      PP[i]=(double *)comMalloc(sizeof(double)*n);
      df[i]=(double *)comMalloc(sizeof(double)*n);
      dx[i]=(double *)comMalloc(sizeof(double)*n);
      Idf[i]=(double *)comMalloc(sizeof(double)*n);
      Z[i]=(double *)comMalloc(sizeof(double)*n);
      ZZ[i]=(double *)comMalloc(sizeof(double)*n);

      B[i]=(double *)comMalloc(sizeof(double)*m);
    }

    R = (double **)comMalloc(sizeof(double *)*m);
    for(i=0;i<m;i++) {
      R[i]=(double *)comMalloc(sizeof(double)*m);
    }


    /* read file */
    rewind(fp);
    data=(double **)comMalloc(sizeof(double *)*(2*n+m));
    for(i=0;i<2*n+m;i++) {
      data[i]=(double *)comMalloc(sizeof(double)*(n+m));
    }
    i=0;
    k=0;
    while(fgets(record,MAX_REC,fp)) {
      if(record[0] == '#') continue;
      if(record[0] == '$') continue;

      pc=strtok(record,", ");
      j=0;
      while(pc) {
        data[i][j]=atof(pc);
        j++;
        pc=strtok(NULL,", ");
      }
      if(k < j) k=j;
      i++;
    }
    fclose(fp);

    /* 行列の設定 */
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        A[i][j]=data[i][j];
      }
    }
    for(i=n;i<n+n;i++) {
      for(j=0;j<n;j++) {
        Q[i-n][j]=data[i][j];
      }
    }
    for(i=n+n;i<n+n+m;i++) {
      for(j=0;j<n;j++) {
        B[j][i-(n+n)]=data[i][j];
      }
    }
    for(i=n+n;i<n+n+m;i++) {
      for(j=n;j<n+m;j++) {
        R[i-(n+n)][j-n]=data[i][j];
      }
    }

    for(i=0;i<2*n+m;i++) {
      free(data[i]);
    }
    free(data);
    
    /* riccati Equation */
#ifdef LBFG
    /****************/
    /* 初期値の設定 */
    /****************/  
    /* 領域の確保 */
    parm=(struct riccatiPrm *)malloc(sizeof(struct riccatiPrm));
    memset(parm,'\0',sizeof(struct riccatiPrm));

    x=lbfgs_malloc(n*n);
    xold=(double *)malloc(sizeof(double)*n*n);

    /******************
      初期値の設定
    *******************/
    for(i=0;i<n*n;i++) x[i]=0;
    //for(i=0;i<n;i++) x[i]=0.5;

    /*********************
      計算可能かのCheck 
    /*********************/
    parm->A=A;
    parm->B=B;
    parm->Q=Q;
    parm->R=R;
    parm->n=n;
    parm->m=m;

    sum=xFunc(parm,n,x);
    if(fabs(sum - -DBL_MAX) < EPS) {
      fprintf(stderr,"cannot execute riccati equation\n");
      goto freeBlock;
    }


    lbfgs_parameter_init(&param);

    /*****************
      SUMT Loop
    ******************/
    maxLoop=100;


    for(loop=1;loop<maxLoop;loop++) {
 	  parm->iter=loop;

      /******************/
      /* 準ニュートン法 */
      /******************/
      ret=lbfgs(n*n,x,&r0,evaluate,progress,parm,&param);
#if 0
	  sum=1/EPS;
      for(i=0;i<n*n;i++) {
        //fprintf(fw,"%lf,",x[i]);
	    if(loop > 1) sum=pow(xold[i]-x[i],2.0);
      }
	  //fprintf(fw,"\n");
#else
     sum=xFunc(parm,n*n,x);
#endif
	  if(loop > 20 && sum < EPS) break;

	  for(i=0;i<n*n;i++) xold[i]=x[i];
    }
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        P[i][j]=x[i*n+j];
      }
    }

    /* 領域の開放 */
    if(x) lbfgs_free(x);
    if(xold) free(xold);

#else
    k=0;
    while(k < 10000) {
      r0=riccati(n,m,A,B,Q,R,P,Z);
      if(fabs(r0) < EPS) break;

      fprintf(stderr,"no=%5d r0=%lf\n",k,r0);
      for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
          PP[i][j] = P[i][j] + 0.00001;
          riccati(n,m,A,B,Q,R,PP,ZZ);
          df[i][j]=(ZZ[i][j] - Z[i][j])/0.00001;
        }
      }
      mxRevGJ(df,n,Idf);
      mxSub(ZZ,Z,n,n,ZZ);
      mxMult(ZZ,n,n,n,Idf,dx);
      mxSub(P,dx,n,n,P);
      k++;
    }
#endif
    fprintf(fw,"iter=%d sum=%lf\n",loop,sum);
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        fprintf(fw,"%lf",P[i][j]);
        if(j < n-1) fprintf(fw,",");
        else        fprintf(fw,"\n");
      }
    }
    fclose(fw);

freeBlock:

    /* 領域の開放 */
    for(i=0;i<n;i++) {
      free(A[i]);
      free(B[i]);
      free(Q[i]);
      free(P[i]);
      free(PP[i]);
      free(df[i]);
      free(dx[i]);
      free(Idf[i]);
      free(Z[i]);
      free(ZZ[i]);
    }
    free(A);
    free(B);
    free(Q);
    free(P);
    free(PP);
    free(df);
    free(dx);
    free(Idf);
    free(Z);
    free(ZZ);

    for(i=0;i<m;i++) {
      free(R[i]);
    }
    free(R);

    for(i=0;i<15;i++) {
      for(j=0;j<N;j++) {
        free(g_wk[i][j]);
      }
      free(g_wk[i]);
    }
    free(g_wk);
    for(i=0;i<n;i++) {
      free(g_P[i]);
      free(g_Z[i]);
    }
    free(g_P);
    free(g_Z);
}
/******************
  riccati Function

  A'*P+P*A+Q-P*B*(IR)*B'*P
*******************/
double riccati(n,m,A,B,Q,R,P,Z)
int n,m;
double **A; /* n*n */
double **B; /* n*m */
double **Q; /* n*n */
double **R; /* m*m */
double **P; /* n*n */
double **Z; /* n*n out */
{

    int i,j;
    //int j,N;
    //double ***wk;
    double **IR;
    double ans;

    IR = (double **)comMalloc(sizeof(double *)*m);
    for(i=0;i<m;i++) {
      IR[i]=(double *)comMalloc(sizeof(double)*m);
    }


    mxTrns(A,n,n,g_wk[1]);             /* wk1 n*n */
    mxMult(g_wk[1],n,n,n,P,g_wk[2]);   /* wk2 n*n */
    mxMult(P,n,n,n,A,g_wk[3]);         /* wk3 n*n */
    mxAdd(g_wk[2],g_wk[3],n,n,g_wk[4]);/* wk4 n*n */
    mxAdd(g_wk[4],Q,n,n,g_wk[5]);      /* wk5 n*n */

    mxRevGJ(R,m,IR);                   /* IR  m*m */
    mxTrns(B,n,m,g_wk[6]);             /* wk6 m*n */ 
    mxMult(g_wk[6],m,n,n,P,g_wk[7]);   /* wk7 m*n */
    mxMult(IR,m,m,n,g_wk[7],g_wk[8]);  /* wk8 m*n */ 
    mxMult(B,n,m,n,g_wk[8],g_wk[9]);   /* wk9 n*n */
    mxMult(P,n,n,n,g_wk[9],g_wk[10]);  /* wk10 n*n */ 
    mxSub(g_wk[5],g_wk[10],n,n,Z);     /* Z   n*n */

#if 0
    ans=mxDet(Z,n);
#else
    ans=0;
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        ans += fabs(Z[i][j]);
      }
    }
#endif
    for(i=0;i<m;i++) {
      free(IR[i]);
    }
    free(IR);

    return(ans);

}

#if 0
/******************
  Riccati Equation
*******************/
double A[2][2] = {{-2, 0}, {2, 1}};   //行列A
double B[2] = {-1, 1};                //行列B
double Q[2][2] = {{200, 0}, {0, 200}};  //重み行列Q
double R = 10;                          //重み行列R
double f[2] = {-4, -6};               //設計ゲインF
double Af[2][2];                      //行列AF(i)
double P[2][2];                       //リアプノフ方程式の解P(i)
double oldP[2][2] = {0};              //収束条件計算用P(i - 1)

/* 行列表示関数 */
void disp_array(char *str, double *a, int n, int m)
{
    int i;  //行
    int j;  //列
    
    puts(str);
    for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            printf("%f ",*(a + i * m + j));
        }
        printf("\n");
    }
}

/* 行列AF(i)を求める関数 */
void setAf(void)
{
    int i;  //行
    int j;  //列
    
    for(i = 0; i < 2; i++){
        for(j = 0; j < 2; j++){
            Af[i][j] = A[i][j] + B[i] * f[j];
        }
    }
}

/* 逆行列を求める関数 */
void inv(double a[3][3])
{    
    double inv_a[3][3];  //ここに逆行列が入る
    double buf;          //一時的なデータを蓄える
    int i;
    int j;
    int k;
    
    /* 単位行列を作る */
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            if(i == j){
                inv_a[i][j] = 1.0;
            }else{
                inv_a[i][j] = 0.0;
            }
        }
    }
    
    disp_array("inv_a:", (double *)inv_a, 3, 3);
    
    /* 掃き出し法 */
    for(i = 0; i < 3; i++){
        buf = 1 / a[i][i];
        for(j = 0; j < 3; j++){
            a[i][j] *= buf;
            inv_a[i][j] *= buf;
        }
        for(j = 0; j < 3; j++){
            if(i != j){
                buf = a[j][i];
                for(k = 0; k < 3; k++){
                    a[j][k] -= a[i][k] * buf;
                    inv_a[j][k] -= inv_a[i][k] * buf;
                }
            }
        }
    }
        
    /* 求めた逆行列を代入 */
    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
            a[i][j] = inv_a[i][j];
}

/* リアプノフ方程式の解P(i)を求める関数 */
void makeP(void)
{
    double a[3][3] = {0};
    double c[3];
    double p[3] = {0};
    int i;
    int j;
    
    a[0][0] = 2 * Af[0][0];
    a[0][1] = 2 * Af[1][0];
    a[1][0] = Af[0][1];
    a[1][1] = Af[0][0] + Af[1][1];
    a[1][2] = Af[1][0];
    a[2][1] = 2 * Af[0][1];
    a[2][2] = 2 * Af[1][1];
    
    disp_array("a:", (double *)a, 3, 3);     //debug
    
    /* 逆行列を計算 */
    inv(a);
    
    disp_array("inv_a:", (double *)a, 3, 3); //debug
    
    c[0] = R * f[0] * f[0] - Q[0][0];
    c[1] = R * f[0] * f[1];
    c[2] = R * f[1] * f[1] - Q[1][1];
    
    /* 行列P(i)を計算する */
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            p[i] += (a[i][j] * c[j]);
        }
    }
    
    P[0][0] = p[0];
    P[0][1] = p[1];
    P[1][0] = p[1];
    P[1][1] = p[2];
}

/* 収束条件計算関数 */
double e_range(void)
{
    int i;           //行
    int j;           //列
    double e[2][2];  //P(i)とP(i - 1)の差を格納する行列
    double ret;      //行列式
    
    // 差を求める
    for(i = 0; i < 2; i++){
        for(j = 0; j < 2; j++){
            e[i][j] = P[i][j] - oldP[i][j];
        }
    }
    
    //行列式
    ret = e[0][0] * e[1][1] - e[0][1] * e[1][0];
    
    //oldP(P(i - 1))へコピー
    for(i = 0; i < 2; i++){
        for(j = 0; j < 2; j++){
            oldP[i][j] = P[i][j];
        }
    }
    
    if(ret < 0)
        return -1 * ret;
    return ret;
}

/* ゲインF(i + 1)を求める関数 */
void makenewf(void)
{
    f[0] = -1 / R * (B[0] * P[0][0] + B[1] * P[0][1]);
    f[1] = -1 / R * (B[0] * P[0][1] + B[1] * P[1][1]);
}

/* 行列AF(i + 1)を求める関数 */
void makenewAf(void)
{
    Af[0][0] = A[0][0] + B[0] * f[0];
    Af[0][1] = A[0][1] + B[0] * f[1];
    Af[1][0] = A[1][0] + B[1] * f[0];
    Af[0][0] = A[1][1] + B[1] * f[1];
}

/* ゲインFを逐次表示する関数 */
void dispf(int i,double e)
{
    printf("%d回目f1 = %f f2 = %f e = %f\n", i, f[0], f[1], e);
}

/* メイン関数 */
int main(void)
{
    int i = 0;
    //int ch;
    double e;
    
    setAf();  //AF0を求める。
    disp_array("Af:", (double *)Af, 2, 2);
    
    /* 収束条件を満たすまで繰り返す */
    while(1){
        makeP();
        e = e_range();
        dispf(i, e);
        if(e < 1.0e-6)
            break;
        i++;
        makenewf();
        makenewAf();
    }
    return 0;
}
#endif
