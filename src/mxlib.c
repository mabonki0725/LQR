#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

char *comMalloc(int);

/* 行列式 */
typedef struct{
  double *val;
  int dim;
} Matrix;
static Matrix *newMatrix(int);
static double getDeterminant(Matrix *);
static void freeMatrix(Matrix *);

#define MXVALUE(m,r,c) ((m)->val[(r)*(m)->dim+(c)])
/********************/
/*  行列＊ベクトル演算  */
/********************/
int mxVec(mx,n,c,ar,ans)
double *mx[];  /* n x c */
double ar[];   /* c */
int n;
int c;
double *ans[];  /* n * 1 */
{
    int j,k;
    double sum;
    for(j=0;j<n;j++) {
      sum=0;
      for(k=0;k<c;k++) {
        sum += mx[j][k]*ar[k];
      }
      ans[j][0]=sum;
    }
    return(n);
}
/********************/
/*  行列＊ベクトル演算  */
/********************/
int mxVecR(mx,n,c,ar,vec)
double *mx[];  /* n x c */
double ar[];   /* c */
int n;
int c;
double *vec;  /* n * */
{
    int j,k;
    double sum;
    for(j=0;j<n;j++) {
      sum=0;
      for(k=0;k<c;k++) {
        sum += mx[j][k]*ar[k];
      }
      vec[j]=sum;
    }
    return(n);
}
/********************/
/*  行列＊行列演算  */
/********************/
int mxMult(mx,n,c,m,ar,ans)
double *mx[];  /* n x c */
double *ar[];  /* c x m */
int n;
int c;
int m;
double *ans[]; /* n x m */
{
    int i,j,k;
    double sum;
    for(i=0;i<m;i++) {
      for(j=0;j<n;j++) {
         sum=0;
         for(k=0;k<c;k++) {
           sum += mx[j][k]*ar[k][i];
         }
         ans[j][i]=sum;
      }
    }
    return(n);
}
/******************/
/* 逆行列演算     */
/* Gauss_Jordan法 */
/******************/
int mxRevGJ(mx,n,ans)
double *mx[];
int n;
double *ans[];
{
    int ipv,i,j;
    double inv_pivot,temp;
    double big;
    int pivot_row;
    int *row;

    double **wk;
 
    if(n == 0) {
      return(1);
    }
    else if(n == 1) {
      temp= mx[0][0];
      if(temp) {
        ans[0][0]=1/mx[0][0];
        return(0);
      }
      else return(1);
    }

    /* 作業領域の確保 */
    wk=(double **)comMalloc(sizeof(double *)*n);
    for(i=0;i<n;i++) wk[i]=(double *)comMalloc(sizeof(double)*n);
    row=(int *)comMalloc(sizeof(int)*n);
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) wk[i][j]=mx[i][j];
    }



    /* 単位行列作成 */
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        if(i == j) ans[i][j]=1.0;
        else       ans[i][j]=0.0;
      }
    }

    for(ipv=0;ipv<n;ipv++) {
      /* 最大値探索 */
      big=0.0;
      for(i=ipv;i<n;i++) {
        if(fabs(wk[i][ipv]) > big) {
          big=fabs(wk[i][ipv]);
          pivot_row = i;
        }
      }
      if(big == 0.0) {  /* 計算不能 */
        free(row);
        for(i=0;i<n;i++) free(wk[i]);
        free(wk);

        return(0);
      }

      row[ipv]=pivot_row;
     
      /* 行の入替え */
      if(ipv != pivot_row) {
        for(i=0;i<n;i++) {
          temp=wk[ipv][i];
          wk[ipv][i]=wk[pivot_row][i];
          wk[pivot_row][i]=temp;

          temp=ans[ipv][i];
          ans[ipv][i]=ans[pivot_row][i];
          ans[pivot_row][i]=temp;
        }
      }

      /* 対角成分=1化 ピボット行処理 */
      inv_pivot=1.0/wk[ipv][ipv];
      for(j=0;j<n;j++) {
        wk[ipv][j]  *= inv_pivot;
        ans[ipv][j] *= inv_pivot;
      }

      /* 非ピポット行処理 */
      for(i=0;i<n;i++) {
        if(i != ipv) {
          temp = wk[i][ipv];
          for(j=0;j<n;j++) {
            wk[i][j]  -= temp*wk[ipv][j];
            ans[i][j] -= temp*ans[ipv][j];
          }
        }
      }
    }
    /* 領域の開放　*/
    free(row);
    for(i=0;i<n;i++) free(wk[i]);
    free(wk);

    return(1);
}
/**********************
  行列の和算
***********************/
int mxAdd(mx1,mx2,n,m,ans)
double **mx1;
double **mx2;
double **ans;
int n,m;
{
    int i,j;

    for(i=0;i<n;i++) {
      for(j=0;j<m;j++) {
        ans[i][j]=mx1[i][j]+mx2[i][j];
      }
    }
    return(i);
}
int mxSub(mx1,mx2,n,m,ans)
double **mx1;
double **mx2;
double **ans;
int n,m;
{
    int i,j;

    for(i=0;i<n;i++) {
      for(j=0;j<m;j++) {
        ans[i][j]=mx1[i][j]+mx2[i][j];
      }
    }
    return(i);
}
/**********************
  行列のスカラー倍
***********************/
int mxScr(mx,factor,n,m,ans)
double **mx;
double factor;
double **ans;
int n,m;
{
    int i,j;

    for(i=0;i<n;i++) {
      for(j=0;j<m;j++) {
        ans[i][j]=factor*mx[i][j];
      }
    }
    return(i);
}
/**********************
  行列の転値
***********************/
int mxTrns(mx,n,m,ans)
double **mx;
double **ans;
int n,m;
{
    int i,j;

    for(i=0;i<n;i++) {
      for(j=0;j<m;j++) {
        ans[j][i]=mx[i][j];
      }
    }
    return(i);
}
/******************
  SVD分解
  注）
  データはm行n列（普通と逆)
  a,w,vはFORTRANの配列
  引用先)
  Numerical Recipe
*******************/
#define IMIN(ix,iy) ix > iy ? iy : ix
#define FMAX(x, y ) x  >  y ?  x :  y
#define SQR(x) pow(x,2.0)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double pythag(double,double);

int svdcmp(a,m,n,w,v)
double **a; /* 分析対象行列 */
int n;  /* 列数 */
int m;  /* 行数 */
double *w; /* 特異値対角成分 m個 */
double **v; /* m行m列 */
{
   int flag,i,its,j,jj,k,l,nm;
   double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

   rv1=(double *)malloc(sizeof(double)*(n+1));
   memset(rv1,'\0',sizeof(double)*(n+1));

   g=scale=anorm=0.0;
   /* HouseHolder法の２重対角の形に直す */
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=pythag(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30) {
			  fprintf(stderr,"No convergence in 30 SVDCMP iterations\n");
			  return(-its);
			}
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	return(m);
}

/* (a**2 + b**2)**0.5の計算　オーバフローを起こしにくい */
static double pythag(a,b)
double a;
double b;
{
   double absa,absb;
   absa = fabs(a);
   absb = fabs(b);
   if(absa > absb) return(absa*sqrt(1.0+SQR(absb/absa)));

   return(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
/*********************
  行列式
*********************/
double mxDet(mx,n) 
double **mx;
int n;
{
   Matrix *mtx;
   int i,j;
   double ans;

   mtx = newMatrix(n);
   for(i=0;i<n;i++) {
     for(j=0;j<n;j++) {
       //mtx->val[i*n+j] = mx[i][j];
       MXVALUE(mtx,i,j)=mx[i][j];
     }
   }
   /* 行列値 */
   ans = getDeterminant(mtx);
   /* 開放 */
   freeMatrix(mtx);

   return(ans);

}
static Matrix *newMatrix(int dim) {
  Matrix *matrix;
  matrix = (Matrix*)comMalloc(sizeof(Matrix));
  matrix->dim = dim;
  matrix->val = (double*)comMalloc(sizeof(double)*dim*dim);
  return(matrix);
}
static void freeMatrix(Matrix *matrix) {
  free(matrix->val);
  free(matrix);
}

static double getDeterminant(Matrix *mx) 
{
   double value;
   int i, r, c;
   Matrix *sub;

   double wrk,val;

   if (mx->dim == 1) {
     value = MXVALUE(mx,0,0);
   }
   else {
     value = 0;
     for (i=0 ; i<mx->dim ; i++) {
       sub = newMatrix(mx->dim-1);
       if (sub) {
         for (r=0 ; r<i ; r++) {
           for (c=1 ; c<mx->dim ; c++) {
             val = MXVALUE(mx,r,c);
             MXVALUE(sub,r,c-1) = val;
           }
         }
         for (r=i+1 ; r<mx->dim ; r++) {
           for (c=1 ; c<mx->dim ; c++) {
             val = MXVALUE(mx,r,c);
             MXVALUE(sub,r-1,c-1) = val;
           }
         }
         if (i%2 == 0) {
           wrk = getDeterminant(sub);
           if(wrk == DBL_MAX) return(DBL_MAX);
           val = MXVALUE(mx,i,0);
           value += (val * wrk);
         }
         else {
           wrk = getDeterminant(sub);
           if(wrk == DBL_MAX) return(DBL_MAX);
           val = MXVALUE(mx,i,0);
           value -= (val * wrk);
         }
         /* 領域の開放 */
         free(sub->val);
         free(sub);
       }
       else {
         fprintf(stderr,"cannot get memory\n");    
         return(DBL_MAX);     
       }
     }

   }
   return value;
}