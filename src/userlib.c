#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "usrlib.h"

/**********************
  共通メモリー割り当て
***********************/
char *comMalloc(size)
int size;
{
   char *p;
   p=(char *)malloc(size);
   if(p) memset(p,'\0',size); 
   else  printf("malloc0 領域は確保出来ません!!\n");
   return(p);
}
/********************
  共通文字割り当て
*********************/
char *comAssign(adr,str)
char **adr;
char *str;
{
    int len;

    len=strlen(str);
    if(!len) return(NULL);

    *adr=(char *)malloc(len+1);
    memset(*adr,'\0',len+1);
    strcpy(*adr,str);
    return(*adr);
}
/*********************
  文字列のリスト
**********************/
int comCmember(memb,num,pc)
char **memb;
int num;
char *pc;
{
   int i;
   for(i=0;i<num;i++) {
     if(!strcmp(memb[i],pc)) break;
   }
   if(i == num) return(-1);
   return(i);
}
/*********************
  文字列のリスト
**********************/
int comNmember(memb,num,id)
int *memb;
int num;
int id;
{
   int i;
   for(i=0;i<num;i++) {
     if(memb[i] == id) break;
   }
   if(i == num) return(-1);
   return(i);
}
/**************
数値のソート
***************/
int comDsort(v,order,n)
double v[];
int  order;
int  n;
{
   int gap,i,j;
   double dwrk;

   for(gap = n/2;gap > 0;gap /= 2) {
     for(i = gap;i < n; i++) { 
       for(j = i-gap; j >= 0; j -= gap) {
    	 if((order == FLOW && v[j] <= v[j+gap]) ||
	        (order != FLOW && v[j] >= v[j+gap])) break;
	       dwrk=v[j];
	       v[j]=v[j+gap];
	       v[j+gap]=dwrk;
       } 
     }
   }
   return(n);
}
/**************
数値のソート順番
***************/
int comDsortJun(org,order,n,jun)
double org[];
int  order;
int  n;
int  jun[];
{
   int gap,i,j;
   double dwrk;
   double *v;
   int jwrk;

   v=(double *)comMalloc(sizeof(double)*n);
   for(i=0;i<n;i++) v[i]=org[i];

   for(gap = n/2;gap > 0;gap /= 2) {
     for(i = gap;i < n; i++) { 
       for(j = i-gap; j >= 0; j -= gap) {
    	 if((order == FLOW && v[j] <= v[j+gap]) ||
	        (order != FLOW && v[j] >= v[j+gap])) break;
	       dwrk=v[j];
	       v[j]=v[j+gap];
	       v[j+gap]=dwrk;

           jwrk=jun[j];
           jun[j]=jun[j+gap];
           jun[j+gap]=jwrk;
       } 
     }
   }
   free(v);
   return(n);
}