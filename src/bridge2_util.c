#include "bridge2_util.h"



double code_miss=-9999999;
double c=9;


/****************************************************************************************************/
/*                                                                   stdd                                                                                                  */
/*  Purpose:   Return  the standard deviation of a vector                                                                                       */
/* Argument description:                                                                                                                                       */
/* vector: The sample vector                                                                                                                                  */
/* length: the length of the vector                                                                                                                          */
/* is_finite : (output) The number of finite values in vector                                                                                   */
/****************************************************************************************************/ 


double stdd(double *vector,int *length, int *is_finite)
{
  /* Compute the standard deviation of a vector */

  int i,count=0;
  double sum=0;
  double x_bar;
  double result;
  
  
  x_bar=mean_vec(vector,length);
  if(x_bar==code_miss)
    return(code_miss);
  else
    {
      for(i=0;i<*length;i++)
	{
	  if(vector[i]!=code_miss)
	    {
	      count=count+1;
	      sum=sum+(vector[i]-x_bar)*(vector[i]-x_bar);
	    }
	}
      *is_finite=count;
      if(count>1)
	{
	  result=(sqrt(sum/((double)count-1)));
	  return(result);
	}
      else
	{
	  return(code_miss);
	}
    }
}

/****************************************************************************************************/
/*                                                                   vec_mat                                                                                            */
/*  Purpose:   Coerce a vector into a matrix                                                                                                           */
/* Argument description:                                                                                                                                       */
/* array_vec: The vector to coerce                                                                                                                          */
/* nb_row: The number of row for the matrix                                                                                                       */
/* nb_col: The number of column for the matrix                                                                                                    */
/* array :(output) The two dimmensional array                                                                                                     */
/****************************************************************************************************/ 
void vec_mat(double *array_vec,int* nb_row,int *nb_col,double **array)
{
  int i,j;

  for(i=0;i<*nb_row;i++)
    for(j=0;j<*nb_col;j++)
      array[i][j]=array_vec[i**nb_col+j];
}

/****************************************************************************************************/
/*                                                                   mat_vec                                                                                            */
/*  Purpose:   Coerce a matrix into a vector                                                                                                           */
/* Argument description:                                                                                                                                       */
/* array_vec: (outpout)The vector                                                                                                                         */
/* nb_row: The number of row for the matrix                                                                                                       */
/* nb_col: The number of column for the matrix                                                                                                    */
/* array : The two dimmensional array to coerce                                                                                                    */
/****************************************************************************************************/ 

void mat_vec(double *array_vec,int* nb_row,int *nb_col,double **array)
{
  int i,j;

  for(i=0;i<*nb_row;i++)
    for(j=0;j<*nb_col;j++)
      array_vec[i**nb_col+j]=array[i][j];
}

/****************************************************************************************************/
/*                                                                   init_ivector                                                                                       */
/*  Purpose:   Initialize a vector of type int to zero                                                                                                */
/* Argument description:                                                                                                                                       */
/* vector: The vector to initialize                                                                                                                            */
/* length: the length of the vector                                                                                                                          */
/****************************************************************************************************/ 

void init_ivector(int *vector,int *length, int  value)
{
  int i;

  for(i=0;i<*length;i++)
    vector[i]=value;
}

/****************************************************************************************************/
/*                                                                   dmatrix                                                                                            */
/*  Purpose:  Allocate the memory for a matrix of type double                                                                             */
/* Argument description:                                                                                                                                       */
/* nb_row: The number of row for the matrix                                                                                                       */
/* nb_col: The number of column for the matrix                                                                                                    */
/****************************************************************************************************/ 
double **dmatrix(int nb_row,int nb_col)
{
  double **array;
  int i,j;

  /* Allocate the memory */
  array=Calloc(nb_row,double*);

  for(i=0;i<nb_row;i++)
    array[i]=Calloc(nb_col, double);

  /* Initialize to zero*/
  for(i=0;i<nb_row;i++)
    for(j=0;j<nb_col;j++)
      array[i][j]=0;

  return(array);
}
/****************************************************************************************************/
/*                                                                   dvector                                                                                             */
/*  Purpose:  Allocate the memory for a vector of type double                                                                              */
/* Argument description:                                                                                                                                       */
/* length: The length of the vector                                                                                                                         */
/* init: The value to initialize the vector                                                                                                                 */
/****************************************************************************************************/ 

double *dvector(int length, int init)
{
  int i;
  double *vector;

  /* Allocate the memory */
    vector=Calloc(length, double);
  
  /* Initialize the memory */
  for(i=0;i<length;i++)
    vector[i]=init;

  return(vector);
}

/****************************************************************************************************/
/*                                                                    ivector                                                                                             */
/*  Purpose:  Allocate the memory for a vector of type integer                                                                              */
/* Argument description:                                                                                                                                       */
/* length: The length of the vector                                                                                                                         */
/* init: The value to initialize the vector                                                                                                                 */
/****************************************************************************************************/ 


int *ivector(int length, int init)
{
  int i;
  int *vector;

  /* Allocate the memory */
  vector=Calloc(length, int);

  /* Initialize the memory */
  for(i=0;i<length;i++)
    vector[i]=init;

  return(vector);
}

/****************************************************************************************************/
/*                                                                mean_vec                                                                                             */
/*  Purpose:  Return the mean of vector (remove the missing values)                                                                   */
/* Argument description:                                                                                                                                       */
/* vector : The sample vector                                                                                                                                 */
/* length: The length of the vector                                                                                                                         */
/****************************************************************************************************/ 


double  mean_vec(double *vector,int *length)
{
  int i,count=0;
  double sum=0;

  for(i=0;i<*length;i++)
    {
      if(vector[i]!=code_miss)
	{
	  count=count+1;
	  sum=sum+vector[i];
	}
    }
  if (count>0)
    {
      return(sum/(double)count);
    }
  else
    {
      return(code_miss);
    }
  
}

/****************************************************************************************************/
/*                                                                free_dmatrix                                                                                        */
/*  Purpose:  Free the memory of a matrix of type double                                                                                     */
/* Argument description:                                                                                                                                       */
/* array: the two dimmensional array to free                                                                                                         */
/* nb_row : its number of row                                                                                                                               */
/****************************************************************************************************/ 


void free_dmatrix(double **array, int nb_row)
{

  int i;
  for(i=0;i<nb_row;i++)
    Free(array[i]);

  Free(array);
}
/****************************************************************************************************/
/*                                                                   quicksort2                                                                                        */
/*  Purpose:   Sort a vector using the quicksort algorithm  and move another vector at the same time                */
/* Argument description:                                                                                                                                       */
/* a: The vector to sort  from p to r                                                                                                                       */
/* b: The second vector to sort  from p to r                                                                                                            */
/* p: The first index                                                                                                                                                */
/* r: The last index                                                                                                                                                  */
/*****************************************************************************************************/ 


void quicksort2(double *a, int *b, int *p, int *r)
{
  int q;
  int q_p;

  if (*p<*r)
    {
    q=rand_part2(a,b, *p, *r);
    quicksort2(a,b,p,&q);
    q_p=q+1;
    quicksort2(a,b,&q_p,r);
    }
}
int partition2(double *a, int *b, int p, int r)
{
  double x=a[p];
  int i=p-1;
  int j=r+1;
  double temp;
  int itemp;

  for(;;)
    {
      do
	{
	  j--;
	}while(a[j]>x);
      do
	{
	  i++;
	}while(a[i]<x);
      if(i<j)
	{
	  temp=a[i];
	  a[i]=a[j];
	  a[j]=temp;
	  itemp=b[i];
	  b[i]=b[j];
	  b[j]=itemp;
	}
      else
	{
	  return(j);
	}
    }

}

int rand_part2(double *a, int *b,  int p, int r)
{
  int i;
  double temp;
  int itemp;

  i=uni_rand(p,r);
  temp=a[p];
  a[p]=a[i];
  a[i]=temp;
  itemp=b[p];
  b[p]=b[i];
  b[i]=itemp;

  return(partition2(a,b,p,r));

}
/****************************************************************************************************/
/*                                                                   idquicksort2                                                                                       */
/*  Purpose:   Sort a vector using the quicksort algorithm  and move another vector at the same time                */
/* Argument description:                                                                                                                                       */
/* a: The vector to sort  from p to r                                                                                                                       */
/* b: The second vector to sort  from p to r                                                                                                            */
/* p: The first index                                                                                                                                                */
/* r: The last index                                                                                                                                                  */
/*****************************************************************************************************/ 


void idquicksort2(int *a, double *b, int *p, int *r)
{
  int q;
  int q_p;

  if (*p<*r)
    {
    q=idrand_part2(a,b, *p, *r);
    idquicksort2(a,b,p,&q);
    q_p=q+1;
    idquicksort2(a,b,&q_p,r);
    }
}
int idpartition2(int *a, double *b, int p, int r)
{
  int x=a[p];
  int i=p-1;
  int j=r+1;
  double temp;
  int itemp;

  for(;;)
    {
      do
	{
	  j--;
	}while(a[j]>x);
      do
	{
	  i++;
	}while(a[i]<x);
      if(i<j)
	{
	  itemp=a[i];
	  a[i]=a[j];
	  a[j]=itemp;
	  temp=b[i];
	  b[i]=b[j];
	  b[j]=temp;
	}
      else
	{
	  return(j);
	}
    }

}

int idrand_part2(int *a, double *b,  int p, int r)
{
  int i;
  double temp;
  int itemp;

  i=uni_rand(p,r);
  itemp=a[p];
  a[p]=a[i];
  a[i]=itemp;
  temp=b[p];
  b[p]=b[i];
  b[i]=temp;

  return(idpartition2(a,b,p,r));

}

/****************************************************************************************************/
/*                                                                  uni_rand                                                                                            */
/*  Purpose:  Return a random number between min and max                                                                               */
/****************************************************************************************************/ 

int uni_rand(int min,int max)
{
  int rand_nb;


  GetRNGstate();
  rand_nb=(int)(unif_rand()*(max+1-min)+min);
  PutRNGstate();
  
  return(rand_nb);

}


/****************************************************************************************************/
/*                                                              init_dvector                                                                                           */
/*  Purpose:  Initialize a vector of type double                                                                                                       */
/* Argument description:                                                                                                                                       */
/* vector : The vector                                                                                                                                             */
/* length: The length of the vector                                                                                                                         */
/* value: The value to initialize it with                                                                                                                    */
/****************************************************************************************************/ 

void init_dvector(double *vector, int *length, int value)
{
  int i;

 
  /* Initialize to value */
  for(i=0;i<*length;i++)
    vector[i]=value;
}

/****************************************************************************************************/
/*                                                                sum _vec                                                                                             */
/*  Purpose:  Return the sum  of vector (remove the missing values)                                                                   */
/* Argument description:                                                                                                                                       */
/* vector : The sample vector                                                                                                                                 */
/* length: The length of the vector                                                                                                                         */
/****************************************************************************************************/ 

double  sum_vec(double *vector,int *length)
{
  int i,count=0;
  double sum=0;

  for(i=0;i<*length;i++)
    {
      if(vector[i]!=code_miss)
	{
	  count=count+1;
	  sum=sum+vector[i];
	}
    }
  if (count>0)
    {
      return(sum);
    }
  else
    {
      return(code_miss);
    }
}

/****************************************************************************************************/
/*                                                                   quicksort                                                                                         */
/*  Purpose:   Sort a vector using the quicksort algorithm                                                                                     */
/* Argument description:                                                                                                                                       */
/* a: The vector to sort  from p to r                                                                                                                       */
/* p: The first index                                                                                                                                                */
/* r: The last index                                                                                                                                                  */
/*****************************************************************************************************/ 

void quicksort(double *a, int *p, int *r)
{
  int q;
  int q_p;

  if (*p<*r)
    {
    q=rand_part(a, *p, *r);
    quicksort(a,p,&q);
    q_p=q+1;
    quicksort(a,&q_p,r);
    }
}


int partition(double *a, int p, int r)
{
  double x=a[p];
  int i=p-1;
  int j=r+1;
  double temp;

  for(;;)
    {
      do
	{
	  j--;
	}while(a[j]>x);
      do
	{
	  i++;
	}while(a[i]<x);
      if(i<j)
	{
	  temp=a[i];
	  a[i]=a[j];
	  a[j]=temp;
	}
      else
	return(j);
    }

}

int rand_part(double *a, int p, int r)
{
  int i;
  double temp;

  i=uni_rand(p,r);
  temp=a[p];
  a[p]=a[i];
  a[i]=temp;
  return(partition(a,p,r));

}
/****************************************************************************************************/
/*                                                                     dabs                                                                                                */
/*  Purpose:  Return the abs value of a number of type double                                                                              */
/****************************************************************************************************/ 
double dabs(double a)
{
  if (a<0)
    return(-a);
  else
    return(a);
}
/****************************************************************************************************/
/*                                                                     dmax                                                                                               */
/*  Purpose:  Return the maximum of 2 doubles                                                                                                    */
/****************************************************************************************************/ 
double dmax(double a, double b)
{
  if(a<b)
    return(b);
  else
    return(a); 
}


void qr_solve(double **x, int *n1, double ** y, double **coef)
/* Translation of the R function qr.solve into pure C
   NB We have to transpose the matrices since the ordering of an array is different in Fortran
   NB2 We have to copy x to avoid it being overwritten.
*/
{
    int i,j, info = 0, rank, *pivot, n, p;
    char *vmax;
    double tol = 1.0E-7, *qraux, *work;
    double * xt, *yt, *coeft;

    qraux = dvector(*n1,0);
    pivot = ivector(*n1,0);
    work  = dvector(2*(*n1),0);
    
    for(i = 0; i < *n1; i++)
        pivot[i] = i+1;

    /** Copy the matrix by column **/
    xt = dvector((*n1)*(*n1),0);
    for(i=0;i<*n1;i++)
      for(j=0;j<*n1;j++)
	xt[i*(*n1)+j]=x[j][i];


    n = *n1;
    p = *n1;

    F77_CALL(dqrdc2)(xt, &n, &n, &p, &tol, &rank,qraux, pivot, work);

    if (rank != p)
        error("Singular matrix in qr_solve\n");

 
    coeft=dvector((*n1)*(*n1),0);
    
    /** Copy the matrix by column **/
    yt = dvector((*n1)*(*n1),0);
    for(i=0;i<*n1;i++)
      for(j=0;j<*n1;j++)
	yt[i*(*n1)+j]=y[j][i];


    F77_CALL(dqrcf)(xt, &n, &rank, qraux,yt, &n, coeft, &info);

    /** Put back into a matrix **/
    for(i=0;i<*n1;i++)
      for(j=0;j<*n1;j++)
	coef[j][i]=coeft[i*(*n1)+j];

    
    
    Free(qraux);
    Free(pivot);
    Free(work);	
    Free(xt);
    Free(yt);
    Free(coeft);
}



void inverse(double **mat1, int *n ,double **res)
{
  /** QR decomposition solve mat1*x=I **/
  /** res contains the result **/
  int i,j;
  double **iden;
  iden=dmatrix(*n,*n);
  
  for(i=0;i<*n;i++)
    iden[i][i]=1;

  qr_solve(mat1, n, iden, res);
  
  free_dmatrix(iden,*n);

}
	  
	 
void product_matrix(double **mat1, int *n1, int *n2, double **mat2, int *m1, int *m2, double **res)
{
  
  /** res has to be n1*n2 **/
  int i,j,k;
  double sum;
  
  for(i=0;i<*n1;i++)
    for(j=0;j<*m2;j++)
      {
	sum=0;
	for(k=0;k<*n2;k++)
	  sum+=mat1[i][k]*mat2[k][j];
	res[i][j]=sum;
      }
}

double product_vec_vec(double *vec1, int *n1, double *vec2)
{

  int i;
  double sum=0;
  
  for(i=0;i<*n1;i++)
    sum=vec1[i]*vec2[i]+sum;
  
  return(sum);
}

void product_mat_vec(double **mat, int *n1, int *n2, double *vec, double *res)
{
  
  /** res has to be n1*n2 **/
  int i,j,k;
  double sum;
  
  for(i=0;i<*n1;i++)
    {
      sum=0;
      for(j=0;j<*n2;j++)
	sum+=mat[i][j]*vec[j];
      res[i]=sum;
    }
}

double ldet(double ** x, int *n1)
/* Log determinant of square matrix */
{
  int i,j, rank, *pivot, n, p;
  double ll, tol = 1.0E-7, *qraux, *work;
  double *xtmp;
    
    
    qraux=dvector(*n1,0);
    pivot=ivector(*n1,0);
    work=dvector(2*(*n1),0);
    /** Copy the matrix by column **/
    xtmp=dvector((*n1)*(*n1),0);
    
    for(i=0;i<*n1;i++)
      for(j=0;j<*n1;j++)
	xtmp[i*(*n1)+j]=x[j][i];

    n = *n1;
    p = *n1;
    
    
    for(i = 0; i < *n1; i++)
      pivot[i] = i+1;


    F77_CALL(dqrdc2)(xtmp, &n, &n, &p, &tol, &rank,qraux, pivot, work);

    if (rank != p)
        error("Singular matrix in ldet\n");

    for (i = 0, ll=0.0; i < rank; i++) {
      ll += log(fabs(xtmp[i*(*n1)+i]));
    }

    Free(xtmp);
    Free(qraux);
    Free(pivot);
    Free(work);
    return ll;
}


double rexp_trunc(double lambda, double min, double max)
     /* generate a random deviate from a truncated exponential */
{
  
  return(-1/lambda*log(runif(exp(-lambda*max),exp(-lambda*min))));
  
}

double dexp_trunc(double x, double lambda, double min, double max)
     /* generate a random deviate from a truncated exponential */
{
  if(x<min || x>max)
    return(0);
  else
    return(lambda*exp(-lambda*x)/(exp(-lambda*min)-exp(-lambda*max)));
}

double log2(double x)
{
  return(log(x)/log(2.));
}

double log_f_ab(double sum_log_lambda, double sum_lambda, double a, double b, int n)
{
  
  double log_f;
  
  log_f=(a*a/b-1.)*sum_log_lambda-a/b*sum_lambda+n*a*a/b*log(a/b)-n*lgammafn(a*a/b)+dunif(a,0,10000.,1.)+dunif(b,0,10000.,1.);
    
  return(log_f);
   
}

double slice_sampling_b(double b0, double w, int p, double sum_log_lambda, double sum_lambda, double a, int n)
{

  double L, R, z, b1;
  int J,K;
  double log_f_R=0,log_f_L=0;
  

  z=log_f_ab(sum_log_lambda, sum_lambda,a,b0,n)-rgamma(1.,1.);
  /** lower bound **/
  L=b0-w*runif(0,1);
  R=L+w;
  J=(int)(p*runif(0,1));
  K=(p-1)-J;
  
  
  /** Find the interval **/
  log_f_R=log_f_ab(sum_log_lambda, sum_lambda,a,R,n);
  log_f_L=log_f_ab(sum_log_lambda, sum_lambda,a,L,n);

  while(J>0 & z<log_f_L)
    {
      L=L-w;
      log_f_L=log_f_ab(sum_log_lambda, sum_lambda,a,L,n);
      
      J--;
    }
  
  while(K>0 & z<log_f_R)
    {
      R=R+w;
      log_f_R=log_f_ab(sum_log_lambda, sum_lambda,a,R,n);
      
      K--;
    }
  
  /** Make sure it is in the range **/  
  L=fmax2(0.,L);
  R=fmin2(10000.,R);
  
  b1=runif(L,R);
  while(log_f_ab(sum_log_lambda, sum_lambda, a, b1, n)<z)
    {
      /** Shrink the interval **/
      if(b1<b0)
      	L=b1;
      else
      	R=b1;
      b1=runif(L,R);
    }
  
  return(b1);
  
}

double slice_sampling_a(double a0, double w, int p, double sum_log_lambda, double sum_lambda, double b, int n)
{

  double L, R, z, a1;
  int J,K;
  double log_f_R=0,log_f_L=0;
  

  z=log_f_ab(sum_log_lambda, sum_lambda,a0,b,n)-rgamma(1.,1.);
  /** lower bound **/
  L=a0-w*runif(0,1);
  R=L+w;
  J=(int)(p*runif(0,1));
  K=(p-1)-J;
  
  
  /** Find the interval **/
  log_f_R=log_f_ab(sum_log_lambda, sum_lambda,R,b,n);
  log_f_L=log_f_ab(sum_log_lambda, sum_lambda,L,b,n);

  while(J>0 & z<log_f_L)
    {
      L=L-w;
      log_f_L=log_f_ab(sum_log_lambda, sum_lambda,L,b,n);
      
      J--;
    }
  
  while(K>0 & z<log_f_R)
    {
      R=R+w;
      log_f_R=log_f_ab(sum_log_lambda, sum_lambda,R,b,n);
      
      K--;
    }
  
  /** Make sure it is in the range **/  
  L=fmax2(0.,L);
  R=fmin2(10000.,R);
  
  a1=runif(L,R);
  while(log_f_ab(sum_log_lambda, sum_lambda, a1, b, n)<z)
    {
      /** Shrink the interval **/
      if(a1<a0)
      	L=a1;
      else
      	R=a1;
      a1=runif(L,R);
    }
  
  return(a1);
  
}

double slice_sampling_rho(double rho, double w, int p, double SSR1, double SSR2, double SS12, int n)
{

  double L, R, z, rho_new;
  int J,K;
  double log_f_R=0,log_f_L=0;


  
  z=log_f_rho(SSR1, SSR2, SS12, rho, n)-rgamma(1,1);

  /** lower bound **/
  L=rho-w*runif(0,1);
  R=L+w;
  J=(int)(p*runif(0,1));
  K=(p-1)-J;
  

  /** Find the interval **/
  log_f_L=log_f_rho(SSR1, SSR2, SS12, L, n);
  log_f_R=log_f_rho(SSR1, SSR2, SS12, R, n);
  
  while(J>0 & z<log_f_L)
    {
      L=L-w;
      log_f_L=log_f_rho(SSR1, SSR2, SS12,L,n);
      
      J--;
      
    }
  
  while(K>0 & z<log_f_R)
    {
      R=R+w;
      log_f_R=log_f_rho(SSR1, SSR2, SS12,R,n);
      
      K--;
    }
  

  /** Make sure it is in the range **/
  L=fmax2(-1.,L);
  R=fmin2(1.,R);


  rho_new=runif(L,R);
  
  while(log_f_rho(SSR1, SSR2, SS12, rho_new, n)<z)
    {
      /** Shrink the interval **/
      if(rho_new<rho)
      	L=rho_new;
      else
      	R=rho_new;
      rho_new=runif(L,R);
    }
  
  return(rho_new);
  
}


double log_f_rho(double SSR1, double SSR2, double SS12, double rho, int n)
{
  
  double log_f;
  
  log_f=-n/2.*log(1-rho*rho)-1./(2.*(1-rho*rho))*(SSR1-2*rho*SS12+SSR2)+dunif(rho,-1,1,1);
  
  return(log_f);
   
}


double log_f_p(int sum_gamma_equ, int sum_gamma_diff, double prob)
{
  
  double log_f;
  
  log_f=sum_gamma_equ*log(1.-prob)+sum_gamma_diff*log(prob)+dunif(prob,0.,1.,1);
    
  return(log_f);
   
}


double slice_sampling_p(double prob0, double w, int p, int sum_gamma_equ, int sum_gamma_diff)
{

  double L, R, z, prob1;
  int J,K;
  double log_f_R=0,log_f_L=0;


  z=log_f_p(sum_gamma_equ, sum_gamma_diff, prob0)-rgamma(1,1);
  
  /** lower bound **/
  L=prob0-w*runif(0,1);
  R=L+w;
  J=(int)(p*runif(0,1));
  K=(p-1)-J;
  

  /** Find the interval **/
  log_f_R=log_f_p(sum_gamma_equ, sum_gamma_diff, R);
  log_f_L=log_f_p(sum_gamma_equ, sum_gamma_diff, L);

  /** Find the interval **/
  while(J>0 & z<log_f_L)
    {
      L=L-w;
      log_f_L=log_f_p(sum_gamma_equ, sum_gamma_diff, L);
      
      J--;
    }
  
  while(K>0 & z<log_f_R)
    {
      R=R+w;
      log_f_R=log_f_p(sum_gamma_equ, sum_gamma_diff, R);
      
      K--;
    }  


  /** Make sure it is in the range **/  
  L=fmax2(0.,L);
  R=fmin2(1.,R);
  
  
  prob1=runif(L,R);
  while(log_f_p(sum_gamma_equ, sum_gamma_diff, prob1)<z)
    {
      if(prob1<prob0)
      	L=prob1;
      else
      	R=prob1;
      prob1=runif(L,R);
    }
  
  return(prob1);
  
}
