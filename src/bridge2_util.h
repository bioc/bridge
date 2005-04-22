#include <R.h>
#include <Rmath.h>
//#include "/usr/local/lib/R/include/R_ext/Applic.h"	/* Fortran routines */
#include <float.h>


double stdd(double *vector,int *length, int *is_finite);
double **dmatrix(int nb_row,int nb_col);
void mat_vec(double *array_vec,int* nb_row,int *nb_col,double **array);
void vec_mat(double *array_vec,int* nb_row,int *nb_col,double **array);
double *dvector(int length, int init);
int *ivector(int length, int init);
double  mean_vec(double *vector,int *length);
void free_dmatrix(double **array, int nb_row);
void init_ivector(int *vector, int *length, int  value);
void quicksort2(double *a, int *b, int *p, int *r);
int partition2(double *a, int *b, int p, int r);
int rand_part2(double *a, int *b,  int p, int r);
void idquicksort2(int *a, double *b, int *p, int *r);
int idpartition2(int *a, double *b, int p, int r);
int idrand_part2(int *a, double *b,  int p, int r);
int uni_rand(int min,int max);
void init_dvector(double *vector, int *length, int value);
double  sum_vec(double *vector,int *length);
void quicksort(double *a, int *p, int *r);
int partition(double *a, int p, int r);
int rand_part(double *a, int p, int r);
double dabs(double a);
double dmax(double a, double b);
void qr_solve(double **x, int *n1, double **y, double ** coef);
void inverse(double **mat1, int *n ,double **res);
void product_matrix(double **mat1, int *n1, int *n2, double **mat2, int *m1, int *m2, double **res);
void product_mat_vec(double **mat1, int *n1, int *n2, double *vec2, double *res);
double product_vec_vec(double *vec1, int *n1, double *vec2);
void product_mat_vec(double **mat, int *n1, int *n2, double *vec, double *res);
double ldet(double ** x, int *n1);
double rexp_trunc(double lambda, double min, double max);
double dexp_trunc(double x, double lambda, double min, double max);
double log2(double x);
double log_f_ab(double sum_log_lambda, double sum_lambda, double a, double b, int n);
double slice_sampling_rho(double rho, double w, int p, double SSR1, double SSR2, double SS12, int n);
double log_f_rho(double SSR1, double SSR2, double SS12, double rho, int n);
double slice_sampling_b(double b0, double w, int p, double sum_log_lambda, double sum_lambda, double a, int n);
double slice_sampling_a(double a0, double w, int p, double sum_log_lambda, double sum_lambda, double b, int n);
double slice_sampling_p(double prob0, double w, int p, int sum_gamma_equ, int sum_gamma_diff);
double log_f_p(int sum_gamma_equ, int sum_gamma_diff, double prob);
