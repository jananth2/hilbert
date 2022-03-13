/*
 * CS371 Computing Project 2 -- Due 3/12/22
 * Joshua Anantharaj and Nick Drake
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Prints a n x n matrix to stdout
void print_matrix(double **M, const size_t n){
  for(size_t i = 0; i < n; i++){
    for(size_t j = 0; j < n; j++){
      printf("%f ", M[i][j]);
    }
    printf("\n");
  }
}

void generate_h(double **H, const size_t n){
  for(size_t i = 0; i < n; i++){
    for(size_t j = 0; j < n; j++){
      H[i][j] = 1.0/(i+j+1.0);
    }
  }
}

void print_vector(double *x, const size_t n){
  for(size_t i = 0; i < n; i++){
    printf("%f ", x[i]);
  }
  printf("\n");
}

//Poplulates ALLOCATED 2d matrix of doubles L
void choleski(double **A, double **L, const size_t n){
  for(size_t i = 0; i < n; i++){
    for(size_t j = 0; j <= i; j++){
      double s = A[i][j];
      for(size_t k = 0; k < j; k++){
        s -= L[j][k]*L[i][k];
      }
      if(j<i){
        L[i][j] = L[j][i] = s/L[j][j];
      } else {
        L[i][j] = sqrt(s);
      }
    }
  }
}

void fwd_sub(double **L, double *b, size_t n, double *y){
  for(size_t i = 0; i < n; i++){
    y[i] = b[i];
    for(size_t j = 0; j < i; j++){
      y[i] -= L[i][j]*y[j];
    }
    y[i] /= L[i][i];
  }
}

void back_sub(double **U, double *y, size_t n, double *x){
  for(int i = (n-1); i >= 0; i--){
    x[i] = y[i];
    for(int j = i+1; j < n; j++){
      x[i] -= U[i][j]*x[j];
    }
    x[i] /= U[i][i];
  }
}

double max_norm(double **H, size_t n){
  double m = 0;
  for(size_t i = 0; i < n; i++){
    double ar = 0;
    for(size_t j = 0; j < n; j++){
      ar += fabs(H[i][j]);
    }
    m = fmax(m, ar);
  }
  return m;
}

void generate_inverse_h(double **H, double **I, size_t n){
  for(int i = 1; i <= n; i++){
    for(int j = 1; j <= n; j++){
      double n1 = pow(-1, i+j);
      double n2 = 1, n3 = 1, n4 = 1, n5 = 1, n6 = 1, n7 = 1;
      double n8 = i+j-1;
      for(int x = 1; x <= i+n-1
                     || x <= j+n-1
                     || x <= i-1
                     || x <= j-1
                     || x <= n-j
                     || x <= n-i; x++){
        if(x <= i+n-1) n2 *= x;
        if(x <= j+n-1) n3 *= x;
        if(x <= i-1)   n4 *= x;
        if(x <= j-1)   n5 *= x;
        if(x <= n-j)   n6 *= x;
        if(x <= n-i)   n7 *= x;
      }
      I[i-1][j-1] = (n1*n2*n3)/(n4*n5*n6*n7*n8);
    }
  }
}

void calculate_final(double **X, double **I, const size_t n, double **O){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      O[i][j] = I[i][j]-X[i][j];
    }
  }
}

//O(n) allocation time vs O(1) for 1d
//2d indexing makes worth?
//Returns a double pointer to memory allocated for n x n matrix
double** alloc_matrix(const size_t n, const char type){
  int row_len = n;
  double** ret = (double**)malloc(8*n);
  for(size_t i = 0; i < n; i++){
    if(type == 'L') row_len = i+1;
    if(type == 'U') row_len = n-i;
    ret[i] = (double*)malloc(8*row_len);
  }
  return ret;
}

int main(){
  const size_t n = 8;

  double **H = alloc_matrix(n, 'S');
  double **L = alloc_matrix(n, 'S');
  //double **U = alloc_matrix(n, 'U');
  double **I = alloc_matrix(n, 'S');
  double **X = alloc_matrix(n, 'S');
  double **O = alloc_matrix(n, 'S');

  printf("Creating Hilbert Matrix H:\n");
  generate_h(H, n);
  print_matrix(H, n);

  double h_max = max_norm(H, n);
  printf("H has Max Norm of %f\n", h_max);

  printf("\nComputing Inverse Hilbert Matrix I:\n");
  generate_inverse_h(H, I, n);
  print_matrix(I, n);
  double i_max = max_norm(I, n);
  printf("I has Max Norm of %f\n", i_max);

  double cond_num = i_max*h_max;
  printf("Actual condition number: %f\n\n", cond_num);

  printf("Computing A = LL^T Choleski Factorization of H, L:\n");
  choleski(H, L, n);
  print_matrix(L, n);

  for(size_t i = 0; i < n; i++){
    double *e = (double *)calloc(n, 8);
    double *tmp = (double *)malloc(8*n);
    e[i] = 1;
    fwd_sub(L, e, n, tmp);
    back_sub(L, tmp, n, X[i]);
  }
  printf("\nCalculated X:\n");
  print_matrix(X, n);
  calculate_final(X, I, n, O);
  printf("\nDifference in Norm: %f\n", max_norm(O, n));

  return 0;
}
