#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <iostream>
#include "matrix.h"

double residual1(double *a, double *egn, int n, double eps, double norm)
{
  double ret = 0;
  for (int i = 0; i < n; i++)
    ret += a[i * n + i] - egn[i];
  ret = fabs(ret);
  if (fabs (norm) > eps * norm)
    return ret / norm;
  return ret;
}

double residual2(double *a, double *egn, int n, double eps, double norm)
{
  double sqrt1 = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      sqrt1 += a[i * n + j] * a[j * n + i];
    }
  }
  if (sqrt1 < 0)
    return 0;
  sqrt1 = sqrt(sqrt1);
  double sqrt2 = 0;
  for (int i = 0; i < n; i++)
    sqrt2 += egn[i] * egn[i];
  if (sqrt2 < 0)
    return 0;
  sqrt2 = sqrt(sqrt2);
  if (fabs(norm) > eps * norm)
    return fabs(sqrt1 - sqrt2) / norm;
  return fabs(sqrt1 - sqrt2);
}

int main(int argc, char *argv[])
{
    int n, m, k;
    double eps;
    char * filename = nullptr;
    if (!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n) == 1 && 
        sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%lf", &eps) == 1 && sscanf(argv[4], "%d", &k) == 1)) 
    {
        printf("Usage %s n m eps k (file) \n", argv[0]);
        return 1;
    }
    if (n <= 0)
    {
        printf("n should be at less 1\n");
        return 1;
    }
    
    if (argc == 6)
      filename = argv[5];
    //Исходная матрица, которую не будем преобразовывать, чтобы сверить ответ
    double *a = new double[n*n];
    if (!a)
    {
        printf("Not enough memory\n");
        return 2;
    }
    if (filename)
    {
        int ret = read_matrix(a,n,filename);
        if (ret != SUCCESS)
        {
            switch (ret)
            {
                case ERROR_OPEN:
                    printf("Cannot open %s\n", filename);
                    break;
                case ERROR_READ:
                    printf("Cannot read %s\n", filename);
                    break;
                default:
                    printf("Unknown error %d in file %s\n", ret, filename);
            }
            delete [] a;
            return 3;
        }
    }
    else
        init_matrix(a,n,k);
    
    double *egn = new double[n];
    print_matrix(a, n, n, m);
    double norm = norma_matrix(a, n, n);
    double t1 = 0;
    double t2 = 0;
    int its = 0;
    int ret = find_eigenvalues(egn, a, n, eps, norm, t1, t2, its);
    if (ret != SUCCESS)
    {
      printf("Can't find eigenvalues\n");
      delete [] egn;
      delete [] a;
      return 3;
    }
    printf("RESULT : \n");
    print_matrix(egn, 1, n, m);
    
    if (filename)
    {
        int ret = read_matrix(a,n,filename);
        if (ret != SUCCESS)
        {
            switch (ret)
            {
                case ERROR_OPEN:
                    printf("Cannot open %s\n", filename);
                    break;
                case ERROR_READ:
                    printf("Cannot read %s\n", filename);
                    break;
                default:
                    printf("Unknown error %d in file %s\n", ret, filename);
            }
            delete [] egn;
            delete [] a;
            return 4;
        }
    }
    else
        init_matrix(a,n,k);
    
    double res1 = residual1(a, egn, n, eps, norm);
    double res2 = residual2(a, egn, n, eps, norm);
    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, its, its / n, t1, t2);
    
    delete [] a;
    delete [] egn;
    return 0;
}
