#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "matrix.h"

int read_matrix(double *a, int n, char *filename)
{
    FILE *fp;
    int i,j;
    if (!(fp = fopen(filename,"r")))
        return ERROR_OPEN;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (fscanf(fp, "%lf", a+i*n+j) != 1)
            {
                fclose(fp);
                return ERROR_READ;
            }
    fclose(fp);
    return SUCCESS;
}

double f(int s, int n, int i, int j)
{
    i++; j++;
    if (s == 1)
        return n - (i > j ? i : j) + 1;
    if (s == 2)
    {
        if (i == j)
          return 2;
        if (abs(i - j) == 1)
          return -1;
        return 0;
    }
    if (s == 3)
    {
        if (i == j && j < n)
          return 1;
        if (j == n)
          return i;
        if (i == n)
          return j;
        return 0;
    }
    if (s == 4)
    {
        return (double)1/(i+j-1); 
    }
    return 1;
}

void init_matrix(double *a, int n, int s)
{
    int i,j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            a[i*n+j] = f(s,n,i,j);
}

void print_matrix(double *a, int m, int n, int r)
{
    double n_max = (n > r ? r : n);
    double m_max = (m > r ? r : m);
    for (int i = 0; i < m_max; i++)
    {
        for (int j = 0; j < n_max; j++)
        {
            printf(" %10.3e", a[i*n+j]);
        }
        printf("\n");
    }
}

double norma_vector(double *b, int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
    	double check = fabs(b[i]);
        if (check > res)
        	res = check;
    }
    return res;
}


//horizontal norm
double norma_matrix(double *a, int m, int n)
{
  double res = 0;
  for (int i = 0; i < m; i++)
  {
    double check = 0;
    for (int j = 0; j < n; j++)
    {
      check += fabs(a[i*n + j]);
    }
    if (check > res)
    {
      res = check;
    }
  }
  return res;
}


static void multiply_U_and_a(double *x, double *a, int n, int start)
{
  for (int j = start; j < n; j++)
  {
    double scal = 0;
    for (int i = start + 1; i < n; i++)
      scal += a[i * n + j] * x[i - start - 1];
    for (int i = start + 1; i < n; i++)
      a[i * n + j] -= 2 * x[i - start - 1] * scal;
  }
}

static void multiply_a_and_U(double *a, double *x, int n, int start)
{
  for (int i = 0; i < n; i++)
  {
    double scal = 0;
    for (int j = start + 1; j < n; j++)
      scal += x[j - start - 1] * a[i * n + j];
    for (int j = start + 1; j < n; j++)
      a[i * n + j] -= 2 * x[j - start - 1] * scal;
  }
}


int to_partly_triangular_reflection(double *a, int n, double eps, double norm)
{
  double *x = new double[n];
  for (int k = 0; k < n - 2; k++)
  {
    if (fabs(a[(k + 2) * n + k]) < eps * norm)
    {
      bool check = true;
      for (int i = k + 3; i < n && check; i++)
        if (fabs(a[i * n + k]) > eps * norm)
          check = false;
      if (check)
        continue;
    }
    double s_k = 0;
    for (int j = k + 2; j < n; j++)
      s_k += a[j * n + k] * a[j * n + k];
    double a1_norm = a[(k + 1) * n + k] * a[(k + 1) * n + k] + s_k;
    if (a1_norm < 0)
    {
      delete [] x;
      return NEGATIVE_SQUARE;
    }
    a1_norm = sqrt(a1_norm);
    x[0] = a[(k + 1) * n + k] - a1_norm;
    //CHECK
    for (int j = 1; j < n - k - 1; j++)
      x[j] = a[(k + 1 + j) * n + k];
    double x_norm = x[0] * x[0] + s_k;
    if (x_norm < 0)
    {
      delete [] x;
      return NEGATIVE_SQUARE;
    }
    x_norm = sqrt(x_norm);
    if (x_norm < eps * norm)
    {
      bool check = true;
      for (int i = k + 1; i < n && check; i++)
      {
        for (int j = 0; j < i && check; j++)
        {
          if (fabs(a[i * n + j]) > eps * norm)
            check = false;
        }
      }
      delete [] x;
      if (check)
        return SUCCESS;
      else
        return DIVISION_ZERO;
    }
    for (int j = 0; j < n - k - 1; j++)
      x[j] /= x_norm;

    multiply_U_and_a(x, a, n, k);
    multiply_a_and_U(a, x, n, k);
  }
  delete [] x;
  return SUCCESS;
}


int qr_rotate(double *a, double *cos, double *sin, int n, int size, double eps, double norm)
{
  for (int k = 0; k < size - 1; k++)
  {
    double x = a[k * n + k];
    double y = a[(k + 1) * n + k];
    if (fabs (y) < eps * norm)
    {
      a[(k + 1) * n + k] = 0;
      cos[k] = 1;
      sin[k] = 0;
      continue;
    }
    double value = x * x + y * y;
    if (value < 0)
      return NEGATIVE_SQUARE;
    double square_root = sqrt(value);
    if (square_root < eps * norm)
      return DIVISION_ZERO;
    double cos_k = x / square_root;
    double sin_k = -y / square_root;
    cos[k] = cos_k;
    sin[k] = sin_k;
    a[k * n + k] = a[k * n + k] * cos_k - a[(k + 1) * n + k] * sin_k;
    a[(k + 1) * n + k] = 0;
    for (int j = k + 1; j < size; j++)
    {
      double up = a[k * n + j];
      double down = a[(k + 1) * n + j];
      a[k * n + j] = up * cos_k - down * sin_k;
      a[(k + 1) * n + j] = up * sin_k + down * cos_k;
    }
  }
  return SUCCESS;
}


int find_eigenvalues(double *egn, double *a, int n, double eps, double norm, double & t1, double & t2, int & it)
{
  t1 = clock();
  int ret = to_partly_triangular_reflection(a, n, eps, norm);
  t1 = (clock() - t1) / CLOCKS_PER_SEC;
  if (ret != SUCCESS)
    return ret;
  double *cos = new double[n - 1];
  double *sin = new double[n - 1];
  
  it = 0;
  t2 = clock();
  for (int size = n; size > 2; size--)
  {
    while (fabs (a[(size - 1) * n + size - 2]) > eps * norm)
    {
      double shift = a[(size - 1) * n + size - 1] - 0.5 * a[(size - 1) * n + size - 2];
      for (int i = 0; i < size; i++)
      {
        a[i * n + i] -= shift;
      }
      ret = qr_rotate(a, cos, sin, n, size, eps, norm);
      if (ret != SUCCESS)
      {
        delete [] cos;
        delete [] sin;
        return ret;
      }
      //В a лежит Q_k, теперь надо а умножить сначала на T_0^t потом T_1^t и т.д.  
      for (int j = 0; j < size - 1; j++)
      {
        double cosa = cos[j];
        double sina = sin[j];
        int ind = j;
        //A
        for (int i = 0; i <= ind; i++)
        {
          double a_i1 = a[i * n + ind];
          double a_i2 = a[i * n + ind + 1];
          a[i * n + ind] = a_i1 * cosa - a_i2 * sina;
          a[i * n + ind + 1] = a_i1 * sina + a_i2 * cosa;
        }
        //B
        double a_i2 = a[(ind + 1) * n + ind + 1];
        a[(ind + 1) * n + ind] = - a_i2 * sina;
        a[(ind + 1) * n + ind + 1] = a_i2 * cosa;
      }
      for (int i = 0; i < size; i++)
      {
        a[i * n + i] += shift;
      }
      it++;
    }
    egn[size - 1] = a[(size - 1) * n + size - 1];
  }
  double A = a[0];
  double B = a[1];
  double C = a[n];
  double D = a[n + 1];
  double discr = (A + D) * (A + D) - 4 * (A * D - B * C);
  if (discr < 0)
  {
    delete [] cos;
    delete [] sin;
    return NEGATIVE_SQUARE;
  }
  double sqrt_discr = sqrt(discr);
  egn[0] = (A + D + sqrt_discr) / 2;
  egn[1] = (A + D - sqrt_discr) / 2;
  t2 = (clock() - t2) / CLOCKS_PER_SEC;
  delete [] cos;
  delete [] sin;
  return SUCCESS;
}



