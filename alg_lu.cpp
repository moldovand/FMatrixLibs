/***********************************************************
	alg_lu.c -- LUʬ��
***********************************************************/
#include <stdlib.h>
#include <math.h>
#include "alg_matutil.h"  /* �������ξ�ƻ�� */

double _alg_lu(int n, matrix a, int *ip)
{
  int i, j, k, ii, ik;
  double t, u, det;
  vector weight;
  weight = _alg_new_vector(n);    /* weight[0..n-1] �ε����ΰ���� */
  det = 0;                   /* ���� */
  for (k = 0; k < n; k++) {  /* �ƹԤˤĤ��� */
    ip[k] = k;             /* �Ը򴹾���ν���� */
    u = 0;                 /* ���ιԤ������ͺ�������Ǥ���� */
    for (j = 0; j < n; j++) {
      t = fabs(a[k][j]);  if (t > u) u = t;
    }
    if (u == 0) goto EXIT; /* 0 �ʤ�����LUʬ��Ǥ��ʤ� */
    weight[k] = 1 / u;     /* ���������ͤεտ� */
  }
  det = 1;                   /* ���󼰤ν���� */
  for (k = 0; k < n; k++) {  /* �ƹԤˤĤ��� */
    u = -1;
    for (i = k; i < n; i++) {  /* ��겼�γƹԤˤĤ��� */
      ii = ip[i];        /* �Ťߡ������� ������ιԤ򸫤Ĥ��� */
      t = fabs(a[ii][k]) * weight[ii];
      if (t > u) {  u = t;  j = i;  }
    }
    ik = ip[j];
    if (j != k) {
      ip[j] = ip[k];  ip[k] = ik;  /* ���ֹ��� */
      det = -det;  /* �Ԥ�򴹤���й��󼰤���椬�Ѥ�� */
    }
    u = a[ik][k];  det *= u;  /* �г���ʬ */
    if (u == 0) goto EXIT;    /* 0 �ʤ�����LUʬ��Ǥ��ʤ� */
    for (i = k + 1; i < n; i++) {  /* Gauss�õ�ˡ */
      ii = ip[i];
      t = (a[ii][k] /= u);
      for (j = k + 1; j < n; j++)
	a[ii][j] -= t * a[ik][j];
    }
  }
EXIT:
  _alg_free_vector(weight);  /* �����ΰ����� */
  return det;                /* ����ͤϹ��� */
}

void _alg_solve(int n, matrix a, vector b, int *ip, vector x)
{
  int i, j, ii;
  double t;
  for (i = 0; i < n; i++) {       /* Gauss�õ�ˡ�λĤ� */
    ii = ip[i];  t = b[ii];
    for (j = 0; j < i; j++) t -= a[ii][j] * x[j];
    x[i] = t;
  }
  for (i = n - 1; i >= 0; i--) {  /* �������� */
    t = x[i];  ii = ip[i];
    for (j = i + 1; j < n; j++) t -= a[ii][j] * x[j];
    x[i] = t / a[ii][i];
  }
}

double _alg_matinv(int n, matrix a, matrix a_inv)
{
  int i, j, k, ii;
  double t, det;
  int *ip;   /* �Ը򴹤ξ��� */
  ip = (int *)malloc(sizeof(int) * n);
  if (ip == NULL) _alg_error("�����ΰ���­");
  det = _alg_lu(n, a, ip);
  if (det != 0)
    for (k = 0; k < n; k++) {
      for (i = 0; i < n; i++) {
	ii = ip[i];  t = (ii == k);
	for (j = 0; j < i; j++)
	  t -= a[ii][j] * a_inv[j][k];
	a_inv[i][k] = t;
      }
      for (i = n - 1; i >= 0; i--) {
	t = a_inv[i][k];  ii = ip[i];
	for (j = i + 1; j < n; j++)
	  t -= a[ii][j] * a_inv[j][k];
	a_inv[i][k] = t / a[ii][i];
      }
    }
  free(ip);
  return det;
}

double _alg_gauss(int n, matrix a, vector b, vector x)
{
  double det;  /* ���� */
  int *ip;     /* �Ը򴹤ξ��� */
  ip = (int *)malloc(sizeof(int) * n);       /* �����ΰ���� */
  if (ip == NULL) _alg_error("�����ΰ���­");
  det = _alg_lu(n, a, ip);                   /* LUʬ�� */
  if (det != 0) _alg_solve(n, a, b, ip, x);  /* LUʬ��η�̤�Ȥä�ϢΩ��������� */
  free(ip);                                  /* �����ΰ�β��� */
  return det;                                /* ����ͤϹ��� */
}
double _alg_det(int n, matrix a)
{
  double det;  /* ���� */
  int *ip;     /* �Ը򴹤ξ��� */
  ip = (int *)malloc(sizeof(int) * n);       /* �����ΰ���� */
  if (ip == NULL) _alg_error("�����ΰ���­");
  det = _alg_lu(n, a, ip);                   /* LUʬ�� */
  free(ip);                                  /* �����ΰ�β��� */
  return det;                                /* ���� */
}
