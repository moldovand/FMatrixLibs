#include "Vector.h"

const Vector Zv(0.0, 0.0, 0.0);
const Vector Iv(1.0, 0.0, 0.0);
const Vector Jv(0.0, 1.0, 0.0);
const Vector Kv(0.0, 0.0, 1.0);

//
//     Class Vector
//

// constructor
Vector::Vector() : VectorN(3)
{
  // dummy
}

Vector::Vector(double m1, double m2, double m3) : VectorN(3)
{ 
  v[0] = m1; v[1] = m2; v[2] = m3;
}

Vector::Vector(double *a) : VectorN(a, 3)
{
  //  for (int i=0; i<3; ++i) v[i] = a[i];
}

// copy constructor
Vector::Vector(const Vector& N) : VectorN(N.v, 3)
{ 
  //  for (int i=0; i<3; ++i) v[i] = N.v[i];
}

// destructor
Vector::~Vector()
{
  // dummy
}

// assingment operator
Vector& Vector::operator=(const Vector& N)
{
  if (this != &N) 
    for (int i=0; i<3; ++i) v[i] = N.v[i];

  return *this;
}

// �٥��ȥ������
//double& Vector::operator[](int x)
//{ 
//  return v[x];
//}

// ����
int Vector::Print(FILE *fp, char *str) const 
{ 
  fprintf(fp, "%s(%13.7g,%13.7g,%13.7g )\n", str, v[0],v[1],v[2]);

  return true;
}

int Vector::LongPrint(FILE *fp, char *str) const 
{ 
  fprintf(fp, "%s(%22.16g,%22.16g,%22.16g )\n", str, v[0],v[1],v[2]);

  return true;
}

// ���� 
int Vector::Scan(FILE *fp)
{
  int    i;

  for(i=0; i<3; ++i)
    if (ScanForScan(fp) == DIGIT) fscanf(fp, "%lf", &v[i]);
    else {
      fprintf(stderr, "Scan error, in Vector.\n");
      return false;
    }
  
  return true;
}

Vector& Vector::operator+=(const Vector& N)
{
  for (int i=0; i<3; ++i) v[i] += N.v[i];

  return *this;
}
  
Vector& Vector::operator-=(const Vector& N)
{
  for (int i=0; i<3; ++i) v[i] -= N.v[i];

  return *this;
}
  
Vector& Vector::operator*=(double x)
{
  for (int i=0; i<3; ++i) v[i] *= x;

  return *this;
}
  
Vector& Vector::operator/=(double x)
{
  if (fabs(x) < eps)
    fprintf(stderr, "Division by zero, Vector/=.\n");
  else
    for (int i=0; i<3; ++i) v[i] /= x;

  return *this;
}

// ������1�������䤹
VectorN operator+(const Vector& N, double x)
{
  int i;
  VectorN ret(4);

  for (i=0; i<3; ++i) ret[i] = N.v[i];

  ret[i] = x;

  return ret;
}

VectorN operator+(double x, const Vector& N)
{
  VectorN ret(4);

  ret[0] = x;
  for (int i=1; i<4; ++i) ret[i] = N.v[i-1];

  return ret;
}

// �٥��ȥ��ñ��ޥ��ʥ�
Vector operator-(const Vector& N)
{
  return Vector(-N.v[0], -N.v[1], -N.v[2]);
}

// �٥��ȥ�Ʊ�Τκ�
Vector operator-(const Vector& N1, const Vector& N2)
{ 
  return Vector(N1.v[0] - N2.v[0], N1.v[1] - N2.v[1], N1.v[2] - N2.v[2]);
}

// �٥��ȥ�Ʊ�Τ���
Vector operator+(const Vector& N1, const Vector& N2)
{ 
  return Vector(N1.v[0] + N2.v[0], N1.v[1] + N2.v[1], N1.v[2] + N2.v[2]);
}

// �٥��ȥ�������
Vector operator*(double r, const Vector& N)
{
  return Vector(N.v[0] * r, N.v[1] * r, N.v[2] * r);
}

// �٥��ȥ�������
Vector operator*(const Vector& N, double r)
{ 
  return Vector(N.v[0] * r, N.v[1] * r, N.v[2] * r);
}

// �٥��ȥ�γ���
Vector operator%(const Vector& N1, const Vector& N2)
{
  return Vector(N1.v[1] * N2.v[2] - N1.v[2] * N2.v[1],
		N1.v[2] * N2.v[0] - N1.v[0] * N2.v[2],
		N1.v[0] * N2.v[1] - N1.v[1] * N2.v[0]);
}

// �٥��ȥ�ξ�
Vector operator/(const Vector& N, double r)
{
  if (fabs(r) < eps) {
    fprintf(stderr, "Error : division by zero, Vector/\n");
    return Zv;
  } 
  else return Vector(N.v[0] / r, N.v[1] / r, N.v[2] / r);
}

// �٥��ȥ�Ʊ�Τ����
int operator==(const Vector& N1, const Vector& N2)
{ 
  if (fabs(Norm(N1 - N2)) < eps) return true;
  else return false;
}

// �٥��ȥ�Ʊ�Τ����
int operator!=(const Vector& N1, const Vector& N2)
{ 
  if (N1 == N2) return false;
  else return true;
}

// ����
double operator,(const Vector& N1, const Vector& N2)
{
  double x=0.0;

  for(int i=0; i<3; ++i) x += N1.v[i] * N2.v[i];

  return x;
}

// �����飳����
double ScalarTP(const Vector& N1, const Vector& N2, const Vector& N3)
{
  return ((N1 % N2), N3);
}

// �٥��ȥ�ΥΥ��
double Norm(const Vector& N)
{
  return sqrt(sqr(N.v[0]) + sqr(N.v[1]) + sqr(N.v[2]));
}

// �٥��ȥ��������
Vector Normalize(const Vector& N)
{ 
  double d;

  d = Norm(N);

  if (fabs(d) < eps) return Zv;
  else
    return Vector(N.v[0] / d, N.v[1] / d, N.v[2] / d);
}

// Ʃ�������
Vector Proj(const Vector& N)
{
  if (fabs(N.v[2]) > eps)
    return Vector(N.v[0] / N.v[2], N.v[1] / N.v[2], 1.0);
  else {
    fprintf(stderr, "Can't projection.\n");
    return N;
  }
}
