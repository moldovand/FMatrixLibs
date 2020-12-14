#include <stdlib.h>
#include "Matrix.h"
#include "alg.h"

//
//     Class Matrix (3-Dimension Matrix)
//

const Matrix Im(1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0);

const Matrix Zm(0.0, 0.0, 0.0,
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0);

const Matrix Pk(1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 0.0);

// constructor
Matrix::Matrix() : MatrixMN(3, 3)
{
  // dummy
}

Matrix::Matrix(const Vector& N1, const Vector& N2, const Vector& N3) : MatrixMN(3, 3)
{
  m[0][0] = N1[0]; m[0][1] = N2[0]; m[0][2] = N3[0]; 
  m[1][0] = N1[1]; m[1][1] = N2[1]; m[1][2] = N3[1];
  m[2][0] = N1[2]; m[2][1] = N2[2]; m[2][2] = N3[2];
}

// ��A�β���angle�٤β�ž����
Matrix::Matrix(const Vector& A, double angle) : MatrixMN(3, 3)
{ 
  double c, s;
  Vector B;

  B = Normalize(A);
  
  angle = angle / 180.0 * M_PI;
  c = cos(angle); s = sin(angle);
  m[0][0] = c + B[0] * B[0] * (1.0 - c); 
  m[0][1] = B[0] * B[1] * (1.0 - c) - B[2] * s;
  m[0][2] = B[0] * B[2] * (1.0 - c) + B[1] * s;
  m[1][0] = B[1] * B[0] * (1.0 - c) + B[2] * s; 
  m[1][1] = c + B[1] * B[1] * (1.0 - c);
  m[1][2] = B[1] * B[2] * (1.0 - c) - B[0] * s;
  m[2][0] = B[2] * B[0] * (1.0 - c) - B[1] * s;
  m[2][1] = B[2] * B[1] * (1.0 - c) + B[0] * s;
  m[2][2] = c + B[2] * B[2] * (1.0 - c); 
}

Matrix::Matrix(double* a) : MatrixMN(a, 3, 3)
{
  // dummy
}

Matrix::Matrix(double** a) : MatrixMN(a, 3, 3)
{
  // dummy
}

Matrix::Matrix(double a11, double a12, double a13,
	       double a21, double a22, double a23,
	       double a31, double a32, double a33) : MatrixMN(3, 3)
{
  m[0][0] = a11; m[0][1] = a12; m[0][2] = a13;
  m[1][0] = a21; m[1][1] = a22; m[1][2] = a23;
  m[2][0] = a31; m[2][1] = a32; m[2][2] = a33;
}

Matrix::Matrix(const Vector& N1, const Vector& N2) : MatrixMN(3, 3)
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) m[i][j] = N1[i] * N2[j];
}

// �гѹ�����г�����
Matrix::Matrix(double d1, double d2, double d3) : MatrixMN(3, 3)
{
  m[0][0] = d1; m[1][1] = d2; m[2][2] = d3;
}

// copy constructor
Matrix::Matrix(const Matrix& M) : MatrixMN(3, 3)
{
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j) m[i][j] = M.m[i][j];
}

// destructor
Matrix::~Matrix()
{
  // dummy
}

// assignment operator
Matrix& Matrix::operator=(const Matrix& M)
{
  if (this != &M) 
    for (int i=0; i<3; ++i) 
      for (int j=0; j<3; ++j) m[i][j] = M.m[i][j];
  
  return *this;
}

// i�Ԥιԥ٥��ȥ�
double* Matrix::operator[](int i)
{
  if ((0 <= i) && (i <= 2)) return m[i];
  else{
    fprintf(stderr, "Domain error in Matrix[].\n");
    return NULL;
  }
}

const double* Matrix::operator[](int i) const
{
  if ((0 <= i) && (i <= 2)) return m[i];
  else{
    fprintf(stderr, "Domain error in Matrix[].\n");
    return NULL;
  }
}

// j�����٥��ȥ�
Vector Matrix::operator()(int j) const
{
  if ((0 <= j) && (j <= 2)) return Vector(m[0][j], m[1][j], m[2][j]);
  else{
    fprintf(stderr, "Domain error in Matrix().\n");
    return Zv;
  }
}

// Matrix�����
int Matrix::Print(FILE *fp, char *str) const 
{ 
  char sp[20];
  int  i, len;
  
  len = strlen(str);
  for (i = 0; i < len; ++i) sp[i] = ' ';
  sp[len] = '\0';
  
  for (i=0; i<3; ++i) {
    fprintf(fp, "%s|", (i == 1) ? str : sp);
    for(int j = 0; j < 3; ++j) fprintf(fp, "%13.7g ", m[i][j]);
    fprintf(fp, "|\n");
  }

  return true;
}

// Matrix����� (Ĺ����)
int Matrix::LongPrint(FILE *fp, char *str) const 
{ 
  char sp[20];
  int  i, len;
  
  len = strlen(str);
  for (i = 0; i < len; ++i) sp[i] = ' ';
  sp[len] = '\0';
  
  for (i=0; i<3; ++i) {
    fprintf(fp, "%s|", (i == 1) ? str : sp);
    for(int j = 0; j < 3; ++j) fprintf(fp, "%22.16g ", m[i][j]);
    fprintf(fp, "|\n");
  }

  return true;
}

// ���������
int Matrix::Scan(FILE *fp)
{
  int    i, j;
  
  for (i=0; i<3; ++i) 
    for (j=0; j<3; ++j)
      if (ScanForScan(fp) == DIGIT) fscanf(fp, "%lf", &m[i][j]);
      else {
	fprintf(stderr, "Scan error, in Matrix.\n");
	return false;
      }

  return true;
}

// �оι��󤫤ɤ���
int Matrix::Symmetric() const
{
  int ret = true;

  for (int i=0; i<2; ++i)
    for (int j=i+1; j<3; ++j)
      if (fabs(m[i][j] - m[j][i]) > eps) {
        ret = false;
        break;
      }
  
  return ret;
}

// �����
int Matrix::Init()
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) m[i][j] = 0.0;

  return true;
}

// ñ�̲�
int Matrix::Unit()
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      if (i == j) m[i][j] = 1.0;
      else        m[i][j] = 0.0;

  return true;
}

// ���
int operator==(const Matrix& M1, const Matrix& M2)
{
  if (fabs(Norm(M1 - M2)) <= eps) return true;
  else                            return false;
}

int operator!=(const Matrix& M1, const Matrix& M2)
{
  if (M1 == M2) return false;
  else          return true;
}

Matrix& Matrix::operator+=(const Matrix& M)
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) m[i][j] += M.m[i][j];

  return *this;
}

Matrix& Matrix::operator-=(const Matrix& M)
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) m[i][j] -= M.m[i][j];

  return *this;
}

Matrix& Matrix::operator*=(double x)
{
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) m[i][j] *= x;

  return *this;
}

Matrix& Matrix::operator/=(double x)
{
  if (fabs(x) < eps)
    fprintf(stderr, "Division by ZERO, Matrix/=.\n");
  else
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j) m[i][j] /= x;

  return *this;
}

// ����ȥ٥��ȥ����
MatrixMN operator+(const Matrix& M, double x)
{
  MatrixMN ret(4, 4);
  
  for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
      if ((i == 3) || (j == 3)) ret[i][j] = 0.0;
      else    	                ret[i][j] = M.m[i][j];

  ret[3][3] = x;

  return ret;
}
  
MatrixMN operator+(double x, const Matrix& M)
{
  MatrixMN ret(4, 4);
  
  for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
      if ((i == 0) || (j == 0))	ret[i][j] = 0.0;
      else     	                ret[i][j] = M.m[i - 1][j - 1];

  ret[0][0] = x;

  return ret;
}
  
// ����ȥ٥��ȥ����
Vector operator*(const Matrix& M, const Vector& N)
{ 
  Vector R;
  
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) R[i] += M.m[i][j] * N[j];
 
  return R;
}

// ����ȹ������
Matrix operator*(const Matrix& M1, const Matrix& M2)
{
  Matrix R;

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      for (int k=0; k<3; ++k) R.m[i][j] += M1.m[i][k] * M2.m[k][j];

  return R;
}

// �٥��ȥ�ȹ������
Matrix operator%(const Vector& N, const Matrix& M)
{
  return Matrix(N % M(0), N % M(1), N % M(2));
}

// ����ȥ٥��ȥ�γ���
Matrix operator%(const Matrix& M, const Vector& N)
{
  return Trans(N % Trans(M));
}

// ����ȹ���γ���
Matrix operator%(const Matrix& M1, const Matrix& M2)
{
  Matrix ret;

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      for (int k=0; k<3; ++k)
	for (int l=0; l<3; ++l)
	  for (int mm=0; mm<3; ++mm)
	    for (int n=0; n<3; ++n) 
	      ret.m[i][j] += EddingtonEpsilon(i, k, l) *
   		             EddingtonEpsilon(j, mm, n) *
			     M1.m[k][mm] * M2.m[l][n];

  return ret;
}

// ����ȥ����顼����
Matrix operator*(double r, const Matrix& M)
{
  Matrix R;
  
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) R.m[i][j] = r * M.m[i][j];
 
  return R;
}

Matrix operator*(const Matrix& M, double r)
{
  return r * M;
}

// ����Υ�����ʬ��1
Matrix operator/(const Matrix& M, double r)
{
  if (fabs(r) > eps) return (1 / r) * M;
  else {
    fprintf(stderr, "Error : division by zero in Matrix/\n");
    return Zm;
  }
}

// �������
Matrix operator+(const Matrix& M1, const Matrix& M2)
{
  Matrix R;
  
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) R.m[i][j] = M1.m[i][j] + M2.m[i][j];
 
  return R;
}

// �����ñ��ޥ��ʥ�
Matrix operator-(const Matrix& M)
{
  Matrix R;
  
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) R.m[i][j] = -M.m[i][j];
 
  return R;
}

Matrix operator-(const Matrix& M1, const Matrix& M2)
{
  return M1 + (-M2);
}

// �ͱƹ���
Matrix Projection(const Vector& a)
{
  return Im - Matrix(Normalize(a), Normalize(a));
}

// ���������
double operator,(const Matrix& M1, const Matrix& M2)
{
  double x = 0.0;

  for (int i=0; i<3; ++i) 
    for (int j=0; j<3; ++j) x += M1.m[i][j] * M2.m[i][j];

  return x;
}

// ����ΥΥ��
double Norm(const Matrix& M)
{
  double x = 0.0;

  for (int i=0; i<3; ++i) 
    for (int j=0; j<3; ++j) x += sqr(M.m[i][j]);

  return sqrt(x);
}

// ������
Matrix Normalize(const Matrix& M)
{
  return M / Norm(M);
}

// ���
int Rank(const Matrix& M)
{
  Matrix U;
  Vector L;
  int    i, ret = 0;

  (M*Trans(M)).Householder(L,U);

  for (i=0; i<3; ++i) 
    if (L[i] > eps) ++ret;

  return ret;
}

// �����ž��
Matrix Trans(const Matrix& M)
{
  Matrix R;
  
  for (int i=0; i<3; ++i) 
    for (int j=0; j<3; ++j) R.m[i][j] = M.m[j][i];

  return R;
}

// �ȥ졼��
double Trace(const Matrix& M)
{
  double x = 0.0;
  
  for (int i=0; i<3; ++i)  x += M.m[i][i];

  return x;
}

// �չ���
Matrix Inverse(const Matrix& M)
{
  double det;
  int    i, j;

  double *mm[3], mm0[3], mm1[3], mm2[3];
  mm[0] = mm0; mm[1] = mm1; mm[2] = mm2;

  double *inv[3], inv0[3], inv1[3], inv2[3];
  inv[0] = inv0; inv[1] = inv1; inv[2] = inv2;

  Matrix  R;

  for(i=0; i<3; ++i)
    for(j=0; j<3; ++j)
       mm[i][j] = M.m[i][j];

  det = _alg_matinv(M.dim_i,mm,inv);

  for(i=0; i<3; ++i)
    for(j=0; j<3; ++j)
      R.m[i][j] = inv[i][j];

  return R;
}

// ���̵չ���
Matrix GeneralInverse(const Matrix& M, int times)
{
  Matrix V, R;
  Vector D;
  int    t;

  if (times == 0) t = Rank(M);
  else            t = times;

//  if (fabs(Det(M)) <= eps)
    if (M.Householder(D, V)) {
      for (int i=0; i<t; ++i)
	if (fabs(D[i]) >= eps)
	  R += 1.0 / D[i] * Matrix(V(i), V(i));
      
      return R;
    }
    else {
      fprintf(stderr, "Matrix is not regular.\n");
      return Zm;
    }
//  else
//    return Inverse(M);
}

// ����
double Det(const Matrix& M)
{
  int    i, j;
  double ret;

  double *mm[3], mm0[3], mm1[3], mm2[3];
  mm[0] = mm0; mm[1] = mm1; mm[2] = mm2;

  for(i=0; i<3; ++i)
    for(j=0; j<3; ++j)
       mm[i][j] = M.m[i][j];

  ret = _alg_det(M.dim_i,mm);

  return ret;
}

// ;���ҹ���
Matrix Cofactor(const Matrix& M)
{
  return Matrix(M.m[1][1] * M.m[2][2] - M.m[2][1] * M.m[1][2],
		M.m[2][1] * M.m[0][2] - M.m[0][1] * M.m[2][2],
		M.m[0][1] * M.m[1][2] - M.m[1][1] * M.m[0][2],
		M.m[1][2] * M.m[2][0] - M.m[2][2] * M.m[1][0],
		M.m[2][2] * M.m[0][0] - M.m[0][2] * M.m[2][0],
		M.m[0][2] * M.m[1][0] - M.m[1][2] * M.m[0][0],
		M.m[1][0] * M.m[2][1] - M.m[2][0] * M.m[1][1],
		M.m[2][0] * M.m[0][1] - M.m[0][0] * M.m[2][1],
		M.m[0][0] * M.m[1][1] - M.m[1][0] * M.m[0][1]);
}
  
// ��ͭ�͡���ͭ�٥��ȥ�
int Matrix::Householder(Vector &D, Matrix &V) const
{
  int ret = 0;
  
//  if (Symmetric()) {
    int    i, j;
    double *mm[3], mm0[3], mm1[3], mm2[3];
    double d[3];

    mm[0] = mm0; mm[1] = mm1; mm[2] = mm2;

    for(i=0; i<3; ++i)
      for(j=0; j<3; ++j)
	 mm[i][j] = m[i][j];

    ret = _alg_eigen(3,mm,d);

    D = Vector(d);
    V = Matrix(mm[0][0],mm[1][0],mm[2][0],
               mm[0][1],mm[1][1],mm[2][1],
               mm[0][2],mm[1][2],mm[2][2]);

//  }
  
  if(ret == EXIT_SUCCESS) {
    ret = 1; }
  else {
    ret = 0; }
  return ret;
}

// // �ð���ʬ�� �ˡ�֣̣�t
// int Matrix::Svdecomp(Matrix& V, Matrix& L, Matrix& U) const 
// {
//   double *v = new double[9], *l = new double[9], *u = new double[9],
//          *mm = new double[9];
//   int     ret;
// 
//   for (int i=0; i<3; ++i)
//     for (int j=0; j<3; ++j) mm[i * 3 + j] = m[i][j];
// 
//   ret = MTXsvdecomp3(mm, v, l, u);
// 
//   V = Matrix(v);
//   L = Matrix(l);
//   U = Matrix(u);
// 
//   delete [] v; delete [] l; delete [] u; delete [] mm;
//   
//   return ret;
// }

// ��ž������ž��(��)�Ȳ�ž����ʬ�򤹤�
int Matrix::Rtrans(Vector& l, double& theta) const 
{
  theta = 0.5 * (m[0][0] + m[1][1] + m[2][2] - 1.0);
  if (fabs(theta) >= 1.0) {
    theta = 0.0;
    l = Zv;
  }
  else {
    theta = acos(theta);
    l[0] = 0.5 * (m[2][1] - m[1][2]) / sin(theta);
    l[1] = 0.5 * (m[0][2] - m[2][0]) / sin(theta);
    l[2] = 0.5 * (m[1][0] - m[0][1]) / sin(theta);
  }
  theta *= 180.0 / M_PI;

  return 1;
}



