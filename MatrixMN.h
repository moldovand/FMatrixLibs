#ifndef _MatrixMN_h_
#define _MatrixMN_h_

#include "VectorN.h"
#include "Vector.h"
//#include "Matrix.h"

//
//     Class  MatrixMN (MN-Matrix)
//

class MatrixMN     
{ 
private:
protected:
  double **m;
  int      dim_i, dim_j;

public:
  // constructor
  MatrixMN(int di = 3, int dj = 3);
  MatrixMN(const VectorN&, const VectorN&, const VectorN&);
  MatrixMN(double*, int, int);
  MatrixMN(double**, int, int);
  MatrixMN(const VectorN&, const VectorN&);
//  MatrixMN(const Matrix&, const Vector&, const Vector&, double);

  // copy constructor
  MatrixMN(const MatrixMN&);

  // destructor
  ~MatrixMN();
  
  // assignment operator
  MatrixMN& operator=(const MatrixMN&);
  
  double* operator[](int i);                       // ��i��
  const double* operator[](int i) const;
  VectorN operator()(int j) const;                 // ��j��
  int     DimI() const { return dim_i; }
  int     DimJ() const { return dim_j; }
  
  int     Print(FILE *fp = stdout, char *str = "") const;
  int     LongPrint(FILE *fp = stdout, char *str = "") const;
  int     Scan(FILE *fp = stdin);

  // �оι���?
  int Symmetric() const;
  // �����
  int     Init(int i = 0, int j = 0);
  // ñ�̲�
  int     Unit();
  
  // ���
  friend int operator==(const MatrixMN&, const MatrixMN&);
  friend int operator!=(const MatrixMN&, const MatrixMN&);
  
  MatrixMN& operator+=(const MatrixMN&);
  MatrixMN& operator-=(const MatrixMN&);
  MatrixMN& operator*=(double);
  MatrixMN& operator/=(double);
  
  // ������dim_i + 1, dim_j + 1 ���������䤹
  friend MatrixMN  operator+(const MatrixMN&, double);
  friend MatrixMN  operator+(double, const MatrixMN&);
  
  // ���ѱ黻
  friend VectorN   operator*(const MatrixMN&, const VectorN&);
  friend MatrixMN  operator*(const MatrixMN&, const MatrixMN&);
  friend MatrixMN  operator*(double, const MatrixMN&);
  friend MatrixMN  operator*(const MatrixMN&, double);
  friend MatrixMN  operator%(const VectorN&, const MatrixMN&);
  friend MatrixMN  operator/(const MatrixMN&, double);
  friend MatrixMN  operator+(const MatrixMN&, const MatrixMN&);
  friend MatrixMN  operator-(const MatrixMN&);         
  friend MatrixMN  operator-(const MatrixMN&, const MatrixMN&);

  // �ͱƹ���
  friend MatrixMN  Projection(const VectorN&);
  // ����
  friend double    operator,(const MatrixMN&, const MatrixMN&);      
  // �Υ��
  friend double    Norm(const MatrixMN&);
  // ������
  friend MatrixMN  Normalize(const MatrixMN&);
  // ���
  friend int       Rank(const MatrixMN&);
  // ž��
  friend MatrixMN  Trans(const MatrixMN&);
  // �ȥ졼��
  friend double    Trace(const MatrixMN&);
  // �չ���
  friend MatrixMN  Inverse(const MatrixMN&);
  // ���̵չ���
  friend MatrixMN  GeneralInverse(const MatrixMN&, int times = 0);
  // ����
  friend double    Det(const MatrixMN&);
  // ;���ҹ���
//  friend MatrixMN  Cofactor(const MatrixMN&);
  // ��ͭ�͡���ͭ�٥��ȥ�
  int    Householder(VectorN&, MatrixMN&) const;
  // �ð���ʬ�� �ˡ�֣̣�t
//  int    Svdecomp(MatrixMN&, MatrixMN&, MatrixMN&) const;
};

int EddingtonEpsilon(int, int, int);
int KroneckerDelta(int, int);
     
#endif  // _MatrixMN_h_


