#ifndef _Matrix_h_
#define _Matrix_h_

#include "MatrixMN.h"

//
//     Class  Matrix (3-Dimensional Square Matrix)
//

class Matrix : public MatrixMN
{ 
private:

public:
  // constructor
  Matrix();
  Matrix(const Vector&, const Vector&, const Vector&);
  Matrix(const Vector&, double);
  Matrix(double*);
  Matrix(double**);
  Matrix(double, double, double,
	 double, double, double,
	 double, double, double);
  Matrix(const Vector&, const Vector&);
  Matrix(double, double, double);

  // copy constructor
  Matrix(const Matrix&);

  // destructor
  ~Matrix();

  // assignment operator
  Matrix& operator=(const Matrix&);
  
  double* operator[](int i);                       // ��i��
  const double* operator[](int i) const;
  Vector  operator()(int j) const;                 // ��j��

  int     Print(FILE *fp = stdout, char *str = "") const;
  int     LongPrint(FILE *fp = stdout, char *str = "") const;
  int     Scan(FILE *fp = stdin);

  // �оι���?
  int Symmetric() const;
  // �����
  int     Init();
  // ñ�̲�
  int     Unit();

  // ���
  friend int operator==(const Matrix&, const Matrix&);
  friend int operator!=(const Matrix&, const Matrix&);

  Matrix& operator+=(const Matrix&);
  Matrix& operator-=(const Matrix&);
  Matrix& operator*=(double);
  Matrix& operator/=(double);

  // ������4, 4 ���������䤹
  friend MatrixMN operator+(const Matrix&, double);
  friend MatrixMN operator+(double, const Matrix&);
  
  // ���ѱ黻
  friend Vector  operator*(const Matrix&, const Vector&);
  friend Matrix  operator*(const Matrix&, const Matrix&);
  friend Matrix  operator*(double, const Matrix&);
  friend Matrix  operator*(const Matrix&, double);
  friend Matrix  operator%(const Vector&, const Matrix&);  // ����
  friend Matrix  operator%(const Matrix&, const Vector&);  // ����
  friend Matrix  operator%(const Matrix&, const Matrix&);  // ����
  friend Matrix  operator/(const Matrix&, double);
  friend Matrix  operator+(const Matrix&, const Matrix&);
  friend Matrix  operator-(const Matrix&);         
  friend Matrix  operator-(const Matrix&, const Matrix&);

  // �ͱƹ���
  friend Matrix  Projection(const Vector&);
  // ����
  friend double  operator,(const Matrix&, const Matrix&);      
  // �Υ��
  friend double  Norm(const Matrix&);
  // ������
  friend Matrix  Normalize(const Matrix&);
  // ���
  friend int     Rank(const Matrix&);
  // ž��
  friend Matrix  Trans(const Matrix&);
  // �ȥ졼��
  friend double  Trace(const Matrix&);
  // �չ���
  friend Matrix  Inverse(const Matrix&);
  // ���̵չ���
  friend Matrix  GeneralInverse(const Matrix&, int times = 0);
  // ����
  friend double  Det(const Matrix&);
  // ;���ҹ���
  friend Matrix  Cofactor(const Matrix&);
  // ��ͭ�͡���ͭ�٥��ȥ�
  int    Householder(Vector&, Matrix&) const;
  // �ð���ʬ�� �ˡ�֣̣�t
//  int    Svdecomp(Matrix&, Matrix&, Matrix&) const;
  // ��ž���Ȳ�ž��(��)��ʬ�򤹤�
  int    Rtrans(Vector&, double&) const;
};

extern const Matrix Im;
extern const Matrix Zm;
extern const Matrix Pk;

#endif  // _Matrix_h_
