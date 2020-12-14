#ifndef _VectorN_h_
#define _VectorN_h_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "FileScan.h"
#include "defVM.h"

//
//     Class  VectorN (N-Dimension Vector)
//

class VectorN     
{ 
private:
protected:
  double *v;
  int     dim;

public:
  // constructor
  VectorN(int d = 3);                                        
  VectorN(double*, int d);

  // copy constructor
  VectorN(const VectorN&);

  // destructor
  ~VectorN();

  // assingment operator
  VectorN& operator=(const VectorN&);
  
  double&       operator[](int i)       { return v[i]; }
  const double& operator[](int i) const { return v[i]; }
  int     Dim() const { return dim; }
  
  int     Print(FILE *fp = stdout, char *str = "") const;
  int     LongPrint(FILE *fp = stdout, char *str = "") const;
  int     Scan(FILE *fp = stdin);

  // �����
  int     Init(int d = 0);
  
  // ���
  friend int operator==(const VectorN&, const VectorN&);     
  friend int operator!=(const VectorN&, const VectorN&);

  VectorN& operator+=(const VectorN&);
  VectorN& operator-=(const VectorN&);
  VectorN& operator*=(double);
  VectorN& operator/=(double);

  // ������1�������䤹
  friend VectorN  operator+(const VectorN&, double);
  friend VectorN  operator+(double, const VectorN&);
  
  // ���ѱ黻
  friend VectorN  operator-(const VectorN&);
  friend VectorN  operator-(const VectorN&, const VectorN&);
  friend VectorN  operator+(const VectorN&, const VectorN&);
  friend VectorN  operator*(double, const VectorN&);
  friend VectorN  operator*(const VectorN&, double);
  friend VectorN  operator%(const VectorN&, const VectorN&); // ����
  friend VectorN  operator/(const VectorN&, double);
  // ����
  friend double   operator,(const VectorN&, const VectorN&);
  // �Υ��
  friend double   Norm(const VectorN&);
  // ������
  friend VectorN  Normalize(const VectorN&);                       
  // Ʃ�������
  friend VectorN  Proj(const VectorN&);
};

double sqr(double);
int    sgn(double);
double cbrt(double);

#endif  // _VectorN_h_


