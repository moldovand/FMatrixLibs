#ifndef _Vector_h_
#define _Vector_h_

#include "VectorN.h"

//
//     Class  Vector (3-Dimension Vector)
//

class Vector : public VectorN 
{ 
private:
  
public:
  // constructor
  Vector();                                        
  Vector(double, double, double);
  Vector(double*);

  // copy constructor
  Vector(const Vector&);

  // destructor
  ~Vector();

  // assingment operator
  Vector& operator=(const Vector&);
  
  double&       operator[](int i)       { return v[i]; }
  const double& operator[](int i) const { return v[i]; }

  int     Print(FILE *fp = stdout, char *str = "") const;
  int     LongPrint(FILE *fp = stdout, char *str = "") const;
  int     Scan(FILE *fp = stdin);

  // ���
  friend int operator==(const Vector&, const Vector&);     
  friend int operator!=(const Vector&, const Vector&);

  Vector& operator+=(const Vector&);
  Vector& operator-=(const Vector&);
  Vector& operator*=(const double);
  Vector& operator/=(const double);

  // ������1�������䤹
  friend VectorN operator+(const Vector&, double);
  friend VectorN operator+(double, const Vector&);

  // ���ѱ黻
  friend Vector  operator-(const Vector&);
  friend Vector  operator-(const Vector&, const Vector&);
  friend Vector  operator+(const Vector&, const Vector&);
  friend Vector  operator*(double, const Vector&);
  friend Vector  operator*(const Vector&, double);
  friend Vector  operator%(const Vector&, const Vector&); // ����
  friend Vector  operator/(const Vector&, double);
  // ����
  friend double  operator,(const Vector&, const Vector&);
  // �����飳����
  friend double  ScalarTP(const Vector&, const Vector&, const Vector&);
  // �Υ��
  friend double  Norm(const Vector&);
  // ������
  friend Vector  Normalize(const Vector&);                       
  // Ʃ�������
  friend Vector  Proj(const Vector&);
};

extern const Vector Zv;
extern const Vector Iv;
extern const Vector Jv;
extern const Vector Kv;

#endif  // _Vector_h_

