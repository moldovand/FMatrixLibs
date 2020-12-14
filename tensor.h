#ifndef _TENSOR_
#define _TENSOR_

//ijkl�ƥ󥽥륯�饹
class Tensor3333
{
  int in, jn, kn, ln;       //���줾������ǿ�
  double**** t;             //����

public:
  //1st constructor
  Tensor3333(int i = 3, int j = 3, int k = 3, int l = 3);

  //copy constructor
  Tensor3333(const Tensor3333&);

  //destructor
  ~Tensor3333();

  //���Ǽ��Ф�
  double&       operator()(int, int, int, int);
  const double& operator()(int, int, int, int) const;


  //�����黻  
  const Tensor3333& operator=(const Tensor3333&);

  //���ѱ黻
  Tensor3333 operator+(const Tensor3333&) const;  //­����
  Tensor3333 operator-(const Tensor3333&) const;  //������
  Tensor3333 operator*(double) const;             //������
  friend Tensor3333 operator*(double, const Tensor3333&);
  //3333�ƥ󥽥�*33����
  Matrix operator*(const Matrix&);
  Tensor3333 operator/(double);  //�任  
  //��洹��
  friend Tensor3333 operator-(const Tensor3333&);
  //�ƥ󥽥���(�����Σ��Ĥι���ϣ��ߣ��˸¤�)
  friend Tensor3333 TensorProduct(const Matrix&, const Matrix&);

  //���줾��μ���
  int DimI() const { return in; }
  int DimJ() const { return jn; }
  int DimK() const { return kn; }
  int DimL() const { return ln; }

  //�����
  int Init();

  //ʣ���ˤ��3333�ƥ󥽥����������
  void Fill();

  //ʣ���ˤ��(33)(33)�ƥ󥽥����������
  void FillSym();

  //(33)(33)�ƥ󥽥�����Ƥθ�ͭ�ͤ����Ƥθ�ͭ��������
  int EigenSym(VectorN&, Matrix*);

  //3333�ƥ󥽥�����Ƥθ�ͭ�ͤ����Ƥθ�ͭ��������
  int Eigen(VectorN&, Matrix*);
                                                 
  //(33)(33)�ƥ󥽥�ΰ��̵դ����
  friend Tensor3333 GeneralInverseSym(const Tensor3333&, int times = 0);

  //3333�ƥ󥽥�ΰ��̵դ����
  friend Tensor3333 GeneralInverse(const Tensor3333&, int times = 0);
  
  //���ǽ���
  void Print();
  void Fprint(FILE*);
};

#endif //_TENSOR_
