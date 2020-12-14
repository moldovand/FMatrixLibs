#ifndef _TRANSFORM_
#define _TRANSFORM_

//3333�ƥ󥽥�->99����
MatrixMN Type99(const Tensor3333&);

//(33)(33)�ƥ󥽥�->66����
MatrixMN Type66(const Tensor3333&);

//6�٥��ȥ�->(33)����
Matrix Type33sym(const VectorN&);

//66����->(33)(33)�ƥ󥽥�
Tensor3333 Type3333sym(const MatrixMN&);

//9�٥��ȥ�->33����
Matrix Type33(const VectorN&);

//33����->9�٥��ȥ�
VectorN Type9(const Matrix&);

//99����->3333�ƥ󥽥�
Tensor3333 Type3333(const MatrixMN&);

#endif //_TRANSFORM_
