//---------------------------------------------------------------------------

#ifndef matrixH
#define matrixH
#include <vcl.h>
#include <vector>
//---------------------------------------------------------------------------
//��������� ��������� ������� A 3�3 �� ������ r
void matr_X_vect(double A[3][3], double r[3], double r1[3]);
//��������� ��������� ������� A 3�3 �� ������� 3�3 B
void matr_X_matr(double A[3][3], double B[3][3], double C[3][3]);
//��������� ���������������� ������� � 3�3
void matr_T(const double A[3][3], double B[3][3]);
//��������� �������������� ������� � ����������
void MatrToQ(const double a[3][3], double q[4]);
//��������� �������������� ����������� � �������
void QToMatr(const double q[4], double a[3][3]);
//����������� ������� �������� ������ ����� �� ���� �� �� �������� ����
void RotMat(int i, double a, double A[3][3]);
 //������� � ������� 6�6
//---------------------------------------------------------------------------
//��������� �������� 6-�� ������� �������
void vect_zero(double X[6]);
//---------------------------------------------------------------------------
//��������� ������ ������� � � ������ Y
void vect_copy(double X[6], double Y[6]);
//---------------------------------------------------------------------------
//��������� �� 6-�� ������� ������� 6-�� ������� �������
void vect_minus(double X[6], double Y[6], double Z[6]);
//---------------------------------------------------------------------------
//��������� �� 6-�� ������� ������� 6-�� ������� �������
void vect_minus(double X[6], double Y[6]);
//�������� 6-�� ������� ������� � 6-�� ������ ��������
void vect_add(double A[6], double B[6], double Z[6]);
//---------------------------------------------------------------------------
//�������� 6-�� ������� ������� � 6-�� ������ ��������
void vect_add(double A[6], double B[6]);
//---------------------------------------------------------------------------
//��������� ������� �� �����
void vect_x_n(double X[6], double n, double Y[6]);
//---------------------------------------------------------------------------
//��������� ������� �� �����
void vect_x_n(double X[6], double n);
//---------------------------------------------------------------------------
//��������� ��������� ������� 6�6
void matr_zero(double A[6][6]);
//---------------------------------------------------------------------------
//��������� ������ ������� � 6�6 � ������� � 6�6
void matr_copy(double A[6][6], double B[6][6]);
//---------------------------------------------------------------------------
//��������� ������� 6�6 �� ������� 6�6
void matr_minus(double A[6][6], double B[6][6], double C[6][6]);
//---------------------------------------------------------------------------
//��������� ������� 6�6 �� ������� 6�6
void matr_minus(double A[6][6], double B[6][6]);
//---------------------------------------------------------------------------
//�������� ������� 6�6 � �������� 6�6
void matr_add(double A[6][6], double B[6][6], double C[6][6]);
//---------------------------------------------------------------------------
//�������� ������� 6�6 � �������� 6�6
void matr_add(double A[6][6], double B[6][6]);
//---------------------------------------------------------------------------
//��������� ������� 6�6 �� �����
void matr_X_n(double A[6][6], double n, double B[6][6]);
//---------------------------------------------------------------------------
//��������� ������� 6�6 �� �����
void matr_X_n(double A[6][6], double n);
//---------------------------------------------------------------------------
//��������� ��������� ������� A 6�6 �� ������ X 6�1
void matr_X_vect(double A[6][6], double X[6], double Y[6]);
//---------------------------------------------------------------------------
//��������� ��������� ������ X 1�6 �� ������� A 6�6
void vect_X_matr(double X[6], double A[6][6], double Y[6]);
//---------------------------------------------------------------------------
//��������� �������� ������� A 6�6 �� ������� B 6�6
void matr_X_matr(double A[6][6], double B[6][6], double C[6][6]);
//---------------------------------------------------------------------------
//��������� �������� ������� A 6�6 �� ������� B 6�6
void matr_X_matr(double A[6][6], double B[6][6]);
//---------------------------------------------------------------------------
//��������� ��������������� ������� A 6�6
void matr_T(double A[6][6], double B[6][6]) ;
//---------------------------------------------------------------------------
//��������� ��������������� ������� A 6�6
void matr_T(double A[6][6]);
//---------------------------------------------------------------------------
//��������� ������� 6x6
void inversion(double A[6][6]);
//---------------------------------------------------------------------------


class matrix{
	public:

	//������ �������
	std::vector<double> X;
	//������ ������
	std::vector<double> Y;

	//���������� ����� n
	int n;
	//���������� �������� m
	int m;
	//������� nxm
	std::vector< std::vector<double> > A;

	//---------------------------------------------------------------------------
	//��������� ������������� ������� nxm
	void matrix::create_A(int n_, int m_);
	//---------------------------------------------------------------------------
	//��������� ������������� ������� �������
	void matrix::create_X(int n_);
	//---------------------------------------------------------------------------
	//��������� ������������� ������� ������
	void matrix::create_Y(int m_);
	//---------------------------------------------------------------------------
	//��������� ��������� ������� nxm
	void matrix::zero_A();
	//---------------------------------------------------------------------------
	//��������� ��������� ������� nxm
	void matrix::zero_X();
	//---------------------------------------------------------------------------
	//��������� ��������� ������� nxm
	void matrix::zero_A();
	//---------------------------------------------------------------------------
	//��������� ������������� ������������ �������
	void diagonal_matrix();
	//---------------------------------------------------------------------------
	//��������� ������������� ��������� �������
	void identity_matrix();
    //��������� �������� ������ A=A+B
	void add_A(matrix B);
	//---------------------------------------------------------------------------
	//��������� �������� �������� �������� X=X+B.X
	void add_X(matrix B);
	//---------------------------------------------------------------------------
	//��������� �������� �������� �������� Y=Y+B.Y
	void add_Y(matrix B);
    //---------------------------------------------------------------------------
	//��������� ��������� ������ A=A-B.A
	void matrix::minus_A(matrix B);
	//---------------------------------------------------------------------------
	//��������� ��������� �������� �������� X=X-B.X
	void matrix::minus_X(matrix B);
	//---------------------------------------------------------------------------
	//��������� ��������� �������� �������� Y=Y-B.Y
	void matrix::minus_Y(matrix B);
    //---------------------------------------------------------------------------
	//��������� ��������� ������� �� ����� A=k*A
	void matrix::multiplication_A(double k);
	//---------------------------------------------------------------------------
	//��������� ��������� �������� �������� X=k*X
	void matrix::multiplication_X(double k);
	//---------------------------------------------------------------------------
	//��������� ��������� �������� �������� Y=k*Y
	void matrix::multiplication_Y(double k);
    //---------------------------------------------------------------------------
	//��������� ���������������� ������� A=A^T
	void matrix::transposition_A();
	//---------------------------------------------------------------------------
	//��������� ��������� �������� �������� Y=X^T
	void matrix::transposition_X();
	//---------------------------------------------------------------------------
	//��������� ��������� �������� �������� X=Y^T
	void matrix::transposition_Y();





	//privete
	protected:


}:




#endif
