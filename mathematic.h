//---------------------------------------------------------------------------

#ifndef mathematicH
#define mathematicH
//---------------------------------------------------------------------------
//������� ���������� ������ ������� ����������� ����� ������������
double norm(double x, double y, double z);
//������� ���������� ������ �������
double norm(double r[5]);
//������� ���������� ������ �������� �������� r1 � r2
double norm(double r1[3], double r2[3]);
//��������� ������������ �������� �1 � �2
double Skalyar(const double x1[3], const double x2[3]);
//���� ����� ���������
double angle_between_vectors(double r1[3], double r2[3]);
//��������� ������������ �������� Vec1 � Vec2
void VectProizv (const double Vec1[3], const double Vec2[3], double proizv[3]);
//������� ������� dV1 ������ ��� ���������� ������ Ort �� ���� �� � ��������
void Povorot_Vektora(const double Ort[], const double dV1[], const double hi,double dV2[]);

#endif
