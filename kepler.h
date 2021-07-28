//---------------------------------------------------------------------------

#ifndef keplerH
#define keplerH
//---------------------------------------------------------------------------
//��������� �������� ��������� ������ � ������ ���������
void ElemToCoord(double Ar1[], double Ar2[]);
//---------------------------------------------------------------------------
//��������� �������� ��������� ������ � ������ ���������
void ElemToCoord(double Ar1[], double r1[], double v1[]);
//---------------------------------------------------------------------------
//��������� �������� ������� ��������� � �������� ������
void CoordToElem(double Ar1[], double Ar2[]);
//---------------------------------------------------------------------------
//��������� �������� ������� ��������� � �������� ������
void CoordToElem(double r1[], double v1[], double Ar2[]);
//---------------------------------------------------------------------------
//��������� ���������� u (��������� ������ � ��������) � ����������� �� ������� � �������� (t=0 - ���������)
double PoiskU(const double &t, const double &a, const double &e, const double &omega);
//---------------------------------------------------------------------------
//��������� ���������� �������(� ��� �� ����������) �� �������� �������� � ���������� ������
double poisk_vrema(const double &teta, const double &a, const double &e);
//---------------------------------------------------------------------------
//��������� ����������� ��������� �������� ��� ������������� ���������� ��������
void Kepler (double r[], double v[], double t);
//---------------------------------------------------------------------------
//��������� ����������� ��������� �������� ��� ������������� ���������� ��������
void Kepler (double Elem[6], double t, double r[], double v[]);
#endif
