//---------------------------------------------------------------------------

#ifndef coordinate_systemH
#define coordinate_systemH
//---------------------------------------------------------------------------
//�������������� ����������� ��������� � ���������
void philamTOxyz (double flr[3], double r[3]);
//�������������� ����������� ��������� � ���������
void philamTOxyz (const double &phi,const double &lam,const double &r,double &x,double &y,double &z);
// �������������� ���������� ��������� � �����������
void xyzTOphilam (double r[3], double flr[3]);
void xyzTOphilam (double r[3], double flr[3], double DS_DX[9]);
// �������������� ���������� ��������� � �����������
void xyzTOphilam (const double &x,const double &y,const double &z,double &phi,double &lam,double &r);
//������� ���������� ��������� � ���������������
void DecToEll(const double Dec[3], double Ell[3]);
//������� ��������������� ��������� � ���������
void EllToDec(const double Ell[3], double Dec[3]);
//�������������� ������������ �� ������ ��������������(������� ������ �������)
void GEOtoGEI(const double Ar1[],const double t, double Ar2[]);
//�������������� ������ �������������� � ������������(������� ������ �������)
void GEItoGEO(const double Ar1[],const double t, double Ar2[]);
//��������� ���������� �������� ���������� �������
double GMST(const double JD);
//��������� ���������� ��������� ��������� ������� � ������� �������� � �������
void precessionVSnutation(const double t,double &s, double NP[3][3]);
//�������������� ������ �������������� � �������� ������������
void GEItoGEO(double r[3], double v[3], double t, double r1[3], double v1[3]);
//�������������� �������� ������������ �� ������ ��������������
void GEOtoGEI(double r[3], double v[3], double t, double r1[3], double v1[3]);
//�������������� ��������� � ������������ (�������. ������ ������)
void MAGtoGEO(const double Ar1[], double Ar2[]);
//�������������� ������������ � ��������� (�������. ������ ������)
void GEOtoMAG(const double Ar1[], double Ar2[]);
//�������������� J2000 � ��������� (�������. ������ ������)
void GEItoMAG(const double Ar1[],const double t, double Ar2[]);
//�������������� ��������� � J2000 (�������. ������ ������)
void MAGtoGEI(const double Ar1[],const double t, double Ar2[]);
//��������� �������� �� J2000 � ������������� �� � ������� � ����� L2 � ���� � ������������ �� ������
void J2000toL2(double r[3], double t, double r3[3]);
void L2toJ2000 ( double r0[3],double v[3], double t, double r1[3],double v1[3]);
//��������� �������� �� J2000 � ��������� ��������� (X ��������� �� ������)
void J2000toEkl ( double r[3],double v[3], double t, double r1[3],double v1[3]);  //!!!������� ��������
//��������� �������� �� ��������� ��������� (X ��������� �� ������) � J2000
void EkltoJ2000 ( double r[3],double v[3], double t, double r1[3],double v1[3]); //!!!������� ��������
//��������� �������� �� J2000 � ���������������� ��������� �������
void J2000toSSB(double r_J2000[3], double v_J2000[3], double t, double r_SSB[3], double v_SSB[3]);
//��������� �������� �� ���������������� ��������� ������� � J2000
void SSBtoJ2000(double r_SSB[3], double v_SSB[3], double t, double r_J2000[3], double v_J2000[3]);
//---------------------------------------------------------------------
//��������� �������� �� J2000 � ����������������� ��
void J2000toSG(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]);
//---------------------------------------------------------------------
//��������� �������� �� ����������������� �� � J2000
void SGtoJ2000(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]);
//---------------------------------------------------------------------
//��������� �������� �� J2000 � ������������������ ��
void J2000toSC(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]);
//---------------------------------------------------------------------
//��������� �������� �� ������������������ �� � J2000
void SCtoJ2000(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]);
//---------------------------------------------------------------------
//��������� �������� �� ������������������ �� � ����������������� ��
void SCtoSG(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]);
//---------------------------------------------------------------------
//��������� �������� �� ������������������ �� � ����������������� ��
 void SCtoSG(double R_in[6],double JD, double R_out[6]);
 //---------------------------------------------------------------------
//��������� �������� �� ����������������� �� � ������������������ ��
void SGtoSC(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]);
//---------------------------------------------------------------------
//��������� �������� �� ����������������� �� � ������������������ ��
 void SGtoSC(double R_in[6],double JD, double R_out[6]);

#endif
