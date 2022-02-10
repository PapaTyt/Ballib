//---------------------------------------------------------------------------

#ifndef L2H
#define L2H
#include "integration.h"
//---------------------------------------------------------------------------
class correctionL2{
	public:
	void setNU(double r[3], double v[3], double t);
	void corrl2();
	void corrl2_1();
	//private:
	protected:
	double t0;
	double r0[3];
	double v0[3];
	std::vector<chi::Vect> RV_L2[7];        //������ ���������� � L2
	std::vector<chi::Vect> RV_J2[7];        //������ ���������� � J2
	std::vector<chi::Vect> RV_M87[7];		//������ ������� ���������� ���������� ���������� �87
	std::vector<chi::Vect> RV_SGR_A[7];    //������ ������� ���������� ���������� ���������� SGR-A


    /*������� ������ ��� ������ �����������*/
	void clearfile();
    /*������ �������� ���������*/
	void print(int L);
    /*������ ������� ����� �������*/
	void print();

	/*���������� ������ ���������� � ������ ���� � ��������������� ����������*/
	void basa(chi::integration O, int l);
	/*��������� ����������� ������ ������� ���������� ��������� ������� 1/2/3 - M87/SGR-A/L2*/
	void maxt(int type, int &L);


};








#endif
