//---------------------------------------------------------------------------


#pragma hdrstop

#include "L2.h"
#include "global.h"
#include "integration.h"
//---------------------------------------------------------------------------

#pragma package(smart_init)
void corrl2(){
	double s[7][3]={0, 0, 0,
					1, 0, 0,
					0, 1, 0,
					0, 0, 1,
				   -1, 0, 0,
					0,-1, 0,
					0, 0,-1};

double RA = 3.276086532;    //M87
double DEC = 0.2162659;

double RA1 = 4.649850637;    //SGR_A
double DEC1 = -0.506282045;

double basa_M87;
double basa_SGR_A;
double aM=165*dtr;
double am=90*dtr;
int NT=11, NC=3;
chi::integration O;



/*�������� ���� ����������� �� ���������*/
philamTOxyz(DEC, RA, 1, r_M87[0], r_M87[1], r_M87[2]);
philamTOxyz(DEC1, RA1,  1, r_SGR_A[0], r_SGR_A[1], r_SGR_A[2]);

//������ ��� ��������������
K=11;
O.setTypeCalculation(K);
O.setParametrs();
O.setParametrs(300, 3600);
//����������� ��������� ���������
O.SetNU(r0,v0,t0);
//�����������
O.ABM8();
for(int j=0; j<O.rv_trace.size(); j++){

	/*����������� ����������� �� ������*/
	de405.calculateR(NT, NC, O.rv_trace[j].t, rs);
	/*���������� ������� ����������� ��-������*/
	for(int i=0; i<3; i++){
		r1[i]=rs[i]-O.rv_trace[j].r[i];
	}
	/*���������� ���� ����� ������������� �� ������ � �� �87*/
	a=angle_between_vectors(r1,r_M87);
	/*�������� ����������� ���������� (����������� �� ������)*/
	if(a>am && a<aM){
		/*���������� ���� ����� ������������ �� �� � �� �87*/
		a_M87=angle_between_vectors(O.rv_trace[j].r, r_M87);
		/*���������� ���� �87*/
		basa_M87=norm(trace)*sin(a_M87);
        /**/
		if(basa_M87<minM87){
			minM87=basa_M87;
			tm87=O.trase_RV[j].t;
		}
		if( basa_M87>basa_min && basa_M87<basa_max){
			t_M87[0]+=step/3600;
		}
	}

}







}

