//---------------------------------------------------------------------------


#pragma hdrstop

#include "L2.h"
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
chi::integration O;


//задаем тип интегрирования
K=11;
O.setTypeCalculation(K);
O.setParametrs();
O.setParametrs(300, 3600);
//закладываем начальные параметры
O.SetNU(r0,v0,t0);
//интегрируем
O.ABM8();
for(int j=0; j<O.rv_trace.size(); j++){



}







}

