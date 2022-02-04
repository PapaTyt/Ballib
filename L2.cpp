//---------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>

#pragma hdrstop

#include "L2.h"
#include "constants.h"
#include "global.h"
#include "integration.h"
#include "mathematic.h"
#include "time_convert.h"
#include "coordinate_system.h"
//---------------------------------------------------------------------------

#pragma package(smart_init)


void correctoinL2::setNU(double r[3], double v[3], double t){
for(int i=0; i<3; i++){
	r0[i]=r[i];
	v0[i]=v[i];

}
t0=t;
}


void correctoinL2::corrl2(){
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
int K;
chi::Vect temp;
std::vector<chi::Vect> RV_M87;
std::vector<chi::Vect> RV_SGR_A;
chi::integration O;
double  r_M87[3],r_SGR_A[3];
double rs[3];
double r1[3], a, a_M87, a_SGR_A, basa_min=3000, basa_max=50000;

double step=3600;

/*получаем орты направлени€ на источники*/
philamTOxyz(DEC, RA, 1, r_M87[0], r_M87[1], r_M87[2]);
philamTOxyz(DEC1, RA1,  1, r_SGR_A[0], r_SGR_A[1], r_SGR_A[2]);

//задаем тип интегрировани€
K=11;
O.setTypeCalculation(K);
O.setParametrs();
O.setParametrs(300, step);
//закладываем начальные параметры
O.setNU(r0,v0,t0);
//интегрируем
O.ABM8();
for(int j=0; j<O.rv_trace.size(); j++){

	/*определение направлени€ на —олнце*/
	de405.calculateR(NT, NC, O.rv_trace[j].t, rs);
	/*вычисление вектора направлени€  ј-—олнце*/
	for(int i=0; i<3; i++){
		r1[i]=rs[i]-O.rv_trace[j].r[i];
	}
	/*вычисление угла между направлени€ми на —олнце и на ћ87*/
	a=angle_between_vectors(r1,r_M87);
	/*проверка возможности наблюдени€ (ограничение по —олнцу)*/
	if(a>am && a<aM){
		/*вычисление угла между направлением на  ј и на ћ87*/
		a_M87=angle_between_vectors(O.rv_trace[j].r, r_M87);
		/*вычисление базы ћ87*/
		basa_M87=norm(O.rv_trace[j].r)*sin(a_M87);
		/*считаем врем€ наблюдени€ ћ87*/
		if( basa_M87>basa_min && basa_M87<basa_max){
			for(int ii=0; ii<3; ii++){
				temp.r[ii]=O.rv_trace_L2[j].r[ii];
				temp.v[ii]=O.rv_trace_L2[j].v[ii];
			}
			temp.t=O.rv_trace[j].t;
			RV_M87.push_back(temp);
		}
	}

	/*вычисление угла между направлени€ми на —олнце и на SGR-A*/
	a=angle_between_vectors(r1,r_SGR_A);
	/*проверка возможности наблюдени€ (ограничение по —олнцу)*/
	if(a>am && a<aM){
		/*вычисление угла между направлением на  ј и на SGR-A*/
		a_SGR_A=angle_between_vectors(O.rv_trace[j].r, r_SGR_A);
		/*вычисление базы SGR-A*/
		basa_SGR_A=norm(O.rv_trace[j].r)*sin(a_SGR_A);
		/*считаем врем€ наблюдени€ SGR-A*/
		if( basa_SGR_A>basa_min && basa_SGR_A<basa_max){
			for(int ii=0; ii<3; ii++){
				temp.r[ii]=O.rv_trace_L2[j].r[ii];
				temp.v[ii]=O.rv_trace_L2[j].v[ii];
			}
			temp.t=O.rv_trace[j].t;
			RV_SGR_A.push_back(temp);
		}
	}
}





FILE *ff;

ff=fopen("tr_L2.txt", "w");
for(int j=0; j<O.rv_trace_L2.size(); j++){
	fprintf(ff, "%s ", JDToStr(O.rv_trace_L2[j].t, 1));
	for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", O.rv_trace_L2[j].r[ii]);
	for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", O.rv_trace_L2[j].v[ii]);
	for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", O.rv_trace[j].r[ii]);
	for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", O.rv_trace[j].v[ii]);
	fprintf(ff, "\n");
}
fclose(ff);



ff=fopen("tr_M87.txt", "w");
if(RV_M87.size()>0){
	for(int j=0; j<RV_M87.size(); j++){
		fprintf(ff, "%s ", JDToStr(RV_M87[j].t, 1));
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_M87[j].r[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_M87[j].v[ii]);
		fprintf(ff, "\n");
	}
}
fclose(ff);

ff=fopen("tr_SGR_A.txt", "w");
if(RV_SGR_A.size()>0){
	for(int j=0; j<RV_SGR_A.size(); j++){
		fprintf(ff, "%s ", JDToStr(RV_SGR_A[j].t, 1));
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_SGR_A[j].r[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_SGR_A[j].v[ii]);
		fprintf(ff, "\n");
	}
	fprintf(ff, "%5.2f", (RV_SGR_A[RV_SGR_A.size()-1].t-RV_SGR_A[0].t)*24);
}
fclose(ff);









}

