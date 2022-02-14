//---------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>

#pragma hdrstop

#include "L2.h"
#include "constants.h"
#include "global.h"

#include "mathematic.h"
#include "time_convert.h"
#include "coordinate_system.h"
//---------------------------------------------------------------------------

#pragma package(smart_init)


void correctionL2::setNU(double r[3], double v[3], double t){
for(int i=0; i<3; i++){
	r0[i]=r[i];
	v0[i]=v[i];

}
t0=t;
}


void correctionL2::corrl2(){
FILE *ff;
ff=fopen("calc.txt", "w");
fclose(ff);
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

double basa_M87; 	//текущее значение базы для М87
double basa_SGR_A;  //текущее значение базы для SGR-A
double aM=165*dtr;  //максимальный угол между осью Х и направлением на Солнце
double am=90*dtr;   //минимальный угол между осью Х и направлением на Солнце
int NT=11,			//Солнце
	NC=3;           //Земля
int K=11;           //тип расчета при прогнозирование движения (запись трассы в пределах ограничений)
chi::Vect temp;
chi::Vect zero;

chi::integration O;
double  r_M87[3],	//орт направления на М87
		r_SGR_A[3]; //орт направления на SGR-A
double 	rs[3],       //вектор Земля-Солнце
		rkas[3],      //вектор КА-Солнце
		a,          //угол между направлением на источник и на Солнце
		a_M87,      //угол между направлением на КА и на М87
		a_SGR_A,    //угол между направлением на КА и на SGR-A
		basa_min=3000,    //минимальная база
		basa_max=50000;   //максимальная база

double r1[3], v1[3], t1, r[3], v[3], t, dv[3]={0}, DV=0.0001;

double 	T_L2[7],		//массив времен существования в окрестности Л2 при поиске коррекции
		T_M87[7],       //массив времен возможности наблюдения М87 при поиске коррекции
		T_SGR_A[7],     //массив времен возможности наблюдения SGR-A при поиске коррекции
		T_L2_max,       //максимальное время существования в окрестности Л2 при поиске коррекции
		T_M87_max,      //максимальное время возможности наблюдения М87 при поиске коррекции
		T_SGR_A_max;    //максимальное время возможности наблюдения SGR-A при поиске коррекции
int q=0, L;


double  interval=300, //интервал расчета
		step=3600;	//шаг расчета

/*получаем орты направления на источники*/
philamTOxyz(DEC, RA, 1, r_M87[0], r_M87[1], r_M87[2]);
philamTOxyz(DEC1, RA1,  1, r_SGR_A[0], r_SGR_A[1], r_SGR_A[2]);

//задаем тип интегрирования
K=11;
O.setTypeCalculation(K);
O.setParametrs();
O.setParametrs(50, step);
//закладываем начальные параметры
O.setNU(r0,v0,t0);
//интегрируем
O.ABM8();
O.getNU(r0,v0,t0);
O.setParametrs(300, step);



zero.t=2451544.5;
for(int i=0; i<3; i++){
	v1[i]=v0[i];
	r1[i]=r0[i];
}
t1=t0;
bool end=false;

while(!end){
	q++;
	for(int l=0; l<7; l++){
		//формируем ну на очередной прострел
		for(int i=0; i<3; i++){
			v[i]=v1[i]+dv[i]+DV*s[l][i];
			r[i]=r1[i];
		}
		t=t1;
		//закладываем начальные параметры
		O.setNU(r,v,t);
		//интегрируем
		O.ABM8();
		//сохраняем трассы в J2 и L2
		RV_J2[l]=O.rv_trace;
		RV_L2[l]=O.rv_trace_L2;
		//отчищаем трассы для источников
		RV_M87[l].clear();
		RV_SGR_A[l].clear();
		for(int j=0; j<O.rv_trace.size(); j++){
		

	
			/*определение направления на Солнце*/
			de405.calculateR(NT, NC, O.rv_trace[j].t, rs);
			/*вычисление вектора направления КА-Солнце*/
			for(int i=0; i<3; i++){
				rkas[i]=rs[i]-O.rv_trace[j].r[i];
			}
			/*вычисление угла между направлениями на Солнце и на М87*/
			a=angle_between_vectors(rkas,r_M87);
			/*проверка возможности наблюдения (ограничение по Солнцу)*/
			if(a>am && a<aM){
				/*вычисление угла между направлением на КА и на М87*/
				a_M87=angle_between_vectors(O.rv_trace[j].r, r_M87);
				/*вычисление базы М87*/
				basa_M87=norm(O.rv_trace[j].r)*sin(a_M87);
				/*считаем время наблюдения М87*/
				if( basa_M87>basa_min && basa_M87<basa_max){
					for(int ii=0; ii<3; ii++){
						temp.r[ii]=O.rv_trace_L2[j].r[ii];
						temp.v[ii]=O.rv_trace_L2[j].v[ii];
					}
					temp.t=O.rv_trace[j].t;
					RV_M87[l].push_back(temp);
				}
			}

			/*вычисление угла между направлениями на Солнце и на SGR-A*/
			a=angle_between_vectors(rkas,r_SGR_A);
			/*проверка возможности наблюдения (ограничение по Солнцу)*/
			if(a>am && a<aM){
				/*вычисление угла между направлением на КА и на SGR-A*/
				a_SGR_A=angle_between_vectors(O.rv_trace[j].r, r_SGR_A);
				/*вычисление базы SGR-A*/
				basa_SGR_A=norm(O.rv_trace[j].r)*sin(a_SGR_A);
				/*считаем время наблюдения SGR-A*/
				if( basa_SGR_A>basa_min && basa_SGR_A<basa_max){
					for(int ii=0; ii<3; ii++){
						temp.r[ii]=O.rv_trace_L2[j].r[ii];
						temp.v[ii]=O.rv_trace_L2[j].v[ii];
					}
					temp.t=O.rv_trace[j].t;
					RV_SGR_A[l].push_back(temp);
				}
			}
		}
	}
	ff=fopen("calc.txt", "a");
	//вычисляем времена
	for(int l=0; l<7; l++){
		T_L2[l]=RV_L2[l][RV_L2[l].size()-1].t-RV_L2[l][0].t;
		if (RV_M87[l].size()==0) RV_M87[l].push_back(zero);
		T_M87[l]=RV_M87[l][RV_M87[l].size()-1].t-RV_M87[l][0].t;
		if (RV_SGR_A[l].size()==0) RV_SGR_A[l].push_back(zero);
		T_SGR_A[l]=RV_SGR_A[l][RV_SGR_A[l].size()-1].t-RV_SGR_A[l][0].t;
	}

	T_L2_max=0;
	T_M87_max=0.0001;
	T_SGR_A_max=0.0001;
	L=-1;
	for(int l=0; l<7; l++){
		fprintf(ff, "%03i %1i %7.3f %7.3f %7.3f ", q, l, T_L2[l], T_M87[l], T_SGR_A[l]);
		if(T_M87[l]+T_SGR_A[l]>=T_M87_max+T_SGR_A_max && T_L2[l]>180){
			if(T_M87[l]>=T_M87_max && T_SGR_A[l]>=T_SGR_A_max){
				L=l;
				T_M87_max =T_M87[l];
				T_SGR_A_max=T_SGR_A[l];
				fprintf(ff, "M87+SGR_A ");
				if(T_L2[l]>T_L2_max) {
					T_L2_max=T_L2[l];
					fprintf(ff, "T_L2 ");
				}
			}
			else{
				if(T_M87[l]-T_M87_max>0 && T_M87[l]-T_M87_max>fabs(T_SGR_A[l]-T_SGR_A_max)&& T_L2[l]>180){
					L=l;
					T_M87_max =T_M87[l];
					T_SGR_A_max=T_SGR_A[l];
					fprintf(ff, "M87 ");
					if(T_L2[l]>T_L2_max) {
						T_L2_max=T_L2[l];
						fprintf(ff, "T_L2 ");
					}
				}
				if(T_SGR_A[l]-T_SGR_A_max>0 && T_SGR_A[l]-T_SGR_A_max>fabs(T_M87[l]-T_M87_max)&& T_L2[l]>180){
					L=l;
					T_M87_max =T_M87[l];
					T_SGR_A_max=T_SGR_A[l];
					fprintf(ff, "SGR_A ");
					if(T_L2[l]>T_L2_max) {
						T_L2_max=T_L2[l];
						fprintf(ff, "T_L2 ");
					}
				}
			}
		}
		fprintf(ff, "\n");
	}

	if(L>=0){
	   for(int i=0; i<3; i++) dv[i]+=DV*s[L][i];
	   fprintf(ff, "------------------------\n%03i %1i %7.3f %7.3f %7.3f \n", q, L, T_L2[L], T_M87[L], T_SGR_A[L]);
	   fprintf(ff, "%13.10f %13.10f %13.10f\n------------------------\n", dv[0], dv[1], dv[2]);
	}
	else{
		DV*=0.7;
	   fprintf(ff, "------------------------\n без улучшения\n------------------------\n");
	}
	fclose(ff);
	if(DV<0.000001 ||  (T_M87_max>2  && T_SGR_A_max>2)) end=true;





	ff=fopen("tr_L2.txt", "w");
	for(int j=0; j<RV_L2[L].size(); j++){
		fprintf(ff, "%s ", JDToStr(RV_L2[L][j].t, 1));
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_L2[L][j].r[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_L2[L][j].v[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_L2[L][j].r[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_L2[L][j].v[ii]);
		fprintf(ff, "\n");
	}
	fclose(ff);



	ff=fopen("tr_M87.txt", "w");
	if(RV_M87[L].size()>0){
		for(int j=0; j<RV_M87[L].size(); j++){
			fprintf(ff, "%s ", JDToStr(RV_M87[L][j].t, 1));
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_M87[L][j].r[ii]);
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_M87[L][j].v[ii]);
			fprintf(ff, "\n");
		}
	}
	fclose(ff);

	ff=fopen("tr_SGR_A.txt", "w");
	if(RV_SGR_A[L].size()>0){
		for(int j=0; j<RV_SGR_A[L].size(); j++){
			fprintf(ff, "%s ", JDToStr(RV_SGR_A[L][j].t, 1));
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_SGR_A[L][j].r[ii]);
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_SGR_A[L][j].v[ii]);
			fprintf(ff, "\n");
		}
		fprintf(ff, "%5.2f", (RV_SGR_A[L][RV_SGR_A[L].size()-1].t-RV_SGR_A[L][0].t)*24);
	}
	fclose(ff);











}


//ff=fopen("tr_L2.txt", "w");
//for(int j=0; j<O.rv_trace_L2.size(); j++){
//	fprintf(ff, "%s ", JDToStr(O.rv_trace_L2[j].t, 1));
//	for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", O.rv_trace_L2[j].r[ii]);
//	for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", O.rv_trace_L2[j].v[ii]);
//	for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", O.rv_trace[j].r[ii]);
//	for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", O.rv_trace[j].v[ii]);
//	fprintf(ff, "\n");
//}
//fclose(ff);
//
//
//
//ff=fopen("tr_M87.txt", "w");
//if(RV_M87.size()>0){
//	for(int j=0; j<RV_M87.size(); j++){
//		fprintf(ff, "%s ", JDToStr(RV_M87[j].t, 1));
//		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_M87[j].r[ii]);
//		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_M87[j].v[ii]);
//		fprintf(ff, "\n");
//	}
//}
//fclose(ff);
//
//ff=fopen("tr_SGR_A.txt", "w");
//if(RV_SGR_A.size()>0){
//	for(int j=0; j<RV_SGR_A.size(); j++){
//		fprintf(ff, "%s ", JDToStr(RV_SGR_A[j].t, 1));
//		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_SGR_A[j].r[ii]);
//		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_SGR_A[j].v[ii]);
//		fprintf(ff, "\n");
//	}
//	fprintf(ff, "%5.2f", (RV_SGR_A[RV_SGR_A.size()-1].t-RV_SGR_A[0].t)*24);
//}
//fclose(ff);

}





/*процедура определение номера итераци достижения максимума времени 1/2/3 - M87/SGR-A/L2*/
void correctionL2::maxt(int type, int &L){
FILE *ff;
chi::Vect zero;
double 	T_L2[7],		//массив времен существования в окрестности L2 при поиске коррекции
		T_M87[7],       //массив времен возможности наблюдения М87 при поиске коррекции
		T_SGR_A[7],     //массив времен возможности наблюдения SGR-A при поиске коррекции
		T_L2_max,       //максимальное время существования в окрестности L2 при поиске коррекции
		T_M87_max,      //максимальное время возможности наблюдения М87 при поиске коррекции
		T_SGR_A_max;    //максимальное время возможности наблюдения SGR-A при поиске коррекции
zero.t=2451544.5;

  ff=fopen("calc.txt", "a");
	//вычисляем времена
	for(int l=0; l<7; l++){
		T_L2[l]=RV_L2[l][RV_L2[l].size()-1].t-RV_L2[l][0].t;
		if (RV_M87[l].size()==0) RV_M87[l].push_back(zero);
		T_M87[l]=RV_M87[l][RV_M87[l].size()-1].t-RV_M87[l][0].t;
		if (RV_SGR_A[l].size()==0) RV_SGR_A[l].push_back(zero);
		T_SGR_A[l]=RV_SGR_A[l][RV_SGR_A[l].size()-1].t-RV_SGR_A[l][0].t;
	}

	T_L2_max=0;
	T_M87_max=0;
	T_SGR_A_max=0;
	L=-1;



switch(type){
	case 1: for(int l=0; l<7; l++){
				fprintf(ff, "%1i %7.3f %7.3f %7.3f ", l, T_L2[l], T_M87[l], T_SGR_A[l]);
				if(T_M87[l]>T_M87_max && T_L2[l]>180){
					L=l;
					T_M87_max =T_M87[l];
					T_SGR_A_max=T_SGR_A[l];
					T_L2_max=T_L2[l];
				}
				fprintf(ff, "\n");
			}
			break;

	case 2: for(int l=0; l<7; l++){
				fprintf(ff, "%1i %7.3f %7.3f %7.3f ",l, T_L2[l], T_M87[l], T_SGR_A[l]);
				if(T_SGR_A[l]>T_SGR_A_max && T_L2[l]>160){
					L=l;
					T_M87_max =T_M87[l];
					T_SGR_A_max=T_SGR_A[l];
					T_L2_max=T_L2[l];
				}
				fprintf(ff, "\n");
			}
			break;

	case 3: for(int l=0; l<7; l++){
				fprintf(ff, "%1i %7.3f %7.3f %7.3f ",l, T_L2[l], T_M87[l], T_SGR_A[l]);
				if(T_L2[l]>T_L2_max){
					L=l;
					T_M87_max =T_M87[l];
					T_SGR_A_max=T_SGR_A[l];
					T_L2_max=T_L2[l];
				}
				fprintf(ff, "\n");
			}
			break;
}

	if(L>=0){
	   fprintf(ff, "------------------------\n%1i %7.3f %7.3f %7.3f \n", L, T_L2[L], T_M87[L], T_SGR_A[L]);
	}
	else{
	   fprintf(ff, "------------------------\n без улучшения\n");
	}
	fclose(ff);
}





/*печать текущего прострела*/
void correctionL2::print(int L){
FILE *ff;


	ff=fopen("tr_L2.txt", "w");
	for(int j=0; j<RV_L2[L].size(); j++){
		fprintf(ff, "%s ", JDToStr(RV_L2[L][j].t, 1));
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_L2[L][j].r[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_L2[L][j].v[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_J2[L][j].r[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_J2[L][j].v[ii]);
		fprintf(ff, "\n");
	}
	fclose(ff);



	ff=fopen("tr_M87.txt", "w");
	if(RV_M87[L].size()>0){
		for(int j=0; j<RV_M87[L].size(); j++){
			fprintf(ff, "%s ", JDToStr(RV_M87[L][j].t, 1));
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_M87[L][j].r[ii]);
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_M87[L][j].v[ii]);
			fprintf(ff, "\n");
		}
	}
	fclose(ff);

	ff=fopen("tr_SGR_A.txt", "w");
	if(RV_SGR_A[L].size()>0){
		for(int j=0; j<RV_SGR_A[L].size(); j++){
			fprintf(ff, "%s ", JDToStr(RV_SGR_A[L][j].t, 1));
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_SGR_A[L][j].r[ii]);
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_SGR_A[L][j].v[ii]);
			fprintf(ff, "\n");
		}
		fprintf(ff, "%5.2f", (RV_SGR_A[L][RV_SGR_A[L].size()-1].t-RV_SGR_A[L][0].t)*24);
	}
	fclose(ff);
}
/*очистка файлов для записи результатов*/
void correctionL2::clearfile(){
FILE *ff;
ff=fopen("calc.txt", "w");
fclose(ff);
ff=fopen("tr_full.txt", "w");
fclose(ff);
ff=fopen("tr_M87_full.txt", "w");
fclose(ff);
ff=fopen("tr_SGR_A_full.txt", "w");
fclose(ff);
}

/*печать полного файла расчета*/
void correctionL2::print(){
FILE *ff;


	ff=fopen("tr_full.txt", "a");
	for(int j=0; j<RV_L2[0].size(); j++){
		fprintf(ff, "%s ", JDToStr(RV_L2[0][j].t, 1));
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_L2[0][j].r[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_L2[0][j].v[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_J2[0][j].r[ii]);
		for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_J2[0][j].v[ii]);
		fprintf(ff, "\n");
	}
	fclose(ff);



	ff=fopen("tr_M87_full.txt", "a");
	if(RV_M87[0].size()>0){
		for(int j=0; j<RV_M87[0].size(); j++){
			fprintf(ff, "%s ", JDToStr(RV_M87[0][j].t, 1));
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_M87[0][j].r[ii]);
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_M87[0][j].v[ii]);
			fprintf(ff, "\n");
		}
	}
	fclose(ff);

	ff=fopen("tr_SGR_A_full.txt", "a");
	if(RV_SGR_A[0].size()>0){
		for(int j=0; j<RV_SGR_A[0].size(); j++){
			fprintf(ff, "%s ", JDToStr(RV_SGR_A[0][j].t, 1));
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.8f ", RV_SGR_A[0][j].r[ii]);
			for(int ii=0; ii<3; ii++) fprintf(ff, "%15.12f ", RV_SGR_A[0][j].v[ii]);
			fprintf(ff, "\n");
		}
		fprintf(ff, "%5.2f", (RV_SGR_A[0][RV_SGR_A[0].size()-1].t-RV_SGR_A[0][0].t)*24);
	}
	fclose(ff);
}


/*вычисление трассы наблюдений и запись трас в соответствующие переменные*/
void correctionL2::basa(chi::integration O, int l){

chi::Vect temp;
double RA = 3.276086532;    //M87
double DEC = 0.2162659;

double RA1 = 4.649850637;    //SGR_A
double DEC1 = -0.506282045;

double basa_M87; 	//текущее значение базы для М87
double basa_SGR_A;  //текущее значение базы для SGR-A
double aM=165*dtr;  //максимальный угол между осью Х и направлением на Солнце
double am=90*dtr;   //минимальный угол между осью Х и направлением на Солнце
int NT=11,			//Солнце
	NC=3;           //Земля
double  r_M87[3],	//орт направления на М87
		r_SGR_A[3]; //орт направления на SGR-A
double 	rs[3],       //вектор Земля-Солнце
		rkas[3],      //вектор КА-Солнце
		a,          //угол между направлением на источник и на Солнце
		a_M87,      //угол между направлением на КА и на М87
		a_SGR_A,    //угол между направлением на КА и на SGR-A
		basa_min=3000,    //минимальная база
		basa_max=50000;   //максимальная база

/*получаем орты направления на источники*/
philamTOxyz(DEC, RA, 1, r_M87[0], r_M87[1], r_M87[2]);
philamTOxyz(DEC1, RA1,  1, r_SGR_A[0], r_SGR_A[1], r_SGR_A[2]);



		//сохраняем трассы в J2 и L2
		RV_J2[l]=O.rv_trace;
		RV_L2[l]=O.rv_trace_L2;
		//отчищаем трассы для источников
		RV_M87[l].clear();
		RV_SGR_A[l].clear();
		for(int j=0; j<O.rv_trace.size(); j++){



			/*определение направления на Солнце*/
			de405.calculateR(NT, NC, O.rv_trace[j].t, rs);
			/*вычисление вектора направления КА-Солнце*/
			for(int i=0; i<3; i++){
				rkas[i]=rs[i]-O.rv_trace[j].r[i];
			}
			/*вычисление угла между направлениями на Солнце и на М87*/
			a=angle_between_vectors(rkas,r_M87);
			/*проверка возможности наблюдения (ограничение по Солнцу)*/
			if(a>am && a<aM){
				/*вычисление угла между направлением на КА и на М87*/
				a_M87=angle_between_vectors(O.rv_trace[j].r, r_M87);
				/*вычисление базы М87*/
				basa_M87=norm(O.rv_trace[j].r)*sin(a_M87);
				/*считаем время наблюдения М87*/
				if( basa_M87>basa_min && basa_M87<basa_max){
					for(int ii=0; ii<3; ii++){
						temp.r[ii]=O.rv_trace_L2[j].r[ii];
						temp.v[ii]=O.rv_trace_L2[j].v[ii];
					}
					temp.t=O.rv_trace[j].t;
					RV_M87[l].push_back(temp);
				}
			}

			/*вычисление угла между направлениями на Солнце и на SGR-A*/
			a=angle_between_vectors(rkas,r_SGR_A);
			/*проверка возможности наблюдения (ограничение по Солнцу)*/
			if(a>am && a<aM){
				/*вычисление угла между направлением на КА и на SGR-A*/
				a_SGR_A=angle_between_vectors(O.rv_trace[j].r, r_SGR_A);
				/*вычисление базы SGR-A*/
				basa_SGR_A=norm(O.rv_trace[j].r)*sin(a_SGR_A);
				/*считаем время наблюдения SGR-A*/
				if( basa_SGR_A>basa_min && basa_SGR_A<basa_max){
					for(int ii=0; ii<3; ii++){
						temp.r[ii]=O.rv_trace_L2[j].r[ii];
						temp.v[ii]=O.rv_trace_L2[j].v[ii];
					}
					temp.t=O.rv_trace[j].t;
					RV_SGR_A[l].push_back(temp);
				}
			}
		}
}



void correctionL2::corrl2_1(){
FILE *ff;

	double s[7][3]={0, 0, 0,
					1, 0, 0,
					0, 1, 0,
					0, 0, 1,
				   -1, 0, 0,
					0,-1, 0,
					0, 0,-1};




int K=11;           //тип расчета при прогнозирование движения (запись трассы в пределах ограничений)

int type_max[10]={1, 2, 3, 3, 1, 2, 1, 2, 1, 2};
chi::integration O;


double 	r1[3],
		v1[3],
		t1,
		r[3],
		v[3],
		t,
		dv[3]={0},
		DV=0.01;


int q=0, L;


double  dT=90,        //интервал прогноза
		interval=300, //интервал прострела при поиске
		step=3600;	//шаг расчета





clearfile();

//задаем тип интегрирования
K=11;
O.setTypeCalculation(K);
O.setParametrs();
O.setParametrs(300, step);



for(int i=0; i<3; i++){
	v1[i]=v0[i];
	r1[i]=r0[i];
}
t1=t0;
bool end=false;
int qr=-1;

for(double T=0; T<1000; T+=dT){
qr++;
ff=fopen("calc.txt", "a");
	fprintf(ff, "%03i\n", qr);
	fclose(ff);
DV=0.001;
end=false;
O.setParametrs(300, step);
while(!end){
	q++;
	ff=fopen("calc.txt", "a");
	fprintf(ff, "%03i\n", q);
	fclose(ff);
	for(int l=0; l<7; l++){
		//формируем ну на очередной прострел
		for(int i=0; i<3; i++){
			v[i]=v1[i]+dv[i]+DV*s[l][i];
			r[i]=r1[i];
		}
		t=t1;
		//закладываем начальные параметры
		O.setNU(r,v,t);
		//интегрируем
		O.ABM8();
		//считаем базы и возможные интервалы наблюдений
		basa(O, l);
	}
	maxt(3/*type_max[qr]*/, L);
	ff=fopen("calc.txt", "a");
	if(L>0){
	   for(int i=0; i<3; i++) dv[i]+=DV*s[L][i];
	   fprintf(ff, "%13.10f %13.10f %13.10f\n------------------------\n", dv[0], dv[1], dv[2]);
	   print(L);

	   if(	RV_M87[L][RV_M87[L].size()-1].t-RV_M87[L][0].t>2  &&
			RV_SGR_A[L][RV_SGR_A[L].size()-1].t-RV_SGR_A[L][0].t>2) end=true;

	}
	else{
		DV*=0.7;
		fprintf(ff, "%13.10f\n------------------------\n", DV);
	}
	fclose(ff);
	if(DV<0.0001) end=true;


}
ff=fopen("calc.txt", "a");
fprintf(ff, "-------------------Прогноз на 90 дней------------------------\n");
fprintf(ff, "%13.10f %13.10f %13.10f\n", dv[0], dv[1], dv[2]);
fprintf(ff, "%s %13.10f %13.10f %13.10f %13.10f %13.10f %13.10f\n", JDToStr(t1, 1), r1[0], r1[1], r1[2], v1[0], v1[1], v1[2]);

//формируем ну на очередной прострел
for(int i=0; i<3; i++) v1[i]=v1[i]+dv[i];
fprintf(ff, "%s %13.10f %13.10f %13.10f %13.10f %13.10f %13.10f\n", JDToStr(t1, 1), r1[0], r1[1], r1[2], v1[0], v1[1], v1[2]);
fclose(ff);
O.setParametrs(dT, step);
//закладываем начальные параметры
O.setNU(r1,v1,t1);
//интегрируем
O.ABM8();
basa(O, 0);
print();
O.getNU(r1, v1, t1);

}
}
