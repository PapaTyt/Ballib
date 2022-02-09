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

double basa_M87; 	//������� �������� ���� ��� �87
double basa_SGR_A;  //������� �������� ���� ��� SGR-A
double aM=165*dtr;  //������������ ���� ����� ���� � � ������������ �� ������
double am=90*dtr;   //����������� ���� ����� ���� � � ������������ �� ������
int NT=11,			//������
	NC=3;           //�����
int K=11;           //��� ������� ��� ��������������� �������� (������ ������ � �������� �����������)
chi::Vect temp;
chi::Vect zero;
std::vector<chi::Vect> RV_L2[7];        //������ ���������� � L2
std::vector<chi::Vect> RV_J2[7];        //������ ���������� � J2
std::vector<chi::Vect> RV_M87[7];		//������ ������� ���������� ���������� ���������� �87
std::vector<chi::Vect> RV_SGR_A[7];    //������ ������� ���������� ���������� ���������� SGR-A
chi::integration O;
double  r_M87[3],	//��� ����������� �� �87
		r_SGR_A[3]; //��� ����������� �� SGR-A
double 	rs[3],       //������ �����-������
		rkas[3],      //������ ��-������
		a,          //���� ����� ������������ �� �������� � �� ������
		a_M87,      //���� ����� ������������ �� �� � �� �87
		a_SGR_A,    //���� ����� ������������ �� �� � �� SGR-A
		basa_min=3000,    //����������� ����
		basa_max=50000;   //������������ ����

double r1[3], v1[3], t1, r[3], v[3], t, dv[3]={0}, DV=0.0001;

double T_L2[7], T_M87[7], T_SGR_A[7], T_L2_max, T_M87_max, T_SGR_A_max;
int q=0, L;


double  interval=300, //�������� �������
		step=3600;	//��� �������

/*�������� ���� ����������� �� ���������*/
philamTOxyz(DEC, RA, 1, r_M87[0], r_M87[1], r_M87[2]);
philamTOxyz(DEC1, RA1,  1, r_SGR_A[0], r_SGR_A[1], r_SGR_A[2]);

//������ ��� ��������������
K=11;
O.setTypeCalculation(K);
O.setParametrs();
O.setParametrs(50, step);
//����������� ��������� ���������
O.setNU(r,v,t);
//�����������
O.ABM8();
�.getNU(r0,v0,t0);




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
		//��������� �� �� ��������� ��������
		for(int i=0; i<3; i++){
			v[i]=v1[i]+dv[i]+DV*s[l][i];
			r[i]=r1[i];
		}
		t=t1;
		//����������� ��������� ���������
		O.setNU(r,v,t);
		//�����������
		O.ABM8();
		//��������� ������ � J2 � L2
		RV_J2[l]=O.rv_trace;
		RV_L2[l]=O.rv_trace_L2;
		//�������� ������ ��� ����������
		RV_M87[l].clear();
		RV_SGR_A[l].clear();
		for(int j=0; j<O.rv_trace.size(); j++){
		

	
			/*����������� ����������� �� ������*/
			de405.calculateR(NT, NC, O.rv_trace[j].t, rs);
			/*���������� ������� ����������� ��-������*/
			for(int i=0; i<3; i++){
				rkas[i]=rs[i]-O.rv_trace[j].r[i];
			}
			/*���������� ���� ����� ������������� �� ������ � �� �87*/
			a=angle_between_vectors(rkas,r_M87);
			/*�������� ����������� ���������� (����������� �� ������)*/
			if(a>am && a<aM){
				/*���������� ���� ����� ������������ �� �� � �� �87*/
				a_M87=angle_between_vectors(O.rv_trace[j].r, r_M87);
				/*���������� ���� �87*/
				basa_M87=norm(O.rv_trace[j].r)*sin(a_M87);
				/*������� ����� ���������� �87*/
				if( basa_M87>basa_min && basa_M87<basa_max){
					for(int ii=0; ii<3; ii++){
						temp.r[ii]=O.rv_trace_L2[j].r[ii];
						temp.v[ii]=O.rv_trace_L2[j].v[ii];
					}
					temp.t=O.rv_trace[j].t;
					RV_M87[l].push_back(temp);
				}
			}

			/*���������� ���� ����� ������������� �� ������ � �� SGR-A*/
			a=angle_between_vectors(rkas,r_SGR_A);
			/*�������� ����������� ���������� (����������� �� ������)*/
			if(a>am && a<aM){
				/*���������� ���� ����� ������������ �� �� � �� SGR-A*/
				a_SGR_A=angle_between_vectors(O.rv_trace[j].r, r_SGR_A);
				/*���������� ���� SGR-A*/
				basa_SGR_A=norm(O.rv_trace[j].r)*sin(a_SGR_A);
				/*������� ����� ���������� SGR-A*/
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
	//��������� �������
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
	   fprintf(ff, "------------------------\n ��� ���������\n------------------------\n");
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






//void correctionL2::maxt()
