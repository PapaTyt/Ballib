//---------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>

#pragma hdrstop


#include "global.h"
#include "matrix.h" 
#include "report.h"
#include "constants.h"
#include "mathematic.h"
#include "integration.h"
#include "time_convert.h"
#include "coordinate_system.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)


/*������ ��*/
void chi::integration::set_NU(double r[3], double v[3], double t){
	/*
	 *���������:
	 * 	r[3]	- ������ ������ �� [��]
	 *	v[3] 	- ������ �������� �� [��/�]
	 *	t	 	- ������ ������� �� ������� �������� ������ ��������� � �������
				  ��������� ���� (UTC)
	 */
for(int i=0; i<3; i++){
	rv_nu.r[i]=r[i];
	rv_nu.v[i]=v[i];
}
rv_nu.t=0;
t_nu=t;
}
void chi::integration::set_NU(double rv[6],double t){
	/*
	 *���������:
	 *	�v[6] 	- ������ ��������� �� [��, ��/�]
	 *	t	 	- ������ ������� �� ������� �������� ������ ��������� � �������
				  ��������� ���� (UTC)
	 */
for(int i=0; i<3; i++){
	rv_nu.r[i]=rv[i];
	rv_nu.v[i]=rv[i+3];
}
t_nu=t;
}
/*������ ���������� ��������������*/
void chi::integration::setParametrs(double interval_, double step_){
	/*
	 *���������:
	 * 	interval_ 	- �������� �������� � ����
	 *	step_   	- ��� ������ ���������� � ��������
	 ����������:
		- ���� �� ��������� ����� ������������� ����������� �� ������� ��������
		  ��� ������ ���������
	 */
interval=interval_;
step=step_;
}
/*������ ���������� ��������������*/
void chi::integration::setParametrs(){
rp[0]=true;
rp[1]=true;
rp[2]=true;
rp[3]=false;
rp[4]=false;
centralBody=B_EARTH;
harmonicType=H32;
harmonicOrder=8;

planet[0]=true;
planet[1]=true;
planet[2]=false;
planet[3]=true;
planet[4]=true;
planet[5]=true;
planet[6]=true;
planet[7]=true;
planet[8]=false;
planet[9]=true;
planet[10]=true;
calculeteMatrix=0;




Lx=500000;
Ly-1000000;
Lz=400000;

}
























	/* [���������� ������ ������] */

/*���������� ������ ������ ��*/
void chi::integration::rightPart(VECTOR &rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 

	 *  ���������� �������� ����������
	 *
	 *   calculeteMatrix:
	 *		0 - ������ ������ ����������� ���������,
	 *		1 - ������ ������ ����������� ��������� � 
				������� ������� �����������
	 *
	 *   centralBody:
	 *		B_EARTH - 0 � �����;
	 *		B_MOON  - 1 � ����;
	 *		B_SUN   - 2 � ������;
	 */


/*��������� �������� ���������*/
for(int i=0; i<3; i++){
	a_central_field[i]=0;
	a_off_central_field[i]=0;
	a_celestial_bodies[i]=0;
	a_solar_radiation[i]=0;
	a_atmosphere[i]=0;
	a_traction[i]=0;
}


/*���������� ���������*/
switch(centralBody){
	case B_EARTH:	if(rp[0]) central_field(rv);      	//����������� ����
					if(rp[1]) off_central_field(rv);    //��������������� 
					if(rp[2]) celestial_bodies(rv);     //�������� ����
					if(rp[3]) solar_radiation(rv);      //��������� ��������� 
					if(rp[4]) atmosphere(rv);           //���������
					if(rp[5]) traction(rv);             //���� ��
					break;


	case B_MOON: 	if(rp[0]) central_field_moon(rv);       //����������� ����
					if(rp[1]) off_central_field_moon(rv);   //��������������� 
					if(rp[2]) celestial_bodies_moon(rv);    //�������� ����
					break;


	case B_SUN:     break;

}

/*������������ ���������� �������� ����������� ���������*/
for (int i=0; i<3; i++){
	rv.f[i] = a_central_field[i]+
			  a_off_central_field[i]+
			  a_celestial_bodies[i]+
			  a_solar_radiation[i]+
			  a_atmosphere[i]+
			  a_traction[i];
}

/*������ ������� ������� �����������*/
if(calculeteMatrix){

	/*��������� �������*/
	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
			rv.dfdx[i][j]=0;
			
	/*������� ��������� ���������*/
	rv.dfdx[0][3]=1;
	rv.dfdx[1][4]=1;
	rv.dfdx[2][5]=1;

	/*������� ��������� ������� �� ������� ���������� �����������*/
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)  
			rv.dfdx[i+3][j]	= df_central_field[i][j]+
							  df_off_central_field[i][j]+
							  df_celestial_bodies[i][j];

	matr_X_matr(rv.dfdx, rv.F, rv.F_);
}
}

/*���������� ������������ ���������, �������������� �����������
  ������������� ����� ����� � ������� ������� ����������� ����� �������*/
void chi::integration::central_field(VectSost rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 

	 *  ���������� �������� ����������
	 *
	 *   calculeteMatrix:
	 *		0 - ������ ������ ����������� ���������,
	 *		1 - ������ ������ ����������� ��������� � 
				������� ������� �����������
	 *

	 */

/*���������� ������ ������-�������, ������ � ������ ������� */
double R=norm(rv.r);
double R2=R*R;
double R3=R2*R;

/*���������� ������������ ���������, �������������� �����������
  ������������� ����� �����*/
for (int k=0; k<3; k++)
	a_central_field[k] =  - mu*rv.r[k]/R3;

/*���������� ������� ������� ����������� ������� �������������� ���������� 
  ������������ ���� �� ������� ���������*/
if(calculeteMatrix){  
	for(int i=0;i<3;i++){
		for(int j=0; j<3;j++)
			df_central_field[i][j]=rv.r[i]*rv.r[j]/R2;
		df_central_field[i][i]-=1./3.;
	}

	for(int i=0; i<3;i++){
		for(int j=0; j<3; j++)
			df_central_field[i][j]*=3*mu/R3;
	}
}
}

/*���������� ������������ ���������, �������������� �����������
  ������������� ����� ���� � ������� ������� ����������� ����� �������*/
void chi::integration::central_field_moon(VectSost rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 

	 *  ���������� �������� ����������
	 *
	 *   calculeteMatrix:
	 *		0 - ������ ������ ����������� ���������,
	 *		1 - ������ ������ ����������� ��������� � 
				������� ������� �����������
	 *

	 */

/*���������� ������ ������-�������, ������ � ������ ������� */
double R=norm(rv.r);
double R2=R*R;
double R3=R2*R;

/*���������� ������������ ���������, �������������� �����������
  ������������� ����� �����*/
for (int k=0; k<3; k++)
	a_central_field[k] =  - mum*rv.r[k]/R3;

/*���������� ������� ������� ����������� ������� �������������� ���������� 
  ������������ ���� �� ������� ���������*/
if(calculeteMatrix){  
	for(int i=0;i<3;i++){
		for(int j=0; j<3;j++)
			df_central_field[i][j]=rv.r[i]*rv.r[j]/R2;
	}
	for(int i=0; i<3; i++)  df_central_field[i][i]-=1./3.;
	for(int i=0; i<3;i++){
		for(int j=0; j<3; j++)
			df_central_field[i][j]*=3*mum/R3;
	}
}
}

/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� ����� � ������� ������� ����������� ����� �������*/
void chi::integration::off_central_field(VectSost rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 

	 *  ���������� �������� ����������
	 *
	 *   harmonicType:
	 *		C20 - 0 - ���������� �� �������� ������ ��������� ���������
	 *		C40 - 1 - ���������� �� �������� ��������� ��������� ���������
	 *		H32	- 2 - ���������� ���������� �� ������� 32�32    
	 *
	 *   calculeteMatrix:
	 *		0 - ������ ������ ����������� ���������,
	 *		1 - ������ ������ ����������� ��������� � 
				������� ������� �����������
	 *

	 */

double df[3]={0};

switch(harmonicType){
	case C20: off_central_field_C20(rv);
			  break;
	case C40: off_central_field_C40(rv);
			  break;
	case H32: off_central_field_32(rv, a_off_central_field);
			  /*���������� ������� ������� ����������� ������� �������������� 
				���������� ���������������� ��������������� ���� ���� �� 
				������� ���������*/
			  if(calculeteMatrix){
				  double delta=0.001;
				  VECTOR temp;  
				  for(int i=0; i<3; i++) {
					  temp=rv;
					  temp.r[i]+=delta;
					  off_central_field_32(rv, df);
					  for(int j=0; j<3; j++) df_off_central_field[i][j]=df[j];
					  temp=rv;
					  temp.r[i]-=delta;
					  off_central_field_32(rv, df);
					  for(int j=0; j<3; j++) df_off_central_field[i][j]-=df[j];
					  for(int j=0; j<3; j++) df_off_central_field[i][j]/=2*delta;
				  }

			  }
			  break;

}
}

/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� ���� � ������� ������� ����������� ����� �������*/
void chi::integration::off_central_field_moon(VectSost rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 

	 *  ���������� �������� ����������
	 *
	 *   calculeteMatrix:
	 *		0 - ������ ������ ����������� ���������,
	 *		1 - ������ ������ ����������� ��������� � 
	 *			������� ������� �����������
	 *

	 */

double df[3]={0};

off_central_field_75_moon(rv, df);
/*���������� ������� ������� ����������� ������� �������������� ���������� 
  ���������������� ��������������� ���� ���� �� ������� ���������*/
if(calculeteMatrix){
	double delta=0.001;
	VECTOR temp;  
	for(int i=0; i<3; i++) {
		  temp=rv;
		  temp.r[i]+=delta;
		  off_central_field_75_moon(rv, df);
		  for(int j=0; j<3; j++) df_off_central_field[i][j]=df[j];
		  temp=rv;
		  temp.r[i]-=delta;
		  off_central_field_32(rv, df);
		  for(int j=0; j<3; j++) df_off_central_field[i][j]-=df[j];
		  for(int j=0; j<3; j++) df_off_central_field[i][j]/=2*delta;
	}
}
}


/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� �����(������ ��������� ���������� - �20)*/
void chi::integration::off_central_field_C20(VectSost rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 
	 *  ���������� �������� ����������
	 *
	 *   calculeteMatrix:
	 *		0 - ������ ������ ����������� ���������,
	 *		1 - ������ ������ ����������� ��������� � 
	 *			������� ������� �����������
	 *
	 * 	������������ ��������� (CONSTBNO.h):
	 *		double R_Earth = 6378.136;       ���. ������ ����� (��)
	 *		double b2      = 17599253992.788798202843556774502;
	 */

/*���������� ������ ������-������� ��� �������� � ����� �������*/
double R=norm(rv.r);
double R2=R*R;
double R5=R2*R2*R;
/*������������ �����������*/
double b2=17599253992.788798202843556774502;

/*���������� ��������� ������������ ������ ��������� ���������*/
a_off_central_field[0]=3*b2*rv.r[0]*(5*rv.r[2]*rv.r[2]/R2-1)/2/R5;
a_off_central_field[1]=3*b2*rv.r[1]*(5*rv.r[2]*rv.r[2]/R2-1)/2/R5;
a_off_central_field[2]=3*b2*rv.r[2]*(5*rv.r[2]*rv.r[2]/R2-3)/2/R5;

/*���������� ������� ������� ����������� ������� �������������� ����������
  ���������������� ��������������� ���� ���� �� ������� ���������*/
if(calculeteMatrix){
	double e[3];
	for(int i=0; i<3; i++) e[i]=rv.r[i]/R;
	double H=15*J20*mu*R_Earth*R_Earth/2/R5;
	double ez2=e[2]*e[2];
	df_off_central_field[0][0]=H*(e[0]*e[0]*(7*ez2-1)-(ez2-0.2));
	df_off_central_field[0][1]=H*(e[0]*e[1]*(7*ez2-1));
	df_off_central_field[0][2]=H*(e[0]*e[2]*(7*ez2-3));

	df_off_central_field[1][0]=H*(e[0]*e[1]*(7*ez2-1));
	df_off_central_field[1][1]=H*(e[1]*e[1]*(7*ez2-1)-(ez2-0.2));
	df_off_central_field[1][2]=H*(e[1]*e[2]*(7*ez2-3));

	df_off_central_field[2][0]=H*(e[0]*e[2]*(7*ez2-3));
	df_off_central_field[2][1]=H*(e[1]*e[2]*(7*ez2-3));
	df_off_central_field[2][2]=H*(ez2*(7*ez2-6)+3./5.);
}



}
	
/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� �����(������ ��������� ���������� - �40)*/
 void chi::integration::off_central_field_C40(VectSost rv){
 /*

	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 
	 *
	 * 	������������ ��������� (CONSTBNO.h):
	 *		double R_Earth = 6378.136;       ���. ������ ����� (��)
	 *		double b0 = 3.9860044e5; ��������������� ��������� ��������������� ���� ����� � ������ ���������
	 * 		double b2 =-1.755513e10; ���������� ��� ����� ������ ��������� ���������
	 *		double b4 = 1.564027e15; ���������� ��� ����� ���������� ��������� ���������
	 *
*/

static double b0; /*��������������� ��������� ��������������� ���� ����� � ������ ���������*/
static double b2; /*���������� ��� ����� ������ ��������� ���������*/
static double b4; /*���������� ��� ����� ���������� ��������� ���������*/

static double R;  /*������ ������-������� �� �� ������� �������� ���������� �������������					*/
static double R3; /*������ ������� ������ ������-������� �� �� ������� �������� ���������� �������������  */
static double R5; /*����� ������� ������ ������-������� �� �� ������� �������� ���������� �������������   */
static double R7; /*������� ������� ������ ������-������� �� �� ������� �������� ���������� �������������   */

static double S;  /*������������ ������� �� ��� z*/
static double S2; /*������ ������� ������������� �������� �� ��� z*/
static double S4; /*��������� ������� ������������� �������� �� ��� z*/

static double q;  /*������������*/
static double dq; /*�������� ������������� �� ��� z*/

/*���������� �������� ��������������� ��������*/
b0 = 3.9860044e5;
b2 =-1.755513e10;
b4 = 1.564027e15;

/*���������� ������ ������-������� �� � ��� 3, 5, 7 ��������*/
R = norm(rv.r);
if(R<R_Earth) report(100);
R3 = R*R*R;
R5 = R3*R*R;
R7 = R5*R*R;

/*���������� ������������� �������� �� ��� z � ��� 2 � 4 ��������*/
S  = rv.r[2]/R;
S2 = S*S;
S4 = S2*S2;

/*���������� �������������*/
dq = 3*b2/R5 + 5*b4*(7*S2 - 3)/(2*R7);
q  = -b0/R3 - 3*b2*(5*S2 - 1)/(2*R5) - 15*b4*(21*S4 - 14*S2 + 1)/(8*R7);

/*���������� ���������*/
a_off_central_field[0] =  rv.r[0]*q;
a_off_central_field[1] =  rv.r[1]*q;
a_off_central_field[2] =  rv.r[2]*(q + dq);

 }


/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� ����� � ������ �������� �� 32�32 � ������� �������
  ����������� ����� �������*/
void chi::integration::off_central_field_32(VectSost rv, double df[3]){
/*

	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *   df - ��������� �������������� ���������������� ���������������
	 *  ����������:
	 *    � 
	 *
	 *  ���������� �������� ����������
	 *
	 *   harmonicOrder: �������� [2,32] ������
	 *		
	 *
	 * 	������������ ��������� (CONSTBNO.h):
	 *		double R_Moon       ���. ������ ���� (��)
	 *
*/

int nm;
double U[34][34]={0};
double V[34][34]={0};
double dU[33][33][3]={0};
double dV[33][33][3]={0};
double a_G_gr[3];
double r[3]={0};
double R=norm(rv.r);
double R2=R*R;
double t=t_nu+rv.t/86400.;
double RZi;


/*���������� ������� ��������*/
if(harmonicOrder<=32 && harmonicOrder>=2) nm=harmonicOrder+1;
else nm=33;

// ������� � ��������
GEItoGEO(rv.r,t,r);

/*������������ ������� ��������*/
U[0][0] = 1./R;
V[0][0] = 0.;

/*������������ ��������*/
for (int i=0; i<nm; i++){
	U[i+1][i+1] = (2*i + 1)*(r[0]*U[i][i] - r[1]*V[i][i])/R2;
	V[i+1][i+1] = (2*i + 1)*(r[0]*V[i][i] + r[1]*U[i][i])/R2;
}

/*��������� ��������*/
for (int i=0; i<nm; i++){
	for (int j=0; j<nm; j++){
		if (j>i){
			U[i][j] = 0.;
			V[i][j] = 0.;
		}
		else{
			if (i==0){
				U[i+1][j] = (2*i + 1)*r[2]*U[i][j]/(i - j + 1)/R2;
				V[i+1][j] = (2*i + 1)*r[2]*V[i][j]/(i - j + 1)/R2;
			}
			else{
				U[i+1][j] = ((2*i + 1)*r[2]*U[i][j] - (i + j)*U[i-1][j])/(i - j + 1)/R2;
				V[i+1][j] = ((2*i + 1)*r[2]*V[i][j] - (i + j)*V[i-1][j])/(i - j + 1)/R2;
			}
		}
	}
}
/*���������� ����������� ����������� �������*/
for (int i=0; i<nm; i++){
	/*������� ������� � i-�� ������ ����������� ����������� ������� �� �*/
	dU[i][0][0] = -U[i+1][1];
	dV[i][0][0] = 0.;
	/*������� ������� � i-�� ������ ����������� ����������� ������� �� Y*/
	dU[i][0][1] = -V[i+1][1];
	dV[i][0][1] = 0.;
	/*������� ������� � i-�� ������ ����������� ����������� ������� �� Z*/
	dU[i][0][2] = -(i + 1.)*U[i+1][0];
	dV[i][0][2] = -(i + 1.)*V[i+1][0];
	for (int j=1; j<nm; j++){
		/*��������� �������� ����������� ����������� ������� �� �*/
		dU[i][j][0] = -0.5*U[i+1][j+1] + 0.5*(i - j + 2)*(i - j + 1)*U[i+1][j-1];
		dV[i][j][0] = -0.5*V[i+1][j+1] + 0.5*(i - j + 2)*(i - j + 1)*V[i+1][j-1];
		/*��������� �������� ����������� ����������� ������� �� �*/
		dU[i][j][1] = -0.5*V[i+1][j+1] - 0.5*(i - j + 2)*(i - j + 1)*V[i+1][j-1];
		dV[i][j][1] =  0.5*U[i+1][j+1] + 0.5*(i - j + 2)*(i - j + 1)*U[i+1][j-1];
		/*��������� �������� ����������� ����������� ������� �� Z*/
		dU[i][j][2] = -(i - j + 1.)*U[i+1][j];
		dV[i][j][2] = -(i - j + 1.)*V[i+1][j];
	}
}

/*���������� ���������*/
for (int k=0; k<3; k++) {
	/*���������� ������� ����� � ������� �������*/
	RZi=1;
	/*��������� ������� �������� ������� ���������*/
	a_G_gr[k]=0;
	for (int i=0; i<nm; i++){
		for (int j=0; j<nm; j++){
			/*���������� ������� �������� ������� ���������*/
			a_G_gr[k] += mu*RZi*(Cnn[i][j]*dU[i][j][k] + Snn[i][j]*dV[i][j][k]);
		}
		/*���������� ������� ����� � ��������� ������� �������*/
		RZi*=R_Earth;
	}

}
//������� ������� � 2���
GEOtoGEI(a_G_gr,t,df);
}

/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� ���� � ������ �������� �� 75�75 � ������� �������
  ����������� ����� �������*/
void chi::integration::off_central_field_75_moon(VectSost rv, double df[3]){
/*

	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *   df - ��������� �������������� ���������������� ���������������
	 *  ����������:
	 *    � 
	 *
	 *  ���������� �������� ����������
	 *
	 *   harmonicOrder: �������� [2,75] 
	 *		
	 *
	 * 	������������ ��������� (CONSTBNO.h):
	 *		double R_Moon       ���. ������ ���� (��)
	 *
*/
int nm;
double U[34][34]={0};
double V[34][34]={0};
double dU[33][33][3]={0};
double dV[33][33][3]={0};
double r[3]={0};
double v[3]={0};
double df_G_gr[3]={0};
double R=norm(rv.r);        /*���������� ������ ������-�������*/
double R2=R*R;              /*���������� �������� ������ ������-������� */
double t=t_nu+rv.t/86400.;  /*������������ ������� � ������� UTC*/
double RZi;
double zero[3]={0};


/*���������� ������� ��������*/
if(harmonicOrder<=32 && harmonicOrder>=2) nm=harmonicOrder+1;
else nm=33;

/*������� ������� ��������� � ����������������� ��*/
SCtoSG(rv.r, rv.v, t,r, v);

/*������������ ������� ��������*/
U[0][0] = 1./R;
V[0][0] = 0.;

/*������������ ��������*/
for (int i=0; i<nm; i++){
	U[i+1][i+1] = (2*i + 1)*(r[0]*U[i][i] - r[1]*V[i][i])/R2;
	V[i+1][i+1] = (2*i + 1)*(r[0]*V[i][i] + r[1]*U[i][i])/R2;
}

/*��������� ��������*/
for (int i=0; i<nm; i++){
	for (int j=0; j<nm; j++){
		if (j>i){
			U[i][j] = 0.;
			V[i][j] = 0.;
		}
		else{
			if (i==0){
				U[i+1][j] = (2*i + 1)*r[2]*U[i][j]/(i - j + 1)/R2;
				V[i+1][j] = (2*i + 1)*r[2]*V[i][j]/(i - j + 1)/R2;
			}
			else{
				U[i+1][j] = ((2*i + 1)*r[2]*U[i][j] - (i + j)*U[i-1][j])/(i - j + 1)/R2;
				V[i+1][j] = ((2*i + 1)*r[2]*V[i][j] - (i + j)*V[i-1][j])/(i - j + 1)/R2;
			}
		}
	}
}
/*���������� ����������� ����������� �������*/
for (int i=0; i<nm; i++){
	/*������� ������� � i-�� ������ ����������� ����������� ������� �� �*/
	dU[i][0][0] = -U[i+1][1];
	dV[i][0][0] = 0.;
	/*������� ������� � i-�� ������ ����������� ����������� ������� �� Y*/
	dU[i][0][1] = -V[i+1][1];
	dV[i][0][1] = 0.;
	/*������� ������� � i-�� ������ ����������� ����������� ������� �� Z*/
	dU[i][0][2] = -(i + 1.)*U[i+1][0];
	dV[i][0][2] = -(i + 1.)*V[i+1][0];
	for (int j=1; j<nm; j++){
		/*��������� �������� ����������� ����������� ������� �� �*/
		dU[i][j][0] = -0.5*U[i+1][j+1] + 0.5*(i - j + 2)*(i - j + 1)*U[i+1][j-1];
		dV[i][j][0] = -0.5*V[i+1][j+1] + 0.5*(i - j + 2)*(i - j + 1)*V[i+1][j-1];
		/*��������� �������� ����������� ����������� ������� �� �*/
		dU[i][j][1] = -0.5*V[i+1][j+1] - 0.5*(i - j + 2)*(i - j + 1)*V[i+1][j-1];
		dV[i][j][1] =  0.5*U[i+1][j+1] + 0.5*(i - j + 2)*(i - j + 1)*U[i+1][j-1];
		/*��������� �������� ����������� ����������� ������� �� Z*/
		dU[i][j][2] = -(i - j + 1.)*U[i+1][j];
		dV[i][j][2] = -(i - j + 1.)*V[i+1][j];
	}
}

/*���������� ���������*/
for (int k=0; k<3; k++) {
	/*���������� ������� ����� � ������� �������*/
	RZi=1;
	/*��������� ������� �������� ������� ���������*/
	df_G_gr[k]=0;
	for (int i=0; i<nm; i++){
		for (int j=0; j<nm; j++){
			/*���������� ������� �������� ������� ���������*/
			df_G_gr[k] += mu*RZi*(Cnn[i][j]*dU[i][j][k] + Snn[i][j]*dV[i][j][k]);
		}
		/*���������� ������� ����� � ��������� ������� �������*/
		RZi*=R_Moon;
	}

}
//������� ������� � ������������������ ��
SGtoSC(df_G_gr, zero, t, df, zero);
}

/*���������� ��������� ������������� ��������� �������� ���
  (���������� ��� 10 �������� ���, ����������� ���� �����) � ������� �������
  ����������� ����� �������*/
void chi::integration::celestial_bodies(VectSost rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 
	 *  ���������� �������� ����������
	 *   planet[j] -  ������� ������ ������
	 *   calculeteMatrix:
	 *		0 - ������ ������ ����������� ���������,
	 *		1 - ������ ������ ����������� ��������� � 
	 *			������� ������� �����������
	 *

	 */

unsigned centerBody=dph::B_EARTH; 	/*���������� ����������� ���� ����� - 3 */
unsigned targetBody;              	/*������� �������� ����*/
double t=t_nu+rv.t/86400.;  		/*������������ ������� � ������� UTC*/
double r_planet[3];
double r_ka_planet[3];
double r[3];
double R_planet;
double R_planet3;
double R_ka_planet;
double R_ka_planet2;
double R_ka_planet3;
double a_planet[11][3]={0};
double df_planet[11][3][3]={0};


for(int j=0; j<11; j++){
	if(planet[j]){
		/*����������� �������� ��������� ����*/
		targetBody=j+1;
		/*��������� ������ ������������ ��������� ����*/ 
		if(targetBody==centerBody) {j++; targetBody++;}  
		/*����������� ��������� �������� ��������� ����*/	
		de405.calculateR(targetBody, centerBody, t, r_planet);
		/*����������� ���������� �� �������� ��������� ���� �� ������������*/
		R_planet=norm(r_planet); 
		/*����������� ���������� �� �������� ��������� ���� �� ������������*/
		R_planet3=R_planet*R_planet*R_planet;  
		/*����������� ������� � �������� �� ������� �������� ����*/
		for(int i=0; i<3; i++) 
			r_ka_planet[i]=r_planet[i]-rv.r[i];
		/*����������� ���������� �� �������� �� �������� ��������� ����	*/
		R_ka_planet=norm(r_ka_planet); 
		/*����������� �������� ���������� �� �������� �� �������� ��������� ����*/
		R_ka_planet2=R_ka_planet*R_ka_planet;   
		/*����������� ���� ���������� �� �������� �� �������� ��������� ����*/
		R_ka_planet3=R_ka_planet2*R_ka_planet;    
		//���������� ����� �������� ��������� ������������� �������� �������� ���
		for(int i=0; i<3; i++)
			a_planet[j][i]=mu_planet[j]*(r_ka_planet[i]/R_ka_planet3-r_planet[i]/R_planet3);

		/*���������� ������� ������� ����������� ������� �������������� 
		  ���������� �������� ��� �� ������� ���������*/
		if(calculeteMatrix){
			
			/*���������� ������� �������-��*/
			for(int i=0; i<3; i++) r[i]=-r_ka_planet[i];
			/*���������� ������� ���������� �������*/
			for(int i=0;i<3;i++){
				for(int k=0; k<3;k++)
					df_planet[j][i][k]=r[i]*r[k]/R_ka_planet2;
			}
			/*���������� ������� ���������� �������*/
			for(int i=0; i<3; i++)  df_planet[j][i][i]-=1./3.;
			/*����������� ��������� ������� �� �����������*/
			for(int i=0; i<3;i++){
				for(int k=0; k<3; k++)
					df_planet[j][i][k]*=3*mu_planet[j]/R_ka_planet3;
			}
		}
	}
}

/*���������� �������������� ��������� �� ���� �������� ��� ����� ������������*/
for(int i=0; i<3; i++){
	a_celestial_bodies[i]=0;
	for(int j=0; j<11; j++)
		a_celestial_bodies[i]+=a_planet[j][i];
}

/*���������� ������� ������� ����������� ������� �������������� 
  ���������� �������� ��� �� ������� ���������*/
if(calculeteMatrix){
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			df_celestial_bodies[i][j]=0;
			for(int k=0; k<11; k++){
				df_celestial_bodies[i][j]+=df_planet[k][i][j];
			}
		}
	}
}
}

/*���������� ��������� ������������� ��������� �������� ���
  (���������� ��� 10 �������� ���, ����������� ���� ����) � ������� �������
  ����������� ����� �������*/
void chi::integration::celestial_bodies_moon(VectSost rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 
	 *  ���������� �������� ����������
	 *   planet[j] -  ������� ������ ������
	 *   calculeteMatrix:
	 *		0 - ������ ������ ����������� ���������,
	 *		1 - ������ ������ ����������� ��������� � 
	 *			������� ������� �����������
	 *

	 */
unsigned centerBody=dph::B_MOON; 	/*���������� ����������� ���� ����� - 3 */
unsigned targetBody;              	/*������� �������� ����*/
double t=t_nu+rv.t/86400.;  		/*������������ ������� � ������� UTC*/
double r_planet[3];
double r_ka_planet[3];
double r[3];
double R_planet;
double R_planet3;
double R_ka_planet;
double R_ka_planet2;
double R_ka_planet3;
double a_planet[11][3];
double df_planet[11][3][3];


for(int j=0; j<11; j++){
	if(planet[j]){
		/*����������� �������� ��������� ����*/
		targetBody=j+1;
		/*��������� ������ ������������ ��������� ����*/ 
		if(targetBody==centerBody) {j++; targetBody++;}  
		/*����������� ��������� �������� ��������� ����*/	
		de405.calculateR(targetBody, centerBody, t, r_planet);
		/*����������� ���������� �� �������� ��������� ���� �� ������������*/
		R_planet=norm(r_planet); 
		/*����������� ���������� �� �������� ��������� ���� �� ������������*/
		R_planet3=R_planet*R_planet*R_planet;  
		/*����������� ������� � �������� �� ������� �������� ����*/
		for(int i=0; i<3; i++) 
			r_ka_planet[i]=r_planet[i]-rv.r[i];
		/*����������� ���������� �� �������� �� �������� ��������� ����	*/
		R_ka_planet=norm(r_ka_planet); 
		/*����������� �������� ���������� �� �������� �� �������� ��������� ����*/
		R_ka_planet2=R_ka_planet*R_ka_planet;   
		/*����������� ���� ���������� �� �������� �� �������� ��������� ����*/
		R_ka_planet3=R_ka_planet2*R_ka_planet;    
		//���������� ����� �������� ��������� ������������� �������� �������� ���
		for(int i=0; i<3; i++)
			a_planet[j][i]=mu_planet[j]*(r_ka_planet[i]/R_ka_planet3-r_planet[i]/R_planet3);

		/*���������� ������� ������� ����������� ������� �������������� 
		  ���������� �������� ��� �� ������� ���������*/
		if(calculeteMatrix){
			
			/*���������� ������� �������-��*/
			for(int i=0; i<3; i++) r[i]=-r_ka_planet[i];
			/*���������� ������� ���������� �������*/
			for(int i=0;i<3;i++){
				for(int k=0; k<3;k++)
					df_planet[j][i][k]=r[i]*r[k]/R_ka_planet2;
			}
			/*���������� ������� ���������� �������*/
			for(int i=0; i<3; i++)  df_planet[j][i][i]-=1./3.;
			/*����������� ��������� ������� �� �����������*/
			for(int i=0; i<3;i++){
				for(int k=0; k<3; k++)
					df_planet[j][i][k]*=3*mu_planet[j]/R_ka_planet3;
			}
		}
	}
}

/*���������� �������������� ��������� �� ���� �������� ��� ����� ������������*/
for(int i=0; i<3; i++){
	a_celestial_bodies[i]=0;
	for(int j=0; j<11; j++)
		a_celestial_bodies[i]+=a_planet[j][i];
}

/*���������� ������� ������� ����������� ������� �������������� 
  ���������� �������� ��� �� ������� ���������*/
if(calculeteMatrix){
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			df_celestial_bodies[i][j]=0;
			for(int k=0; k<11; k++){
				df_celestial_bodies[i][j]+=df_planet[k][i][j];
			}
		}
	}
}
}

	/*���������� ��������� �������������� ��������� ���������� ��������� �
	  ������� ������� ����������� ����� �������*/
	void chi::integration::solar_radiation(VectSost rv){
	}

	/*���������� ��������� �������������� ������������ ���� �������������
	  ��������� � ������������ � ���� � 25645.166-2004*/
	void chi::integration::atmosphere(VectSost rv){
	}

	/*���������� ��������� �������������� ������������ ���� �������������
	  ��������� � ������������ � ���� � 25645.166-2004*/
	void chi::integration::atmosphereGOST2004(VectSost rv){
	}

	/*���������� ��������� �������������� ������� ������������ ���������*/
	void chi::integration::traction(VectSost rv){
	}




void chi::integration::RK4_ODE1()
{

double dt=1;      /*��� �������������� �� ����� ����� �������	*/
bool end=false; //���� ��������� �������
double tk=interval*86400.; //�������� ��������������

/*������������ ��������������*/ 
double 	k1[3]={0},
		k2[3]={0},
		k3[3]={0},
		k4[3]={0};
/*������������ ��������������*/		
double  v1[3]={0},
		v2[3]={0},
		v3[3]={0},
		v4[3]={0};
VectSost temp;
/*��������������� ������� �� ��������� ��*/
VectSost rv=rv_nu;


while(rv.t+dt<=tk){
/*���� ��������������*/

	/*���������� k1*/
	rightPart(rv);
	for(int i=0; i<3; i++) k1[i]=rv.f[i]*dt;
	for(int i=0; i<3; i++) v1[i]=rv.v[i]*dt;


	/*���������� k2*/
	temp.t=rv.t+dt/2.;
	for (int i=0; i<3; i++) {
		temp.r[i] = rv.r[i]+ v1[i]/2.;
		temp.v[i]= rv.v[i]+k1[i]/2.;
	}
	rightPart(temp);
	for (int i=0;i<3;i++) k2[i]=temp.f[i]*dt;
	for (int i=0;i<3;i++) v2[i]=temp.v[i]*dt;


	/*���������� k3*/
	temp.t=rv.t+dt/2.;
	for (int i=0; i<3; i++) {
		temp.r[i] = rv.r[i]+ v2[i]/2.;
		temp.v[i]= rv.v[i]+k2[i]/2.;
	}
	rightPart(temp);
	for (int i=0;i<3;i++) k3[i]=temp.f[i]*dt;
	for (int i=0;i<3;i++) v3[i]=temp.v[i]*dt;

	/*���������� k4*/
	temp.t=rv.t+dt;
	for (int i=0; i<3; i++) {
		temp.r[i] = rv.r[i]+ v3[i];
		temp.v[i]= rv.v[i]+k3[i];
	}
	rightPart(temp);
	for (int i=0;i<3;i++) k4[i]=temp.f[i]*dt;
	for (int i=0;i<3;i++) v4[i]=temp.v[i]*dt;
	/*���������� ��������� �����*/
	for (int i=0;i<3;i++){
		temp.r[i] =  rv.r[i]+((v1[i]+2.*v2[i]+2.*v3[i]+v4[i]))/6.;
		temp.v[i] = rv.v[i]+((k1[i]+2.*k2[i]+2.*k3[i]+k4[i]))/6.;
	}
	temp.t=rv.t+dt;
	rv=temp; 


	//���� ��� ����. ���� ������ �� ��������, �� ����. ��� ������ ����� �� tk
	if ((rv.t+dt>tk) && !end){
		dt=tk-rv.t;
		end = true; //���� ��������� ��������������
	}
}



}

/*����� �����-����� 4-�� ������� ��� ������������ ���������������� ���������
  2-�� �������*/
void chi::integration::RK4_ODE2()
{

double dt=1;      /*��� �������������� �� ����� ����� �������	*/
bool end=false; //���� ��������� �������
double tk=interval*86400.;  /*���������� ��������� ��������������*/
VectSost temp;
VectSost rv;

/*������������ ��������������*/
double 	k1[3]={0},
		k2[3]={0},
		k3[3]={0},
		k4[3]={0};


/*��������������� ������� �� ��������� ��*/
rv=rv_nu;


//call_function();

while(rv.t+dt<tk)
{
/*���� ��������������*/

	/*���������� k1*/
	rightPart(rv);
	for(int i=0; i<3; i++) k1[i]=rv.f[i]*dt;

	/*���������� k2*/
	temp.t=rv.t+dt/2.;
	for (int i=0; i<3; i++) {
		temp.r[i] = rv.r[i]+ rv.v[i]*dt/2.;
		temp.v[i]= rv.v[i]+k1[i]/2.;
	}
	rightPart(temp);
	for (int i=0;i<3;i++) k2[i]=temp.f[i]*dt;

	/*���������� k3*/
	temp.t=rv.t+dt/2.;
	for (int i=0; i<3; i++) {
		temp.r[i] = rv.r[i]+ rv.v[i]*dt/2.+k1[i]*dt/4.;
		temp.v[i]= rv.v[i]+k2[i]/2.;
	}
	rightPart(temp);
	for (int i=0;i<3;i++) k3[i]=temp.f[i]*dt;

	/*���������� k4*/
	temp.t=rv.t+dt;
	for (int i=0; i<3; i++) {
		temp.r[i] = rv.r[i]+ rv.v[i]*dt+k2[i]*dt/2.;
		temp.v[i]= rv.v[i]+k3[i];
	}
	rightPart(temp);
	for (int i=0;i<3;i++) k4[i]=temp.f[i]*dt;

	/*���������� ��������� �����*/
	for (int i=0;i<3;i++){
		temp.r[i] =  rv.r[i]+rv.v[i]*dt+((k1[i]+k2[i]+k3[i])*dt)/6.;
		temp.v[i] = rv.v[i]+((k1[i]+2.*k2[i]+2.*k3[i]+k4[i]))/6.;
	}
	temp.t=rv.t+dt;
	rv=temp;

	//���� ��� ����. ���� ������ �� ��������, �� ����. ��� ������ ����� �� tk
	if ((rv.t+dt>tk) && !end){
		dt=tk-rv.t;
		end = true; //���� ��������� ��������������
	}
	//call_function();

 }



}

/*����� �������-������ 7-�� �������*/
void chi::integration::DP7(){

bool end=false;
bool end_step=false;
double t=0;
double tk=interval*86400.;
double t_step=step;
double DT=0;
double sum;
double TE, TEmax, TOL=1e-8, TOLmax=1e-10;
VectSost rv[7];
VectSost temp;
VectSost rv_next;
VectSost rv_next_cap;
rv[0]=rv_nu;
rv[0].t=0;

while(t+dt<=tk){


	if(t+dt>=t_step || (t+dt>=tk)) { DT=dt; dt=t_step-t;end_step=true;}
	//������� �� 12 ���� f[13][6] ��� x y z Vx Vy Vz
	rightPart(rv[0]);

	for(int i=1; i<=6; i++){
		for(int j=0; j<3; j++){
			sum=0;
			for(int m=0; m<=i-1; m++)
				sum+=(dp5_beta[i][m]*rv[m].v[j]);
			rv[i].r[j] = rv[0].r[j]+dt*sum;

			sum=0;
			for(int m=0; m<=i-1; m++)
				sum+=(dp5_beta[i][m]*rv[m].f[j]);
			rv[i].v[j] = rv[0].v[j]+dt*sum;


		}
		rv[i].t=rv[0].t+dp5_alpha[i]*dt;
		rightPart(rv[i]);
	}

	//������� x y z Vx Vy Vz �� ���� ���� � �� �� � ������

	for(int i=0; i<3; i++){
		sum = 0;
		for (int j=0; j<=6; j++)
			sum += (dp5_ce[j]*rv[j].v[i]);
		rv_next.r[i]=rv[0].r[i] + dt*sum;

		sum = 0;
		for (int j=0; j<=6; j++)
			sum += (dp5_ce[j]*rv[j].f[i]);
		rv_next.v[i]=rv[0].v[i] + dt*sum;
	}

	for(int i=0; i<3; i++){
		sum = 0;
		for (int j=0; j<=6; j++)
			sum += (dp5_ce_cap[j]*rv[j].v[i]);
		rv_next_cap.r[i]=rv[0].r[i] + dt*sum;

		sum = 0;
		for (int j=0; j<=6; j++)
			sum += (dp5_ce_cap[j]*rv[j].f[i]);
		rv_next_cap.v[i]=rv[0].v[i] + dt*sum;
	}
	//���������� ����������� ���� - ������� ���� ����������� �� ����� ����
	TE = norm(rv_next.r, rv_next_cap.r);
	//���� ��� ����������� ������� �� ��������, �� ������ ������ ����, ����� ������ ������
	if (TE<TOLmax && !end && !end_step){
		dt *= 1.4;
	}
	else {
		if (TE>TOL && !end && !end_step){
			dt *= 0.7;
		}
		else{
			t+=dt;
			rv_next.t=t;
			rv[0]=rv_next;
			if(end_step) {
				t_step+=step;
				dt=DT;
				end_step=false;
			}
			//������ �� �������
		}
	}
	//���� ��� ����. ���� ������ �� ��������, �� ����. ��� ������ ����� �� tk
	if ((t+dt>tk) && !end){
		dt=tk-t;
		end = true; //���� ��������� ��������������
	}

}
}

/*��������� ���������� �������������� ������� ������-��������-�������
		8-�� ������� � ���������� ����� */
void chi::integration::ABM8 (){

double t0=t_nu; //��������� ������ ������� � ��������� ������ � UTC
double t=0; 	//������ ������ ������� ���� ������ ����������
double max;     //������������ ������� ����� �������������� � �������������� ���������
bool end=false; //���� ��������� �������
double e1=10e-8, e2=10e-14; //������� �������� � �������� ������� �������� ���
double INTERVAL=interval*86400.; //�������� ��������������
int INC=0, INC_DEC=0;    //�������� ���������� ���������� � ���������� ����
VectSost tmp;
//------------------
//�������� ���� ������ �� ������ ��������������
//CleanFlag();
RV[0]=rv_nu;
RV[0].t=0;

/*��������� ���������� ����������� ���� ��������������.
�� ��������� ����� 1� � ���������� �������� ��� �������� ��� � �����������
�� ��������. ����� ���������� ����������� ���� � 1� �������������� ��� ���
���������� ���� ������������ �������� ����� ����� � ��� ���������� ����������
�������� ����� ����� �����-������.*/
dt=1.;
//�������� �������� ���������
if(interval==0)return;
//�������� ����������� �������������� ������(+)/�����(-)
if(interval<0) {
	dt=-dt;
	step=-step;
}

//��������� �������� ��������������
if(fabs(INTERVAL)<=fabs(10*dt)) dt=interval/20.;
//��������� ��������� ������
Iteracii_ABM8();


while(!end){
	//��������� ��������� �����
	ProgKor();
	max=-1;
	//��������� �������� ����������� ��������� �����
	for(int i=0; i<3; i++)
		if(fabs(rv_prog.r[i]-rv_corr.r[i])>max)
			max=fabs(rv_prog.r[i]-rv_corr.r[i]);

	if(max>e1){
		//��������� ��� ��������������
		Decrease_dt_ABM8 ();
		INC_DEC--;
		INC++;
	}
	else{
		if(max<e2 && abs(INC-INC_DEC)<5){
			//����������� ��� ��������������
			Increase_dt_ABM8 ();
			INC++;
			INC_DEC++;
		}
		else{
			//��������� ��� ���� ������ ������� ������������ � �������� 0 � 2 �������
			while(fabs(RV[0].t)<=fabs(t) && fabs(t)<=fabs(RV[1].t) && !end){
				//������������� �� �������� ������ �������
				rv_time=RV[0];
				dt_step=t-RV[0].t;
				Extrapolation(rv_time);
				//��������� ����������� �������� �� ����
				if(!stepCalculation()) end=true;
				//��������� �� ����� ���������
				if(t==INTERVAL) end=true;
				//��������� �� �������� ���������� ����
				else if(fabs(t+step)>fabs(INTERVAL)) step=INTERVAL-t;
				t+=step;
			}
			//������ ����� �������� � �������� �������
			for(int i=0; i<7; i++) RV[i]=RV[i+1];
			RV[7]=rv_prog;

			//���������� ��������� ����������/���������� ����
			INC=0;
			INC_DEC=0;
			//�������� ������ �� ��������
			if (fabs(INTERVAL)<fabs(RV[0].t))
				end =true;
		}
	}
}
//�������� ������ ��� �� ��������� �����
//CloseFiles();
//��������� ����� ����� ��������� �������
//if(RP[4]) Mass=Mass-rv_time.t*P/Pyd;
//��������������� ��������� �������


}
/*��������� ����������� ���������� ������� ����������������� �������
		��� �������������� ������� ������-��������-������� 8-�� ������� */
void chi::integration::Iteracii_ABM8(){
VectSost temp;
dt_step=dt;
//��������� �������
rightPart(RV[0]);
//���������� �������� �� ��������� ����
	for (int k=0; k<7; k++){
		temp=RV[k];
		Extrapolation(temp);
		rightPart(temp);
		RV[k+1]=temp;
	}

}  

/*��������� ����������� �������� � ��������� �� ������
		������-��������-������� 8-�� �������*/
void chi::integration::ProgKor (){
//---------------------------------------------------------------------------


	double 	Fr_AB[3]={0},
			Fv_AB[3]={0},
			Fr_AM[3]={0},
			Fv_AM[3]={0},
			Fdf_AB[6][6]={0},
			Fdf_AM[6][6]={0};
	double T;
	rv_prog.t=RV[7].t+dt;
	rv_corr.t=rv_prog.t;
	T=t_nu+rv_prog.t/86400.;
	//�������
	for (int i=0; i<3; i++){
		Fv_AB[i]=0;
		Fr_AB[i]=0;
		for (int j=0; j<8; j++){
			Fv_AB[i]+=abm8_AB[j]*RV[j].f[i];
			Fr_AB[i]+=abm8_AB[j]*RV[j].v[i];
		}
		rv_prog.v[i]=RV[7].v[i]+dt*Fv_AB[i];
		rv_prog.r[i]=RV[7].r[i]+dt*Fr_AB[i];
	}
	for(int i=0; i<6; i++)
		for(int ii=0; ii<6; ii++){
			Fdf_AB[i][ii]=0;
			for (int j=0; j<8; j++) Fdf_AB[i][ii]+=abm8_AB[j]*RV[j].F_[i][ii];
			rv_prog.F[i][ii]=RV[7].F[i][ii]+dt*Fdf_AB[i][ii];
		}

	rightPart(rv_prog);
	//���������
	for (int i=0; i<3; i++){
		Fv_AM[i]=abm8_AM[7]*rv_prog.f[i];
		Fr_AM[i]=abm8_AM[7]*rv_prog.v[i];
		for (int j=1; j<8; j++){
			Fv_AM[i]+=abm8_AM[j-1]*RV[j].f[i];
			Fr_AM[i]+=abm8_AM[j-1]*RV[j].v[i];
		}
		rv_corr.v[i]=RV[7].v[i]+dt*Fv_AM[i];
		rv_corr.r[i]=RV[7].r[i]+dt*Fr_AM[i];
	}

    for(int i=0; i<6; i++)
		for(int ii=0; ii<6; ii++){
			Fdf_AM[i][ii]=abm8_AM[7]*rv_prog.F_[i][ii];
			for (int j=0; j<8; j++) Fdf_AM[i][ii]+=abm8_AM[j-1]*RV[j].F_[i][ii];
			rv_corr.F[i][ii]=RV[7].F[i][ii]+dt*Fdf_AM[i][ii];
		}
	rightPart(rv_corr);

}
/*��������� ��������� ���� �������������� ��� ������
		������-��������-������� 8-�� ������� */
void chi::integration::Decrease_dt_ABM8 (){

dt/=2.;
Iteracii_ABM8();
//<-4.4
}
//---------------------------------------------------------------------------
/*��������� ���������� ���� �������������� ��� ������
		������-��������-������� 8-�� ������� */
void chi::integration::Increase_dt_ABM8 (){

VectSost RV_tmp[4];
for(int i=0; i<4; i++){
	 RV_tmp[i]=RV[2*i];

}

	//����������� ��������� 3 ��������
	for (int k=0; k<7; k++){
		ProgKor();
		//������ ����� �������� � �������� �������
		for(int i=0; i<7; i++) RV[i]=RV[i+1];
		RV[7]=rv_prog;
	}

	RV[7]=RV[7];
	RV[6]=RV[5];
	RV[5]=RV[3];
	RV[4]=RV[1];
	RV[3]=RV_tmp[3];
	RV[2]=RV_tmp[2];
	RV[1]=RV_tmp[1];
	RV[0]=RV_tmp[0];

	//����������� ��� � ��� ����
	dt*=2;
//<-4.5
}


/* 4.7-> ��������� ������������� �����-�������-����*/
void chi::integration::Extrapolation(VECTOR &rv0){
//---------------------------------------------------------------------------
/* 4.7-> ��������� ������������� �����-�������-����
=============================================================================
���������� ���������:
		  rv - ��������� ���� VECTOR
*/


int n[30];
bool end=false, END=false, fl=false;
VectSost rv[3], N[30][30];
double  eps=1e-12, t, t0, H=0.1, NN, VV, n_n,n1,
		h[30]={0};
//������� ���������� �������� �������
t=t0=rv0.t;
//��������� ������ ��� ����� ������������������ ����� �����
n[1]=2;
n[2]=4;
n[3]=6;
//��������� �� ����������� ��������� ��� H=1�. �������� ��������������
///////////////////////////////////
//if(STEP<H) {H=STEP; fl=true;}
if(dt_step==0) return;
if(dt_step<0) H=-H;
if(fabs(dt_step)<=fabs(H)) H=dt_step/10.;
////////////////////////////////
//�������� ���� �� �������� END - �������������� � �������� ��������� STEP
while (!END){
	//�������� �� ��������� ��� ��������������
	if(fl)END=true;
	//���� �� �������� end - ������������ ����������������� ������� ��� �������� ���� �
	for(int k=1; !end; k++){
		//����������� ����� ������ ��� ���������� ������ ����� ��������� �� ������� ��� ������ �� ����� �����������
		if(k>3) n[k]=2*n[k-2];
		//���������� �������� ��������� ��� ������ ������ ����������������� ������� �� ������ ���� �
		h[k]=H/n[k];
		//��������������� ��������� ������
		rv[1]=rv0;
		//��������� ��������� ��� ����������� ������� r[1] �� ������� ������ �������  t
		rightPart(rv[1]);
		//���������� ������ ��������� �� ��������� ���������
		for(int i=0; i<3; i++){
			rv[2].v[i]=rv[1].v[i]+h[k]*rv[1].f[i];
			rv[2].r[i]=rv[1].r[i]+h[k]*rv[1].v[i];
		}
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				rv[2].F[i][j]=rv[1].F[i][j]+h[k]*rv[1].F_[i][j];


		//���������� ������ ��������� �� ����������� ����������
		for(int i=1; i<=n[k]; i++){
			//��������������� ������� ���������
			rv[0]=rv[1];
			rv[1]=rv[2];
			//��������� ��������� ��� ����������� ������� r[1] �� ������� ������ ������� t+i*h[k]
			rv[1].t=t+i*h[k];
			rightPart(rv[1]);
			//���������� ������ ��������� �� ��������� ���� ���������
			for(int j=0; j<3; j++){
				rv[2].v[j]=rv[0].v[j]+2.*h[k]*rv[1].f[j];
				rv[2].r[j]=rv[0].r[j]+2.*h[k]*rv[1].v[j];
			}
			for(int i=0; i<6; i++)
				for(int j=0; j<6; j++)
					rv[2].F[i][j]=rv[1].F[i][j]+h[k]*rv[1].F_[i][j];

		}
		//��������� ������ ���������������� ������� ����������������� ������� � ������� ������
		for(int i=0; i<3; i++){
			 N[k][1].r[i]= rv[0].r[i]/4.+ rv[1].r[i]/2.+rv[2].r[i]/4.;
			 N[k][1].v[i]= rv[0].v[i]/4.+ rv[1].v[i]/2.+rv[2].v[i]/4.;
		}
		for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
			N[k][1].F[i][j]=rv[0].F[i][j]/4.+rv[1].F[i][j]/2.+rv[2].F[i][j]/4.;

		//��������� ����������� ���������������� �������� ����������������� ������� � ������� ������
		for(int j=1; j<k; j++){
			n1=n[k-j];
			n_n=n[k]*n[k]/n1/n1;
			//if(n_n==1) ShowMessage("");
			for(int i=0; i<3; i++){
			  N[k][j+1].r[i]=N[k][j].r[i]+(N[k][j].r[i]-N[k-1][j].r[i])/(n_n-1);
			  N[k][j+1].v[i]=N[k][j].v[i]+(N[k][j].v[i]-N[k-1][j].v[i])/(n_n-1);
			}
			for(int i=0; i<6; i++)
			for(int ii=0; ii<6; ii++)
				N[k][j+1].F[i][ii]=N[k][j].F[i][ii]+(N[k][j].F[i][ii]-N[k-1][j].F[i][ii])/(n_n-1);
		}
		//���������� ����� ���������� ���� ��������� ��������� ����������������� ������� � ������� ������
		NN=norm(N[k][k].r,N[k][k-1].r);
		VV=norm(N[k][k].v,N[k][k-1].v);
		//����������� ������� �� ��������
		if(NN<eps && VV<eps){
				//������������ �������� ��������� ���������� ����� ����������������� �������
				end=true;
				//����������� ������� �� ���
				if(fabs(VV)>fabs(NN)) NN=VV;
				//������� �� ������ ���� ����� ����� ����
				if(NN==0) NN=1e-20;
				//��������������� ��������� ������
				rv[1]=N[k][k];
				rv0=rv[1];
				//��������������� �����
				t=t+H;
				rv0.t=t;
				//��������� ����� �������� ��� ���������� ����
				H=0.94*H*pow((0.65*eps/NN),1./(2.*k-1));
				//��������� ������� �� ������ �� �������� �������������� �� ��������� ����
				if(fabs(t0+dt_step-t)<=fabs(H)) {
					//��������� ��� �� ����� ��������� ���������������
					H=t0+dt_step-t;
					//��������� �������� �� ������� ���������� ���� �� ��������� STEP
					fl=true;
				}
			}
			if(k==29) {end=true; fl=true;}
			}
			//�������������� ������� �� ������������ ����� ����������������� �������
			end=false;
	}
//��������������� �����!
rv0.t=t;

}



























/*������ ���� ���������� �� ������ ���� ��������������*/
void chi::integration::setTypeCalculation(typeCalculation_){
 typeCalculation=typeCalculation_;
}



bool chi::integration::stepCalculation(){

bool call;

switch(typeCalculation){
	case 0: call=printSteate();
			break;













	/*������� ��� ���*/
	case 11: call=max_t();
			 break;















}
return call;
}

bool chi::integration::printSteate(){
FILE *ff;
AnsiString str=JDToStr(t_nu+rv_time.t/86400, 1);

ff=fopen("1.txt", "a");
 
fprintf(ff, "%10.1f %s %15.7f %15.7f %15.7f %15.7f %15.7f %15.7f\n",
		rv_time.t,
		str.c_str(),
		rv_time.r[0], rv_time.r[1], rv_time.r[2],
		rv_time.v[0], rv_time.v[1], rv_time.v[2]);

fclose(ff);
}


bool chi::integration::trajectory(){
Vect temp;

if(rv_time.t==0) rv_trace.clear();

for(int i=0; i<3; i++){
	temp.r[i]=rv_time.r[i];
	temp.v[i]=rv_time.v[i];
}
temp.t=rv_time.t;
rv_trace.push_back(temp);
if(norm(rv_time.r)<2000000) return 1;
else return 0;
}


bool chi::integration::max_t(){

Vect temp;
Vect temp_L2;
double r[3];
double t=t_nu+rv_time.t/86400.;

if(rv_time.t==0) rv_trace.clear();

J2000toL2(rv_time.r,  t, temp_L2.r);

for(int i=0; i<3; i++){
	temp.r[i]=rv_time.r[i];
	temp.v[i]=rv_time.v[i];
}
temp.t=t;
temp_L2.t=t;
rv_trace.push_back(temp);
rv_trace_L2.push_back(temp_L2);

if(fabs(temp_L2.r[0])<Lx && fabs(temp_L2.r[1])<Ly && fabs(temp_L2.r[2])<Lz) {
	return 0;
}
else return 1;



}

