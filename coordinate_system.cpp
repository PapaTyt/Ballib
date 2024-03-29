//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include "coordinate_system.h"
#include "mathematic.h"
#include "global.h"
#include "constants.h"
#include "time_convert.h"
#include "matrix.h"


//---------------------------------------------------------------------------

#pragma package(smart_init)
//----------------------------------------------------------------------------
//�������������� ����������� ��������� � ���������
void philamTOxyz (double flr[3], double r[3])
				 //flr[0] - ������ � ��������
				 //flr[1] - ������� � ��������
				 //flr[2] - ������ � ����������
				 //r[0] - �������� ������� �� ��� �
				 //r[1] - �������� ������� �� ��� y
				 //r[3] - �������� ������� �� ��� z
{
	r[0] = flr[2]*cos(flr[0])*cos(flr[1]);
	r[1] = flr[2]*cos(flr[0])*sin(flr[1]);
	r[2] = flr[2]*sin(flr[0]);
}
//----------------------------------------------------------------------------
//�������������� ����������� ��������� � ���������
void philamTOxyz (	const double &phi,	//������ � ��������
					const double &lam,  //������� � ��������
					const double &r,    //������ � ����������
					double &x,          //�������� ������� �� ��� �
					double &y,          //�������� ������� �� ��� y
					double &z)          //�������� ������� �� ��� z
{
	x = r*cos(phi)*cos(lam);
	y = r*cos(phi)*sin(lam);
	z = r*sin(phi);
}
//-------------------------------------------------------------------------
// �������������� ���������� ��������� � �����������
void xyzTOphilam (double r[3], double flr[3])
				 //flr[0] - ������ � ��������
				 //flr[1] - ������� � ��������
				 //flr[2] - ������ � ����������
				 //r[0] - �������� ������� �� ��� �
				 //r[1] - �������� ������� �� ��� y
				 //r[3] - �������� ������� �� ��� z
{
	flr[0] = atan(r[2]/sqrt(r[0]*r[0] + r[1]*r[1]));
	flr[1] = atan2(r[1], r[0]);
	flr[2] = norm(r);
}
//-------------------------------------------------------------------------
// �������������� ���������� ��������� � �����������
void xyzTOphilam (double r[3], double flr[3], double DS_DX[9])
				 //flr[0] - ������ � ��������
				 //flr[1] - ������� � ��������
				 //flr[2] - ������ � ����������
				 //r[0] - �������� ������� �� ��� �
				 //r[1] - �������� ������� �� ��� y
				 //r[3] - �������� ������� �� ��� z
				 //DS_DX - ������� ������� ����������� ����������� ��������� flr �� ����������� ������� r
{
	flr[0] = atan(r[2]/sqrt(r[0]*r[0] + r[1]*r[1]));
	flr[1] = atan2(r[1], r[0]);
	flr[2] = norm(r);


double R=norm(r);
double xy=sqrt(r[0]*r[0]+r[1]*r[1]);

   //������
   DS_DX[0]=r[0]/R;
   DS_DX[1]=r[1]/R;
   DS_DX[2]=r[2]/R;
   //�������
   DS_DX[3]=-r[1]/(xy*xy);
   DS_DX[4]=r[0]/(xy*xy);
   DS_DX[5]=0.0;
   //������
   DS_DX[6]=-(r[2]*r[0])/(R*R*xy);
   DS_DX[7]=-(r[2]*r[1])/(R*R*xy);
   DS_DX[8]=xy/(R*R);








}
// �������������� ���������� ��������� � �����������
void xyzTOphilam (	const double &x,   //�������� ������� �� ��� �
					const double &y,   //�������� ������� �� ��� y
					const double &z,   //�������� ������� �� ��� z
					double &phi,       //������ � ��������
					double &lam,       //������� � ��������
					double &r)         //������ � ����������
{
	phi = atan(z/sqrt(x*x + y*y));
	lam = atan2(y, x);
	r = norm(x,y,z);
}
//---------------------------------------------------------------------------
/*������� ���������� ��������� � ���������������
������� ���������:
		  Dec - ������ ���������� ��������� (x, y, z) (���������)
�������� ���������:
			Ell - ������ ��������������� ��������� (fi, lambda, h)
				  fi     - ������  (�������)
				  lambda - ������� (�������)
				  h      - ������  (���������)
������������ ������:
		<math.h>
*/
void DecToEll(const double Dec[3], double Ell[3])
{
	double x  = Dec[0]*1000.;
	double y  = Dec[1]*1000.;
	double z  = Dec[2]*1000.;

	double a = 6378137.;
	double f = 1./298.257223563;
	double e = sqrt(2*f-f*f);

	double lambda;
	if( x == 0.)
		{
			if( y > 0.) lambda = M_PI/2.;
			if( y < 0.) lambda = M_PI*3./2.;
		}
	else
		{
			lambda = atan2( y, x );
		}

	double fi = atan(z/sqrt(x*x+y*y)); //��������� �������� ��� �����
	double fi1, N, h;
	fi1=fi+1;
   //	for( ;  fabs(fi-fi1) > 10e-10  ; )
   while(fabs(fi-fi1)>10e-8)
		{
			fi1=fi;
			N = a/sqrt(1. - e*e*sin(fi1)*sin(fi1));
			h = ( sqrt(x*x+y*y) / cos(fi1) ) - N;
			fi = atan(z/sqrt(x*x+y*y) / (1.-e*e* N/(N+h)));
		}

	Ell[0] = fi;
	Ell[1] = lambda;
	Ell[2] = h/1000.;
}
//---------------------------------------------------------------------------
/*������� ��������������� ��������� � ���������
������� ���������:
			Ell - ������ ��������������� ��������� (fi, lambda, h)
				  fi     - ������  (�������)
				  lambda - ������� (�������)
				  h      - ������  (���������)
�������� ���������:
		  Dec - ������ ���������� ��������� (x, y, z) (���������)
������������ ������:
		<math.h>
*/
void EllToDec(const double Ell[3], double Dec[3])
{
	double fi      = Ell[0];
	double lambda  = Ell[1];
	double h       = Ell[2]*1000.;

	double a = 6378137.;
	double f = 1./298.257223563;
	double e = sqrt(2*f-f*f);

	double N       = a/sqrt(1.-e*e*sin(fi)*sin(fi));

	double x       = (N + h) * cos(fi) * cos(lambda);
	double y       = (N + h) * cos(fi) * sin(lambda);
	double z       = ( N*( 1.-e*e ) + h ) * sin(fi) ;

	Dec[0] = x/1000.;
	Dec[1] = y/1000.;
	Dec[2] = z/1000.;
}
//---------------------------------------------------------------------------
/*�������������� ������������ �� ������ ��������������
(������� ������ �������)
������� ���������:
		Ar1[3] - ���������� ������ (x, y, z) � ����������
		t - ����� � ��������� ������
�������� ���������:
		Ar2[3] - ���������� ������ (x, y, z) � ����������
������������ ���������:
		precessionVSnutation - ���� ��������� �������
������������ ������:
		<math.h>
*/
void GEOtoGEI(const double Ar1[],const double t, double Ar2[]){
	 double cs, ss, s,
			S[3][3], NP[3][3],
			PNS[3][3], SNP[3][3];
precessionVSnutation(t, s, NP);
cs = cos(s);
ss = sin(s);
S[0][0]=cs;		S[0][1]=-ss; 	S[0][2]=0;
S[1][0]=ss;		S[1][1]=cs; 	S[1][2]=0;
S[2][0]=0;		S[2][1]=0; 	 	S[2][2]=1;
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++)   {
		PNS[i][c]=0;
		for(int j=0; j<3; j++)
			PNS[i][c]+=NP[j][i]*S[j][c];
	}
for(int i=0; i<3; i++){
	Ar2[i]=0;
	for(int j=0; j<3; j++)
		Ar2[i]+=PNS[i][j]*Ar1[j];
}
}
//-------------------------------------------------------------------
/*�������������� ������ �������������� � ������������
(������� ������ �������)
������� ���������:
		Ar1[3] - ���������� ������ (x, y, z) � ����������
		t - ����� � ��������� ������
�������� ���������:
		Ar2[3] - ���������� ������ (x, y, z) � ����������
������������ ���������:
		precessionVSnutation - ���� ��������� �������
������������ ������:
		<math.h>
*/
void GEItoGEO(const double Ar1[],const double t, double Ar2[]){
	double cs, ss, s,
		   S[3][3], NP[3][3],
		   SNP[3][3];
precessionVSnutation(t, s, NP);
cs = cos(s);
ss = sin(s);
S[0][0]=cs;		S[0][1]=ss; 	S[0][2]=0;
S[1][0]=-ss;	S[1][1]=cs; 	S[1][2]=0;
S[2][0]=0;		S[2][1]=0; 	 	S[2][2]=1;
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++)  {
		SNP[i][c]=0;
		for(int j=0; j<3; j++)
			SNP[i][c]+=S[i][j]*NP[j][c];
	}
for(int i=0; i<3; i++){
	Ar2[i]=0;
	for(int j=0; j<3; j++)
		Ar2[i]+=SNP[i][j]*Ar1[j];
}
}
//---------------------------------------------------------------------------
/*��������� ���������� �������� ���������� �������
������� ���������:
		JD - ����� � ��������� ������
�������� ���������:
		�������� ����� � ��������
������������ ���������:
		modf
		pow
������������ ������:
		<math.h>
		"constants.h"
*/
double GMST(const double JD)
{
   double st,JD0,Tut,n;
   modf(JD - 0.5,&n);
   JD0 = n + 0.5;
   Tut = (JD0 - 2451545.)/36525.;
   st = 1.7533685592 + 628.3319706889*Tut + 6.7707139*pow(10,-6)*Tut*Tut
   - 4.50876*pow(10,-10)*pow(Tut,3) + 1.002737909350795*2*M_PI*(JD - JD0);
   modf(st/(2*M_PI),&n);
   return st-n*2*M_PI; // � ���������� �������� �� 0 �� 2*PI
}

//----------------------------------------------------------------------------
/*��������� ���������� ��������� ��������� ������� � ������� �������� � �������
������� ���������:
		t - ����� � ��������� ������
�������� ���������:
		s - �������� ����� � ��������
		NP - ������� 3�3 ��������� ������������ ������� ������� �� ������� ���������
������������ ���������:
		UTCtoTDB - ������� UTC � ������������ �����
		GMST - ��������� ���������� �������� ��������� �������
������������ ������:
		<math.h>
		#include "time.h"
*/
void precessionVSnutation(const double t,double &s, double NP[3][3]){
double ksi, zeta, teta, psi, eps, eps0;
double 	c_ksi, s_ksi,
		c_zeta, s_zeta,
		c_teta, s_teta,
		c_psi, s_psi,
		c_eps, s_eps,
		c_eps0, s_eps0;
double T, T2, T3;
double Ml,Ms,ul,Ds,Uzl;
double P[3][3], N[3][3];
T   = (UTCtoTDB(t) - 2451545.0)/36525.0;
T2=T*T;
T3=T2*T;
ksi  = 0.011180860*T + 1.464E-6*T2 + 8.70E-8*T3;
zeta = 0.011180860*T + 5.308E-6*T2 + 8.90E-8*T3;
teta = 0.009717173*T - 2.068E-6*T2 - 2.02E-7*T3;
eps0 = 0.4090928042 - 0.2269655E-3*T - 2.86E-9*T2 + 8.80E-9*T3;

Ml   = 2.355548394 + (1325*pi_2 + 3.470890873)*T + 1.517952E-4*T2 + 3.103E-7*T3;
Ms   = 6.240035940 +   (99*pi_2 + 6.266610600)*T - 2.797400E-6*T2 - 5.820E-8*T3;
ul   = 1.627901934 + (1342*pi_2 + 1.431476084)*T - 6.427170E-5*T2 + 5.340E-8*T3;
Ds   = 5.198469514 + (1236*pi_2 + 5.360106500)*T - 3.340860E-5*T2 + 9.220E-8*T3;
Uzl  = 2.182438624 -    (5*pi_2 + 2.341119397)*T + 3.614290E-5*T2 + 3.880E-8*T3;

eps = 0; psi = 0;
for(int i=0;i<=105;i++){
	psi +=   (coef_A[i] + coef_B[i]*T)*sin(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
	eps +=   (coef_C[i] + coef_D[i]*T)*cos(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
}
psi  = psi*1E-4/3600./180.*pi;
eps  = eps*1E-4/3600./180.*pi;
eps += eps0;

c_ksi=cos(ksi);
s_ksi=sin(ksi);
c_zeta=cos(zeta);
s_zeta=sin(zeta);
c_teta=cos(teta);
s_teta=sin(teta);
c_psi=cos(psi);
s_psi=sin(psi);
c_eps=cos(eps);
s_eps=sin(eps);
c_eps0=cos(eps0);
s_eps0=sin(eps0);

s=GMST(t)+psi*c_eps;

P[0][0] =  c_ksi*c_zeta*c_teta - s_ksi*s_zeta;
P[0][1] = -s_ksi*c_zeta*c_teta - c_ksi*s_zeta;
P[0][2] = -c_zeta*s_teta;

P[1][0] =  c_ksi*s_zeta*c_teta + s_ksi*c_zeta;
P[1][1] = -s_ksi*s_zeta*c_teta + c_ksi*c_zeta;
P[1][2] = -s_zeta*s_teta;

P[2][0] =  c_ksi*s_teta;
P[2][1] = -s_ksi*s_teta;
P[2][2] =  c_teta;

N[0][0] =  c_psi;
N[0][1] = -s_psi*c_eps0;
N[0][2] = -s_psi*s_eps0;

N[1][0] =  s_psi*c_eps;
N[1][1] =  c_psi*c_eps*c_eps0+s_eps*s_eps0;
N[1][2] =  c_psi*c_eps*s_eps0-s_eps*c_eps0;

N[2][0] =  s_psi*s_eps;
N[2][1] =  c_psi*s_eps*c_eps0-c_eps*s_eps0;
N[2][2] =  c_psi*s_eps*s_eps0+c_eps*c_eps0;

for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		NP[i][c]=0;
		for(int j=0; j<3; j++)
			NP[i][c]+=N[i][j]*P[j][c];
	}

}
//---------------------------------------------------------------------------
/*�������������� ������ �������������� � �������� ������������
������� ���������:
		r[3] - ������-������ (x, y, z) [��]
		v[3] - ������ �������� (vx, vy, vz) [��/�]
		t - ����� � ��������� ������
�������� ���������:
		r1[3] - ������-������ (x, y, z) [��]
		v1[3] - ������ �������� (vx, vy, vz) [��/�]
������������ ���������:
		UTCtoTDB - ������� UTC � ������������ �����
		precessionVSnutation - ���� ��������� �������
������������ ������:
		<math.h>
		"constants.h"
		#include "time.h"
*/
void GEItoGEO(double r[3], double v[3], double t, double r1[3], double v1[3]){
double ksi, zeta, teta, psi, eps, eps0;
double 	c_ksi, s_ksi,
		c_zeta, s_zeta,
		c_teta, s_teta,
		c_psi, s_psi,
		c_eps, s_eps,
		c_eps0, s_eps0,
		c_l, s_l;
double T, T2, T3;
double xp, yp;
double Ml,Ms,ul,Ds,Uzl,s;
double w=7.292115e-5;
double A[3][3], B[3][3], B1[3][3], P[3][3], N[3][3], a[3],
		AB[3][3]={0}, AB1[3][3]={0},
		ABN[3][3]={0}, AB1N[3][3]={0},
		ABNP[3][3]={0}, AB1NP[3][3]={0};
T   = (UTCtoTDB(t) - 2451545.0)/36525.0;
T2=T*T;
T3=T2*T;
ksi  = 0.011180860*T + 1.464E-6*T2 + 8.70E-8*T3;
zeta = 0.011180860*T + 5.308E-6*T2 + 8.90E-8*T3;
teta = 0.009717173*T - 2.068E-6*T2 - 2.02E-7*T3;
eps0 = 0.4090928042 - 0.2269655E-3*T - 2.86E-9*T2 + 8.80E-9*T3;

Ml   = 2.355548394 + (1325*pi_2 + 3.470890873)*T + 1.517952E-4*T2 + 3.103E-7*T3;
Ms   = 6.240035940 +   (99*pi_2 + 6.266610600)*T - 2.797400E-6*T2 - 5.820E-8*T3;
ul   = 1.627901934 + (1342*pi_2 + 1.431476084)*T - 6.427170E-5*T2 + 5.340E-8*T3;
Ds   = 5.198469514 + (1236*pi_2 + 5.360106500)*T - 3.340860E-5*T2 + 9.220E-8*T3;
Uzl  = 2.182438624 -    (5*pi_2 + 2.341119397)*T + 3.614290E-5*T2 + 3.880E-8*T3;

eps = 0; psi = 0;
for(int i=0;i<=105;i++){
	psi +=   (coef_A[i] + coef_B[i]*T)*sin(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
	eps +=   (coef_C[i] + coef_D[i]*T)*cos(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
}
psi  = psi*1E-4/3600./180.*pi;
eps  = eps*1E-4/3600./180.*pi;
eps += eps0;

c_ksi=cos(ksi);
s_ksi=sin(ksi);
c_zeta=cos(zeta);
s_zeta=sin(zeta);
c_teta=cos(teta);
s_teta=sin(teta);
c_psi=cos(psi);
s_psi=sin(psi);
c_eps=cos(eps);
s_eps=sin(eps);
c_eps0=cos(eps0);
s_eps0=sin(eps0);



P[0][0] =  c_ksi*c_zeta*c_teta - s_ksi*s_zeta;
P[0][1] = -s_ksi*c_zeta*c_teta - c_ksi*s_zeta;
P[0][2] = -c_zeta*s_teta;

P[1][0] =  c_ksi*s_zeta*c_teta + s_ksi*c_zeta;
P[1][1] = -s_ksi*s_zeta*c_teta + c_ksi*c_zeta;
P[1][2] = -s_zeta*s_teta;

P[2][0] =  c_ksi*s_teta;
P[2][1] = -s_ksi*s_teta;
P[2][2] =  c_teta;

N[0][0] =  c_psi;
N[0][1] = -s_psi*c_eps0;
N[0][2] = -s_psi*s_eps0;

N[1][0] =  s_psi*c_eps;
N[1][1] =  c_psi*c_eps*c_eps0+s_eps*s_eps0;
N[1][2] =  c_psi*c_eps*s_eps0-s_eps*c_eps0;

N[2][0] =  s_psi*s_eps;
N[2][1] =  c_psi*s_eps*c_eps0-c_eps*s_eps0;
N[2][2] =  c_psi*s_eps*s_eps0+c_eps*c_eps0;

s=GMST(t)+psi*c_eps;
c_l=cos(s);
s_l=sin(s);
B[0][0]=c_l;  B[0][1]=s_l; B[0][2]=0;
B[1][0]=-s_l; B[1][1]=c_l; B[1][2]=0;
B[2][0]=0;    B[2][1]=0;   B[2][2]=1;


B1[0][0]=-w*s_l; B1[0][1]=w*c_l;  B1[0][2]=0;
B1[1][0]=-w*c_l; B1[1][1]=-w*s_l; B1[1][2]=0;
B1[2][0]=0;      B1[2][1]=0;      B1[2][2]=0;

xp=0;//0.10157/3600*M_PI/180;
yp=0;//0.49738/3600*M_PI/180;

A[0][0]=1;   A[0][1]=0;  A[0][2]=xp;
A[1][0]=0;   A[1][1]=1;  A[1][2]=-yp;
A[2][0]=-xp; A[2][1]=yp; A[2][2]=1;



for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB[i][c]=0;
		for(int j=0; j<3; j++)
			AB[i][c]+=A[i][j]*B[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		ABN[i][c]=0;
		for(int j=0; j<3; j++)
			ABN[i][c]+=AB[i][j]*N[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		ABNP[i][c]=0;
		for(int j=0; j<3; j++)
			ABNP[i][c]+=ABN[i][j]*P[j][c];
	}
/////////////////////////////////
	for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1[i][c]=0;
		for(int j=0; j<3; j++)
			AB1[i][c]+=A[i][j]*B1[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1N[i][c]=0;
		for(int j=0; j<3; j++)
			AB1N[i][c]+=AB1[i][j]*N[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1NP[i][c]=0;
		for(int j=0; j<3; j++)
			AB1NP[i][c]+=AB1N[i][j]*P[j][c];
	}
//////////////////////////////////////////////
for(int i=0; i<3; i++){
	r1[i]=0;
	for(int j=0; j<3; j++)
		r1[i]+=ABNP[i][j]*r[j];
}

for(int i=0; i<3; i++){
	v1[i]=0;
	for(int j=0; j<3; j++)
		v1[i]+=ABNP[i][j]*v[j]+AB1NP[i][j]*r[j];
}
}
//--------------------------------------------------------------------------
/*�������������� �������� ������������ �� ������ ��������������
������� ���������:
		r[3] - ������-������ (x, y, z) [��]
		v[3] - ������ �������� (vx, vy, vz) [��/�]
		t - ����� � ��������� ������
�������� ���������:
		r1[3] - ������-������ (x, y, z) [��]
		v1[3] - ������ �������� (vx, vy, vz) [��/�]
������������ ���������:
		UTCtoTDB - ������� UTC � ������������ �����
		precessionVSnutation - ���� ��������� �������
������������ ������:
		<math.h>
		"constants.h"
		#include "time.h"
*/
void GEOtoGEI(double r[3], double v[3], double t, double r1[3], double v1[3]){
double ksi, zeta, teta, psi, eps, eps0;
double 	c_ksi, s_ksi,
		c_zeta, s_zeta,
		c_teta, s_teta,
		c_psi, s_psi,
		c_eps, s_eps,
		c_eps0, s_eps0,
		c_l, s_l;
double T, T2, T3;
double xp, yp;
double Ml,Ms,ul,Ds,Uzl,s;
double w=7.292115e-5;
double A[3][3], B[3][3], B1[3][3], P[3][3], N[3][3], a[3],
		AB[3][3]={0}, AB1[3][3]={0},
		ABN[3][3]={0}, AB1N[3][3]={0},
		ABNP[3][3]={0}, AB1NP[3][3]={0};
T   = (UTCtoTDB(t) - 2451545.0)/36525.0;
T2=T*T;
T3=T2*T;
ksi  = 0.011180860*T + 1.464E-6*T2 + 8.70E-8*T3;
zeta = 0.011180860*T + 5.308E-6*T2 + 8.90E-8*T3;
teta = 0.009717173*T - 2.068E-6*T2 - 2.02E-7*T3;
eps0 = 0.4090928042 - 0.2269655E-3*T - 2.86E-9*T2 + 8.80E-9*T3;

Ml   = 2.355548394 + (1325*pi_2 + 3.470890873)*T + 1.517952E-4*T2 + 3.103E-7*T3;
Ms   = 6.240035940 +   (99*pi_2 + 6.266610600)*T - 2.797400E-6*T2 - 5.820E-8*T3;
ul   = 1.627901934 + (1342*pi_2 + 1.431476084)*T - 6.427170E-5*T2 + 5.340E-8*T3;
Ds   = 5.198469514 + (1236*pi_2 + 5.360106500)*T - 3.340860E-5*T2 + 9.220E-8*T3;
Uzl  = 2.182438624 -    (5*pi_2 + 2.341119397)*T + 3.614290E-5*T2 + 3.880E-8*T3;

eps = 0; psi = 0;
for(int i=0;i<=105;i++){
	psi +=   (coef_A[i] + coef_B[i]*T)*sin(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
	eps +=   (coef_C[i] + coef_D[i]*T)*cos(coef_a1[i]*Ml + coef_a2[i]*Ms +coef_a3[i]*ul +coef_a4[i]*Ds +coef_a5[i]*Uzl);
}
psi  = psi*1E-4/3600./180.*pi;
eps  = eps*1E-4/3600./180.*pi;
eps += eps0;

c_ksi=cos(ksi);
s_ksi=sin(ksi);
c_zeta=cos(zeta);
s_zeta=sin(zeta);
c_teta=cos(teta);
s_teta=sin(teta);
c_psi=cos(psi);
s_psi=sin(psi);
c_eps=cos(eps);
s_eps=sin(eps);
c_eps0=cos(eps0);
s_eps0=sin(eps0);



P[0][0] =  c_ksi*c_zeta*c_teta - s_ksi*s_zeta;
P[0][1] = -s_ksi*c_zeta*c_teta - c_ksi*s_zeta;
P[0][2] = -c_zeta*s_teta;

P[1][0] =  c_ksi*s_zeta*c_teta + s_ksi*c_zeta;
P[1][1] = -s_ksi*s_zeta*c_teta + c_ksi*c_zeta;
P[1][2] = -s_zeta*s_teta;

P[2][0] =  c_ksi*s_teta;
P[2][1] = -s_ksi*s_teta;
P[2][2] =  c_teta;

N[0][0] =  c_psi;
N[0][1] = -s_psi*c_eps0;
N[0][2] = -s_psi*s_eps0;

N[1][0] =  s_psi*c_eps;
N[1][1] =  c_psi*c_eps*c_eps0+s_eps*s_eps0;
N[1][2] =  c_psi*c_eps*s_eps0-s_eps*c_eps0;

N[2][0] =  s_psi*s_eps;
N[2][1] =  c_psi*s_eps*c_eps0-c_eps*s_eps0;
N[2][2] =  c_psi*s_eps*s_eps0+c_eps*c_eps0;

s=GMST(t)+psi*c_eps;
c_l=cos(s);
s_l=sin(s);
B[0][0]=c_l;  B[0][1]=s_l; B[0][2]=0;
B[1][0]=-s_l; B[1][1]=c_l; B[1][2]=0;
B[2][0]=0;    B[2][1]=0;   B[2][2]=1;


B1[0][0]=-w*s_l; B1[0][1]=w*c_l;  B1[0][2]=0;
B1[1][0]=-w*c_l; B1[1][1]=-w*s_l; B1[1][2]=0;
B1[2][0]=0;      B1[2][1]=0;      B1[2][2]=0;

xp=0;
yp=0;

A[0][0]=1;   A[0][1]=0;  A[0][2]=xp;
A[1][0]=0;   A[1][1]=1;  A[1][2]=-yp;
A[2][0]=-xp; A[2][1]=yp; A[2][2]=1;

for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB[i][c]=0;
		for(int j=0; j<3; j++)
			AB[i][c]+=A[i][j]*B[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		ABN[i][c]=0;
		for(int j=0; j<3; j++)
			ABN[i][c]+=AB[i][j]*N[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		ABNP[i][c]=0;
		for(int j=0; j<3; j++)
			ABNP[i][c]+=ABN[i][j]*P[j][c];
	}
/////////////////////////////////
	for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1[i][c]=0;
		for(int j=0; j<3; j++)
			AB1[i][c]+=A[i][j]*B1[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1N[i][c]=0;
		for(int j=0; j<3; j++)
			AB1N[i][c]+=AB1[i][j]*N[j][c];
	}
for(int c=0; c<3; c++)
	for(int i=0; i<3; i++){
		AB1NP[i][c]=0;
		for(int j=0; j<3; j++)
			AB1NP[i][c]+=AB1N[i][j]*P[j][c];
	}
//////////////////////////////////////////////
for(int i=0; i<3; i++){
	r1[i]=0;
	for(int j=0; j<3; j++)
		r1[i]+=ABNP[j][i]*r[j];
}
for(int i=0; i<3; i++){
	v1[i]=0;
	for(int j=0; j<3; j++)
		v1[i]+=ABNP[j][i]*v[j]+AB1NP[j][i]*r[j];
}

}
//---------------------------------------------------------------------
/*�������������� ��������� � ������������ (�������. ������ ������)
������� ���������:
		Ar1[][3] - ������-������ (x, y, z) [��]
�������� ���������:
		Ar2[3] - ������-������ (x, y, z) [��]
������������ ���������:
		DegToRad
������������ ������:
		<math.h>
*/
void MAGtoGEO(const double Ar1[], double Ar2[]){
     const double phiM = 11.283;
     const double lambdaM = 289.371;
     double fi, lam;
	 fi = DegToRad(phiM);
     lam = DegToRad(lambdaM);
     Ar2[0] = cos(fi)*cos(lam)*Ar1[0] - sin(lam)*Ar1[1] + sin(fi)*cos(lam)*Ar1[2];
     Ar2[1] = cos(fi)*sin(lam)*Ar1[0] + cos(lam)*Ar1[1] + sin(fi)*sin(lam)*Ar1[2];
     Ar2[2] = - sin(fi)*Ar1[0] + cos(fi)*Ar1[2];
}
//---------------------------------------------------------------------------
/*�������������� ������������ � ��������� (�������. ������ ������)
������� ���������:
		Ar1[][3] - ������-������ (x, y, z) [��]
�������� ���������:
		Ar2[3] - ������-������ (x, y, z) [��]
������������ ���������:
		DegToRad
������������ ������:
		<math.h>
*/
void GEOtoMAG(const double Ar1[], double Ar2[]){
     const double phiM = 11.283;
     const double lambdaM = 289.371;
     double fi, lam;
     fi = DegToRad(phiM);
     lam = DegToRad(lambdaM);
     Ar2[0] = cos(fi)*cos(lam)*Ar1[0] + cos(fi)*sin(lam)*Ar1[1] - sin(fi)*Ar1[2];
     Ar2[1] = - sin(lam)*Ar1[0] + cos(lam)*Ar1[1];
     Ar2[2] =  sin(fi)*cos(lam)*Ar1[0] + sin(fi)*sin(lam)*Ar1[1] + cos(fi)*Ar1[2];
}
//---------------------------------------------------------------------------
/*�������������� J2000 � ��������� (�������. ������ ������)
������� ���������:
		Ar1[][3] - ������-������ (x, y, z) [��]
		t - ����� � ��������� ������
�������� ���������:
		Ar2[3] - ������-������ (x, y, z) [��]
������������ ���������:
		GEItoGEO
		GEOtoMAG
������������ ������:
		<math.h>
*/
void GEItoMAG(const double Ar1[],const double t, double Ar2[]){
  double Ar[3]={0};
  GEItoGEO(Ar1, t, Ar);
  GEOtoMAG(Ar, Ar2);
}
//---------------------------------------------------------------------------
/*�������������� ��������� � J2000 (�������. ������ ������)
������� ���������:
		Ar1[][3] - ������-������ (x, y, z) [��]
		t - ����� � ��������� ������
�������� ���������:
		Ar2[3] - ������-������ (x, y, z) [��]
������������ ���������:
		GEItoGEO
		GEOtoMAG
������������ ������:
		<math.h>
*/
void MAGtoGEI(const double Ar1[],const double t, double Ar2[]){
  double Ar[3]={0};
  MAGtoGEO(Ar1, Ar);
  GEOtoGEI(Ar, t, Ar2);
}
//---------------------------------------------------------------------------
/*��������� �������� �� J2000 �
������������� �� � ������� � ����� L2 � ���� � ������������ �� ������
������� ���������:
		r[3] - ������-������ (x, y, z) [��]
		t - ����� � ��������� ������
�������� ���������:
		r3[3] - ������-������ (x, y, z) [��]
������������ ���������:
		DegToRad
		PLEPH_R
		VectProizv
������������ ������:
		<math.h>
		"jpl_eph.h"
		"mathematics.h"
*/
void J2000toL2(double r[3], double t, double r3[3]){
int NCENTR, NTARG;
double  R_sol, N, S, mu1,
		EKL=DegToRad(23.4390346928518), //���������� ��� �������� ����� � ��������� ���������
		RL=1501000,   //���������� �� ����� L2
	   r_sol[3], v_sol[3], n[3], r2[3]={0},r1[3]={0}, rl[3]={0}, rs[3]={0},
	   AM[3][3]={0};

	   AM[0][0]=1; AM[1][1]=AM[2][2]=cos(EKL); AM[1][2]=sin(EKL); AM[2][1]=-AM[1][2];
NCENTR=3;    //���������� ����������� ���� ����� - 3
NTARG=11;    //���������� �������� ���� ������ - 11

de405.calculateR(NTARG,NCENTR, t, r_sol); //���������� ���������� ������
R_sol=norm(r_sol);

//��������� ������ �� � ������ � ������������� ��������������� ��
for(int i=0; i<3; i++) {
	r1[i]=0;
	rs[i]=0;
	for(int j=0; j<3; j++){
		r1[i]+=AM[i][j]*r[j];
		rs[i]+=AM[i][j]*r_sol[j];
	}
}
for(int i=0; i<3; i++) rl[i]=-rs[i]/R_sol*RL;  //���������� ��������� ����� L2 � ������������� ��
for(int i=0; i<3; i++) r2[i]=r1[i]-rl[i];      //�������� ����� �� � ����� L2
//��������� ������� �������� ���� ��� � ���� ����������� �� �����
for(int i=0; i<3; i++) AM[0][i]=rs[i]/R_sol;
AM[2][0]=AM[2][1]=0; AM[2][2]=1;
VectProizv(AM[2], AM[0], AM[1]);
//�������� �������������� ������ �� ��������� ���������� ��������.
for(int i=0; i<3; i++){
	r3[i]=0;
	for(int j=0; j<3; j++)
		r3[i]+=AM[i][j]*r2[j];
}
}
/*��������� �������� �� ��������� ��������� (X ��������� �� ������) � J2000
������� ���������:
		r[3] - ������-������ (x, y, z) [��]
		v[3] - ������ �������� (vx, vy, vz) [��/�]
		t - ����� � ��������� ������
�������� ���������:
		r1[3] - ������-������ (x, y, z) [��]
		v1[3] - ������ �������� (vx, vy, vz) [��/�]
������������ ���������:
		DegToRad
		PLEPH_RV
		matr_X_vect
		matr_X_matr
������������ ������:
		<math.h>
		"matrix.h"
		"jpl_eph.h"
		"mathematics.h"
*/
   void L2toJ2000 ( double r0[3],double v[3], double t, double r1[3],double v1[3])
   {
	  double 	EKL,    					//���� ����� ���������� � ���������
				RS,                         //���������� �� ������
				VS,                         //�������� ������
				WS,                         //������� �������� �������� ����� ������������ ������
				Ang,                        //���� ����� ����� � ������������� ������ ��
				r_sol[3],                   //������-������ ������
				v_sol[3],                   //������ �������� ������
				r_sol_E[3],                 //������ ������ � ���������
				v_per[3];					//���� ���������� ������������ ��������
	  double 	MatrixX[3][3]={0},          //������� �������� ������ ��� �
				MatrixZ[3][3]={0},          //������� �������� ������ ��� Z ������������� ������� ���������
				MatrixPOV[3][3]={0},        //�������������� ������� �������� J2000->Ekl
				MatrixPOV_T[3][3]={0},      //�������������� ������� �������� Ekl->J2000
				MatrixPOVperen[3][3]={0},   //������� ����������� ���������� �������� J2000->Ekl
				MatrixPOVperen_T[3][3]={0}; //������� ����������� ���������� �������� Ekl->J2000
	  int 	NT=11,
			NC=3;
	  double r[3];
	  EKL = DegToRad(23.4390346928518);	//���� ������� ������� �������� �
									//��������� ��������� �� 01.01.2000
	  r[0]=r0[0]-1500000;
	  r[1]=r0[1];
	  r[2]=r0[2];
	  de405.calculateRV(NT, NC,t,r_sol,v_sol);

	  //����� ������� �������� ������ (� L2)
	  RS = norm(r_sol);
	  VS = norm(v_sol);
	  WS = VS/RS;

	  //������� ������ ��� X �� ���� EKL
	  MatrixX[0][0] = 1.;	MatrixX[0][1] = 0.;			MatrixX[0][2] = 0.;
	  MatrixX[1][0] = 0.;	MatrixX[1][1] = cos(EKL);	MatrixX[1][2] = sin(EKL);
	  MatrixX[2][0] = 0.;	MatrixX[2][1] = -sin(EKL);	MatrixX[2][2] = cos(EKL);

	  //������ ������ � ��������� ���������
	  matr_X_vect(MatrixX, r_sol, r_sol_E);

	  //���� ����� ������-�������� ������ � ������������ �� �����
	  //��������� ������������� [���.]
	  Ang = atan2(r_sol_E[1]/RS, r_sol_E[0]/RS);

	  //������� ������� ��������� �� � ��������� �� (������� ������ Z)
	  MatrixZ[0][0] =  cos(Ang);	MatrixZ[0][1] = sin(Ang);	MatrixZ[0][2] = 0.;
	  MatrixZ[1][0] = -sin(Ang);	MatrixZ[1][1] = cos(Ang);	MatrixZ[1][2] = 0.;
	  MatrixZ[2][0] = 0.;			MatrixZ[2][1] = 0.;			MatrixZ[2][2] = 1.;


	  matr_X_matr( MatrixZ, MatrixX, MatrixPOV);
	  matr_T(MatrixPOV, MatrixPOV_T);

	  //������ ��������� �� � ��, ����� ������� ��������� � �����,
	  //�������� ��������� - ��������� ���������, ��� X ���������� �� ������
	  matr_X_vect(MatrixPOV_T, r, r1);
	  matr_X_vect(MatrixPOV_T, v, v1);


	  //���� ���������� ��������
	  MatrixPOVperen[0][0] = -sin(Ang);		MatrixPOVperen[0][1] =  cos(Ang)*cos(EKL);   MatrixPOVperen[0][2] =  cos(Ang)*sin(EKL);
	  MatrixPOVperen[1][0] = -cos(Ang);		MatrixPOVperen[1][1] = -sin(Ang)*cos(EKL);   MatrixPOVperen[1][2] = -sin(Ang)*sin(EKL);
	  MatrixPOVperen[2][0] =  0.;			MatrixPOVperen[2][1] =  0.;                  MatrixPOVperen[2][2] =  0.;
	  matr_T(MatrixPOVperen, MatrixPOVperen_T);
	  matr_X_vect(MatrixPOVperen_T, r, v_per);

      ///����������������� ��� ��� ���� ������ ���������� ������ � �������� �������

	  //for(int i=0; i<3; i++)
		for(int i=0; i<3; i++) v1[i]+=v_per[i]*WS;
   }
//---------------------------------------------------------------------------
/*��������� �������� �� J2000 � ��������� ��������� (X ��������� �� ������)
������� ���������:
		r[3] - ������-������ (x, y, z) [��]
		v[3] - ������ �������� (vx, vy, vz) [��/�]
		t - ����� � ��������� ������
�������� ���������:
		r1[3] - ������-������ (x, y, z) [��]
		v1[3] - ������ �������� (vx, vy, vz) [��/�]
������������ ���������:
		DegToRad
		PLEPH_RV
		matr_X_vect
		matr_X_matr
������������ ������:
		<math.h>
		"matrix.h"
		"jpl_eph.h"
		"mathematics.h"
*/
   void J2000toEkl ( double r[3],double v[3], double t, double r1[3],double v1[3])
   {
	  double 	EKL,    					//���� ����� ���������� � ���������
				RS,                         //���������� �� ������
				VS,                         //�������� ������
				WS,                         //������� �������� �������� ����� ������������ ������
				Ang,                        //���� ����� ����� � ������������� ������ ��
				r_sol[3],                   //������-������ ������
				v_sol[3],                   //������ �������� ������
				r_sol_E[3],                 //������ ������ � ���������
				v_per[3];					//���� ���������� ������������ ��������
	  double 	MatrixX[3][3]={0},          //������� �������� ������ ��� �
				MatrixZ[3][3]={0},          //������� �������� ������ ��� Z ������������� ������� ���������
				MatrixPOV[3][3]={0},        //�������������� ������� ��������
				MatrixPOVperen[3][3]={0};   //������� ����������� ���������� ��������
	  int 	NT=11,
			NC=3;

	  EKL = DegToRad(23.4390346928518);	//���� ������� ������� �������� �
									//��������� ��������� �� 01.01.2000
	  de405.calculateRV(NT, NC,t,r_sol,v_sol);
	  //����� ������� �������� ������ (� L2)
	  RS = norm(r_sol);
	  VS = norm(v_sol);
	  WS = VS/RS;

	  //������� ������ ��� X �� ���� EKL
	  MatrixX[0][0] = 1.;	MatrixX[0][1] = 0.;			MatrixX[0][2] = 0.;
	  MatrixX[1][0] = 0.;	MatrixX[1][1] = cos(EKL);	MatrixX[1][2] = sin(EKL);
	  MatrixX[2][0] = 0.;	MatrixX[2][1] = -sin(EKL);	MatrixX[2][2] = cos(EKL);

	  //������ ������ � ��������� ���������
	  matr_X_vect(MatrixX, r_sol, r_sol_E);

	  //���� ����� ������-�������� ������ � ������������ �� �����
	  //��������� ������������� [���.]
	  Ang = atan2(r_sol_E[1]/RS, r_sol_E[0]/RS);

	  //������� ������� ��������� �� � ��������� �� (������� ������ Z)
	  MatrixZ[0][0] =  cos(Ang);	MatrixZ[0][1] = sin(Ang);	MatrixZ[0][2] = 0.;
	  MatrixZ[1][0] = -sin(Ang);	MatrixZ[1][1] = cos(Ang);	MatrixZ[1][2] = 0.;
	  MatrixZ[2][0] = 0.;			MatrixZ[2][1] = 0.;			MatrixZ[2][2] = 1.;


	  matr_X_matr( MatrixZ, MatrixX, MatrixPOV);


	  //������ ��������� �� � ��, ����� ������� ��������� � �����,
	  //�������� ��������� - ��������� ���������, ��� X ���������� �� ������
	  matr_X_vect(MatrixPOV, r, r1);
	  matr_X_vect(MatrixPOV, v, v1);


	  //���� ���������� ��������
	  MatrixPOVperen[0][0] = -sin(Ang);		MatrixPOVperen[0][1] =  cos(Ang)*cos(EKL);   MatrixPOVperen[0][2] =  cos(Ang)*sin(EKL);
	  MatrixPOVperen[1][0] = -cos(Ang);		MatrixPOVperen[1][1] = -sin(Ang)*cos(EKL);   MatrixPOVperen[1][2] = -sin(Ang)*sin(EKL);
	  MatrixPOVperen[2][0] =  0.;			MatrixPOVperen[2][1] =  0.;                  MatrixPOVperen[2][2] =  0.;

	  matr_X_vect(MatrixPOVperen, r, v_per);
	 // matr_X_vect(MatrixPOVperen, r, v_per);        ��� ������ � �������� � ���
	 // for(int i=0; i<3; i++)                        �����!!!!!
	 //	for(int i=0; i<3; i++) v1[i]+=v_per[i]*WS;    ��� ����� ���� ����� ������� ������� � ���������
   }

//---------------------------------------------------------------------------
/*��������� �������� �� ��������� ��������� (X ��������� �� ������) � J2000
������� ���������:
		r[3] - ������-������ (x, y, z) [��]
		v[3] - ������ �������� (vx, vy, vz) [��/�]
		t - ����� � ��������� ������
�������� ���������:
		r1[3] - ������-������ (x, y, z) [��]
		v1[3] - ������ �������� (vx, vy, vz) [��/�]
������������ ���������:
		DegToRad
		PLEPH_RV
		matr_X_vect
		matr_X_matr
������������ ������:
		<math.h>
		"matrix.h"
		"jpl_eph.h"
		"mathematics.h"
*/
   void EkltoJ2000 ( double r[3],double v[3], double t, double r1[3],double v1[3])
   {
	  double 	EKL,    					//���� ����� ���������� � ���������
				RS,                         //���������� �� ������
				VS,                         //�������� ������
				WS,                         //������� �������� �������� ����� ������������ ������
				Ang,                        //���� ����� ����� � ������������� ������ ��
				r_sol[3],                   //������-������ ������
				v_sol[3],                   //������ �������� ������
				r_sol_E[3],                 //������ ������ � ���������
				v_per[3];					//���� ���������� ������������ ��������
	  double 	MatrixX[3][3]={0},          //������� �������� ������ ��� �
				MatrixZ[3][3]={0},          //������� �������� ������ ��� Z ������������� ������� ���������
				MatrixPOV[3][3]={0},        //�������������� ������� �������� J2000->Ekl
				MatrixPOV_T[3][3]={0},      //�������������� ������� �������� Ekl->J2000
				MatrixPOVperen[3][3]={0},   //������� ����������� ���������� �������� J2000->Ekl
				MatrixPOVperen_T[3][3]={0}; //������� ����������� ���������� �������� Ekl->J2000
	  int 	NT=11,
			NC=3;

	  EKL = DegToRad(23.4390346928518);	//���� ������� ������� �������� �
									//��������� ��������� �� 01.01.2000

	  de405.calculateRV(NT, NC,t,r_sol,v_sol);

	  //����� ������� �������� ������ (� L2)
	  RS = norm(r_sol);
	  VS = norm(v_sol);
	  WS = VS/RS;

	  //������� ������ ��� X �� ���� EKL
	  MatrixX[0][0] = 1.;	MatrixX[0][1] = 0.;			MatrixX[0][2] = 0.;
	  MatrixX[1][0] = 0.;	MatrixX[1][1] = cos(EKL);	MatrixX[1][2] = sin(EKL);
	  MatrixX[2][0] = 0.;	MatrixX[2][1] = -sin(EKL);	MatrixX[2][2] = cos(EKL);

	  //������ ������ � ��������� ���������
	  matr_X_vect(MatrixX, r_sol, r_sol_E);

	  //���� ����� ������-�������� ������ � ������������ �� �����
	  //��������� ������������� [���.]
	  Ang = atan2(r_sol_E[1]/RS, r_sol_E[0]/RS);

	  //������� ������� ��������� �� � ��������� �� (������� ������ Z)
	  MatrixZ[0][0] =  cos(Ang);	MatrixZ[0][1] = sin(Ang);	MatrixZ[0][2] = 0.;
	  MatrixZ[1][0] = -sin(Ang);	MatrixZ[1][1] = cos(Ang);	MatrixZ[1][2] = 0.;
	  MatrixZ[2][0] = 0.;			MatrixZ[2][1] = 0.;			MatrixZ[2][2] = 1.;


	  matr_X_matr( MatrixZ, MatrixX, MatrixPOV);
	  matr_T(MatrixPOV, MatrixPOV_T);

	  //������ ��������� �� � ��, ����� ������� ��������� � �����,
	  //�������� ��������� - ��������� ���������, ��� X ���������� �� ������
	  matr_X_vect(MatrixPOV_T, r, r1);
	  matr_X_vect(MatrixPOV_T, v, v1);


	  //���� ���������� ��������
	  MatrixPOVperen[0][0] = -sin(Ang);		MatrixPOVperen[0][1] =  cos(Ang)*cos(EKL);   MatrixPOVperen[0][2] =  cos(Ang)*sin(EKL);
	  MatrixPOVperen[1][0] = -cos(Ang);		MatrixPOVperen[1][1] = -sin(Ang)*cos(EKL);   MatrixPOVperen[1][2] = -sin(Ang)*sin(EKL);
	  MatrixPOVperen[2][0] =  0.;			MatrixPOVperen[2][1] =  0.;                  MatrixPOVperen[2][2] =  0.;
	  matr_T(MatrixPOVperen, MatrixPOVperen_T);
	  matr_X_vect(MatrixPOVperen_T, r, v_per);

      ///����������������� ��� ��� ���� ������ ���������� ������ � �������� �������

	 // for(int i=0; i<3; i++)
		//for(int i=0; i<3; i++) v1[i]+=v_per[i]*WS;
   }
//---------------------------------------------------------------------
//��������� �������� �� J2000 � ���������������� ��������� �������
void J2000toSSB(double r_J2000[3], double v_J2000[3], double t, double r_SSB[3], double v_SSB[3]){

int	NT=3, 	//earth
	NC=12; 	// solar-system barycenter
double r[3], v[3];
	de405.calculateRV(NT, NC,t,r,v);

	for(int i=0; i<3; i++){
		r_SSB[i]=r_J2000[i]+r[i];
		v_SSB[i]=v_J2000[i]+v[i];
	}
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//��������� �������� �� ���������������� ��������� ������� � J2000
void SSBtoJ2000(double r_SSB[3], double v_SSB[3], double t, double r_J2000[3], double v_J2000[3]){

int	NT=3, 	//earth
	NC=12; 	// solar-system barycenter
double r[3], v[3];
	de405.calculateRV(NT, NC,t,r,v);

	for(int i=0; i<3; i++){
		r_J2000[i]=r_SSB[i]-r[i];
		v_J2000[i]=v_SSB[i]-v[i];
	}
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//��������� �������� �� J2000 � ����������������� ��
void J2000toSG(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
double r[3], v[3];
J2000toSC(r_in,v_in,JD, r, v);
SCtoSG(r, v, JD, r_out, v_out);
}
//---------------------------------------------------------------------
//��������� �������� �� ����������������� �� � J2000
void SGtoJ2000(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
double r[3], v[3];

SGtoSC(r_in,v_in,JD, r, v);
SCtoJ2000(r, v, JD, r_out, v_out);
}
//---------------------------------------------------------------------
//��������� �������� �� J2000 � ������������������ ��
void J2000toSC(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
int	NT=10, 	// moon
	NC=3; 	 //earth
double r_moon[3], v_moon[3];
	de405.calculateRV(NT, NC,JD,r_moon,v_moon);
for(int i=0; i<3; i++){
	r_out[i]=r_in[i]-r_moon[i];
	v_out[i]=v_in[i]-v_moon[i];
}

}

//---------------------------------------------------------------------
//��������� �������� �� ������������������ �� � J2000
void SCtoJ2000(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
int	NT=10, 	// moon
	NC=3; 	//earth
double r_moon[3], v_moon[3];
	de405.calculateRV(NT, NC,JD,r_moon,v_moon);

for(int i=0; i<3; i++){
	r_out[i]=r_in[i]+r_moon[i];
	v_out[i]=v_in[i]+v_moon[i];
}

}
//---------------------------------------------------------------------
//��������� �������� �� ������������������ �� � ����������������� ��
void SCtoSG(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
double R_in[6];
double R_out[6];
for(int i=0; i<3; i++){
	R_in[i]=r_in[i];
	R_in[i+3]=v_in[i];
}
SCtoSG(R_in,JD, R_out);
for(int i=0; i<3; i++){
	 r_out[i]=R_out[i];
	 v_out[i]=R_out[i+3];
}
}
//---------------------------------------------------------------------
//��������� �������� �� ������������������ �� � ����������������� ��
 void SCtoSG(double R_in[6],double JD, double R_out[6]){

	double dt = JD + (dUTC+dUT1)/86400. - 2451545.;
	double T = dt/36525.;
	double E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13;
	double sin_E1,sin_E2,sin_E3,sin_E4,sin_E6,sin_E10,sin_E13;

	E1= DegToRad(125.045- 0.0529921*dt);
	E2= DegToRad(250.089- 0.1059842*dt);
	E3= DegToRad(260.008+13.0120009*dt);
	E4= DegToRad(176.625+13.3407154*dt);
	E5= DegToRad(357.529+ 0.9856003*dt);
	E6= DegToRad(311.589+26.4057084*dt);
	E7= DegToRad(134.963+13.0649930*dt);
	E8= DegToRad(276.617+ 0.3287146*dt);
	E9= DegToRad( 34.226+ 1.7484877*dt);
	E10=DegToRad( 15.134- 0.1589763*dt);
	E11=DegToRad(119.743+ 0.0036096*dt);
	E12=DegToRad(239.961+ 0.1643573*dt);
	E13=DegToRad( 25.053+12.9590088*dt);

	sin_E1  = sin(E1);
    sin_E2  = sin(E2);
	sin_E3  = sin(E3);
    sin_E4  = sin(E4);
	sin_E6  = sin(E6);
    sin_E10 = sin(E10);
    sin_E13 = sin(E13);

	double alpha = 269.9949 + 0.0031*T      - 3.8787*sin_E1	    - 0.1204*sin_E2
				            + 0.07*sin_E3	- 0.0172*sin_E4	    + 0.0072*sin_E6
											- 0.0052*sin_E10	+ 0.0043*sin_E13;

    double beta  = 66.5392  + 0.013*T	    + 1.5419*cos(E1)	+ 0.0239*cos(E2)
							- 0.0278*cos(E3)+ 0.0068*cos(E4)    - 0.0029*cos(E6)
				            + 0.0009*cos(E7)+ 0.0008*cos(E10)   - 0.0009*cos(E13);

    double om = 38.3213	+ 13.17635815*dt	- 1.4e-12*dt*dt		+ 3.5610*sin_E1
                        + 0.1208*sin_E2     - 0.0642*sin_E3 	+ 0.0158*sin_E4
						+ 0.0252*sin(E5)	- 0.0066*sin_E6  	- 0.0047*sin(E7)
                        - 0.0046*sin(E8)	+ 0.0028*sin(E9)	+ 0.0052*sin_E10
                        + 0.0040*sin(E11)	+ 0.0019*sin(E12)	- 0.0044*sin_E13;

    alpha=DegToRad(alpha);
	beta=DegToRad(beta);
    om=DegToRad(om);

	double cos_om    = cos(om);
    double cos_alpha = cos(alpha);
	double cos_beta  = cos(beta);
    double sin_om    = sin(om);
    double sin_alpha = sin(alpha);
    double sin_beta  = sin(beta);


	double M[3][3]={0};

        M[0][0] = -sin_alpha*cos_om - cos_alpha*sin_beta*sin_om;
		M[1][0] = cos_alpha*cos_om - sin_alpha*sin_beta*sin_om;
        M[2][0] = sin_om*cos_beta;
        M[0][1] = sin_alpha*sin_om - cos_alpha*sin_beta*cos_om;
		M[1][1] = -cos_alpha*sin_om - sin_alpha*sin_beta*cos_om;
		M[2][1] = cos_om*cos_beta;
		M[0][2] = cos_alpha*cos_beta;
		M[1][2] = sin_alpha*cos_beta;
		M[2][2] = sin_beta;


	for(int i=0;i<6;i++) R_out[i]=0;
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			R_out[i]   += R_in[j]   *M[j][i];
			R_out[i+3] += R_in[j+3] *M[j][i];
		}
	}
}
//---------------------------------------------------------------------
//��������� �������� �� ����������������� �� � ������������������ ��
void SGtoSC(double r_in[3], double v_in[3],double JD, double r_out[3], double v_out[3]){
double R_in[6];
double R_out[6];
for(int i=0; i<3; i++){
	R_in[i]=r_in[i];
	R_in[i+3]=v_in[i];
}
SGtoSC(R_in,JD, R_out);
for(int i=0; i<3; i++){
	 r_out[i]=R_out[i];
	 v_out[i]=R_out[i+3];
}
}
//---------------------------------------------------------------------
//��������� �������� �� ����������������� �� � ������������������ ��
 void SGtoSC(double R_in[6],double JD, double R_out[6]){

	double dt = JD + (dUTC+dUT1)/86400. - 2451545.;
	double T = dt/36525.;
	double E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13;
	double sin_E1,sin_E2,sin_E3,sin_E4,sin_E6,sin_E10,sin_E13;

    E1= DegToRad(125.045- 0.0529921*dt);
	E2= DegToRad(250.089- 0.1059842*dt);
	E3= DegToRad(260.008+13.0120009*dt);
	E4= DegToRad(176.625+13.3407154*dt);
	E5= DegToRad(357.529+ 0.9856003*dt);
	E6= DegToRad(311.589+26.4057084*dt);
	E7= DegToRad(134.963+13.0649930*dt);
	E8= DegToRad(276.617+ 0.3287146*dt);
	E9= DegToRad( 34.226+ 1.7484877*dt);
	E10=DegToRad( 15.134- 0.1589763*dt);
	E11=DegToRad(119.743+ 0.0036096*dt);
	E12=DegToRad(239.961+ 0.1643573*dt);
	E13=DegToRad( 25.053+12.9590088*dt);

	sin_E1  = sin(E1);
    sin_E2  = sin(E2);
    sin_E3  = sin(E3);
    sin_E4  = sin(E4);
	sin_E6  = sin(E6);
    sin_E10 = sin(E10);
    sin_E13 = sin(E13);

	double alpha = 269.9949 + 0.0031*T      - 3.8787*sin_E1	    - 0.1204*sin_E2
				            + 0.07*sin_E3	- 0.0172*sin_E4	    + 0.0072*sin_E6
                                            - 0.0052*sin_E10	+ 0.0043*sin_E13;

    double beta  = 66.5392  + 0.013*T	    + 1.5419*cos(E1)	+ 0.0239*cos(E2)
							- 0.0278*cos(E3)+ 0.0068*cos(E4)    - 0.0029*cos(E6)
				            + 0.0009*cos(E7)+ 0.0008*cos(E10)   - 0.0009*cos(E13);

    double om = 38.3213	+ 13.17635815*dt	- 1.4e-12*dt*dt		+ 3.5610*sin_E1
                        + 0.1208*sin_E2     - 0.0642*sin_E3 	+ 0.0158*sin_E4
						+ 0.0252*sin(E5)	- 0.0066*sin_E6  	- 0.0047*sin(E7)
                        - 0.0046*sin(E8)	+ 0.0028*sin(E9)	+ 0.0052*sin_E10
                        + 0.0040*sin(E11)	+ 0.0019*sin(E12)	- 0.0044*sin_E13;

    alpha=DegToRad(alpha);
	beta=DegToRad(beta);
    om=DegToRad(om);

    double cos_om    = cos(om);
    double cos_alpha = cos(alpha);
	double cos_beta  = cos(beta);
    double sin_om    = sin(om);
    double sin_alpha = sin(alpha);
    double sin_beta  = sin(beta);


    double M[3][3]={0};

        M[0][0] = -sin_alpha*cos_om - cos_alpha*sin_beta*sin_om;
        M[0][1] = cos_alpha*cos_om - sin_alpha*sin_beta*sin_om;
		M[0][2] = sin_om*cos_beta;
        M[1][0] = sin_alpha*sin_om - cos_alpha*sin_beta*cos_om;
        M[1][1] = -cos_alpha*sin_om - sin_alpha*sin_beta*cos_om;
        M[1][2] = cos_om*cos_beta;
        M[2][0] = cos_alpha*cos_beta;
		M[2][1] = sin_alpha*cos_beta;
		M[2][2] = sin_beta;


	for(int i=0;i<6;i++) R_out[i]=0;
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			R_out[i]   += R_in[j]   *M[j][i];
			R_out[i+3] += R_in[j+3] *M[j][i];
		}
	}
}














//---------------------------------------------------------------------

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//��������� ������������ ��� �������������� ��������� ����� ��� ������ ������ IGRF
void GEOtoNEPONYATNAYA(double GEO[], double N[], double &phi, double &lam)
{
   double R, cf,cl, sf,sl;
   xyzTOphilam(GEO[0], GEO[1], GEO[2], phi, lam, R);
   cf=cos(phi);
   cl=cos(lam);
   sf=sin(phi);
   sl=sin(lam);
   N[0] = GEO[0]*cl*cf + GEO[1]*sl - GEO[2]*cl*sf;
   N[1] = -GEO[0]*sl*cf + GEO[1]*cl + GEO[2]*sl*sf;
   N[2] = GEO[0]*sf + GEO[2]*cf;
}

//---------------------------------------------------------------------
//��������� ������������ ��� �������������� ��������� ����� ��� ������ ������ IGRF
void NEPONYATNAYAtoGEO(double N[], const double &phi, const double &lam, double GEO[])
{
   double R, cf,cl, sf,sl;
   cf=cos(phi);
   cl=cos(lam);
   sf=sin(phi);
   sl=sin(lam);
   GEO[0] = -N[0]*cl*sf - N[1]*sl - N[2]*cl*cf;
   GEO[1] = -N[0]*sl*sf + N[1]*cl - N[2]*sl*cf;
   GEO[2] =  N[0]*cf              - N[2]*sf;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//������� �� ���� ������ ��������� � ���������� ������������ � ��������� ���������������� XY (��������������)
void  perexod_xyz_v_nb(double xyz[], double &x_n, double &y_n)
{
  y_n = xyz[2];
  x_n = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
}
//---------------------------------------------------------------------------
//������� �� ���������� ������������ � ��������� ���������������� XY (��������������) �  ���� ������ �����������
void  perexod_nb_v_xyz(const double x_n, const double y_n, double xyz[])
{
  double k;
  k=y_n/xyz[2];
  xyz[2]=y_n;
  xyz[1] = xyz[1]*k;
  xyz[0] = xyz[0]*k;
}