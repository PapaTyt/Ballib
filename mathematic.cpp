
//---------------------------------------------------------------------------
#include <math.h>

#pragma hdrstop

#include "mathematic.h"
#include "report.h"


//---------------------------------------------------------------------------

#pragma package(smart_init)
//������� ���������� ������ ������� ����������� ����� ������������
double norm(double x, double y, double z){
	return sqrt(x*x+y*y+z*z);
}
//������� ���������� ������ �������
double norm(double r[5]){
	return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}
//������� ���������� ������ �������� �������� r1 � r2
double norm(double r1[3], double r2[3]){
return sqrt((r1[0]-r2[0])*(r1[0]-r2[0])+
			(r1[1]-r2[1])*(r1[1]-r2[1])+
			(r1[2]-r2[2])*(r1[2]-r2[2]));
}
//��������� ������������ �������� �1 � �2
double Skalyar(const double x1[3], const double x2[3]){
   return x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
}
//---------------------------------------------------------------------------
//���� ����� ���������
double angle_between_vectors(double r1[3], double r2[3]){
double angle;
double R1=norm(r1);
double R2=norm(r2);
if(R1==0) report(1);
if(R2==0) report(2);
	angle=acos(Skalyar(r1, r2)/R1/R2);
	return angle;
}
//��������� ������������ �������� Vec1 � Vec2
void VectProizv (const double Vec1[3], const double Vec2[3], double proizv[3]){
  proizv[0] = Vec1[1]*Vec2[2]-Vec1[2]*Vec2[1];
  proizv[1] = -Vec1[0]*Vec2[2]+Vec1[2]*Vec2[0];
  proizv[2] = Vec1[0]*Vec2[1]-Vec1[1]*Vec2[0];
}
//---------------------------------------------------------------------------
//������� ������� dV1 ������ ��� ���������� ������ Ort �� ���� �� � ��������
void Povorot_Vektora(const double Ort[], const double dV1[], const double hi,double dV2[])
{
//////////////////////////////////////////////////////////////////////////
//������� ���������:                                                    //
//          Ort - ���� ��� ��������                                     //
//          dV1 - ������ ������� ���� ���������                         //
//          hi - ���� �������� � ��������                               //
//�������� ���������:                                                   //
//        	dV2 - ������ - ��������� �������� ��������� �������         //
//������������ ���������:                                               //
//			VectProizv (mathematics.h)                             		//
//////////////////////////////////////////////////////////////////////////
	double  Var1[3]={0};
	double  Var2[3]={0};
	double  Var3[3]={0};
	double  teta[3]={0};
	double  r;

    teta[0]=2*Ort[0]*tan(hi/2.);
    teta[1]=2*Ort[1]*tan(hi/2.);
    teta[2]=2*Ort[2]*tan(hi/2.);
    VectProizv (teta,dV1,Var1);
    Var2[0] = dV1[0] + Var1[0]/2.;
    Var2[1] = dV1[1] + Var1[1]/2.;
    Var2[2] = dV1[2] + Var1[2]/2.;
    VectProizv (teta,Var2, Var3);
    r = sqrt(teta[0]*teta[0]+teta[1]*teta[1]+teta[2]*teta[2]);
    dV2[0] = dV1[0] + Var3[0]/(1. + r*r/4.);
    dV2[1] = dV1[1] + Var3[1]/(1. + r*r/4.);
    dV2[2] = dV1[2] + Var3[2]/(1. + r*r/4.);
}