//---------------------------------------------------------------------------


#pragma hdrstop

#include "integration.h"

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
	r_nu[i]=r[i];
	v_nu[i]=v[i];
}
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
	r_nu[i]=rv[i];
	v_nu[i]=rv[i+3];
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


	case B_MOON: 	if(rp[0]) central_field_moon(VECTOR rv);       //����������� ����
					if(rp[1]) off_central_field_moon(VECTOR rv);   //��������������� 
					if(rp[2]) celestial_bodies_moon(VECTOR rv);    //�������� ����
					break;


	case B_SUN:     break;

}

/*������������ ���������� �������� ����������� ���������*/
for (int i=0; i<3; i++){
	rv.f[i] = a_central_fild[i]+
			  a_off_central_fild[i]+
			  a_celestial_bodies[i]+
			  a_solar_pressure[i]+
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
			rv.dfdx[i+3][j]=df_central_fi�ld[i][j]+
					   df_off_central_fi�ld[i][j]+
					   df_celestial_bodies[i][j];

	matr_X_matr(rv.dfdx, rv.F, rv.F_);
}
}

/*���������� ������������ ���������, �������������� �����������
  ������������� ����� ����� � ������� ������� ����������� ����� �������*/
void chi::integration::central_field(VECTOR rv){
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
	a_central_fi�ld[k] =  - mu*rv.r[k]/R3;

/*���������� ������� ������� ����������� ������� �������������� ���������� 
  ������������ ���� �� ������� ���������*/
if(calculeteMatrix){  
	for(int i=0;i<3;i++){
		for(int j=0; j<3;j++)
			df_central_field[i][j]=rv.r[i]*rv.r[j]/R2;
	}
	for(int i=0; i<3; i++)  df_central_fild[i][i]-=1./3.;
	for(int i=0; i<3;i++){
		for(int j=0; j<3; j++)
			df_central_field[i][j]*=3*mu/R3;
	}
}
}

/*���������� ������������ ���������, �������������� �����������
  ������������� ����� ���� � ������� ������� ����������� ����� �������*/
void chi::integration::central_field_moon(VECTOR rv){
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
	a_central_fi�ld[k] =  - mum*rv.r[k]/R3;

/*���������� ������� ������� ����������� ������� �������������� ���������� 
  ������������ ���� �� ������� ���������*/
if(calculeteMatrix){  
	for(int i=0;i<3;i++){
		for(int j=0; j<3;j++)
			df_central_field[i][j]=rv.r[i]*rv.r[j]/R2;
	}
	for(int i=0; i<3; i++)  df_central_fild[i][i]-=1./3.;
	for(int i=0; i<3;i++){
		for(int j=0; j<3; j++)
			df_central_field[i][j]*=3*mum/R3;
	}
}
}

/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� ����� � ������� ������� ����������� ����� �������*/
void chi::integration::off_central_field(VECTOR rv){
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
	case H32: off_central_field_32(rv, df);
			  if(calculeteMatrix){
				  double delta=0.001;
				  VECTOR temp;  
				  for(int i=0; i<3; i++) {
					  temp=rv;
					  temp.r[i]+=delta
					  off_central_field_32(rv, df);
					  for(int j=0; j<3; j++) df_off_central_field[i][j]=df[j];
					  temp=rv;
					  temp.r[i]-=delta
					  off_central_field_32(rv, df);
					  for(int j=0; j<3; j++) df_off_central_field[i][j]-=df[j];
					  for(int j=0; j<3; j++) df_off_central_fild[i][j]/=2*delta;
				  }

			  }
			  break;

}
}

/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� ���� � ������� ������� ����������� ����� �������*/
void chi::integration::off_central_field_moon(VECTOR rv){
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

double df[3]={0};

off_central_field_75_moon(rv, df);

if(calculeteMatrix){
	double delta=0.001;
	VECTOR temp;  
	for(int i=0; i<3; i++) {
		  temp=rv;
		  temp.r[i]+=delta
		  off_central_field_75_moon(rv, df);
		  for(int j=0; j<3; j++) df_off_central_field[i][j]=df[j];
		  temp=rv;
		  temp.r[i]-=delta
		  off_central_field_32(rv, df);
		  for(int j=0; j<3; j++) df_off_central_field[i][j]-=df[j];
		  for(int j=0; j<3; j++) df_off_central_fild[i][j]/=2*delta;
	}
}
}


/*���������� ������������ ��������� �������������� ����������������
  ��������������� ���� �����(������ ��������� ���������� - �20)*/
void chi::integration::off_central_field_C20(VECTOR rv){
/*
	 *  ���������:
	 *   rv � ��������� ���� VECTOR (������ ��������� � �������)
	 *
	 *  ����������:
	 *    � 
	 *
	 */

/*���������� ������ ������-������� ��� �������� � ����� �������*/
double R=norm(rv.r);
double R2=R*R;
double R5=R2*R2*R;
/*������������ �����������*/
double b2=17599253992.788798202843556774502;

/*���������� ��������� ������������ ������ ��������� ���������*/
a_off_central_fild[0]=3*b2*rv.r[0]*(5*rv.r[2]*rv.r[2]/R2-1)/2/R5;
a_off_central_fild[1]=3*b2*rv.r[1]*(5*rv.r[2]*rv.r[2]/R2-1)/2/R5;
a_off_central_fild[2]=3*b2*rv.r[2]*(5*rv.r[2]*rv.r[2]/R2-3)/2/R5;






}
	
	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� �����(������ ��������� ���������� - �40)*/
	void chi::integration::off_central_field_C40(VECTOR rv);


	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� ����� � ������ �������� �� 32�32 � ������� �������
	  ����������� ����� �������*/
	void chi::integration::off_central_field_32(VECTOR rv, double df[3]);

	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� ���� � ������ �������� �� 75�75 � ������� �������
	  ����������� ����� �������*/
	void chi::integration::off_central_field_75_moon(VECTOR rv, double df[3]);

	/*���������� ��������� ������������� ��������� �������� ���
	  (���������� ��� 10 �������� ���, ����������� ���� �����) � ������� �������
	  ����������� ����� �������*/
	void chi::integration::celestial_bodies(VECTOR rv);

    /*���������� ��������� ������������� ��������� �������� ���
	  (���������� ��� 10 �������� ���, ����������� ���� ����) � ������� �������
	  ����������� ����� �������*/
	void chi::integration::celestial_bodies_moon(VECTOR rv);

	/*���������� ��������� �������������� ��������� ���������� ��������� �
	  ������� ������� ����������� ����� �������*/
	void chi::integration::solar_radiation(VECTOR rv);

	/*���������� ��������� �������������� ������������ ���� �������������
	  ��������� � ������������ � ���� � 25645.166-2004*/
	void chi::integration::atmosphere(VECTOR rv);

	/*���������� ��������� �������������� ������������ ���� �������������
	  ��������� � ������������ � ���� � 25645.166-2004*/
	void chi::integration::atmosphereGOST2004(VECTOR rv);

	/*���������� ��������� �������������� ������� ������������ ���������*/
	void chi::integration::traction(VECTOR rv);
