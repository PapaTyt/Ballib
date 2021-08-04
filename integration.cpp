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
void rightPart(VECTOR &rv){


for(int i=0; i<3; i++){
	a_central_field[i]=0;
	a_off_central_field[i]=0;
	a_celestial_bodies[i]=0;
	a_solar_radiation[i]=0;
	a_atmosphere[i]=0;
	a_traction[i]=0;
}



switch(centralBody){
	case B_EARTH:	if(rp[0]) central_field(rv);
					if(rp[1]) off_central_field(rv);
					if(rp[2]) celestial_bodies(rv);
					if(rp[3]) solar_radiation(rv);
					if(rp[4]) atmosphere(rv);
					if(rp[5]) traction(rv);
					break;


	case B_MOON: 	if(rp[0]) central_field_moon(VECTOR rv);
					if(rp[1]) off_central_field_moon(VECTOR rv);
					if(rp[2]) celestial_bodies_moon(VECTOR rv);
					break;


	case B_SUN:     break;

}

for (int i=0; i<3; i++){
	rv.f[i] = a_central_fild[i]+
			  a_off_central_fild[i]+
			  a_celestial_bodies[i]+
			  a_solar_pressure[i]+
			  a_atmosphere[i]+
			  a_traction[i];
}

if(calculeteMatrix){

for(int i=0; i<6; i++)
	for(int j=0; j<6; j++)
		rv.dfdx[i][j]=0;

rv.dfdx[0][3]=1;
rv.dfdx[1][4]=1;
rv.dfdx[2][5]=1;

for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
		rv.dfdx[i+3][j]=df_central_fild[i][j]+
				   df_off_central_fild[i][j]+
				   df_celestial_bodies[i][j];


matr_X_matr(rv.dfdx, rv.F, rv.F_);


}








}

	/*���������� ������������ ���������, �������������� �����������
	  ������������� ����� ����� � ������� ������� ����������� ����� �������*/
	void central_field(VECTOR rv);

	/*���������� ������������ ���������, �������������� �����������
	  ������������� ����� ���� � ������� ������� ����������� ����� �������*/
	void central_field_moon(VECTOR rv);

	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� ����� � ������� ������� ����������� ����� �������*/
	void off_central_field(VECTOR rv);

	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� ���� � ������� ������� ����������� ����� �������*/
	void off_central_field_moon(VECTOR rv);


	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� �����(������ ��������� ���������� - �20)*/
	void off_central_field_C20(VECTOR rv);

	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� �����(������ ��������� ���������� - �40)*/
	void off_central_field_C40(VECTOR rv);


	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� ����� � ������ �������� �� 32�32 � ������� �������
	  ����������� ����� �������*/
	void off_central_field_32(VECTOR rv, double df[3]);

	/*���������� ������������ ��������� �������������� ����������������
	  ��������������� ���� ���� � ������ �������� �� 75�75 � ������� �������
	  ����������� ����� �������*/
	void off_central_field_75_moon(VECTOR rv, double df[3]);

	/*���������� ��������� ������������� ��������� �������� ���
	  (���������� ��� 10 �������� ���, ����������� ���� �����) � ������� �������
	  ����������� ����� �������*/
	void celestial_bodies(VECTOR rv);

    /*���������� ��������� ������������� ��������� �������� ���
	  (���������� ��� 10 �������� ���, ����������� ���� ����) � ������� �������
	  ����������� ����� �������*/
	void celestial_bodies_moon(VECTOR rv);

	/*���������� ��������� �������������� ��������� ���������� ��������� �
	  ������� ������� ����������� ����� �������*/
	void solar_radiation(VECTOR rv);

	/*���������� ��������� �������������� ������������ ���� �������������
	  ��������� � ������������ � ���� � 25645.166-2004*/
	void atmosphere(VECTOR rv);

	/*���������� ��������� �������������� ������������ ���� �������������
	  ��������� � ������������ � ���� � 25645.166-2004*/
	void atmosphereGOST2004(VECTOR rv);

	/*���������� ��������� �������������� ������� ������������ ���������*/
	void traction(VECTOR rv);
