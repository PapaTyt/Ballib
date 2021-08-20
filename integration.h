//---------------------------------------------------------------------------

#ifndef integrationH
#define integrationH
//---------------------------------------------------------------------------
namespace chi {

//��������� ������� ��� ���������� ��������������
struct VECTOR{
	double r[3];	//������-������ ��
	double v[3];	//������ �������� ��
	double f[3];	//������ ������ ������ (����������� ���������)
	double t;		//����� � �������� �� ������� ������� T_NU
	double F[6][6];
	double dfdx[6][6];
	double F_[6][6];
};

/* �������� ���������� ��� ������ dph::EphemerisReelase::calculateBody(...).
   ������ �������� ������������� �������� ���, ��� ������� ����� ��������
   ��������� ����������. */
enum Body{
	B_EARTH   = 0,
	B_MOON    = 1,
	B_SUN     = 2
};
/*�������� ��������� ��� ������ chi::integration::off_central_field(VECTOR rv)
  ������ �������� ���������� �������� ���������� ��������������� ����������
  �����*/
enum Harmonics{
	C20 = 0, /*���������� �� �������� ������ ��������� ���������*/
	C40 = 1, /*���������� �� �������� ��������� ��������� ���������*/
	H32	= 2  /*���������� ���������� �� ������� 32�32*/
};



class integration{

public:
	/* [����������� ������] */

	/*������ ��*/
	void set_NU(double r[3], double v[3], double t);
	void set_NU(double rv[6],  double t);

	/*������ ���������� ��������������*/
	void setParametrs(double interval_, double step_);

	/*������ ���� ���������� �� ������ ���� ��������������*/
	void setTypeCalculation(typeCalculation_);




	/* [������ ���������� ��������������] */

	/*������������ ����� �����-����� 4-�� ������� ��� ������������
	  ���������������� ��������� 1-�� �������*/
	void RK4_ODE1();
	/*����� �����-����� 4-�� ������� ��� ������������
	  ���������������� ��������� 2-�� �������*/
	void RK4_ODE2();
	/*����� �����-�����-������� 4-�� �������*/
	/*����� �������-������ 7-�� �������*/
	void DP7();
	/*����� �����-�����-���������� 8-�� �������*/
	void RKF8();
	/*����� ������-��������-������� 8-�� �������*/
	void ABM8();
	/*����� �������������*/
	void Extrapolation();
	/**/


	/**/

//private:
protected:
	/*��������� ������� ��������������*/
	VECTOR rv_nu;	//��������� ������-������ [��]
	double t_nu;	//��������� ������ ������� � ������� ��������� ���� (UTC)


	/*��������� ���������� ��������������*/
	double dt;
	double interval;
	double step;

	/*��������� ���������� ��������������*/
	unsigned harmonicType; 	  /*������� ����� ���� ������������� ��������:
								0 - ���������� ��������������� ����������:
									������ ��������� ��������� �20
								1 - ���������� ��������������� ����������:
									��������� ��������� ��������� �40
								2 - ���������� ��������������� ����������:
									������� ���������� �� 32�32*/
	unsigned harmonicOrder;    /*������� ���������� �� 32�32*/
	bool planet[11];		/*������� ����� �������������� ���������� �������� ���:
						   0 - ��������
						   1 - ������
						   2 - �����
						   3 - ����
						   4 - ������
						   5 - ������
						   6 - ����
						   7 - ������
						   8 - ������
						   9 - ����
						  10 - ������ */
	unsigned	centralBody; /*������� ������������ ����:
								B_EARTH - 0 � �����;
								B_MOON  - 1 � ����;
								B_SUN   - 2 � ������;*/

	unsigned calculeteMatrix; /*������� ������� ������� ������� �����������*/


	unsigned rp[5];	/*������� ����� ����������:
					  rp[0] - ���� ������������ ��������������� ����
					  rp[1] - ���� ��������������� ��������������� ����
					  rp[2] - ���� �������� ��� ��������� �������
					  rp[3] - ����
					  rp[4] - ���� ���������
					  rp[5] - ���� ���� ������������ ���������*/

	/*��������� ���� ���������� �� ������ ���� ��������������*/
	unsigned typeCalculation;
	/* [���������� �� ������ ���� ��������������] */



	/*������� ���������: */
	double a_central_field[3];      /*��������� ������������� ���������
									 ������������ ��������������� ����*/
	double a_off_central_field[3];  /*��������� ������������� ���������
									 ���������������� ��������������� ����*/
	double a_celestial_bodies[3];  /*��������� ������������� ���������
									 �������� ��� ��������� �������*/
	double a_solar_radiation[3];    /*��������� ������������� ���������
									 �������� ���������� ���������*/
	double a_atmosphere[3];        /*��������� ������������� ���������
									 ���������*/
	double a_traction[3];          /*��������� ������������� ���������
									 ������������ ���������*/

	/*������� ������� �����������:*/
	double df_central_field[3][3];     /*������� ������� ����������� �������
										�������������� ���������� ������������
										���� �� ������� ���������*/
	double df_off_central_field[3][3]; /*������� ������� ����������� �������
										�������������� ���������� ��������������
										���� ������������ ���� �� �������
										���������*/
	double df_celestial_bodies[3][3]; /*������� ������� ����������� �������
										�������������� ���������� �������� ���
										�� ������� ���������*/
	double df_solar_radiation[3][3];	  /*������� ������� ����������� �������
										����������, ������������� ���������
										���������� ���������*/


	/* [���������� ������ ������] */

	/*���������� ������ ������ ��*/
	void rightPart(VECTOR &rv);

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




	/*ABM8*/
	VECTOR RV[8];

	VECTOR rv_prog;
	VECTOR rv_corr;
	VECTOR rv_time;
	/*��������� ����������� ���������� ������� ����������������� �������
	  ��� �������������� ������� ������-��������-������� 8-�� ������� */
	void Iteracii_ABM8();

	/*��������� ����������� �������� � ��������� �� ������
	  ������-��������-������� 8-�� �������*/
	void ProgKor ();
	/*��������� ��������� ���� �������������� ��� ������
	  ������-��������-������� 8-�� ������� */
	void Decrease_dt_ABM8 ();
	/*��������� ���������� ���� �������������� ��� ������
	  ������-��������-������� 8-�� ������� */
	void Increase_dt_ABM8 ();
	/*��������� ������������� �����-�������-����*/
	void Extrapolation(VECTOR &rv0);





	bool stepCalculation();

};



};




#endif
