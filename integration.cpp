//---------------------------------------------------------------------------


#pragma hdrstop

#include "integration.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)


/*запись НУ*/
void chi::integration::set_NU(double r[3], double v[3], double t){
	/*
	 *Параметры:
	 * 	r[3]	- радиус вектор КА [км]
	 *	v[3] 	- вектор скорости КА [км/с]
	 *	t	 	- момент времени на который приведен вектор состояния в формате
				  Юлианской даты (UTC)
	 */
for(int i=0; i<3; i++){
	r_nu[i]=r[i];
	v_nu[i]=v[i];
}
t_nu=t;
}
void chi::integration::set_NU(double rv[6],double t){
	/*
	 *Параметры:
	 *	кv[6] 	- вектор состояния КА [км, км/с]
	 *	t	 	- момент времени на который приведен вектор состояния в формате
				  Юлианской даты (UTC)
	 */
for(int i=0; i<3; i++){
	r_nu[i]=rv[i];
	v_nu[i]=rv[i+3];
}
t_nu=t;
}
/*запись параметров интегрирования*/
void chi::integration::setParametrs(double interval_, double step_){
	/*
	 *Параметры:
	 * 	interval_ 	- интервал прогноза в днях
	 *	step_   	- шаг выдачи информации в секундах
	 Примечание:
		- Если не требуется вывод промежуточных результатов то следует задавать
		  шаг равный интервалу
	 */
interval=interval_;
step=step_;
}

























	/* [ВЫЧИСЛЕНИЕ ПРАВЫХ ЧАСТЕЙ] */

/*вычисление правых частей ДУ*/
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

	/*вычисление возмущающего ускорения, обусловленного центральным
	  гавитационным полем Земли и матрицы частных производных этого вектора*/
	void central_field(VECTOR rv);

	/*вычисление возмущающего ускорения, обусловленного центральным
	  гавитационным полем Луны и матрицы частных производных этого вектора*/
	void central_field_moon(VECTOR rv);

	/*вычисление возмущающего ускорения обусловленного нецентральностью
	  гравитационного поля Земли и матрицы частных производных этого вектора*/
	void off_central_field(VECTOR rv);

	/*вычисление возмущающего ускорения обусловленного нецентральностью
	  гравитационного поля Луны и матрицы частных производных этого вектора*/
	void off_central_field_moon(VECTOR rv);


	/*вычисление возмущающего ускорения обусловленного нецентральностью
	  гравитационного поля Земли(второй зональной гармоникой - С20)*/
	void off_central_field_C20(VECTOR rv);

	/*вычисление возмущающего ускорения обусловленного нецентральностью
	  гравитационного поля Земли(второй зональной гармоникой - С40)*/
	void off_central_field_C40(VECTOR rv);


	/*вычисление возмущающего ускорения обусловленного нецентральностью
	  гравитационного поля Земли с учетом гармоник до 32х32 и матрицы частных
	  производных этого вектора*/
	void off_central_field_32(VECTOR rv, double df[3]);

	/*вычисление возмущающего ускорения обусловленного нецентральностью
	  гравитационного поля Луны с учетом гармоник до 75х75 и матрицы частных
	  производных этого вектора*/
	void off_central_field_75_moon(VECTOR rv, double df[3]);

	/*вычисление ускорения обусловленных действием небесных тел
	  (реализация для 10 небесных тел, центральное тело Земля) и матрицы частных
	  производных этого вектора*/
	void celestial_bodies(VECTOR rv);

    /*вычисление ускорения обусловленных действием небесных тел
	  (реализация для 10 небесных тел, центральное тело Луна) и матрицы частных
	  производных этого вектора*/
	void celestial_bodies_moon(VECTOR rv);

	/*вычисление ускорения обусловленного действием солнечного излучения и
	  матрицы частных производных этого вектора*/
	void solar_radiation(VECTOR rv);

	/*вычисление ускорения обусловленного воздействием силы сопротивления
	  атмосферы в соответствии с ГОСТ Р 25645.166-2004*/
	void atmosphere(VECTOR rv);

	/*вычисление ускорения обусловленного воздействием силы сопротивления
	  атмосферы в соответствии с ГОСТ Р 25645.166-2004*/
	void atmosphereGOST2004(VECTOR rv);

	/*вычисление ускорения обусловленного работой двигательной установки*/
	void traction(VECTOR rv);
