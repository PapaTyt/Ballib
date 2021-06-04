
//---------------------------------------------------------------------------
#include <math.h>

#pragma hdrstop

#include "mathematic.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)
//функция вычисления модуля вектора задаваемого тремя координатами
double norm(double x, double y, double z){
	return sqrt(x*x+y*y+z*z);
}
//функция вычисления модуля вектора
double norm(double r[5]){
	return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}
//функция вычисления модуля разности векторов r1 и r2
double norm(double r1[3], double r2[3]){
return sqrt((r1[0]-r2[0])*(r1[0]-r2[0])+
			(r1[1]-r2[1])*(r1[1]-r2[1])+
			(r1[2]-r2[2])*(r1[2]-r2[2]));
}