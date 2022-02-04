//---------------------------------------------------------------------------

#ifndef mathematicH
#define mathematicH
//---------------------------------------------------------------------------
//функция вычисления модуля вектора задаваемого тремя координатами
double norm(double x, double y, double z);
//функция вычисления модуля вектора
double norm(double r[5]);
//функция вычисления модуля разности векторов r1 и r2
double norm(double r1[3], double r2[3]);
//Скалярное произведение векторов х1 и х2
double Skalyar(const double x1[3], const double x2[3]);
//угол между векторами
double angle_between_vectors(double r1[3], double r2[3]);
//Векторное произведение векторов Vec1 и Vec2
void VectProizv (const double Vec1[3], const double Vec2[3], double proizv[3]);
//Поворот вектора dV1 вокруг оси задоваемой ортами Ort на угол хи в радианах
void Povorot_Vektora(const double Ort[], const double dV1[], const double hi,double dV2[]);

#endif
