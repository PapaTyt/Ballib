//---------------------------------------------------------------------------

#ifndef L2H
#define L2H
#include "integration.h"
//---------------------------------------------------------------------------
class correctionL2{
	public:
	void setNU(double r[3], double v[3], double t);
	void corrl2();
	void corrl2_1();
	//private:
	protected:
	double t0;
	double r0[3];
	double v0[3];
	std::vector<chi::Vect> RV_L2[7];        //трасса траектории в L2
	std::vector<chi::Vect> RV_J2[7];        //трасса траектории в J2
	std::vector<chi::Vect> RV_M87[7];		//трасса участка возможного проведения наблюдения М87
	std::vector<chi::Vect> RV_SGR_A[7];    //трасса участка возможного проведения наблюдения SGR-A


    /*очистка файлов для записи результатов*/
	void clearfile();
    /*печать текущего прострела*/
	void print(int L);
    /*печать полного файла расчета*/
	void print();

	/*вычисление трассы наблюдений и запись трас в соответствующие переменные*/
	void basa(chi::integration O, int l);
	/*процедура определение номера итераци достижения максимума времени 1/2/3 - M87/SGR-A/L2*/
	void maxt(int type, int &L);


};








#endif
