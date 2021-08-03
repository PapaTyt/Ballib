//---------------------------------------------------------------------------

#ifndef time_convertH
#define time_convertH
//---------------------------------------------------------------------------
#include <DateUtils.hpp>
#include "global.h"

//Перевод UTC в динамическое время
double UTCtoTDB(double JD);
//Перевод UTC в динамическое время
TDateTime UTCtoTDB(TDateTime TDT);
//Перевод динамического времени в UTC
double TDBtoUTC(double JD);
//Перевод динамического времени в UTC
TDateTime TDBtoUTC(TDateTime TDT);
//Процедура вычисления среднего солнечного времени
double GMST(const double JD);
AnsiString TDTToStr(TDateTime DTM, int p);
AnsiString JDToStr(double JD, int p);
#endif
