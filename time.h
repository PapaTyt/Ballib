//---------------------------------------------------------------------------

#ifndef timeH
#define timeH
#include <DateUtils.hpp>
//Перевод UTC в динамическое время
double UTCtoTDB(double JD);
//Перевод UTC в динамическое время
TDateTime UTCtoTDB(TDateTime TDT);
//Перевод динамического времени в UTC
double TDBtoUTC(double JD);
//Перевод динамического времени в UTC
TDateTime TDBtoUTC(TDateTime TDT);

AnsiString TDTToStr(TDateTime DTM, int p);
AnsiString JDToStr(double JD, int p);


#endif
