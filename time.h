//---------------------------------------------------------------------------

#ifndef timeH
#define timeH
#include <DateUtils.hpp>
//������� UTC � ������������ �����
double UTCtoTDB(double JD);
//������� UTC � ������������ �����
TDateTime UTCtoTDB(TDateTime TDT);
//������� ������������� ������� � UTC
double TDBtoUTC(double JD);
//������� ������������� ������� � UTC
TDateTime TDBtoUTC(TDateTime TDT);

AnsiString TDTToStr(TDateTime DTM, int p);
AnsiString JDToStr(double JD, int p);


#endif
