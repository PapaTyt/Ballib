//---------------------------------------------------------------------------

#ifndef time_convertH
#define time_convertH
//---------------------------------------------------------------------------
#include <DateUtils.hpp>
#include "global.h"

//������� UTC � ������������ �����
double UTCtoTDB(double JD);
//������� UTC � ������������ �����
TDateTime UTCtoTDB(TDateTime TDT);
//������� ������������� ������� � UTC
double TDBtoUTC(double JD);
//������� ������������� ������� � UTC
TDateTime TDBtoUTC(TDateTime TDT);
//��������� ���������� �������� ���������� �������
double GMST(const double JD);
AnsiString TDTToStr(TDateTime DTM, int p);
AnsiString JDToStr(double JD, int p);
#endif
