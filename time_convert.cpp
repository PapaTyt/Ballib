//---------------------------------------------------------------------------


#pragma hdrstop

#include "time_convert.h"
#include <stdio.h>
#include <math.h>
//---------------------------------------------------------------------------

#pragma package(smart_init)
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------

#pragma package(smart_init)
//---------------------------------------------------------------------------
/*Перевод UTC в динамическое время
Входные параметры:
		JD - время в юлиянских сутках
Выходные параметры:
		динамическое время в юлианских сутках
Подключаемые модули:
		"constants.h"
*/
 double UTCtoTDB(double JD){
return JD+(dUTC+dUT1)/86400.;
}
//---------------------------------------------------------------------------
/*Перевод UTC в динамическое время
Входные параметры:
		JD - время в юлиянских сутках
Выходные параметры:
		динамическое время в юлианских сутках
Подключаемые модули:
		"constants.h"
*/
TDateTime UTCtoTDB(TDateTime TDT){
return TDT+(dUTC+dUT1)/86400.;
}
//---------------------------------------------------------------------------
/*Перевод динамического времени в UTC
Входные параметры:
		JD - динамическое время в юлианских сутках
Выходные параметры:
		время в юлианских сутках
Подключаемые модули:
		"constants.h"
*/
double TDBtoUTC(double JD){
return JD-(dUTC+dUT1)/86400.;
}
//---------------------------------------------------------------------------
/*Перевод динамического времени в UTC
Входные параметры:
		JD - динамическое время в юлианских сутках
Выходные параметры:
		время в юлианских сутках
Подключаемые модули:
		"constants.h"
*/
TDateTime TDBtoUTC(TDateTime TDT){
return TDT-(dUTC+dUT1)/86400.;
}
//---------------------------------------------------------------------------
/*Процедура вычисления среднего солнечного времени
Входные параметры:
		JD - время в юлиянских сутках
Выходные параметры:
		Звездное время в радианах
Используемые процедуры:
		modf
		pow
Подключаемые модули:
		<math.h>
		"constants.h"
*/
double GMST(const double JD)
{
   double st,JD0,Tut,n;
   modf(JD - 0.5,&n);
   JD0 = n + 0.5;
   Tut = (JD0 - 2451545.)/36525.;
   st = 1.7533685592 + 628.3319706889*Tut + 6.7707139*pow(10.0,-6)*Tut*Tut
   - 4.50876*pow(10.0,-10)*pow(Tut,3) + 1.002737909350795*2*M_PI*(JD - JD0);
   modf(st/(2*M_PI),&n);
   return st-n*2*M_PI; // в приведённых радианах от 0 до 2*PI
}

//---------------------------------------------------------------------------
AnsiString TDTToStr(TDateTime DTM, int p){
Word 	yy,  	//год
		mm,  	//месяц
		dd,  	//день
		hh, 	//час
		nn, 	//минута
		ss,  	//секунда
		zzz; 	//милисекунда
char ch[25];
AnsiString s="";
DecodeDateTime(DTM, yy, mm, dd, hh, nn, ss, zzz);

switch(p){

	case 0:  sprintf(ch,"%4i.%02i.%02i %02i:%02i:%02i.%03i", yy, mm, dd, hh, nn, ss, zzz);
		break;
	case 1:  sprintf(ch,"%4i.%02i.%02i %02i:%02i:%02i", yy, mm, dd, hh, nn, ss);
		break;
	case 2:  sprintf(ch,"%4i %02i %02i %02i %02i %02i %03i", yy, mm, dd, hh, nn, ss, zzz);
		break;
	case 3:  sprintf(ch,"%4i %02i %02i %02i %02i %02i", yy, mm, dd, hh, nn, ss);
		break;
	case 4:  sprintf(ch,"%4i.%02i.%02i", yy, mm, dd);
		break;
	case 5:  sprintf(ch,"%4i %02i %02i", yy, mm, dd);
		break;
	case 6:  sprintf(ch,"%02i:%02i:%02i.%03i", hh, nn, ss, zzz);
		break;
	case 7:  sprintf(ch,"%02i:%02i:%02i", hh, nn, ss);
		break;
	case 8:  sprintf(ch,"%02i %02i %02i %03i", hh, nn, ss, zzz);
		break;
	case 9:  sprintf(ch,"%02i %02i %02i", hh, nn, ss);
		break;
	case 10:  sprintf(ch,"%4i%02i%02i_%02i%02i%02i", yy, mm, dd, hh, nn, ss);
		break;

	case 11:  sprintf(ch,"%02i.%02i.%04i %02i:%02i:%02i", dd, mm, yy, hh, nn, ss);
		break;
	case 12:  sprintf(ch,"%02i.%02i.%04i %02i:%02i:%02i.%03i", dd, mm, yy, hh, nn, ss, zzz);
		break;
	case 13:  sprintf(ch,"%02i.%02i.%04i", dd, mm, yy);
		break;
}

s=ch;
return s;
}
//---------------------------------------------------------------------------
AnsiString JDToStr(double JD, int p) {
return TDTToStr(JulianDateToDateTime(JD), p);
}
