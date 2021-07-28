//---------------------------------------------------------------------------


#pragma hdrstop

#include "time.h"
//#include "constants.h"
#include "global.h"
#include <stdio.h>
#include <math.h>

//---------------------------------------------------------------------------

#pragma package(smart_init)
//---------------------------------------------------------------------------
/*������� UTC � ������������ �����
������� ���������:
		JD - ����� � ��������� ������
�������� ���������:
		������������ ����� � ��������� ������
������������ ������:
		"constants.h"
*/
 double UTCtoTDB(double JD){
return JD+(dUTC+dUT1)/86400.;
}
//---------------------------------------------------------------------------
/*������� UTC � ������������ �����
������� ���������:
		JD - ����� � ��������� ������
�������� ���������:
		������������ ����� � ��������� ������
������������ ������:
		"constants.h"
*/
TDateTime UTCtoTDB(TDateTime TDT){
return TDT+(dUTC+dUT1)/86400.;
}
//---------------------------------------------------------------------------
/*������� ������������� ������� � UTC
������� ���������:
		JD - ������������ ����� � ��������� ������
�������� ���������:
		����� � ��������� ������
������������ ������:
		"constants.h"
*/
double TDBtoUTC(double JD){
return JD-(dUTC+dUT1)/86400.;
}
//---------------------------------------------------------------------------
/*������� ������������� ������� � UTC
������� ���������:
		JD - ������������ ����� � ��������� ������
�������� ���������:
		����� � ��������� ������
������������ ������:
		"constants.h"
*/
TDateTime TDBtoUTC(TDateTime TDT){
return TDT-(dUTC+dUT1)/86400.;
}

//---------------------------------------------------------------------------
AnsiString TDTToStr(TDateTime DTM, int p){
Word 	yy,  	//���
		mm,  	//�����
		dd,  	//����
		hh, 	//���
		nn, 	//������
		ss,  	//�������
		zzz; 	//�����������
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
