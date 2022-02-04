//---------------------------------------------------------------------------

#include <vcl.h>
#include <stdio.h>
#include <DateUtils.hpp>
#pragma hdrstop

#include "MainUnit.h"
#include "integration.h"
#include "time_convert.h"




#include "L2.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TForm1::BitBtn1Click(TObject *Sender)
{
//TStringList *list=new TStringList;
//AnsiString str, name;
//name="C:\\_D\\NPOL\\Ballib\\Debug\\matr.dat";
////name="C:\\_D\\NPOL\\Ballib\\Debug\\rv-grin.dat";
////name="C:\\_D\\NPOL\\Ballib\\Debug\\rv-J2000.dat";
//list->LoadFromFile(name);
//double M[3][3];
//TDateTime TTT;
//FILE *ff;
//ff=fopen("matr.txt", "w");
//Word yy, mm, dd, hh, nn, ss, ms;
//double MM[8];
//
//for(int i=0; i<list->Count-1; i++){
//	str=list->Strings[i];
//	strtok( str.c_str(), " ");
//	strtok( NULL, " ");
//	yy=atoi(strtok( NULL, " "));
//	mm=atoi(strtok( NULL, " "));
//	dd=atoi(strtok( NULL, " "));
//	hh=atoi(strtok( NULL, " "));
//	nn=atoi(strtok( NULL, " "));
//	ss=atoi(strtok( NULL, "."));
//	ms=atoi(strtok( NULL, " "));
//	TTT=EncodeDateTime(yy, mm, dd, hh, nn, 0, 0)+(ss+ms/100.)/86400.;
//	for(int k=0; k<8; k++)
//	MM[k]=atof(strtok( NULL, " "));
//
//	DecodeDateTime(TTT, yy, mm, dd, hh, nn, ss, ms);
//	fprintf(ff, "%4i.%02i.%02i %02i:%02i:%02i.%03i ", yy, mm, dd, hh, nn, ss, ms);
//	for(int k=0; k<8; k++){
//		fprintf(ff, "%15.11f ", MM[k]);
//	}
//		fprintf(ff, "\n");
//}
//fclose(ff);
//
//
//
//
//return;
//for(int i=0; i<list->Count; i++){
//	str=list->Strings[i];
//	yy=atoi(strtok( str.c_str(), " "));
//	mm=atoi(strtok( NULL, " "));
//	dd=atoi(strtok( NULL, " "));
//	hh=atoi(strtok( NULL, " "));
//	nn=atoi(strtok( NULL, " "));
//	ss=atoi(strtok( NULL, "."));
//	ms=atoi(strtok( NULL, " "));
//	TTT=EncodeDateTime(yy, mm, dd, hh, nn, 0, 0)+(ss+ms/100.)/86400.;
//
//
//	i++;
//	str=list->Strings[i];
//	M[0][0]=atof(strtok( str.c_str(), " "));
//	M[0][1]=atof(strtok( NULL, " "));
//	M[0][2]=atof(strtok( NULL, " "));
//	i++;
//	str=list->Strings[i];
//	M[1][0]=atof(strtok( str.c_str(), " "));
//	M[1][1]=atof(strtok( NULL, " "));
//	M[1][2]=atof(strtok( NULL, " "));
//	i++;
//	str=list->Strings[i];
//	M[2][0]=atof(strtok( str.c_str(), " "));
//	M[2][1]=atof(strtok( NULL, " "));
//	M[2][2]=atof(strtok( NULL, " "));
//	DecodeDateTime(TTT, yy, mm, dd, hh, nn, ss, ms);
//	fprintf(ff, "%4i.%02i.%02i %02i:%02i:%02i.%03i\n", yy, mm, dd, hh, nn, ss, ms);
//	for(int k=0; k<3; k++){
//		fprintf(ff, "%15.11f %15.11f %15.11f\n", M[k][0], M[k][1], M[k][2]);
//	}
//}
//fclose(ff);
//
//return;
//


//double  t0=DateTimeToJulianDate(EncodeDateTime(2022, 12, 22, 5, 35, 0, 0)),
//		r0[3]={-6894.61065462528, 159.43744304995, 940.85721424237},
//		v0[3]={1.03809306208,   0.99332852891,   7.43011822517};
//
////FILE *ff;
////ff=fopen("1.txt", "w");
////fclose(ff);
//chi::integration O;
//O.set_NU(r0, v0, t0);
//O.setParametrs(900/86400., 0.01);
//O.setParametrs();
//O.setTypeCalculation(0);
//O.ABM8();



ShowMessage("");
}
//---------------------------------------------------------------------------

void __fastcall TForm1::BitBtn2Click(TObject *Sender)
{

double  t0=DateTimeToJulianDate(EncodeDateTime(2030, 2, 1, 0, 0, 0, 0)),
		r0[3]={-659183.3490256558,  1514125.5247464842,   614206.3570772199},
		v0[3]={0.0253602582,       -0.0577662914,        0.1258380499};
correctoinL2 spm;

spm.setNU(r0, v0, t0);
spm.corrl2();

ShowMessage("");

}
//---------------------------------------------------------------------------

