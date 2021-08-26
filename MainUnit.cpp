//---------------------------------------------------------------------------

#include <vcl.h>
#include <stdio.h>
#include <DateUtils.hpp>
#pragma hdrstop

#include "MainUnit.h"
#include "integration.h"
#include "time_convert.h"

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
double  t0=DateTimeToJulianDate(EncodeDateTime(2029, 12, 30, 18, 8, 31, 713)),
		r0[3]={315638.388482,	1063147.993512,	534324.266362},
		v0[3]={-0.550641351,	-0.012697326,	0.033769707};
chi::integration O;
O.set_NU(r0, v0, t0);
O.setParametrs(1, 1);
O.setParametrs();
O.setTypeCalculation(0);
O.ABM8();

}
//---------------------------------------------------------------------------

