//---------------------------------------------------------------------------


#pragma hdrstop

#include "integration.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)


/*запись НУ*/
void chi::integration::set_NU(double r[3], double v[3], double t){

for(int i=0; i<3; i++){
	r_nu[i]=r[i];
	v_nu[i]=v[i];
}
t_nu=t;
}
void chi::integration::set_NU(double rv[6],double t){
for(int i=0; i<3; i++){
	r_nu[i]=rv[i];
	v_nu[i]=rv[i+3];
}
t_nu=t;
}