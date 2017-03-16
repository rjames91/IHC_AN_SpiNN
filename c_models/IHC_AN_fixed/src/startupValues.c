#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "IHC_Utils.h"
#include "IHC_AN.h"
#include "startupValues.h"

/*everything here performed at startup using floating point representation
  this is then converted into fixed point integers for later processing*/

startupVars generateStartupVars(void)
{
	startupVars out;
	float Gu0,kt0,IHCVnow,gmax,u0,u1,recips1,recips0,ga,gk,et,ek,rpc,ekp,gamma,beta,mICaCurr,
	gmaxcalsr,gmaxcamsr,gmaxcahsr,eca,tauca;

	/******Calculate Gu0************/
	gmax=6e-9;//(float)Cilia.Gmax/pow(2,30);
	u0=5e-9;//(float)Cilia.u0/pow(2,30);
	u1=1e-9;//(float)Cilia.u1/pow(2,30);
	recips1=1/1e-9;//(float)Cilia.recips1/pow(2,1);
	recips0=1/3e-8;//(float)Cilia.recips0/pow(2,1);
	ga=8e-10;//(float)Cilia.Ga/pow(2,31);
	Gu0=ga+gmax/((1+exp(u0*recips0))*(1+exp(u1*recips1)));

	/******Calculate IHCVnow********/
	gk=2e-8;//(float)Cilia.Gk/pow(2,30);
	et=0.1;//(float)Cilia.Et/pow(2,30);
	ek=-0.08;//(float)Cilia.Ek/pow(2,30);
	rpc=0.04;//(float)Cilia.Rpc/pow(2,30);
	ekp=ek+rpc*et;
	IHCVnow=(gk*ekp*Gu0*et)/(Gu0+gk);

	/******Calculate mICaCurrent****/
	gamma=130;
	beta=400;
	mICaCurr=1/(1+exp(-gamma*IHCVnow)*(1/beta));

	/*****Calculate CaCurrent*******/
	gmaxcalsr=1e-11;
	gmaxcamsr=7e-11;
	gmaxcahsr=5e-10;
	eca=0.066;
	tauca=1e-4;

	out.IHCVnow0=(long)(IHCVnow*pow(2,FRACBITS));
	out.Ekp=(long)(ekp*pow(2,FRACBITS));
	out.mICaCurr0=(long)(mICaCurr*pow(2,FRACBITS));
	out.CaCurrLSR0=(long)gmaxcalsr*pow(out.mICaCurr0,3)*(out.IHCVnow0-eca)*tauca;
	out.CaCurrMSR0=(long)gmaxcamsr*pow(out.mICaCurr0,3)*(out.IHCVnow0-eca)*tauca;
	out.CaCurrHSR0=(long)gmaxcahsr*pow(out.mICaCurr0,3)*(out.IHCVnow0-eca)*tauca;

	return out;
}
