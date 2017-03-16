/*
 ============================================================================
 Name        : IHC_AN.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "IHC_AN.h"
#include "IHC_Utils.h"
//#include "stdfix-exp.h"
#include "startupValues.h"
#include "createCSV.h"


//N.B. This function is a single IHC model, therefore its
// input will be from a single BM channel model output.
//Using a QF=FRACBITS fixed point format
/*Global Variables*/
bool overflowFlag=false;
long Fs;
float dt;

float testInput[segSize];

//fixedPointVars fixedInput[segSize];
long fixedInput[segSize];

int main(void) {
float fixer=pow(2,FRACBITS);//pow(2,40);

//read input cilia displacement wave
readCSV(testInput,"./CSVs/IHCcilia0dB.csv",segSize,1);

//initialise some global parameters
Fs=44100;
dt=(1.0/Fs);

//convert input floating point cilia displacement into 64bit fixed point
int i;

for(i=0; i<segSize; i++)
{
	long fixed=(long)(testInput[i]*fixer);//pow(2,62));
	fixedInput[i]=fixed;
	float inconv=(float)fixedInput[i]/fixer;
}

//printf("dt=%f\n",(float)dt/fixer);

//initialise the cilia params

IHC_ciliaParams Cilia;

Cilia.tc= 0.00012*fixer;
Cilia.C=0.004*fixer;
Cilia.u0=5e-9*fixer;
Cilia.recips0=(1/3e-8)*pow(2,33);//Q30.33
Cilia.u1=1e-9*fixer;
Cilia.recips1=(1/1e-9)*pow(2,33);//Q30.33
Cilia.Gmax=6e-9*fixer;
Cilia.Ga=8e-10*fixer;
Cilia.dtCap=(dt/4e-12)*pow(2,33);//Q30.33
Cilia.Et=0.1*fixer;
Cilia.Gk=2e-8*fixer;
Cilia.Ek=-0.08*fixer;
Cilia.Rpc=0.04*fixer;

/*********Generate Startup Values************/
startupVars startupValues=generateStartupVars();
long IHCVnow0=-0.060768628470039*fixer;//0dB1kHz-0.067251540515944*fixer;//60dB1kHz-0.061034432519758*fixer;//140dB1kHz-0.063355671303642*fixer;//100dB1kHz-0.061517772780083*pow(2,62);//20dB1kHz
//long IHCVnow0=startupValues.IHCVnow0;
long mICaCurr0=0.127411875118889*fixer;//0dB1kHz 0.127027616182966*fixer;//60dB1kHz
//long mICaCurr0=startupValues.mICaCurr0;
long CaCurrLSR0=2.584290738035695e-15*fixer;//0dB1kHz 1.506850345912777e-15*fixer;//60dB1kHz
//long CaCurrLSR0=startupValues.CaCurrLSR0;
long CaCurrMSR0=1.809003516624987e-14*fixer;//0dB1kHz1.054795242138944e-14*fixer;//60dB1kHz
//long CaCurrMSR0=startupValues.CaCurrMSR0;
long CaCurrHSR0=1.292145369017847e-13*fixer;//0dB1kHz7.534251729563886e-14*fixer;//60dB1kHz
//long CaCurrHSR0=startupValues.CaCurrHSR0;

long IHCVnow=IHCVnow0;
long mICaCurr=mICaCurr0;
long CaCurrLSR[NUMLSR]={CaCurrLSR0,CaCurrLSR0};
long CaCurrMSR[NUMMSR]={CaCurrMSR0,CaCurrMSR0};
long CaCurrHSR[NUMHSR]={CaCurrHSR0,CaCurrHSR0};

/*****initialise the pre synapse params**********/
IHC_preSynParams preSyn;

preSyn.recipBeta=2.5e-3*fixer;
preSyn.gamma=130*fixer;
preSyn.dtTauM=(dt/1e-4)*fixer;
preSyn.dtTauCa=(dt/1e-4)*fixer;
preSyn.Eca=0.066*fixer;

preSyn.GmaxCa[0]=1e-11*fixer;
preSyn.GmaxCa[1]=1e-11*fixer;
preSyn.GmaxCa[2]=7e-11*fixer;
preSyn.GmaxCa[3]=7e-11*fixer;
preSyn.GmaxCa[4]=5e-10*fixer;
preSyn.GmaxCa[5]=5e-10*fixer;
preSyn.GmaxCa[6]=5e-10*fixer;
preSyn.GmaxCa[7]=5e-10*fixer;
preSyn.GmaxCa[8]=5e-10*fixer;
preSyn.GmaxCa[9]=5e-10*fixer;

//Run through samples in segments
long ut[segSize];
long Gu[segSize];
long IHC_RP[segSize];
long mICa[segSize];
long ICaLSR[NUMLSR][segSize];
long ICaMSR[NUMMSR][segSize];
long ICaHSR[NUMHSR][segSize];

float ihcconv[segSize];
float Guconv[segSize];
float utconv[segSize];
float micainfconv[segSize];
float micaconv[segSize];
float calsrconv[NUMLSR][segSize];
float camsrconv[NUMMSR][segSize];
float cahsrconv[NUMHSR][segSize];

for(i=0;i<segSize;i++)
{
	/*if(i==29)
	{
		printf("i=%i\n",i);
	}*/
	/**********Apply Scaler****************/
	ut[i]=longMultiply2(fixedInput[i],Cilia.C);
	utconv[i]=(float)ut[i]/pow(2,FRACBITS);

	/********Apical Conductance************/
	//exponentials
 	long ex0=lApicalExp(ut[i],Cilia.u0,Cilia.recips0);
	float exoconv=(float)ex0/pow(2,FRACBITS);

	long ex1=lApicalExp(ut[i],Cilia.u1,Cilia.recips1);
	float ex1conv=(float)ex1/pow(2,FRACBITS);
	long lone=(long)1<<FRACBITS;
	long div;
	long denom;
	denom=longMultiply2(ex0,(ex1+lone));
	//check for maximum overflow
	if (denom==MAXSIGNEDL || denom==MINSIGNEDL)
	{
		div=0;
	}

	else
	{
		denom=denom+lone;
		float denomconv=(float)denom/pow(2,FRACBITS);
		//TODO:precision is lost here, is there another implementation? Reciprocal approximations
		long recip= ((long)1<<62)/denom;
		float recipconv=(float)recip/pow(2,(62-FRACBITS));
		div=longMultiply2(recip<<(FRACBITS-(62-FRACBITS)),Cilia.Gmax);
	}

	float divconv=(float)div/fixer;
	Gu[i]=Cilia.Ga+div;
	Guconv[i]=(float)Gu[i]/pow(2,FRACBITS);

	/*********Receptor Potential*********/
	float IHCVnowconv= (float)IHCVnow/pow(2,FRACBITS);
	long subt=IHCVnow-startupValues.Ekp;
	float subtconv=(float)subt/pow(2,FRACBITS);
	long b=longMultiply2(Cilia.Gk,subt);

	long subt2=IHCVnow-Cilia.Et;
	float subt2conv=(float)subt2/pow(2,FRACBITS);
	long a=longMultiply2(subt2,(-1*Gu[i]));

	long aminusb=(a-b);
	float aminusbconv=(float)aminusb/pow(2,FRACBITS);
	long mult=longMultiply2(aminusb,Cilia.dtCap);
	mult=mult<<FRACBITS-33;
	float multconv=(float)mult/pow(2,FRACBITS);

	//TODO:add overflow protection here?
	IHCVnow+=mult;
	IHC_RP[i]=IHCVnow;
	ihcconv[i]=(float)IHCVnow/pow(2,FRACBITS);

	/***************mICa*****************/
	mult=longMultiply2(((long)-1*preSyn.gamma),IHC_RP[i]);
	multconv=(float)mult/pow(2,FRACBITS);
	long ex3=lExp(mult);

	mult=longMultiply2(ex3,preSyn.recipBeta);
	//add 1
	mult+=(long)1<<FRACBITS;
	multconv=(float)mult/pow(2,FRACBITS);

	long mICaINF=((long)1<<62)/mult;//TODO:verify format here
	mICaINF=mICaINF<<(FRACBITS-(62-FRACBITS));
	micainfconv[i]=(float)mICaINF/pow(2,FRACBITS);

	mult=longMultiply2((mICaINF-mICaCurr),preSyn.dtTauM);
	multconv=(float)mult/pow(2,FRACBITS);
	mICaCurr=mICaCurr+mult;

	mICa[i]=mICaCurr;
	micaconv[i]=mICaCurr/pow(2,FRACBITS);

	/***************ICa******************/

	long LSRfibres[NUMLSR][segSize];
	long MSRfibres[NUMMSR][segSize];
	long HSRfibres[NUMHSR][segSize];

	if(i==22)
	{
		printf("11thsample\n");
	}

	//copy values for each fibre
	int j;
	long ICa;
	long micapow=lpow(mICa[i],3);
	float micapowconv=(float)micapow/fixer;

	/*********LSR Fibres************/
	for (j=0;j<NUMLSR;j++)
	{
		mult=longMultiply2(micapow,preSyn.GmaxCa[j]);
		multconv=(float)mult/pow(2,FRACBITS);

		ICa=longMultiply2(mult,(IHC_RP[i]-preSyn.Eca));//check this doesn't always equal 0
		float ICaconv=(float)ICa/pow(2,FRACBITS);
		mult=longMultiply2(((-1*ICa)-CaCurrLSR[j]),preSyn.dtTauCa);

		CaCurrLSR[j]+=mult;
		ICaLSR[j][i]=CaCurrLSR[j];
		calsrconv[j][i]=(float)ICaLSR[j][i]/fixer;
		//LSRfibres[j][i]=IHC_RP[i];
	}
	/********MSR Fibres*************/
	for (j=0;j<NUMMSR;j++)
	{
		mult=longMultiply2(micapow,preSyn.GmaxCa[j+2]);
		multconv=(float)mult/pow(2,FRACBITS);

		ICa=longMultiply2(mult,(IHC_RP[i]-preSyn.Eca));//check this doesn't always equal 0
		float ICaconv=(float)ICa/pow(2,FRACBITS);
		mult=longMultiply2(((-1*ICa)-CaCurrMSR[j]),preSyn.dtTauCa);

		CaCurrMSR[j]+=mult;
		ICaMSR[j][i]=CaCurrMSR[j];
		camsrconv[j][i]=(float)ICaMSR[j][i]/fixer;
		//MSRfibres[j][i]=IHC_RP[i];
	}
	/********HSR Fibres************/
	for (j=0;j<NUMHSR;j++)
	{
		mult=longMultiply2(micapow,preSyn.GmaxCa[j+4]);
		multconv=(float)mult/pow(2,FRACBITS);

		ICa=longMultiply2(mult,(IHC_RP[i]-preSyn.Eca));//check this doesn't always equal 0

		float ICaconv=(float)ICa/pow(2,FRACBITS);
		mult=longMultiply2(((-1*ICa)-CaCurrHSR[j]),preSyn.dtTauCa);

		CaCurrHSR[j]+=mult;
		ICaHSR[j][i]=CaCurrHSR[j];
		cahsrconv[j][i]=(float)ICaHSR[j][i]/fixer;
		//HSRfibres[j][i]=IHC_RP[i];
	}

	//printf("OK\n");
}

	printf("Processed %i samples\n",segSize);
	createCSV("./CSVs/IHC_UT.csv",utconv,segSize,1);
	createCSV("./CSVs/IHC_AC.csv",Guconv,segSize,1);
	createCSV("./CSVs/IHC_RP.csv",ihcconv,segSize,1);
	createCSV("./CSVs/IHC_MICAINF.csv",micainfconv,segSize,1);
	createCSV("./CSVs/IHC_MICA.csv",micaconv,segSize,1);
	createCSV("./CSVs/IHC_LSR1.csv",calsrconv,segSize,1);
	createCSV("./CSVs/IHC_MSR1.csv",camsrconv,segSize,1);
	createCSV("./CSVs/IHC_HSR1.csv",cahsrconv,segSize,1);

	return EXIT_SUCCESS;
}

