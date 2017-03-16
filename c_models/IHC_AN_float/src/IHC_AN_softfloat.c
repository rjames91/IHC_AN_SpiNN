/*
 ============================================================================
 Name        : IHC_AN_softfloat.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdfix.h>
#include "startupValues.h"
#include "IHC_AN_softfloat.h"
#include "spin1_api.h"
#include "math.h"

REAL Fs,dt;
REAL testInput[segSize];
uint coreID;
uint chipID;

IHC_ciliaParams Cilia;
IHC_preSynParams preSyn;
IHC_synParams Synapse;

startupVars generateStartupVars(void)
{
	startupVars out;
	REAL Gu0,kt0,IHCVnow,gmax,u0,u1,s1,s0,ga,gk,et,ek,rpc,ekp,gamma,beta,mICaCurr,
	gmaxcalsr,gmaxcamsr,gmaxcahsr,eca,tauca,CaCurrLSR,CaCurrMSR,CaCurrHSR,kt0LSR,kt0MSR,kt0HSR,
	ANCleftLSR,ANCleftMSR,ANCleftHSR,ANAvailLSR,ANAvailMSR,ANAvailHSR,ANReproLSR,ANReproMSR,
	ANReproHSR,z,power,y,x,l,r,MLSR,MMSR,MHSR;

	/******Calculate Gu0************/
	gmax=6e-9;//(float)Cilia.Gmax/pow(2,30);
	u0=5e-9;//(float)Cilia.u0/pow(2,30);
	u1=1e-9;//(float)Cilia.u1/pow(2,30);
	s1=1e-9;//(float)Cilia.recips1/pow(2,1);
	s0=3e-8;//(float)Cilia.recips0/pow(2,1);
	ga=8e-10;//(float)Cilia.Ga/pow(2,31);
	Gu0=ga+gmax/(1+(exp(u0/s0)*(1+exp(u1/s1))));
	float test=exp(u0/s0);
	float test1=exp(u1/s1);
	/******Calculate IHCVnow********/
	gk=2e-8;//(float)Cilia.Gk/pow(2,30);
	et=0.1;//(float)Cilia.Et/pow(2,30);
	ek=-0.08;//(float)Cilia.Ek/pow(2,30);
	rpc=0.04;//(float)Cilia.Rpc/pow(2,30);
	ekp=ek+rpc*et;
	IHCVnow=(gk*ekp+Gu0*et)/(Gu0+gk);

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
	CaCurrLSR=gmaxcalsr*pow(mICaCurr,3)*(IHCVnow-eca)*tauca;
	CaCurrMSR=gmaxcamsr*pow(mICaCurr,3)*(IHCVnow-eca)*tauca;
	CaCurrHSR=gmaxcahsr*pow(mICaCurr,3)*(IHCVnow-eca)*tauca;

	/*****Calculate AN********/
	z=1e38;//reduced from 1e40 fit as single float type
	power=3;
	y=10;
	x=30;
	MLSR=13;
	MMSR=14;
	MHSR=8;
	l=2580;
	r=6580;

	kt0LSR=z*pow(CaCurrLSR,power)*100;//added *100
	kt0MSR=z*pow(CaCurrMSR,power)*100;
	kt0HSR=z*pow(CaCurrHSR,power)*100;
	ANCleftLSR=(kt0LSR*y*MLSR)/(y*(l+r)+kt0LSR*l);
	ANCleftMSR=(kt0MSR*y*MMSR)/(y*(l+r)+kt0MSR*l);
	ANCleftHSR=(kt0HSR*y*MHSR)/(y*(l+r)+kt0HSR*l);
	ANAvailLSR=round((ANCleftLSR*(l+r))/kt0LSR);
	ANAvailMSR=round((ANCleftMSR*(l+r))/kt0MSR);
	ANAvailHSR=round((ANCleftHSR*(l+r))/kt0HSR);
	ANReproLSR=(ANCleftLSR*r)/x;
	ANReproMSR=(ANCleftMSR*r)/x;
	ANReproHSR=(ANCleftHSR*r)/x;

	out.IHCVnow0=IHCVnow;
	out.Ekp=ekp;
	out.mICaCurr0=mICaCurr;
	out.CaCurrLSR0=CaCurrLSR;
	out.CaCurrMSR0=CaCurrMSR;
	out.CaCurrHSR0=CaCurrHSR;
	out.ANCleftLSR0=ANCleftLSR;
	out.ANCleftMSR0=ANCleftMSR;
	out.ANCleftHSR0=ANCleftHSR;
	out.ANAvailLSR0=ANAvailLSR;
	out.ANAvailMSR0=ANAvailMSR;
	out.ANAvailHSR0=ANAvailHSR;
	out.ANReproLSR0=ANReproLSR;
	out.ANReproMSR0=ANReproMSR;
	out.ANReproHSR0=ANReproHSR;

	return out;
}
void process_chan(void) {

	Fs=50000;//44100;
	dt=(1.0/Fs);

	REAL refrac_lsr[NUMLSR]={0,0};
	REAL refrac_msr[NUMMSR]={0,0};
	REAL refrac_hsr[NUMHSR]={0,0,0,0,0,0};

	//read input cilia displacement wave


	//initialise the cilia params
	Cilia.tc= 0.00012;
	Cilia.C=0.004;
	Cilia.u0=5e-9;
	Cilia.s0=3e-8;
	Cilia.u1=1e-9;
	Cilia.s1=1e-9;
	Cilia.Gmax=6e-9;
	Cilia.Ga=8e-10;
	Cilia.dtCap=(dt/4e-12);
	Cilia.Et=0.1;
	Cilia.Gk=2e-8;
	Cilia.Ek=-0.08;
	Cilia.Rpc=0.04;

	startupVars startupValues=generateStartupVars();
	/********sample 0 starting values****************/
	REAL IHCVnow0=startupValues.IHCVnow0;
	REAL mICaCurr0=startupValues.mICaCurr0;
	REAL CaCurrLSR0=startupValues.CaCurrLSR0;
	REAL CaCurrMSR0=startupValues.CaCurrMSR0;
	REAL CaCurrHSR0=startupValues.CaCurrHSR0;
	REAL ANCleftLSR0=startupValues.ANCleftLSR0;
	REAL ANCleftMSR0=startupValues.ANCleftMSR0;
	REAL ANCleftHSR0=startupValues.ANCleftHSR0;
	REAL ANAvailLSR0=startupValues.ANAvailLSR0;
	REAL ANAvailMSR0=startupValues.ANAvailMSR0;
	REAL ANAvailHSR0=startupValues.ANAvailHSR0;
	REAL ANReproLSR0=startupValues.ANReproLSR0;
	REAL ANReproMSR0=startupValues.ANReproMSR0;
	REAL ANReproHSR0=startupValues.ANReproHSR0;

	/********sample 8821 starting values*************/
	/*REAL IHCVnow0=-0.060768628470039;
	REAL mICaCurr0=0.127411875118889;
	REAL CaCurrLSR0=2.584290738035695e-15;
	REAL CaCurrMSR0=1.809003516624987e-14;
	REAL CaCurrHSR0=1.292145369017847e-13;
	REAL ANCleftLSR0=1.365543366496455e-06;
	REAL ANCleftMSR0=1.365494792966599e-06;
	REAL ANCleftHSR0=9.180174283732440e-59;
	REAL ANAvailLSR0=2;
	REAL ANAvailMSR0=0;
	REAL ANAvailHSR0=0;
	REAL ANReproLSR0=1.327509936105308;
	REAL ANReproMSR0=0.792575438323633;
	REAL ANReproHSR0=-0.124454148469082;*/

	/********set recurring values*******************/
	REAL IHCVnow=IHCVnow0;
	REAL mICaCurr=mICaCurr0;
	REAL CaCurrLSR[NUMLSR]={CaCurrLSR0,CaCurrLSR0};
	REAL CaCurrMSR[NUMMSR]={CaCurrMSR0,CaCurrMSR0};
	REAL CaCurrHSR[NUMHSR]={CaCurrHSR0,CaCurrHSR0,CaCurrHSR0,CaCurrHSR0,CaCurrHSR0,CaCurrHSR0};
	REAL ANCleftLSR[NUMLSR]={ANCleftLSR0,ANCleftLSR0};
	REAL ANAvailLSR[NUMLSR]={ANAvailLSR0,ANAvailLSR0};
	REAL ANReproLSR[NUMLSR]={ANReproLSR0,ANReproLSR0};
	REAL ANCleftMSR[NUMMSR]={ANCleftMSR0,ANCleftMSR0};
	REAL ANAvailMSR[NUMMSR]={ANAvailMSR0,ANAvailMSR0};
	REAL ANReproMSR[NUMMSR]={ANReproMSR0,ANReproMSR0};
	REAL ANCleftHSR[NUMHSR]={ANCleftHSR0,ANCleftHSR0,ANCleftHSR0,ANCleftHSR0,ANCleftHSR0,ANCleftHSR0};
	REAL ANAvailHSR[NUMHSR]={ANAvailHSR0,ANAvailHSR0,ANAvailHSR0,ANAvailHSR0,ANAvailHSR0,ANAvailHSR0};
	REAL ANReproHSR[NUMHSR]={ANReproHSR0,ANReproHSR0,ANReproHSR0,ANReproHSR0,ANReproHSR0,ANReproHSR0};

	/*****initialise the pre synapse params**********/
	preSyn.recipBeta=2.5e-3;
	preSyn.gamma=130;
	preSyn.dtTauM=(dt/1e-4);
	preSyn.dtTauCa=(dt/1e-4);
	preSyn.Eca=0.066;
	preSyn.power=3;
	preSyn.z=1e38;//reduced from 1e40 to fit within float range (then an additional *100 is added in code)

	preSyn.GmaxCa[0]=1e-11;
	preSyn.GmaxCa[1]=1e-11;
	preSyn.GmaxCa[2]=7e-11;
	preSyn.GmaxCa[3]=7e-11;
	preSyn.GmaxCa[4]=5e-10;
	preSyn.GmaxCa[5]=5e-10;
	preSyn.GmaxCa[6]=5e-10;
	preSyn.GmaxCa[7]=5e-10;
	preSyn.GmaxCa[8]=5e-10;
	preSyn.GmaxCa[9]=5e-10;

	preSyn.CaTh[0]=1e-14;
	preSyn.CaTh[1]=1e-14;
	preSyn.CaTh[2]=1e-13;
	preSyn.CaTh[3]=1e-13;
	preSyn.CaTh[4]=1.43e-13;
	preSyn.CaTh[5]=1.43e-13;
	preSyn.CaTh[6]=1.43e-13;
	preSyn.CaTh[7]=1.43e-13;
	preSyn.CaTh[8]=1.43e-13;
	preSyn.CaTh[9]=1.43e-13;

	/******initialise the synapse params**************/

	Synapse.l=2580;
	Synapse.y=10;
	Synapse.x=30;
	Synapse.r=6580;
	Synapse.refrac_period=7.5e-4;

	Synapse.M[0]=13;
	Synapse.M[1]=13;
	Synapse.M[2]=14;
	Synapse.M[3]=14;
	Synapse.M[4]=8;
	Synapse.M[5]=8;
	Synapse.M[6]=8;
	Synapse.M[7]=8;
	Synapse.M[8]=8;
	Synapse.M[9]=8;

	/***********OUTPUTS******************************/
	REAL ihcconv[segSize];
	REAL Guconv[segSize];
	REAL utconv[segSize];
	REAL micainfconv[segSize];
	REAL micaconv[segSize];
	REAL calsr[NUMLSR][segSize];
	REAL camsr[NUMMSR][segSize];
	REAL cahsr[NUMHSR][segSize];
	REAL vrrlsr[NUMHSR][segSize];
	REAL vrrmsr[NUMHSR][segSize];
	REAL vrrhsr[NUMHSR][segSize];
	REAL anlsr[NUMLSR][segSize];
	REAL anmsr[NUMMSR][segSize];
	REAL anhsr[NUMHSR][segSize];

	int i;//i corresponds to each sample in the input
	for(i=0;i<segSize;i++)
	{
		/**********Apply Scaler****************/
		utconv[i]=testInput[i]*Cilia.C;

		/********Apical Conductance************/
		REAL div= Cilia.Gmax/(1+exp(-(utconv[i]-Cilia.u0)/Cilia.s0)*(1+exp(-(utconv[i]-Cilia.u1)/Cilia.s1)));
		Guconv[i]=Cilia.Ga+div;

		/*********Receptor Potential*********/
		IHCVnow+=((-Guconv[i]*(IHCVnow-Cilia.Et))-(Cilia.Gk*(IHCVnow-startupValues.Ekp)))*Cilia.dtCap;
		ihcconv[i]=IHCVnow;

		/***************mICa*****************/
		REAL mICaINF=1/(1+exp(-preSyn.gamma*ihcconv[i])*preSyn.recipBeta);
		mICaCurr+=(mICaINF-mICaCurr)*preSyn.dtTauM;

		micaconv[i]=mICaCurr;

		/***************ICa******************/
		long LSRfibres[NUMLSR][segSize];
		long MSRfibres[NUMMSR][segSize];
		long HSRfibres[NUMHSR][segSize];

		/*if(i==22)
		{
			printf("11thsample\n");
		}*/

		//copy values for each fibre
		int j;
		REAL ICa;
		float micapowconv=pow(micaconv[i],preSyn.power);

		/*********LSR Fibres************/
		for (j=0;j<NUMLSR;j++)
		{
			/*****Synaptic Ca**********/
			ICa=preSyn.GmaxCa[j]*micapowconv*(ihcconv[i]-preSyn.Eca);
			CaCurrLSR[j]+=(-ICa-CaCurrLSR[j])*preSyn.dtTauCa;
			calsr[j][i]=CaCurrLSR[j];

			/*****Vesicle Release Rate*****/
			REAL compare=max(100*(pow(calsr[j][i],preSyn.power)-pow(preSyn.CaTh[j],preSyn.power)),0);
			vrrlsr[j][i]=preSyn.z*compare;
			//saturate vrr
			if(vrrlsr[j][i]>Fs)
			{
				vrrlsr[j][i]=Fs;
			}

			/*****Release Probability*****/
			REAL releaseProb=vrrlsr[j][i]*dt;
			REAL M_q=Synapse.M[j]-ANAvailLSR[j];
			if(M_q<0)
			{
				M_q=0;
			}

			/*********Ejected************/
			REAL Probability=1-pow(1-releaseProb,ANAvailLSR[j]);
			REAL ejected;

			if (refrac_lsr[j]>0)
			{
				//decrement refrac counter until it's zero
				refrac_lsr[j]--;
			}

			if(Probability>(rand()/(REAL)RAND_MAX) && refrac_lsr[j]==0)
			{
				ejected=1;
				//refrac value to count down from is always between refrac period + random dist between 0 and refrac period
				refrac_lsr[j]=round(Synapse.refrac_period/dt + ((rand()/(REAL)RAND_MAX)*Synapse.refrac_period)/dt);
			}
			else
			{
				ejected=0;
			}
			anlsr[j][i]=ejected;

			/********Reprocessed*********/
			Probability=1-pow(1-(ANReproLSR[j]*Synapse.x*dt),M_q);
			REAL reprocessed;
			if(Probability>(rand()/(REAL)RAND_MAX))
			{
				reprocessed=1;
			}
			else
			{
				reprocessed=0;
			}

			/********Replenish**********/
			Probability=1-pow(1-Synapse.y*dt,M_q);
			REAL replenish;
			if(Probability>(rand()/(REAL)RAND_MAX))
			{
				replenish=1;
			}
			else
			{
				replenish=0;
			}

			/*******Update Variables*******/
			ANAvailLSR[j]=ANAvailLSR[j]+replenish-ejected+reprocessed;
			REAL reuptakeandlost=((Synapse.r*dt)+(Synapse.l*dt))*ANCleftLSR[j];
			REAL reuptake=Synapse.r*dt*ANCleftLSR[j];
			ANCleftLSR[j]= ANCleftLSR[j]+ejected-reuptakeandlost;
			ANReproLSR[j]=ANReproLSR[j]+ reuptake-reprocessed;

		}

		/********MSR Fibres*************/
		for (j=0;j<NUMMSR;j++)
		{
			ICa=preSyn.GmaxCa[j+2]*micapowconv*(ihcconv[i]-preSyn.Eca);
			CaCurrMSR[j]+=(-ICa-CaCurrMSR[j])*preSyn.dtTauCa;
			camsr[j][i]=CaCurrMSR[j];
			REAL compare=max(100*(pow(camsr[j][i],preSyn.power)-pow(preSyn.CaTh[j+2],preSyn.power)),0);
			vrrmsr[j][i]=preSyn.z*compare;
			//saturate vrr
			if(vrrmsr[j][i]>Fs)
			{
				vrrmsr[j][i]=Fs;
			}

			/*****Release Probability*****/
			REAL releaseProb=vrrmsr[j][i]*dt;
			REAL M_q=Synapse.M[j+2]-ANAvailMSR[j];
			if(M_q<0)
			{
				M_q=0;
			}

			/*********Ejected************/
			REAL Probability=1-pow(1-releaseProb,ANAvailMSR[j]);
			REAL ejected;

			if (refrac_msr[j]>0)
			{
				//decrement refrac counter until it's zero
				refrac_msr[j]--;
			}

			if(Probability>(rand()/(REAL)RAND_MAX) && refrac_msr[j]==0)
			{
				ejected=1;
				//refrac value to count down from is always between refrac period + random dist between 0 and refrac period
				refrac_msr[j]=round(Synapse.refrac_period/dt + ((rand()/(REAL)RAND_MAX)*Synapse.refrac_period)/dt);
			}
			else
			{
				ejected=0;
			}
			anmsr[j][i]=ejected;

			/********Reprocessed*********/
			Probability=1-pow(1-(ANReproMSR[j]*Synapse.x*dt),M_q);
			REAL reprocessed;
			if(Probability>(rand()/(REAL)RAND_MAX))
			{
				reprocessed=1;
			}
			else
			{
				reprocessed=0;
			}

			/********Replenish**********/
			Probability=1-pow(1-Synapse.y*dt,M_q);
			REAL replenish;
			if(Probability>(rand()/(REAL)RAND_MAX))
			{
				replenish=1;
			}
			else
			{
				replenish=0;
			}

			/*******Update Variables*******/
			ANAvailMSR[j]=ANAvailMSR[j]+replenish-ejected+reprocessed;
			REAL reuptakeandlost=((Synapse.r*dt)+(Synapse.l*dt))*ANCleftMSR[j];
			REAL reuptake=Synapse.r*dt*ANCleftMSR[j];
			ANCleftMSR[j]= ANCleftMSR[j]+ejected-reuptakeandlost;
			ANReproMSR[j]=ANReproMSR[j]+ reuptake-reprocessed;
		}
		/********HSR Fibres************/
		for (j=0;j<NUMHSR;j++)
		{
			ICa=preSyn.GmaxCa[j+4]*micapowconv*(ihcconv[i]-preSyn.Eca);
			CaCurrHSR[j]+=(-ICa-CaCurrHSR[j])*preSyn.dtTauCa;
			cahsr[j][i]=CaCurrHSR[j];
			REAL compare=max(100*(pow(cahsr[j][i],preSyn.power)-pow(preSyn.CaTh[j+4],preSyn.power)),0);
			vrrhsr[j][i]=preSyn.z*compare;//to allow for single float range
			//saturate vrr
			if(vrrhsr[j][i]>Fs)
			{
				vrrhsr[j][i]=Fs;
			}

			/*****Release Probability*****/
			REAL releaseProb=vrrhsr[j][i]*dt;
			REAL M_q=Synapse.M[j+4]-ANAvailHSR[j];
			if(M_q<0)
			{
				M_q=0;
			}

			/*********Ejected************/
			REAL Probability=1-pow(1-releaseProb,ANAvailHSR[j]);
			REAL ejected;
			if (refrac_hsr[j]>0)
			{
				//decrement refrac counter until it's zero
				refrac_hsr[j]--;
			}
			if(Probability>(rand()/(REAL)RAND_MAX) && refrac_hsr[j]==0)
			{
				ejected=1;
				//refrac value to count down from is always between refrac period + random dist between 0 and refrac period
				refrac_hsr[j]=round(Synapse.refrac_period/dt + ((rand()/(REAL)RAND_MAX)*Synapse.refrac_period)/dt);
			}
			else
			{
				ejected=0;
			}
			anhsr[j][i]=ejected;

			/********Reprocessed*********/
			Probability=1-pow(1-(ANReproHSR[j]*Synapse.x*dt),M_q);
			REAL reprocessed;
			if(Probability>(rand()/(REAL)RAND_MAX))
			{
				reprocessed=1;
			}
			else
			{
				reprocessed=0;
			}

			/********Replenish**********/
			Probability=1-pow(1-Synapse.y*dt,M_q);
			REAL replenish;
			if(Probability>(rand()/(REAL)RAND_MAX))
			{
				replenish=1;
			}
			else
			{
				replenish=0;
			}

			/*******Update Variables*******/
			ANAvailHSR[j]=ANAvailHSR[j]+replenish-ejected+reprocessed;
			REAL reuptakeandlost=((Synapse.r*dt)+(Synapse.l*dt))*ANCleftHSR[j];
			REAL reuptake=Synapse.r*dt*ANCleftHSR[j];
			ANCleftHSR[j]= ANCleftHSR[j]+ejected-reuptakeandlost;
			ANReproHSR[j]=ANReproHSR[j]+ reuptake-reprocessed;

		}

	}	

}

void c_main()
{
  // Get core and chip IDs

  coreID = spin1_get_core_id ();
  chipID = spin1_get_chip_id ();

  spin1_start (SYNC_WAIT);
  process_chan();

  spin1_exit (0);

}
