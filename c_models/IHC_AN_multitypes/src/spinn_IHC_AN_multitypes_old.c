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
#include "IHC_AN_multitypes.h"
#include "spin1_api.h"
#include "math.h"
#include "random.h"
#include "stdfix-exp.h"

#define TIMER_TICK_PERIOD  30000
#define TOTAL_TICKS 173//197       
#define PROFILE
//#define LOOP_PROFILE
//#define PRINT

//=========GLOBAL VARIABLES============//
REAL Fs,dt,max_rate;
uint coreID;
uint chipID;
uint test_DMA;
uint seg_index;
uint read_switch;
uint write_switch;
uint processing;
uint index_x;
uint index_y;
//REAL r_max_recip;
//unsigned accum r_max_recip; 
REAL dt_spikes;
uint resamp_fac;
SEED_TYPE local_seed;

//uint seed_selection[SEED_SEL_SIZE];//TODO:this needs to be moved to SDRAM

int start_count_process;
int end_count_process;
int start_count_read;
int end_count_read;
int start_count_write;
int end_count_write;
uint refrac_lsr[NUMLSR];
uint refrac_msr[NUMMSR];
uint refrac_hsr[NUMHSR];

REAL *dtcm_buffer_a;
REAL *dtcm_buffer_b;
REAL *dtcm_buffer_x;
REAL *dtcm_buffer_y;
REAL *dtcm_profile_buffer;

REAL *sdramin_buffer;
REAL *sdramout_buffer;
REAL *profile_buffer;

IHC_ciliaParams Cilia;
IHC_preSynParams preSyn;
IHC_synParams Synapse;

startupVars startupValues;
//recurring values
REAL IHCVnow;
REAL mICaCurr;
REAL CaCurrLSR[NUMLSR];
REAL CaCurrMSR[NUMMSR];
REAL CaCurrHSR[NUMHSR];
REAL ANCleftLSR[NUMLSR];
REAL ANAvailLSR[NUMLSR];
REAL ANReproLSR[NUMLSR];
REAL ANCleftMSR[NUMMSR];
REAL ANAvailMSR[NUMMSR];
REAL ANReproMSR[NUMMSR];
REAL ANCleftHSR[NUMHSR];
REAL ANAvailHSR[NUMHSR];
REAL ANReproHSR[NUMHSR];

uint32_t null_func()
{
	return 0;
}

startupVars generateStartupVars(void)
{
	startupVars out;
	REAL Gu0,kt0,IHCV,gmax,u0,u1,s1,s0,ga,gk,et,ek,rpc,ekp,gamma,beta,mICaCurr,
	gmaxcalsr,gmaxcamsr,gmaxcahsr,eca,tauca,CaCurrLSR,CaCurrMSR,CaCurrHSR,kt0LSR,kt0MSR,kt0HSR,
	ANCleftLSR,ANCleftMSR,ANCleftHSR,ANAvailLSR,ANAvailMSR,ANAvailHSR,ANReproLSR,ANReproMSR,
	ANReproHSR,z,power,y,x,l,r,MLSR,MMSR,MHSR;

	//======Calculate Gu0===========//
	gmax=6e-9;
	u0=5e-9;
	u1=1e-9;
	s1=1e-9;
	s0=3e-8;
	ga=8e-10;
	Gu0=ga+gmax/(1+(exp(u0/s0)*(1+exp(u1/s1))));
	//======Calculate IHCVnow=========//
	gk=2e-8;
	et=0.1;
	ek=-0.08;
	rpc=0.04;
	ekp=ek+rpc*et;
	IHCV=(gk*ekp+Gu0*et)/(Gu0+gk);

	//=======Calculate mICaCurrent========//
	gamma=130.;
	beta=400.;
	mICaCurr=1./(1.+exp(-gamma*IHCV)*(1./beta));

	//=======Calculate CaCurrent========//
	gmaxcalsr=2.5e-11;
	gmaxcamsr=8e-11;
	gmaxcahsr=5e-10;
	eca=0.066;
	tauca=1e-4;
	CaCurrLSR=((gmaxcalsr*pow(mICaCurr,3))*(IHCV-eca))*tauca;
	CaCurrMSR=((gmaxcamsr*pow(mICaCurr,3))*(IHCV-eca))*tauca;
	CaCurrHSR=((gmaxcahsr*pow(mICaCurr,3))*(IHCV-eca))*tauca;

	//========Calculate AN==========//
	z=1e38;//reduced from 1e40 fit as single float type
	power=3;
	y=20;
	x=30;
	MLSR=15;
	MMSR=13;
	MHSR=11;
	l=900;
	r=100;

	kt0LSR=-z*pow(CaCurrLSR,power)*100;//added *100 for single float z
	kt0MSR=-z*pow(CaCurrMSR,power)*100;
	kt0HSR=-z*pow(CaCurrHSR,power)*100;
	ANCleftLSR=(kt0LSR*y*MLSR)/(y*(l+r)+kt0LSR*l);
	ANCleftMSR=(kt0MSR*y*MMSR)/(y*(l+r)+kt0MSR*l);
	ANCleftHSR=(kt0HSR*y*MHSR)/(y*(l+r)+kt0HSR*l);
	ANAvailLSR=round((ANCleftLSR*(l+r))/kt0LSR);
	ANAvailMSR=round((ANCleftMSR*(l+r))/kt0MSR);
	ANAvailHSR=round((ANCleftHSR*(l+r))/kt0HSR);
	ANReproLSR=(ANCleftLSR*r)/x;
	ANReproMSR=(ANCleftMSR*r)/x;
	ANReproHSR=(ANCleftHSR*r)/x;

	out.IHCVnow0=IHCV;
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

//application initialisation
void app_init(void)
{
	//Fs=REAL_CONST(44100.);
	Fs=REAL_CONST(32000.);
	dt=(1.0/Fs);
	resamp_fac=10;
	max_rate=Fs/(REAL)resamp_fac;
	dt_spikes=(REAL)resamp_fac*dt;
	seg_index=0;
	read_switch=0;
	write_switch=0;
	//r_max_recip=REAL_CONST(1.)/(REAL)RDM_MAX;
	//r_max_recip=1uk/(unsigned accum)RDM_MAX;
	
	/* say hello */
	
	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);
	
	// Allocate buffers somewhere in SDRAM
	
	//output results buffer, 10x fibre outputs
	sdramout_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					 NUMFIBRES * (uint)44100. * sizeof(REAL),
					 coreID,
					 ALLOC_LOCK);	

	sdramin_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					NUMFIBRES * (uint)44100. *sizeof(REAL),	//this shouldn't have NUMFIBRES
					coreID|32,
					ALLOC_LOCK);
	
	profile_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					3 * ((uint)44100./SEGSIZE) *sizeof(REAL),
					coreID|64,
					ALLOC_LOCK);
	
	// and a buffer in DTCM
	
	dtcm_buffer_a = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_b = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_x = (REAL *) sark_alloc (SEGSIZE*NUMFIBRES, sizeof(REAL));
	dtcm_buffer_y = (REAL *) sark_alloc (SEGSIZE*NUMFIBRES, sizeof(REAL));
	dtcm_profile_buffer = (REAL *) sark_alloc (3*TOTAL_TICKS, sizeof(REAL));
	
	if (dtcm_buffer_a == NULL ||dtcm_buffer_b == NULL ||dtcm_buffer_x == NULL ||dtcm_buffer_y == NULL 
			||  sdramout_buffer == NULL || sdramin_buffer == NULL || profile_buffer == NULL || dtcm_profile_buffer == NULL)
	{
		test_DMA = FALSE;
		io_printf (IO_BUF, "[core %d] error - cannot allocate buffer\n", coreID);
	}
	else
	{
		test_DMA = TRUE;
		// initialize sections of DTCM, system RAM and SDRAM
		for (uint i = 0; i < SEGSIZE; i++)
		{
			dtcm_buffer_a[i]   = 0;
			dtcm_buffer_b[i]   = 0;
		}
		for (uint i = 0; i < SEGSIZE*NUMFIBRES; i++)
		{
			dtcm_buffer_x[i]   = 0;
			dtcm_buffer_y[i]   = 0;
		}
	
		for (uint i=0;i<NUMFIBRES * (uint)44100.;i++)
		{
			sdramout_buffer[i]  = 0;
			sdramin_buffer[i]  = 0;
		}
		
		for (uint i=0;i<3 * TOTAL_TICKS;i++)
		{
			dtcm_profile_buffer[i]  = 0;
			profile_buffer[i]  = 0;
		}
		
		io_printf (IO_BUF, "[core %d] dtcm buffer a @ 0x%08x\n", coreID,
				   (uint) dtcm_buffer_a);
		io_printf (IO_BUF, "[core %d] sdram out buffer @ 0x%08x\n", coreID,
				   (uint) sdramout_buffer);
		io_printf (IO_BUF, "[core %d] sdram in buffer @ 0x%08x\n", coreID,
				   (uint) sdramin_buffer);	
		io_printf (IO_BUF, "[core %d] profile buffer @ 0x%08x\n", coreID,
				   (uint) profile_buffer);	
	}
	
	//============MODEL INITIALISATION================//
	
	//calculate startup values
	startupValues=generateStartupVars();
	
	//initialise random number gen
	//io_printf (IO_BUF, "[core %d] seed_selection=",coreID);
	/*for(uint i=0;i<SEED_SEL_SIZE;i++)
	{
		seed_selection[i]=mars_kiss32();
		//io_printf (IO_BUF, "%u ",seed_selection[i]);

	}*/
	//io_printf (IO_BUF,"\n");

	io_printf (IO_BUF, "[core %d][chip %d] local_seeds=",coreID,chipID);
	for(uint i=0;i<4;i++)
	{
		//local_seed[i]=seed_selection[((chipID<<8) | (coreID<<i))];
		//local_seed[i]=seed_selection[(coreID<<i)|(chipID)];
		local_seed[i]= ((chipID<<8) | (coreID<<i));
		io_printf (IO_BUF, "%u ",local_seed[i]);
	}
	io_printf (IO_BUF,"\n");

	validate_mars_kiss64_seed (local_seed);

	//initialise cilia
	Cilia.tc= REAL_CONST(0.00012);
	Cilia.C=REAL_CONST(0.004);
	Cilia.u0=REAL_CONST(5e-9);
	Cilia.recips0=REAL_CONST(1.)/REAL_CONST(3e-8);
	Cilia.u1=REAL_CONST(1e-9);
	Cilia.recips1=REAL_CONST(1.)/REAL_CONST(1e-9);
	Cilia.Gmax=REAL_CONST(6e-9);
	Cilia.Ga=REAL_CONST(8e-10);
	Cilia.dtCap=(dt/4e-12);
	Cilia.Et=REAL_CONST(0.1);
	Cilia.Gk=REAL_CONST(2e-8);
	Cilia.Ek=REAL_CONST(-0.08);
	Cilia.Rpc=REAL_CONST(0.04);
	
	//==========Recurring Values=================//
	IHCVnow=startupValues.IHCVnow0;
	mICaCurr=startupValues.mICaCurr0;;
	for(uint i=0;i<NUMLSR;i++)
	{
		CaCurrLSR[i]=startupValues.CaCurrLSR0;
		ANCleftLSR[i]=startupValues.ANCleftLSR0;
		ANAvailLSR[i]=startupValues.ANAvailLSR0;
		ANReproLSR[i]=startupValues.ANReproLSR0;
		refrac_lsr[i]=0;
	}
	for(uint i=0;i<NUMMSR;i++)
	{
		CaCurrMSR[i]=startupValues.CaCurrMSR0;
		ANCleftMSR[i]=startupValues.ANCleftMSR0;
		ANAvailMSR[i]=startupValues.ANAvailMSR0;
		ANReproMSR[i]=startupValues.ANReproMSR0;
		refrac_msr[i]=0;
	}
	for(uint i=0;i<NUMHSR;i++) 
	{
		CaCurrHSR[i]=startupValues.CaCurrHSR0;
		ANCleftHSR[i]=startupValues.ANCleftHSR0;
		ANAvailHSR[i]=startupValues.ANAvailHSR0;
		ANReproHSR[i]=startupValues.ANReproHSR0;
		refrac_hsr[i]=0;
	}
	
	//=========initialise the pre synapse params========//
	preSyn.recipBeta=2.5e-3;
	preSyn.gamma=REAL_CONST(130.);
	preSyn.dtTauM=(dt/1e-4);
	preSyn.dtTauCa=(dt/1e-4);
	preSyn.Eca=0.066;
	preSyn.power=REAL_CONST(3.);
	preSyn.z=1e38;//reduced from 1e40 to fit within float range (then an additional *100 is added in code)

	preSyn.GmaxCa[0]=2.5e-11;
	preSyn.GmaxCa[1]=2.5e-11;
	preSyn.GmaxCa[2]=8e-11;
	preSyn.GmaxCa[3]=8e-11;
	preSyn.GmaxCa[4]=5e-10;
	preSyn.GmaxCa[5]=5e-10;
	preSyn.GmaxCa[6]=5e-10;
	preSyn.GmaxCa[7]=5e-10;
	preSyn.GmaxCa[8]=5e-10;
	preSyn.GmaxCa[9]=5e-10;

	//the follwing variables are equal to z*CaTh^3 used in the vesicle release rate equations
	preSyn.CaTh[0]=1e-7;
	preSyn.CaTh[1]=1e-7;
	preSyn.CaTh[2]=0.0216;
	preSyn.CaTh[3]=0.0216;
	preSyn.CaTh[4]=0.2924207;
	preSyn.CaTh[5]=0.2924207;
	preSyn.CaTh[6]=0.2924207;
	preSyn.CaTh[7]=0.2924207;
	preSyn.CaTh[8]=0.2924207;
	preSyn.CaTh[9]=0.2924207;

	//=======initialise the synapse params=======//

	Synapse.ldt=REAL_CONST(900.)*dt_spikes;
	Synapse.ydt=REAL_CONST(20.)*dt_spikes;
	if(Synapse.ydt>REAL_CONST(1.))
	{
		Synapse.ydt=REAL_CONST(1.);
	}
	Synapse.xdt=REAL_CONST(30.)*dt_spikes;
	Synapse.rdt=REAL_CONST(100.)*dt_spikes;
	//Synapse.refrac_period=((7.5e-4)/dt_spikes);
	Synapse.refrac_period=(unsigned accum)((7.5e-4)/dt_spikes);
	//io_printf (IO_BUF,"Synapse.refrac_period=%f\n",Synapse.refrac_period);
	/*Synapse.M[0]=REAL_CONST(13.);
	Synapse.M[1]=REAL_CONST(13.);
	Synapse.M[2]=REAL_CONST(14.);
	Synapse.M[3]=REAL_CONST(14.);
	Synapse.M[4]=REAL_CONST(8.);
	Synapse.M[5]=REAL_CONST(8.);
	Synapse.M[6]=REAL_CONST(8.);
	Synapse.M[7]=REAL_CONST(8.);
	Synapse.M[8]=REAL_CONST(8.);
	Synapse.M[9]=REAL_CONST(8.);*/
	
	Synapse.M[0]=REAL_CONST(15.);
	Synapse.M[1]=REAL_CONST(15.);
	Synapse.M[2]=REAL_CONST(13.);
	Synapse.M[3]=REAL_CONST(13.);
	Synapse.M[4]=REAL_CONST(11.);
	Synapse.M[5]=REAL_CONST(11.);
	Synapse.M[6]=REAL_CONST(11.);
	Synapse.M[7]=REAL_CONST(11.);
	Synapse.M[8]=REAL_CONST(11.);
	Synapse.M[9]=REAL_CONST(11.);
	

	
	//===========casting test=========================//
/*	uint test_u;
	long fract test_lf;
	for(uint i=0;i<10;i++)		
	{
		ui_ulf t;
		//test=(uint)(Synapse.refrac_period + ((unsigned accum)((long fract)random_gen())*Synapse.refrac_period)+0.5uk);
		t.ui = random_gen();
		io_printf (IO_BUF,"test uint=%x\n",t.ui);
		io_printf (IO_BUF,"test long fract=%f\n",(unsigned accum)t.ulf); 
		//io_printf (IO_BUF,"sizeof longfract=%d\n",sizeof(test_lf));
	}*/
	  

#ifdef PROFILE
    // configure timer 2 for profiling
    // enabled, free running, interrupt disabled, no pre-scale, 32 bit, free-running mode
    tc[T2_CONTROL] = TIMER2_CONF;
#endif
    
}
void data_write(uint null_a, uint null_b)
{
	REAL *dtcm_buffer_out;
	uint out_index;
	
	if(test_DMA == TRUE)
	{
		if(!write_switch)
		{
			out_index=index_x;
			dtcm_buffer_out=dtcm_buffer_x;
#ifdef PRINT	
			io_printf (IO_BUF, "buff_x write\n");
#endif
		}
		else
		{
			out_index=index_y;
			dtcm_buffer_out=dtcm_buffer_y;
#ifdef PRINT
			io_printf (IO_BUF, "buff_y write\n");
#endif
		}
#ifdef PROFILE
  start_count_write = tc[T2_COUNT];
#endif
		spin1_dma_transfer(DMA_WRITE,&sdramout_buffer[out_index],dtcm_buffer_out,DMA_WRITE,
		  						NUMFIBRES*SEGSIZE*sizeof(REAL));
#ifdef PRINT
		io_printf (IO_BUF, "[core %d] segment %d written to @ 0x%08x - 0x%08x\n", coreID,seg_index,
							  (uint) &sdramout_buffer[out_index],(uint) &sdramout_buffer[out_index+(NUMFIBRES-1)*SEGSIZE+SEGSIZE-1]);
#endif
	}
}

//DMA read
void data_read(uint ticks, uint null)
{
#ifdef PROFILE
  start_count_read = tc[T2_COUNT];
#endif

	REAL *dtcm_buffer_in;

	//read from DMA and copy into DTCM
	if(test_DMA == TRUE)
	{
		//assign recieve buffer
		if(!read_switch)	
		{
			dtcm_buffer_in=dtcm_buffer_a;
			read_switch=1;
#ifdef PRINT
			io_printf (IO_BUF, "buff_a read\n");
#endif
		}
		else
		{
			dtcm_buffer_in=dtcm_buffer_b;
			read_switch=0;
#ifdef PRINT
			io_printf (IO_BUF, "buff_b read\n");
#endif
		}

#ifdef PRINT
		io_printf (IO_BUF, "[core %d] sdram DMA read @ 0x%08x (segment %d)\n", coreID,
					  (uint) &sdramin_buffer[(seg_index)*SEGSIZE],seg_index+1);
#endif
		
		spin1_dma_transfer(DMA_READ,&sdramin_buffer[seg_index*SEGSIZE], dtcm_buffer_in, DMA_READ,
			   SEGSIZE*sizeof(REAL));
	}
	
	// stop if desired number of ticks reached
	if (ticks > TOTAL_TICKS) 
	{
		io_printf (IO_BUF, "spinn_exit\n");
		spin1_exit (0); 
	}
	
}

/*uint seed_scrambled_rgen(void)
{
	//function to scramble the application seed before a PRNG function call
	seed[0]=seed[0] ^ mars_kiss32();
	seed[1]=seed[1] ^ mars_kiss32();
	seed[2]=seed[2] ^ mars_kiss32();
	seed[3]=seed[3] ^ mars_kiss32();	
	
	return mars_kiss64_seed(seed);
}*/

uint process_chan(REAL *out_buffer,REAL *in_buffer) 
{  
	uint segment_offset=NUMFIBRES*SEGSIZE*(seg_index-1);
	uint i,j,k;
	
	REAL utconv;
	accum ex1;
	accum ex2;
	//REAL ex1;
	//REAL ex2;
	
	REAL Guconv;
	REAL mICaINF;
	REAL micapowconv;
	REAL ICa;
	REAL pos_CaCurr;
	REAL CaCurr_pow;
	//REAL releaseProb_pow;
	unsigned long fract releaseProb_pow;
	REAL Repro_rate;
	REAL Synapse_ypow;
	REAL Synapse_xpow;
	REAL compare;
	REAL vrrlsr;
	REAL vrrmsr;
	REAL vrrhsr;	
	REAL releaseProb;
	//unsigned long fract releaseProb;
	REAL M_q;
	REAL Probability;
	//unsigned long fract Probability;
	REAL ejected;
	REAL spikes;
	REAL reprocessed;
	REAL replenish;
	REAL reuptakeandlost;
	REAL reuptake;
	
	ui_ulf rand_gen_union;
		
#ifdef PRINT
	io_printf (IO_BUF, "[core %d] segment %d (offset=%d) starting processing\n", coreID,seg_index,segment_offset);
#endif
	for(i=0;i<SEGSIZE;i++)
	{			  		
		//===========Apply Scaler============//
  	  	utconv=in_buffer[i] * Cilia.C;
  	  	
		//=========Apical Conductance========//	

		ex1=expk((accum)(-(utconv-Cilia.u1)*Cilia.recips1));
		ex2=expk((accum)(-(utconv-Cilia.u0)*Cilia.recips0));
		
		//ex1=exp(-(utconv-Cilia.u1)*Cilia.recips1);
		//ex2=exp(-(utconv-Cilia.u0)*Cilia.recips0);
  	  	
  	  	Guconv=Cilia.Ga + (Cilia.Gmax/(REAL_CONST(1.)+(REAL)ex2*(REAL_CONST(1.)+(REAL)ex1)));
 
		//========Receptor Potential=========//
		IHCVnow+=((-Guconv*(IHCVnow-Cilia.Et))-(Cilia.Gk*(IHCVnow-startupValues.Ekp)))*Cilia.dtCap;
	
		//================mICa===============//
		ex1=expk((accum)-preSyn.gamma*IHCVnow);
		//ex1=exp(-preSyn.gamma*IHCVnow);
		
		mICaINF=REAL_CONST(1.)/(REAL_CONST(1.)+(REAL)ex1*preSyn.recipBeta);
		mICaCurr+=(mICaINF-mICaCurr)*preSyn.dtTauM;
		if(mICaCurr<REAL_CONST(0.))
		{
			mICaCurr=REAL_CONST(0.);
		}

		//================ICa================//
		//calculate values for each fibre	
		micapowconv=mICaCurr;
		for(k=0;k<preSyn.power-1;k++)
		{
			micapowconv=micapowconv*mICaCurr;
		}
		
			

#ifdef LOOP_PROFILE
  int end_apical_count_read = tc[T2_COUNT];
  io_printf (IO_BUF, "first half complete in %d ticks (sample %d, segment %d)\n",start_apical_count_read-end_apical_count_read,i,seg_index);
#endif	
  
#ifdef LOOP_PROFILE
  int start_max_count_read = tc[T2_COUNT];
#endif

/*	if(i==0)
	{
	#ifdef PROFILE
	  start_count_process = tc[T2_COUNT];
	#endif
	}*/

	/	//============LSR Fibres=============//
		for (j=0;j<NUMLSR;j++)
		{
			//======Synaptic Ca========//
			ICa=preSyn.GmaxCa[j]*micapowconv*(IHCVnow-preSyn.Eca);
			CaCurrLSR[j]+=(ICa-CaCurrLSR[j])*preSyn.dtTauCa;
			
			if(CaCurrLSR[j]>REAL_CONST(0.))
			{
				CaCurrLSR[j]=REAL_CONST(0.);
			}
			//invert Ca
			pos_CaCurr=REAL_CONST(-1.)*CaCurrLSR[j];

			if(i%resamp_fac==0)
			{			
				//=====Vesicle Release Rate=====//
				CaCurr_pow=pos_CaCurr;
				for(k=0;k<(uint)preSyn.power-1;k++)
				{			
					CaCurr_pow=CaCurr_pow*pos_CaCurr;
				}
				compare=max((preSyn.z*CaCurr_pow-preSyn.CaTh[j]),REAL_CONST(0.));	
				vrrlsr=compare*REAL_CONST(100.);
				//saturate vrr
				if(vrrlsr>max_rate)
				{
					vrrlsr=max_rate;
				}
		
				//=====Release Probability=======//
				releaseProb=vrrlsr*dt_spikes;
				//releaseProb=(unsigned long fract)(vrrlsr*dt_spikes);
				M_q=Synapse.M[j]-ANAvailLSR[j];
				if(M_q<REAL_CONST(0.))
				{
					M_q=REAL_CONST(0.);
				}
	
				//===========Ejected============//
				/*releaseProb_pow=REAL_CONST(1.);
				for(k=0;k<(uint)ANAvailLSR[j];k++)
				{			
					releaseProb_pow=releaseProb_pow*(REAL_CONST(1.)-releaseProb);
				}
	
				Probability=REAL_CONST(1.)-releaseProb_pow;*/
				
				releaseProb_pow=1ulr;
				for(k=0;k<(uint)ANAvailLSR[j];k++)
				{			
					releaseProb_pow=releaseProb_pow*(1ulr-(unsigned long fract)releaseProb);
				}
	
				//Probability=1ulr-releaseProb_pow;
				Probability=REAL_CONST(1.)-releaseProb_pow;
			
				if (refrac_lsr[j]>0)
				{
					refrac_lsr[j]--;
				}	
								
				/*if(Probability>((REAL)random_gen()*r_max_recip))
				{
					ejected=REAL_CONST(1.);
					if (refrac_lsr[j]<=0)
					{
						spikes=REAL_CONST(1.);
						refrac_lsr[j]=(uint)(Synapse.refrac_period + (((REAL)random_gen()*r_max_recip)*Synapse.refrac_period)+REAL_CONST(0.5));
					}
					else
					{
						spikes=REAL_CONST(0.);
					}
				}*/
				rand_gen_union.ui=random_gen();				
				if((unsigned long fract)Probability>rand_gen_union.ulf)
				{
					ejected=REAL_CONST(1.);
					if (refrac_lsr[j]<=0)
					{
						spikes=REAL_CONST(1.);
						rand_gen_union.ui=random_gen();
						refrac_lsr[j]=(uint)(Synapse.refrac_period + ((unsigned accum)rand_gen_union.ulf * Synapse.refrac_period)+0.5uk);
					}
					else
					{
						spikes=REAL_CONST(0.);
					}
				}
				else
				{
					ejected=REAL_CONST(0.);
					spikes=REAL_CONST(0.);
				}
	
				//=========Reprocessed=========//
				Repro_rate=ANReproLSR[j]*Synapse.xdt;
				if(Repro_rate>REAL_CONST(1.))
				{
					Repro_rate=REAL_CONST(1.);
				}
	
				Synapse_xpow=REAL_CONST(1.);
				Synapse_ypow=REAL_CONST(1.);
				for(k=0;k<(uint)M_q;k++)
				{			
					Synapse_xpow=Synapse_xpow*(REAL_CONST(1.)-Repro_rate);
					Synapse_ypow=Synapse_ypow*(REAL_CONST(1.)-Synapse.ydt);
				}

				/*Probability=REAL_CONST(1.)-Synapse_xpow;			
				if(Probability>((REAL)random_gen()*r_max_recip))
				{
					reprocessed=REAL_CONST(1.);
				}
				else
				{
					reprocessed=REAL_CONST(0.);
				}
	
				//========Replenish==========//	
				Probability=REAL_CONST(1.)-Synapse_ypow;			
				if(Probability>((REAL)random_gen()*r_max_recip))
				{
					replenish=REAL_CONST(1.);
				}
				else
				{
					replenish=REAL_CONST(0.);
				}*/
				
				Probability=REAL_CONST(1.)-Synapse_xpow;		
				//Probability=1ulr-(unsigned long fract)Synapse_xpow;
				rand_gen_union.ui=random_gen();
				if((unsigned long fract)Probability>rand_gen_union.ulf)
				{
					reprocessed=REAL_CONST(1.);
				}
				else
				{
					reprocessed=REAL_CONST(0.);
				}
	
				//========Replenish==========//	
				//Probability=1ulr-(unsigned long fract)Synapse_ypow;
				Probability=REAL_CONST(1.)-Synapse_ypow;			
				rand_gen_union.ui=random_gen();
				if((unsigned long fract)Probability>rand_gen_union.ulf)
				{
					replenish=REAL_CONST(1.);
				}
				else
				{
					replenish=REAL_CONST(0.);
				}
	
				//==========Update Variables=========//
				ANAvailLSR[j]=ANAvailLSR[j]+replenish+reprocessed-ejected;
				if(ANAvailLSR[j]<REAL_CONST(0.))
				{
					ANAvailLSR[j]=REAL_CONST(0.);
				}
				reuptakeandlost=(Synapse.rdt+Synapse.ldt)*ANCleftLSR[j];
				reuptake=Synapse.rdt*ANCleftLSR[j];
				if(reuptake>REAL_CONST(1.))
				{
					reuptake=REAL_CONST(1.);
				}

				ANCleftLSR[j]= ANCleftLSR[j]+ejected-reuptakeandlost;
				if(ANCleftLSR[j]<REAL_CONST(0.))
				{
					ANCleftLSR[j]=REAL_CONST(0.);
				}
				ANReproLSR[j]=ANReproLSR[j]+ reuptake-reprocessed;
				if(ANReproLSR[j]<REAL_CONST(0.))
				{
					ANReproLSR[j]=REAL_CONST(0.);
				}

			}
			else
			{
				spikes=REAL_CONST(0.);
			}	
			//=======write value to SDRAM========//
			//out_buffer[(j*SEGSIZE)+i]  = ANReproLSR[j];//ANAvailLSR[j];//vrrlsr;// compare;//CaCurr_pow;//
			out_buffer[(j*SEGSIZE)+i]  = spikes;
			
		}


		//=======MSR Fibres===========//
		for (j=0;j<NUMMSR;j++)
		{
			ICa=preSyn.GmaxCa[j+2]*micapowconv*(IHCVnow-preSyn.Eca);
			CaCurrMSR[j]+=(ICa-CaCurrMSR[j])*preSyn.dtTauCa;
			
			if(CaCurrMSR[j]>REAL_CONST(0.))
			{
				CaCurrMSR[j]=REAL_CONST(0.);
			}
			
			pos_CaCurr=REAL_CONST(-1.)*CaCurrMSR[j];
			
			if(i%resamp_fac==0)
			{
				CaCurr_pow=pos_CaCurr;
				for(k=0;k<(uint)preSyn.power-1;k++)
				{			
					CaCurr_pow=CaCurr_pow*pos_CaCurr;
				}
				compare=max((preSyn.z*CaCurr_pow-preSyn.CaTh[j+2]),REAL_CONST(0.));
				vrrmsr=compare*REAL_CONST(100.);
				//saturate vrr
				if(vrrmsr>max_rate)
				{
					vrrmsr=max_rate;
				}
				//=====Release Probability=======//
				releaseProb=vrrmsr*dt_spikes;
				//releaseProb=(unsigned long fract)(vrrmsr*dt_spikes);
				M_q=Synapse.M[j+2]-ANAvailMSR[j];
				if(M_q<REAL_CONST(0.))
				{
					M_q=REAL_CONST(0.);
				}
	
				//===========Ejected===========//
			
				/*releaseProb_pow=REAL_CONST(1.);
				for(k=0;k<(uint)ANAvailMSR[j];k++)
				{			
					releaseProb_pow=releaseProb_pow*(REAL_CONST(1.)-releaseProb);
				}
				Probability=REAL_CONST(1.)-releaseProb_pow;			*/
				releaseProb_pow=1ulr;
				for(k=0;k<(uint)ANAvailMSR[j];k++)
				{			
					releaseProb_pow=releaseProb_pow*(1ulr-(unsigned long fract)releaseProb);
				}	
				//Probability=1ulr-releaseProb_pow;
				Probability=REAL_CONST(1.)-(REAL)releaseProb_pow;		
			
				if (refrac_msr[j]>0)
				{
					refrac_msr[j]--;
				}				
				
				/*if(Probability>((REAL)random_gen()*r_max_recip))
				{
					ejected=REAL_CONST(1.);
					if (refrac_msr[j]<=0)
					{
						spikes=REAL_CONST(1.);
						refrac_msr[j]=(uint)(Synapse.refrac_period + (((REAL)random_gen()*r_max_recip)*Synapse.refrac_period)+REAL_CONST(0.5));
					}
					else
					{
						spikes=REAL_CONST(0.);
					}
				}*/
				rand_gen_union.ui=random_gen();
				if((unsigned long fract)Probability>rand_gen_union.ulf)
				{
					ejected=REAL_CONST(1.);
					if (refrac_msr[j]<=0)
					{
						spikes=REAL_CONST(1.);
						rand_gen_union.ui=random_gen();
						refrac_msr[j]=(uint)(Synapse.refrac_period + ((unsigned accum)rand_gen_union.ulf * Synapse.refrac_period)+0.5uk);
					}
					else
					{
						spikes=REAL_CONST(0.);
					}
				}				
				else
				{
					ejected=REAL_CONST(0.);
					spikes=REAL_CONST(0.);
				}
	
				//========Reprocessed==========//		
				Repro_rate=ANReproMSR[j]*Synapse.xdt;
				if(Repro_rate>REAL_CONST(1.))
				{
					Repro_rate=REAL_CONST(1.);
				}
				
				Synapse_xpow=REAL_CONST(1.);
				Synapse_ypow=REAL_CONST(1.);
				for(k=0;k<(uint)M_q;k++)
				{			
					Synapse_xpow=Synapse_xpow*(REAL_CONST(1.)-Repro_rate);
					Synapse_ypow=Synapse_ypow*(REAL_CONST(1.)-(Synapse.ydt));
				}
	
				/*Probability=REAL_CONST(1.)-Synapse_xpow;			
				if(Probability>((REAL)random_gen()*r_max_recip))
				{
					reprocessed=REAL_CONST(1.);
				}
				else
				{
					reprocessed=REAL_CONST(0.);
				}
	
				//=========Replenish============//
				
				Probability=REAL_CONST(1.)-Synapse_ypow;
				if(Probability>((REAL)random_gen()*r_max_recip))
				{
					replenish=REAL_CONST(1.);
				}
				else
				{
					replenish=REAL_CONST(0.);
				}*/
				Probability=REAL_CONST(1.)-Synapse_xpow;	
				//Probability=1ulr-(unsigned long fract)Synapse_xpow;
				rand_gen_union.ui=random_gen();
				if((unsigned long fract)Probability>rand_gen_union.ulf)
				{
					reprocessed=REAL_CONST(1.);
				}
				else
				{
					reprocessed=REAL_CONST(0.);
				}
	
				//=========Replenish============//				
				Probability=REAL_CONST(1.)-Synapse_ypow;
				//Probability=1ulr-(unsigned long fract)Synapse_ypow;
				rand_gen_union.ui=random_gen();
				if((unsigned long fract)Probability>rand_gen_union.ulf)
				{
					replenish=REAL_CONST(1.);
				}
				else
				{
					replenish=REAL_CONST(0.);
				}				
	
				//========Update Variables=======//
				ANAvailMSR[j]=ANAvailMSR[j]+replenish+reprocessed-ejected;
				if(ANAvailMSR[j]<REAL_CONST(0.))
				{
					ANAvailMSR[j]=REAL_CONST(0.);
				}
				reuptakeandlost=(Synapse.rdt+Synapse.ldt)*ANCleftMSR[j];
				reuptake=Synapse.rdt*ANCleftMSR[j];
				if(reuptake>REAL_CONST(1.))
				{
					reuptake=REAL_CONST(1.);
				}				

				ANCleftMSR[j]= ANCleftMSR[j]+ejected-reuptakeandlost;
				if(ANCleftMSR[j]<REAL_CONST(0.))
				{
					ANCleftMSR[j]=REAL_CONST(0.);
				}

				ANReproMSR[j]=ANReproMSR[j]+ reuptake-reprocessed;
				if(ANReproMSR[j]<REAL_CONST(0.))
				{
					ANReproMSR[j]=REAL_CONST(0.);
				}
			
			}
			else
			{
				spikes=REAL_CONST(0.);
			}	
			//=======write value to SDRAM=======//  
			//out_buffer[((j+2)*SEGSIZE)+i] =ANReproMSR[j];//ANAvailMSR[j];//vrrmsr;// compare;//CaCurr_pow;//
			out_buffer[((j+2)*SEGSIZE)+i] =spikes;//ejected;//
			
		}
				
		//==========HSR Fibres============//
		for (j=0;j<NUMHSR;j++)
		{
			ICa=preSyn.GmaxCa[j+4]*micapowconv*(IHCVnow-preSyn.Eca);			
			CaCurrHSR[j]+=(ICa-CaCurrHSR[j])*preSyn.dtTauCa;
			
			if(CaCurrHSR[j]>REAL_CONST(0.))
			{
				CaCurrHSR[j]=REAL_CONST(0.);
			}
			
			pos_CaCurr=REAL_CONST(-1.)*CaCurrHSR[j];			
			
			if(i%resamp_fac==0)
			{
				CaCurr_pow=pos_CaCurr;
				for(k=0;k<(uint)preSyn.power-1;k++)
				{			
					CaCurr_pow=CaCurr_pow*pos_CaCurr;
				}
				compare=max((preSyn.z*CaCurr_pow-preSyn.CaTh[j+4]),REAL_CONST(0.));
				vrrhsr=compare*REAL_CONST(100.);//to allow for single float range
				//saturate vrr
				if(vrrhsr>max_rate)
				{
					vrrhsr=max_rate;
				}
	
				//=======Release Probability======//
				releaseProb=vrrhsr*dt_spikes;
				//releaseProb=(unsigned long fract)(vrrhsr*dt_spikes);
				M_q=Synapse.M[j+4]-ANAvailHSR[j];
				if(M_q<REAL_CONST(0.))
				{
					M_q=REAL_CONST(0.);
				}
	
				//==========Ejected==============//
				//releaseProb_pow=REAL_CONST(1.);
				/*for(k=0;k<(uint)ANAvailHSR[j];k++)
				{			
					releaseProb_pow=releaseProb_pow*(REAL_CONST(1.)-releaseProb);
				}	
				Probability=REAL_CONST(1.)-releaseProb_pow;*/
				
				releaseProb_pow=1ulr;
				for(k=0;k<(uint)ANAvailHSR[j];k++)
				{			
					releaseProb_pow=releaseProb_pow*(1ulr-(unsigned long fract)releaseProb);
				}			
				//Probability=1ulr-releaseProb_pow;	
				Probability=REAL_CONST(1.)-(REAL)releaseProb_pow;


				if (refrac_hsr[j]>0)
				{				
					refrac_hsr[j]--;
				}
				
				/*if(Probability>((REAL)random_gen()*r_max_recip))
				{
					ejected=REAL_CONST(1.);
					if (refrac_hsr[j]<=0)
					{
						spikes=REAL_CONST(1.);
						refrac_hsr[j]=(uint)(Synapse.refrac_period + (((REAL)random_gen()*r_max_recip)*Synapse.refrac_period)+REAL_CONST(0.5));
					}
					else
					{
						spikes=REAL_CONST(0.);
					}
				}*/
				rand_gen_union.ui=random_gen();
				if((unsigned long fract)Probability>rand_gen_union.ulf)
				{
					ejected=REAL_CONST(1.);
					if (refrac_hsr[j]<=0)
					{
						spikes=REAL_CONST(1.);
						rand_gen_union.ui=random_gen();
						refrac_hsr[j]=(uint)(Synapse.refrac_period + ((unsigned accum)rand_gen_union.ulf * Synapse.refrac_period)+0.5uk);
						//io_printf (IO_BUF,"refrac_hsr[j]=%d\n",refrac_hsr[j]);
					}
					else
					{
						spikes=REAL_CONST(0.);
					}
				}				
				else
				{
					ejected=REAL_CONST(0.);
					spikes=REAL_CONST(0.);
				}
	
				//========Reprocessed==========//
				Repro_rate=ANReproHSR[j]*Synapse.xdt;
				if(Repro_rate>REAL_CONST(1.))
				{
					Repro_rate=REAL_CONST(1.);
				}
				
				//(1-Repro_rate)^M_q and (1-Synapse.ydt)^M_q
				Synapse_xpow=REAL_CONST(1.);
				Synapse_ypow=REAL_CONST(1.);				
				for(k=0;k<(uint)M_q;k++)
				{	
					Synapse_ypow=Synapse_ypow*(REAL_CONST(1.)-Synapse.ydt);				
					Synapse_xpow=Synapse_xpow*(REAL_CONST(1.) - Repro_rate);
				}
				/*Probability=REAL_CONST(1.)-Synapse_xpow;					
				/*
				if(Probability>((REAL)random_gen()*r_max_recip))				
				{
					reprocessed=REAL_CONST(1.);
				}
				else
				{
					reprocessed=REAL_CONST(0.);
				}
	
				//===========Replenish=============//
	
				Probability=REAL_CONST(1.)-Synapse_ypow;			
				
				if(Probability>((REAL)random_gen()*r_max_recip))
				{
					replenish=REAL_CONST(1.);
				}
				else
				{
					replenish=REAL_CONST(0.);
				}*/
				
				//Probability=1ulr-(unsigned long fract)Synapse_xpow;	
				Probability=REAL_CONST(1.)-Synapse_xpow;
				rand_gen_union.ui=random_gen();				
				if((unsigned long fract)Probability>rand_gen_union.ulf)				
				{
					reprocessed=REAL_CONST(1.);
				}
				else
				{
					reprocessed=REAL_CONST(0.);
				}
	
				//===========Replenish=============//
	
				//Probability=1ulr-(unsigned long fract)Synapse_ypow;
				Probability=REAL_CONST(1.)-Synapse_ypow;
				rand_gen_union.ui=random_gen();
				if((unsigned long fract)Probability>rand_gen_union.ulf)
				{
					replenish=REAL_CONST(1.);
				}
				else
				{
					replenish=REAL_CONST(0.);
				}				
	
				//=========Update Variables==========//
				ANAvailHSR[j]=ANAvailHSR[j]+replenish+reprocessed-ejected;
				if(ANAvailHSR[j]<REAL_CONST(0.))
				{
					ANAvailHSR[j]=REAL_CONST(0.);
				}
				reuptakeandlost=(Synapse.rdt+Synapse.ldt)*ANCleftHSR[j];
				reuptake=Synapse.rdt*ANCleftHSR[j];
				if(reuptake>REAL_CONST(1.))
				{
					reuptake=REAL_CONST(1.);
				}

				ANCleftHSR[j]= ANCleftHSR[j]+ejected-reuptakeandlost;
				if(ANCleftHSR[j]<REAL_CONST(0.))
				{
					ANCleftHSR[j]=REAL_CONST(0.);
				}

				ANReproHSR[j]=ANReproHSR[j]+ reuptake-reprocessed;
				if(ANReproHSR[j]<REAL_CONST(0.))
				{
					ANReproHSR[j]=REAL_CONST(0.);
				}				
			}
			else
			{
				spikes=REAL_CONST(0.);
			}				
			//=========write value to SDRAM========//  
			//out_buffer[((j+4)*SEGSIZE)+i]  = ANReproHSR[j];//ANAvailHSR[j];//vrrhsr;//compare;//CaCurr_pow;//
			out_buffer[((j+4)*SEGSIZE)+i]  = spikes;
			

		}
		
		/*if(i==0)
		{
			#ifdef PROFILE
			  end_count_process = tc[T2_COUNT];
			  dtcm_profile_buffer[1+((seg_index-1)*3)]=start_count_process-end_count_process;
			#endif
		}*/

#ifdef LOOP_PROFILE
  int end_max_count_read = tc[T2_COUNT];
  io_printf (IO_BUF, "second half complete in %d ticks (sample %d, segment %d)\n",start_max_count_read-end_max_count_read,i,seg_index);
#endif	

	}
#ifdef PRINT
	/*io_printf (IO_BUF, "[core %d] segment %d processed and written to @ 0x%08x\n", coreID,seg_index,
				  (uint) sdramout_buffer);*/
	/*io_printf (IO_BUF, "[core %d] segment %d processed and written to @ 0x%08x - 0x%08x\n", coreID,seg_index,
					  (uint) &sdramout_buffer[segment_offset],(uint) &sdramout_buffer[segment_offset+(NUMFIBRES-1)*SEGSIZE+SEGSIZE-1]);*/
	io_printf (IO_BUF, "[core %d] segment %d processed\n",coreID,seg_index);
#endif			
	return segment_offset;
}

void transfer_handler(uint tid, uint ttag)
{
	if (ttag==DMA_READ)
	{
#ifdef PROFILE
  end_count_read = tc[T2_COUNT];
  dtcm_profile_buffer[seg_index*3]=(REAL)(start_count_read-end_count_read);
#ifdef PRINT
  io_printf (IO_BUF, "read complete in %d ticks\n",start_count_read-end_count_read);
#endif
#endif
		//increment segment index
		seg_index++;
		
		#ifdef PROFILE
		  start_count_process = tc[T2_COUNT];
		#endif
		
		//choose current buffers
		if(!read_switch && !write_switch)
		{
#ifdef PRINT 
io_printf (IO_BUF, "buff_b-->buff_x\n");
#endif
			index_x=process_chan(dtcm_buffer_x,dtcm_buffer_b);
		}
		else if(!read_switch && write_switch)
		{
#ifdef PRINT
io_printf (IO_BUF, "buff_b-->buff_y\n"); 
#endif
			index_y=process_chan(dtcm_buffer_y,dtcm_buffer_b);
		}
		else if(read_switch && !write_switch)
		{
#ifdef PRINT
io_printf (IO_BUF, "buff_a-->buff_x\n");
#endif	
			index_x=process_chan(dtcm_buffer_x,dtcm_buffer_a);

		}
		else
		{
#ifdef PRINT
io_printf (IO_BUF, "buff_a-->buff_y\n");
#endif
			index_y=process_chan(dtcm_buffer_y,dtcm_buffer_a);
		}		  
			
		#ifdef PROFILE
			  end_count_process = tc[T2_COUNT];
			  dtcm_profile_buffer[1+((seg_index-1)*3)]=start_count_process-end_count_process;
		#ifdef PRINT 
			io_printf (IO_BUF, "process complete in %d ticks (segment %d)\n",start_count_process-end_count_process,seg_index);
		#endif	
		#endif		
		
		spin1_trigger_user_event(NULL,NULL);
	}
	else if (ttag==DMA_WRITE)
	{
#ifdef PROFILE
  end_count_write = tc[T2_COUNT];
  dtcm_profile_buffer[2+((seg_index-1)*3)]=start_count_write-end_count_write;
#ifdef PRINT 
  io_printf (IO_BUF, "write complete in %d ticks\n",start_count_write-end_count_write);
#endif
#endif
		//flip write buffers
		write_switch=!write_switch;
	}
	else
	{
		#ifdef PRINT
		io_printf(IO_BUF,"[core %d] invalid %d DMA tag!\n",coreID,ttag);
		#endif
	}

}

void app_done ()
{
  // report simulation time
  io_printf (IO_BUF, "[core %d] simulation lasted %d ticks\n", coreID,
             spin1_get_simulation_time());

  //copy profile data
#ifdef PROFILE
  //io_printf (IO_BUF, "[core %d] saving profile data...\n", coreID);
	for (uint i=0;i<3*TOTAL_TICKS;i++)
	{
		profile_buffer[i]  = dtcm_profile_buffer[i];
	}
#endif
  
  // say goodbye
  io_printf (IO_BUF, "[core %d] stopping simulation\n", coreID);
  io_printf (IO_BUF, "[core %d] -------------------\n", coreID);
}

void c_main()
{
  // Get core and chip IDs
  coreID = spin1_get_core_id ();
  chipID = spin1_get_chip_id ();

  //set timer tick
  spin1_set_timer_tick (TIMER_TICK_PERIOD);

  //setup callbacks
  //process channel once data input has been read to DTCM
  spin1_callback_on (DMA_TRANSFER_DONE,transfer_handler,0);
  //reads from DMA to DTCM every tick
  spin1_callback_on (TIMER_TICK,data_read,-1);
  spin1_callback_on (USER_EVENT,data_write,0);

  app_init();

  spin1_start (SYNC_WAIT);
  
  app_done ();

}

