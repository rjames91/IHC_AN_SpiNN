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
#include "random.h"
#include "stdfix-exp.h"

#define TIMER_TICK_PERIOD (3000 * NUMFIBRES) //TODO: make this dependent on numfibres
#define TOTAL_TICKS 200//173//197       
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
REAL r_max_recip;
REAL dt_spikes;
uint resamp_fac;
uint spike_seg_size;
SEED_TYPE local_seed;

//uint seed_selection[SEED_SEL_SIZE];//TODO:this needs to be moved to SDRAM

int start_count_process;
int end_count_process;
int start_count_read;
int end_count_read;
int start_count_write;
int end_count_write;
uint refrac[NUMFIBRES];

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
REAL CaCurr[NUMFIBRES];
REAL ANCleft[NUMFIBRES];
REAL ANAvail[NUMFIBRES];
REAL ANRepro[NUMFIBRES];

startupVars generateStartupVars(void)
{
	startupVars out;
	REAL Gu0,kt0,IHCV,gmax,u0,u1,s1,s0,ga,gk,et,ek,rpc,ekp,gamma,beta,mICaCurr,
	gmaxcalsr,gmaxcamsr,gmaxcahsr,eca,taucalsr,taucamsr,taucahsr,CaCurrLSR,CaCurrMSR,CaCurrHSR,kt0LSR,kt0MSR,kt0HSR,
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
	gamma=100.;
	beta=400.;
	mICaCurr=1./(1.+exp(-gamma*IHCV)*(1./beta));

	//=======Calculate CaCurrent========//
	gmaxcalsr=14e-9;
	gmaxcamsr=14e-9;
	gmaxcahsr=14e-9;
	eca=0.066;
	taucalsr=30e-6;
	taucamsr=0;
	taucahsr=80e-6;
	CaCurrLSR=((gmaxcalsr*pow(mICaCurr,3))*(IHCV-eca))*taucalsr;
	CaCurrMSR=((gmaxcamsr*pow(mICaCurr,3))*(IHCV-eca))*taucamsr;
	CaCurrHSR=((gmaxcahsr*pow(mICaCurr,3))*(IHCV-eca))*taucahsr;

	//========Calculate AN==========//
	z=1.4142e21;//sqrt of 2e42 to fit as single float type (then multiplied twice)
	power=3;
	y=6;
	x=60;
	MLSR=12;
	MMSR=12;
	MHSR=12;
	l=250;
	r=500;

	kt0LSR=-z*pow(CaCurrLSR,power)*z;//added *z for single float z
	kt0MSR=-z*pow(CaCurrMSR,power)*z;
	kt0HSR=-z*pow(CaCurrHSR,power)*z;
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
	Fs=SAMPLING_FREQUENCY;
	dt=(1.0/Fs);
	resamp_fac=RESAMPLE_FACTOR;
	spike_seg_size=SEGSIZE/resamp_fac;
	max_rate=Fs/(REAL)resamp_fac;
	dt_spikes=(REAL)resamp_fac*dt;
	seg_index=0;
	read_switch=0;
	write_switch=0;
	r_max_recip=REAL_CONST(1.)/(REAL)RDM_MAX;
	
	/* say hello */
	
	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);
	
	// Allocate buffers somewhere in SDRAM
	
	//output results buffer, 10x fibre outputs
	sdramout_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					 MAX_SIGNAL_S*NUMFIBRES * (uint)44100./resamp_fac * sizeof(REAL),
					 coreID,
					 ALLOC_LOCK);	

	sdramin_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					MAX_SIGNAL_S*(uint)44100. *sizeof(REAL),
					coreID|32,
					ALLOC_LOCK);
	
	profile_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					3 * ((uint)44100./SEGSIZE) *sizeof(REAL),
					coreID|64,
					ALLOC_LOCK);
	
	// and a buffer in DTCM
	
	dtcm_buffer_a = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_b = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_x = (REAL *) sark_alloc (spike_seg_size*NUMFIBRES, sizeof(REAL));
	dtcm_buffer_y = (REAL *) sark_alloc (spike_seg_size*NUMFIBRES, sizeof(REAL));
	dtcm_profile_buffer = (REAL *) sark_alloc (3*TOTAL_TICKS, sizeof(REAL));
	
	if (dtcm_buffer_a == NULL ||dtcm_buffer_b == NULL ||dtcm_buffer_x == NULL ||dtcm_buffer_y == NULL 
			||  sdramout_buffer == NULL || sdramin_buffer == NULL || profile_buffer == NULL || dtcm_profile_buffer == NULL)
	/*if (sdramout_buffer == NULL || sdramin_buffer == NULL || dtcm_buffer_y == NULL 
			|| dtcm_buffer_a == NULL || dtcm_buffer_b == NULL || dtcm_profile_buffer == NULL ||dtcm_buffer_x == NULL)*/
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
		for (uint i = 0; i < spike_seg_size*NUMFIBRES; i++)
		{
			dtcm_buffer_x[i]   = 0;
			dtcm_buffer_y[i]   = 0;
		}
	
		for (uint i=0;i<MAX_SIGNAL_S* NUMFIBRES * ((uint)44100./resamp_fac);i++)
		{
			sdramout_buffer[i]  = 0;
		}
		for (uint i=0;i<MAX_SIGNAL_S * (uint)44100.;i++)
		{
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
	Cilia.C=REAL_CONST(0.08);
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
		CaCurr[i]=startupValues.CaCurrLSR0;
		ANCleft[i]=startupValues.ANCleftLSR0;
		ANAvail[i]=startupValues.ANAvailLSR0;
		ANRepro[i]=startupValues.ANReproLSR0;
		refrac[i]=0;
		preSyn.GmaxCa[i]=14e-9;
		preSyn.TauCa[i]=30e-6;
		Synapse.M[i]=REAL_CONST(12.);
	}
	for(uint i=0;i<NUMMSR;i++)
	{
		CaCurr[i+NUMLSR]=startupValues.CaCurrMSR0;
		ANCleft[i+NUMLSR]=startupValues.ANCleftMSR0;
		ANAvail[i+NUMLSR]=startupValues.ANAvailMSR0;
		ANRepro[i+NUMLSR]=startupValues.ANReproMSR0;
		refrac[i+NUMLSR]=0;
		preSyn.GmaxCa[i+NUMLSR]=14e-9;
		preSyn.TauCa[i+NUMLSR]=0;
		Synapse.M[i+NUMLSR]=REAL_CONST(12.);
	}
	for(uint i=0;i<NUMHSR;i++) 
	{
		CaCurr[i+NUMLSR+NUMMSR]=startupValues.CaCurrHSR0;
		ANCleft[i+NUMLSR+NUMMSR]=startupValues.ANCleftHSR0;
		ANAvail[i+NUMLSR+NUMMSR]=startupValues.ANAvailHSR0;
		ANRepro[i+NUMLSR+NUMMSR]=startupValues.ANReproHSR0;
		refrac[i+NUMLSR+NUMMSR]=0;
		preSyn.GmaxCa[i+NUMLSR+NUMMSR]=14e-9;
		preSyn.TauCa[i+NUMLSR+NUMMSR]=80e-6;
		Synapse.M[i+NUMLSR+NUMMSR]=REAL_CONST(12.);
	}
	
	//=========initialise the pre synapse params========//
	preSyn.recipBeta=2.5e-3;
	preSyn.gamma=REAL_CONST(100.);
	preSyn.dtTauM=(dt/5e-5);
	//preSyn.dtTauCa=(dt/1e-4);
	preSyn.Eca=0.066;
	preSyn.power=REAL_CONST(3.);
	preSyn.z=1.4142e21;//sqrt of 2e42 to fit as single float type (then multiplied twice)

	//=======initialise the synapse params=======//

	Synapse.ldt=REAL_CONST(250.)*dt_spikes;
	Synapse.ydt=REAL_CONST(6.)*dt_spikes;
	if(Synapse.ydt>REAL_CONST(1.))
	{
		Synapse.ydt=REAL_CONST(1.);
	}
	Synapse.xdt=REAL_CONST(60.)*dt_spikes;
	Synapse.rdt=REAL_CONST(500.)*dt_spikes;
	Synapse.refrac_period=((7.5e-4)/dt_spikes);

	
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
	uint segment_offset=NUMFIBRES*spike_seg_size*(seg_index-1);
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
	REAL releaseProb_pow;
	REAL Repro_rate;
	REAL Synapse_ypow;
	REAL Synapse_xpow;
	REAL compare;
	REAL vrr;	
	REAL releaseProb;
	REAL M_q;
	REAL Probability;
	REAL ejected;
	REAL spikes;
	REAL reprocessed;
	REAL replenish;
	REAL reuptakeandlost;
	REAL reuptake;
	
	uint si=0;
		
#ifdef PRINT
	io_printf (IO_BUF, "[core %d] segment %d (offset=%d) starting processing\n", coreID,seg_index,segment_offset);
#endif
	for(i=0;i<SEGSIZE;i++)
	{	
		  
	/*		if(i==0)
			{
			#ifdef PROFILE
			  start_count_process = tc[T2_COUNT];
			#endif
			}*/
		
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
		/*if(mICaCurr<REAL_CONST(0.))
		{
			mICaCurr=REAL_CONST(0.);
		}*/

		//================ICa================//
		//calculate values for each fibre	
		micapowconv=mICaCurr;
		for(k=0;k<preSyn.power-1;k++)
		{
			micapowconv=micapowconv*mICaCurr;
		}
		
		
	/*	if(i==0)
		{
			#ifdef PROFILE
			  end_count_process = tc[T2_COUNT];
			  dtcm_profile_buffer[1+((seg_index-1)*3)]=start_count_process-end_count_process;
			#endif
		}*/
		

#ifdef LOOP_PROFILE
  int end_apical_count_read = tc[T2_COUNT];
  io_printf (IO_BUF, "first half complete in %d ticks (sample %d, segment %d)\n",start_apical_count_read-end_apical_count_read,i,seg_index);
#endif	
  
#ifdef LOOP_PROFILE
  int start_max_count_read = tc[T2_COUNT];
#endif

		if(i%resamp_fac==0)
		{
			si++;
		}

		//============Fibres=============//
		for (j=0;j<NUMFIBRES;j++)
		{
			//======Synaptic Ca========//
			ICa=preSyn.GmaxCa[j]*micapowconv*(IHCVnow-preSyn.Eca);
			//CaCurr[j]+=(ICa-CaCurr[j])*preSyn.dtTauCa;
			CaCurr[j]+=ICa*dt - (CaCurr[j]*dt)/preSyn.TauCa[j];
			/*if(CaCurr[j]>REAL_CONST(0.))
			{
				CaCurr[j]=REAL_CONST(0.);
			}*/
			//invert Ca
			pos_CaCurr=REAL_CONST(-1.)*CaCurr[j];

			if(i%resamp_fac==0)
			{	
				//=====Vesicle Release Rate=====//
				//CaCurr_pow=pos_CaCurr;
				CaCurr_pow=pos_CaCurr*preSyn.z;//first sqrt(z) multiply to prevent underflow CaCurr_pow
				for(k=0;k<(uint)preSyn.power-1;k++)
				{			
					CaCurr_pow=CaCurr_pow*pos_CaCurr;
				}
				//compare=max((preSyn.z*CaCurr_pow-preSyn.CaTh[j]),REAL_CONST(0.));				
				vrr = preSyn.z*CaCurr_pow;//second sqrt(z) multiply
				//vrr=compare*preSyn.z;
				//saturate vrr
				/*if(vrr>max_rate)
				{
					vrr=max_rate;
				}*/
		
				//=====Release Probability=======//
				releaseProb=vrr*dt_spikes;
				M_q=Synapse.M[j]-ANAvail[j];
				if(M_q<REAL_CONST(0.))
				{
					M_q=REAL_CONST(0.);
				}
	
				//===========Ejected============//
				releaseProb_pow=REAL_CONST(1.);
				for(k=0;k<(uint)ANAvail[j];k++)
				{			
					releaseProb_pow=releaseProb_pow*(REAL_CONST(1.)-releaseProb);
				}
	
				Probability=REAL_CONST(1.)-releaseProb_pow;
			
				if (refrac[j]>0)
				{
					refrac[j]--;
				}	
								
				if(Probability>((REAL)random_gen()*r_max_recip))
				{
					ejected=REAL_CONST(1.);
					if (refrac[j]<=0)
					{
						spikes=REAL_CONST(1.);
						refrac[j]=(uint)(Synapse.refrac_period + (((REAL)random_gen()*r_max_recip)*Synapse.refrac_period)+REAL_CONST(0.5));
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
				Repro_rate=ANRepro[j]*Synapse.xdt;
				/*if(Repro_rate>REAL_CONST(1.))
				{
					Repro_rate=REAL_CONST(1.);
				}*/
	
				Synapse_xpow=REAL_CONST(1.);
				Synapse_ypow=REAL_CONST(1.);
				for(k=0;k<(uint)M_q;k++)
				{			
					Synapse_xpow=Synapse_xpow*(REAL_CONST(1.)-Repro_rate);
					Synapse_ypow=Synapse_ypow*(REAL_CONST(1.)-Synapse.ydt);
				}

				Probability=REAL_CONST(1.)-Synapse_xpow;			
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
				}
	
				//==========Update Variables=========//
				ANAvail[j]=ANAvail[j]+replenish+reprocessed-ejected;
				/*if(ANAvail[j]<REAL_CONST(0.))
				{
					ANAvail[j]=REAL_CONST(0.);
				}*/
				reuptakeandlost=(Synapse.rdt+Synapse.ldt)*ANCleft[j];
				reuptake=Synapse.rdt*ANCleft[j];
				/*if(reuptake>REAL_CONST(1.))
				{
					reuptake=REAL_CONST(1.);
				}*/

				ANCleft[j]= ANCleft[j]+ejected-reuptakeandlost;
				/*if(ANCleft[j]<REAL_CONST(0.))
				{
					ANCleft[j]=REAL_CONST(0.);
				}*/
				ANRepro[j]=ANRepro[j]+ reuptake-reprocessed;
				/*if(ANRepro[j]<REAL_CONST(0.))
				{
					ANRepro[j]=REAL_CONST(0.);
				}*/
				
				//=======write value to SDRAM========//
				//out_buffer[(j*SEGSIZE)+i]  = ANReproLSR[j];//ANAvailLSR[j];//vrrlsr;// compare;//CaCurr_pow;//
				out_buffer[(j*spike_seg_size)+(si-1)]  = spikes;//pos_CaCurr;//ICa;//in_buffer[i];//
				
				//io_printf (IO_BUF, "[core %d] index=%d\n", coreID,(j*spike_seg_size)+(si-1));
			}
			
		}


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

