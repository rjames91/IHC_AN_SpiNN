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

#define TIMER_TICK_PERIOD  30000//150000//135000//116000//110000
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
REAL r_max;
REAL dt_spikes;
uint resamp_fac;
//REAL ln2_recip;

int start_count_process;
int end_count_process;
int start_count_read;
int end_count_read;
int start_count_write;
int end_count_write;

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
REAL CaCurrMSR[NUMMSR];

REAL ANCleftMSR[NUMMSR];
REAL ANAvailMSR[NUMMSR];
REAL ANReproMSR[NUMMSR];

REAL refrac_msr[NUMMSR];


startupVars generateStartupVars(void)
{
	startupVars out;
	REAL Gu0,kt0,IHCV,gmax,u0,u1,s1,s0,ga,gk,et,ek,rpc,ekp,gamma,beta,mICaCurr,
	gmaxcalsr,gmaxcamsr,gmaxcahsr,eca,tauca,CaCurrLSR,CaCurrMSR,CaCurrHSR,kt0LSR,kt0MSR,kt0HSR,
	ANCleftLSR,ANCleftMSR,ANCleftHSR,ANAvailLSR,ANAvailMSR,ANAvailHSR,ANReproLSR,ANReproMSR,
	ANReproHSR,z,power,y,x,l,r,MLSR,MMSR,MHSR;

	//======Calculate Gu0===========//
	gmax=6e-9;//(float)Cilia.Gmax/pow(2,30);
	u0=5e-9;//(float)Cilia.u0/pow(2,30);
	u1=1e-9;//(float)Cilia.u1/pow(2,30);
	s1=1e-9;//(float)Cilia.recips1/pow(2,1);
	s0=3e-8;//(float)Cilia.recips0/pow(2,1);
	ga=8e-10;//(float)Cilia.Ga/pow(2,31);
	Gu0=ga+gmax/(1+(exp(u0/s0)*(1+exp(u1/s1))));
	//======Calculate IHCVnow=========//
	gk=2e-8;//(float)Cilia.Gk/pow(2,30);
	et=0.1;//(float)Cilia.Et/pow(2,30);
	ek=-0.08;//(float)Cilia.Ek/pow(2,30);
	rpc=0.04;//(float)Cilia.Rpc/pow(2,30);
	ekp=ek+rpc*et;
	IHCV=(gk*ekp+Gu0*et)/(Gu0+gk);

	//=======Calculate mICaCurrent========//
	gamma=130;
	beta=400;
	mICaCurr=1/(1+exp(-gamma*IHCV)*(1/beta));

	//=======Calculate CaCurrent========//
	gmaxcalsr=1e-11;
	gmaxcamsr=7e-11;
	gmaxcahsr=5e-10;
	eca=0.066;
	tauca=1e-4;
	CaCurrLSR=gmaxcalsr*pow(mICaCurr,3)*(IHCV-eca)*tauca;
	CaCurrMSR=gmaxcamsr*pow(mICaCurr,3)*(IHCV-eca)*tauca;
	CaCurrHSR=gmaxcahsr*pow(mICaCurr,3)*(IHCV-eca)*tauca;

	//========Calculate AN==========//
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
	ANCleftMSR=(kt0MSR*y*MMSR)/(y*(l+r)+kt0MSR*l);
	ANAvailMSR=round((ANCleftMSR*(l+r))/kt0MSR);
	ANReproMSR=(ANCleftMSR*r)/x;

	out.IHCVnow0=IHCV;
	out.Ekp=ekp;
	out.mICaCurr0=mICaCurr;
	out.CaCurrMSR0=CaCurrMSR;
	out.ANCleftMSR0=ANCleftMSR;
	out.ANAvailMSR0=ANAvailMSR;
	out.ANReproMSR0=ANReproMSR;

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
	r_max=(REAL)RDM_MAX;
	//ln2_recip=REAL_CONST(1.)/(REAL)log(2);
	
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
					(uint)44100. *sizeof(REAL),
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
			sdramout_buffer[i]  = 0;//((uint)1<<31)-1;//0x5a5a5a5a;
		}
		
		for (uint i=0;i<(uint)44100.;i++)
		{
			sdramin_buffer[i]  = 0;			
		}
		
		for (uint i=0;i<3 * ((uint)44100./SEGSIZE);i++)
		{
			profile_buffer[i]  = 0;
		}
		
		for (uint i=0;i<3 * TOTAL_TICKS;i++)
		{
			dtcm_profile_buffer[i]  = 0;
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
	init_WELL1024a_simp();
	
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

	for(uint i=0;i<NUMMSR;i++)
	{
		CaCurrMSR[i]=startupValues.CaCurrMSR0;
		ANCleftMSR[i]=startupValues.ANCleftMSR0;
		ANAvailMSR[i]=startupValues.ANAvailMSR0;
		ANReproMSR[i]=startupValues.ANReproMSR0;
		refrac_msr[i]=0;
	}
	
	//=========initialise the pre synapse params========//
	preSyn.recipBeta=2.5e-3;
	preSyn.gamma=REAL_CONST(130.);
	preSyn.dtTauM=(dt/1e-4);
	preSyn.dtTauCa=(dt/1e-4);
	preSyn.Eca=0.066;
	preSyn.power=REAL_CONST(3.);
	preSyn.z=1e38;//reduced from 1e40 to fit within float range (then an additional *100 is added in code)

	preSyn.GmaxCa[0]=7e-11;
	preSyn.GmaxCa[1]=7e-11;

	preSyn.CaTh[0]=0.1;//pow(1e-13,preSyn.power);
	preSyn.CaTh[1]=0.1;//pow(1e-13,preSyn.power);

	//=======initialise the synapse params=======//

	Synapse.l=REAL_CONST(2580.);
	Synapse.y=REAL_CONST(10.);
	Synapse.x=REAL_CONST(30.);
	Synapse.r=REAL_CONST(6580.);
	Synapse.refrac_period=7.5e-4;

	Synapse.M[0]=REAL_CONST(14.);
	Synapse.M[1]=REAL_CONST(14.);
	  
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
	if(test_DMA == TRUE)// && ticks<TOTAL_TICKS)
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

uint process_chan(REAL *out_buffer,REAL *in_buffer) 
{  
	uint segment_offset=NUMFIBRES*SEGSIZE*(seg_index-1);
	uint i,j,k;
	
	REAL utconv;
	accum ex1;
	accum ex2;
	
	REAL Guconv;
	REAL mICaINF;
	REAL micapowconv;
	REAL ICa;
	REAL CaCurr_pow;
	REAL releaseProb_pow;
	REAL Synapse_ypow;
	REAL Synapse_xpow;
	REAL compare;
	REAL vrrmsr;
	
	REAL releaseProb;
	REAL M_q;
	REAL Probability;
	REAL ejected;
	REAL reprocessed;
	REAL replenish;
	REAL reuptakeandlost;
	REAL reuptake;
		
#ifdef PRINT
	io_printf (IO_BUF, "[core %d] segment %d (offset=%d) starting processing\n", coreID,seg_index,segment_offset);
#endif
	for(i=0;i<SEGSIZE;i++)
	{		
		//===========Apply Scaler============//
  	  	utconv=in_buffer[i] * Cilia.C;
  	  	
		//=========Apical Conductance========//	
		/*if(i==0)
		{
		#ifdef PROFILE
		  start_count_process = tc[T2_COUNT];
		#endif
		}*/

		ex1=expk((accum)(-(utconv-Cilia.u1)*Cilia.recips1));
		ex2=expk((accum)(-(utconv-Cilia.u0)*Cilia.recips0));
  	  	
  	  	Guconv=Cilia.Ga + (Cilia.Gmax/(REAL_CONST(1.)+(REAL)ex2*(REAL_CONST(1.)+(REAL)ex1)));
		//Guconv=Cilia.Ga + (Cilia.Gmax/(REAL_CONST(1.)+exp(-(utconv-Cilia.u0)*Cilia.recips0)*(REAL_CONST(1.)+exp(-(utconv-Cilia.u1)*Cilia.recips1))));
	/*	if(i==0)
		{
		#ifdef PROFILE
		  end_count_process = tc[T2_COUNT];
		  dtcm_profile_buffer[1+((seg_index-1)*3)]=start_count_process-end_count_process;
		#endif
		}*/  
		//========Receptor Potential=========//
		IHCVnow+=((-Guconv*(IHCVnow-Cilia.Et))-(Cilia.Gk*(IHCVnow-startupValues.Ekp)))*Cilia.dtCap;
	
		//================mICa===============//
		ex1=expk((accum)-preSyn.gamma*IHCVnow);
		mICaINF=REAL_CONST(1.)/(REAL_CONST(1.)+(REAL)ex1*preSyn.recipBeta);
		//mICaINF=REAL_CONST(1.)/(REAL_CONST(1.)+exp(-preSyn.gamma*IHCVnow)*preSyn.recipBeta);
		mICaCurr+=(mICaINF-mICaCurr)*preSyn.dtTauM;

		//================ICa================//
		//calculate values for each fibre	
		//micapowconv=pow(mICaCurr,preSyn.power);
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

		//=======MSR Fibres===========//
		for (j=0;j<NUMMSR;j++)
		{
			ICa=preSyn.GmaxCa[j]*micapowconv*(IHCVnow-preSyn.Eca);
			CaCurrMSR[j]+=(-ICa-CaCurrMSR[j])*preSyn.dtTauCa;
			if(i%resamp_fac==0)
			{
				CaCurr_pow=CaCurrMSR[j];
				for(k=0;k<preSyn.power-1;k++)
				{			
					CaCurr_pow=CaCurr_pow*CaCurrMSR[j];
				}
				compare=max(preSyn.z*CaCurr_pow-preSyn.CaTh[j]),0);
							
				//compare=max(100*(pow(CaCurrMSR[j],preSyn.power)-preSyn.CaTh[j+2]),0);
				vrrmsr=REAL_CONST(100.)*compare;
				//saturate vrr
				if(vrrmsr>max_rate)
				{
					vrrmsr=max_rate;
				}
				//=====Release Probability=======//
				releaseProb=vrrmsr*dt_spikes;
				M_q=Synapse.M[j]-ANAvailMSR[j];
				if(M_q<0)
				{
					M_q=0;
				}
	
				//===========Ejected===========//
			
				releaseProb_pow=REAL_CONST(1.);
				for(k=0;k<ANAvailMSR[j];k++)
				{			
					releaseProb_pow=releaseProb_pow*(REAL_CONST(1.)-releaseProb);
				}
				Probability=REAL_CONST(1.)-releaseProb_pow;			
				//Probability=1-pow(1-releaseProb,ANAvailMSR[j]);
				if (refrac_msr[j]>0)
				{
					//decrement refrac counter until it's zero
					refrac_msr[j]--;
					ejected=REAL_CONST(0.);
				}
	
				else if(Probability>((REAL)random_gen()/r_max))
				{
					ejected=REAL_CONST(1.);
					//refrac_msr[j]=round(Synapse.refrac_period/dt + (((REAL)random_gen()/r_max)*Synapse.refrac_period)/dt);
					refrac_msr[j]=(uint)(Synapse.refrac_period/dt_spikes + (((REAL)random_gen()/r_max)*Synapse.refrac_period)/dt_spikes);
				}
				else
				{
					ejected=REAL_CONST(0.);
				}
	
				//========Reprocessed==========//		
				Synapse_xpow=REAL_CONST(1.);
				Synapse_ypow=REAL_CONST(1.);
				for(k=0;k<M_q;k++)
				{			
					Synapse_xpow=Synapse_xpow*(REAL_CONST(1.)-(ANReproMSR[j]*Synapse.x*dt_spikes));
					Synapse_ypow=Synapse_ypow*(REAL_CONST(1.)-(Synapse.y*dt_spikes));
				}
	
				Probability=REAL_CONST(1.)-Synapse_xpow;			
				//Probability=1-pow(1-(ANReproMSR[j]*Synapse.x*dt),M_q);
				if(Probability>((REAL)random_gen()/r_max))
				{
					reprocessed=REAL_CONST(1.);
				}
				else
				{
					reprocessed=REAL_CONST(0.);
				}
	
				//=========Replenish============//
				
				Probability=REAL_CONST(1.)-Synapse_ypow;
				//Probability=1-pow(1-Synapse.y*dt,M_q);
				if(Probability>((REAL)random_gen()/r_max))
				{
					replenish=REAL_CONST(1.);
				}
				else
				{
					replenish=REAL_CONST(0.);
				}
	
				//========Update Variables=======//
				ANAvailMSR[j]=ANAvailMSR[j]+replenish-ejected+reprocessed;
				reuptakeandlost=((Synapse.r*dt_spikes)+(Synapse.l*dt_spikes))*ANCleftMSR[j];
				reuptake=Synapse.r*dt_spikes*ANCleftMSR[j];
				ANCleftMSR[j]= ANCleftMSR[j]+ejected-reuptakeandlost;
				ANReproMSR[j]=ANReproMSR[j]+ reuptake-reprocessed;
			}
			else
			{
				ejected=REAL_CONST(0.);
			}	
			//=======write value to SDRAM=======//  
			//out_buffer[((j)*SEGSIZE)+i] = vrrmsr;
			out_buffer[(j*SEGSIZE)+i] =ejected;
			
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

