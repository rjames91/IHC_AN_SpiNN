/*
 ============================================================================
 Name        : IHC_AN_softfloat_LSR.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : LSR fibre model for IHC/AN simulation
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
#define TOTAL_TICKS 250      
#define PROFILE
//#define PRINT

//=========GLOBAL VARIABLES============//
REAL Fs,dt;
uint coreID;
uint chipID;
uint test_DMA;
uint seg_index;
uint rx_index;
uint read_switch;
uint write_switch;
uint processing;
uint index_x;
uint index_y;
REAL r_max;
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
REAL refrac_lsr[NUMLSR];
REAL refrac_msr[NUMMSR];
REAL refrac_hsr[NUMHSR];


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
	Fs=REAL_CONST(44100.);
	dt=(1.0/Fs);
	seg_index=0;
	rx_index=0;
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
					 NUMFIBRES * (uint)Fs * sizeof(REAL),
					 coreID,
					 ALLOC_LOCK);	
	
	profile_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					3 * ((uint)Fs/SEGSIZE) *sizeof(REAL),
					coreID|64,
					ALLOC_LOCK);
	
	// and a buffer in DTCM
	dtcm_buffer_a = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_b = (REAL *) sark_alloc (SEGSIZE, sizeof(REAL));
	dtcm_buffer_x = (REAL *) sark_alloc (SEGSIZE*NUMFIBRES, sizeof(REAL));
	dtcm_buffer_y = (REAL *) sark_alloc (SEGSIZE*NUMFIBRES, sizeof(REAL));
	dtcm_profile_buffer = (REAL *) sark_alloc (3*TOTAL_TICKS, sizeof(REAL));
	
	if (dtcm_buffer_x == NULL ||dtcm_buffer_y == NULL || dtcm_buffer_a == NULL ||dtcm_buffer_b == NULL 
			||  sdramout_buffer == NULL || profile_buffer == NULL || dtcm_profile_buffer == NULL)
	{
		test_DMA = FALSE;
		io_printf (IO_BUF, "[core %d] error - cannot allocate buffer\n", coreID);
	}
	else
	{
		test_DMA = TRUE;
		// initialize sections of DTCM, system RAM and SDRAM
		for (uint i = 0; i < SEGSIZE*NUMFIBRES; i++)
		{
			dtcm_buffer_x[i]   = 0;
			dtcm_buffer_y[i]   = 0;
		}
	
		for (uint i=0;i<NUMFIBRES * (uint)Fs;i++)
		{
			sdramout_buffer[i]  = 0;//((uint)1<<31)-1;//0x5a5a5a5a;
		}
		
		for (uint i=0;i<3 * TOTAL_TICKS;i++)
		{
			dtcm_profile_buffer[i]  = 0;
			profile_buffer[i]  = 0;
		}
		
		io_printf (IO_BUF, "[core %d] sdram out buffer @ 0x%08x\n", coreID,
				   (uint) sdramout_buffer);
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

	preSyn.CaTh[0]=pow(1e-14,preSyn.power);
	preSyn.CaTh[1]=pow(1e-14,preSyn.power);
	preSyn.CaTh[2]=pow(1e-13,preSyn.power);
	preSyn.CaTh[3]=pow(1e-13,preSyn.power);
	preSyn.CaTh[4]=pow(1.43e-13,preSyn.power);
	preSyn.CaTh[5]=pow(1.43e-13,preSyn.power);
	preSyn.CaTh[6]=pow(1.43e-13,preSyn.power);
	preSyn.CaTh[7]=pow(1.43e-13,preSyn.power);
	preSyn.CaTh[8]=pow(1.43e-13,preSyn.power);
	preSyn.CaTh[9]=pow(1.43e-13,preSyn.power);

	//=======initialise the synapse params=======//

	Synapse.l=REAL_CONST(2580.);
	Synapse.y=REAL_CONST(10.);
	Synapse.x=REAL_CONST(30.);
	Synapse.r=REAL_CONST(6580.);
	Synapse.refrac_period=7.5e-4;

	Synapse.M[0]=REAL_CONST(13.);
	Synapse.M[1]=REAL_CONST(13.);
	Synapse.M[2]=REAL_CONST(14.);
	Synapse.M[3]=REAL_CONST(14.);
	Synapse.M[4]=REAL_CONST(8.);
	Synapse.M[5]=REAL_CONST(8.);
	Synapse.M[6]=REAL_CONST(8.);
	Synapse.M[7]=REAL_CONST(8.);
	Synapse.M[8]=REAL_CONST(8.);
	Synapse.M[9]=REAL_CONST(8.); 

	  
#ifdef PROFILE
    // configure timer 2 for profiling
    // enabled, free running, interrupt disabled, no pre-scale, 32 bit, free-running mode
    tc[T2_CONTROL] = TIMER2_CONF;
#endif
    
}
void data_write(uint null_a, uint null_b)
{
	REAL *dtcm_buffer_out;
	
	if(test_DMA == TRUE)
	{
		if(!write_switch)
		{
			dtcm_buffer_out=dtcm_buffer_x;
#ifdef PRINT	
			io_printf (IO_BUF, "buff_x write\n");
#endif
		}
		else
		{
			dtcm_buffer_out=dtcm_buffer_y;
#ifdef PRINT
			io_printf (IO_BUF, "buff_y write\n");
#endif
		}
#ifdef PROFILE
  start_count_write = tc[T2_COUNT];
#endif
		spin1_dma_transfer(DMA_WRITE,&sdramout_buffer[SEGSIZE*(seg_index-1)],dtcm_buffer_out,DMA_WRITE,
		  						NUMFIBRES*SEGSIZE*sizeof(REAL));
#ifdef PRINT
		io_printf (IO_BUF, "[core %d] segment %d written to @ 0x%08x - 0x%08x\n", coreID,seg_index,
							  (uint) &sdramout_buffer[SEGSIZE*(seg_index-1)],(uint) &sdramout_buffer[SEGSIZE*(seg_index-1)+(NUMFIBRES-1)*SEGSIZE+SEGSIZE-1]);
#endif
	}
}

void process_chan(REAL *out_buffer,REAL *in_buffer) 
{  
	uint i,j,k;

	REAL ICa;
	REAL CaCurr_pow;
	REAL releaseProb_pow;
	REAL Synapse_ypow;
	REAL Synapse_xpow;
	REAL compare;
	REAL vrrlsr;
	REAL vrrmsr;
	REAL vrrhsr;
	
	REAL releaseProb;
	REAL M_q;
	REAL Probability;
	REAL ejected;
	REAL reprocessed;
	REAL replenish;
	REAL reuptakeandlost;
	REAL reuptake;
		
#ifdef PRINT
	//io_printf (IO_BUF, "[core %d] segment %d starting processing\n", coreID,seg_index);
#endif
	for(i=0;i<SEGSIZE;i++)
	{	
			
			//======Synaptic Ca========//
			ICa=preSyn.GmaxCa[j]*in_buffer[i]*(IHCVnow-preSyn.Eca);
			CaCurrLSR[j]+=(-ICa-CaCurrLSR[j])*preSyn.dtTauCa;
		
			//=====Vesicle Release Rate=====//
			CaCurr_pow=CaCurrLSR[j];
			for(k=0;k<preSyn.power-1;k++)
			{			
				CaCurr_pow=CaCurr_pow*CaCurrLSR[j];
			}
			compare=max(REAL_CONST(100.)*(CaCurr_pow-preSyn.CaTh[j]),0);											
			//compare=max(REAL_CONST(100.)*(pow(CaCurrLSR[j],preSyn.power)-preSyn.CaTh[j]),0);		
			vrrlsr=preSyn.z*compare;
			//saturate vrr
			if(vrrlsr>Fs)
			{
				vrrlsr=Fs;
			}
		
			//=====Release Probability=======//
			releaseProb=vrrlsr*dt;
			M_q=Synapse.M[j]-ANAvailLSR[j];
			if(M_q<0)
			{
				M_q=0;
			}
		
			//===========Ejected============//
			releaseProb_pow=REAL_CONST(1.);
			for(k=0;k<ANAvailLSR[j];k++)
			{			
				releaseProb_pow=releaseProb_pow*(REAL_CONST(1.)-releaseProb);
			}
		
			Probability=REAL_CONST(1.)-releaseProb_pow;
		
			//Probability=REAL_CONST(1.)-pow(REAL_CONST(1.)-releaseProb,ANAvailLSR[j]);
			if (refrac_lsr[j]>0)
			{
				//decrement refrac counter until it's zero
				refrac_lsr[j]--;
				ejected=REAL_CONST(0.);
			}
		
			else if(Probability>((REAL)random_gen()/r_max))
			{
				ejected=REAL_CONST(1.);
				//refrac_lsr[j]=round(Synapse.refrac_period/dt + (((REAL)random_gen()/r_max)*Synapse.refrac_period)/dt);
				refrac_lsr[j]=(uint)(Synapse.refrac_period/dt + (((REAL)random_gen()/r_max)*Synapse.refrac_period)/dt);
			}
			else
			{
				ejected=REAL_CONST(0.);
			}
		
			//=========Reprocessed=========//
		
			Synapse_xpow=REAL_CONST(1.);
			Synapse_ypow=REAL_CONST(1.);
			for(k=0;k<M_q;k++)
			{			
				Synapse_xpow=Synapse_xpow*(REAL_CONST(1.)-(ANReproLSR[j]*Synapse.x*dt));
				Synapse_ypow=Synapse_ypow*(REAL_CONST(1.)-Synapse.y*dt);
			}
		
			Probability=REAL_CONST(1.)-Synapse_xpow;			
			//Probability=1-pow(1-(ANReproLSR[j]*Synapse.x*dt),M_q);			
			if(Probability>((REAL)random_gen()/r_max))
			{
				reprocessed=REAL_CONST(1.);
			}
			else
			{
				reprocessed=REAL_CONST(0.);
			}
		
			//========Replenish==========//
		
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
		
			//==========Update Variables=========//
			ANAvailLSR[j]=ANAvailLSR[j]+replenish-ejected+reprocessed;
			reuptakeandlost=((Synapse.r*dt)+(Synapse.l*dt))*ANCleftLSR[j];
			reuptake=Synapse.r*dt*ANCleftLSR[j];
			ANCleftLSR[j]= ANCleftLSR[j]+ejected-reuptakeandlost;
			ANReproLSR[j]=ANReproLSR[j]+ reuptake-reprocessed;
		
			//=======write value to SDRAM========//
			//out_buffer[(j*SEGSIZE)+i]  = vrrlsr;
			out_buffer[i]  = ejected;		
	}

#ifdef PRINT
	/*io_printf (IO_BUF, "[core %d] segment %d processed and written to @ 0x%08x\n", coreID,seg_index,
				  (uint) sdramout_buffer);*/
	/*io_printf (IO_BUF, "[core %d] segment %d processed and written to @ 0x%08x - 0x%08x\n", coreID,seg_index,
					  (uint) &sdramout_buffer[segment_offset],(uint) &sdramout_buffer[segment_offset+(NUMFIBRES-1)*SEGSIZE+SEGSIZE-1]);*/
	//io_printf (IO_BUF, "[core %d] segment %d processed\n",coreID,seg_index);
#endif			
}

void dma_transfer_handler(uint tid, uint ttag)
{
	if (ttag==DMA_WRITE)
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

void packet_rx(uint key, uint data)
{
	//TODO: check key is from relevent location?
	//increment rx index
	//io_printf (IO_BUF, "MC packet received!\n");
	rx_index++;
	
	#ifdef PROFILE
	  start_count_process = tc[T2_COUNT];
	#endif
	
	dtcm_buffer_a[rx_index-1]=(REAL)data;
  
/*	//choose current buffers
	if(!write_switch)
	{		
		//process_chan(dtcm_buffer_x,(REAL)data);
	}
	
	else
	{
		process_chan(dtcm_buffer_y,(REAL)data);
		//dtcm_buffer_y[rx_index-1]=(REAL)data;
	}		  
*/		
	#ifdef PROFILE
		  end_count_process = tc[T2_COUNT];
		  dtcm_profile_buffer[1+((seg_index-1)*3)]=start_count_process-end_count_process;
	#ifdef PRINT 
		//io_printf (IO_BUF, "process complete in %d ticks (segment %d)\n",start_count_process-end_count_process,seg_index);
	#endif	
	#endif			
	if(rx_index >= SEGSIZE)	
	{		
		seg_index ++;
		rx_index=0;
		
		if(!write_switch)
		{		
			process_chan(dtcm_buffer_x,dtcm_buffer_a);
		}
		
		else
		{
			process_chan(dtcm_buffer_y,dtcm_buffer_a);
		}	
		
		spin1_trigger_user_event(NULL,NULL);
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

void time_check(uint ticks, uint null)
{
	// stop if desired number of ticks reached
	if (ticks > TOTAL_TICKS) 
	{
		io_printf (IO_BUF, "spinn_exit\n");
		spin1_exit (0); 
	}
	
}

void c_main()
{
  // Get core and chip IDs
  coreID = spin1_get_core_id ();
  chipID = spin1_get_chip_id ();

  //set timer tick
  spin1_set_timer_tick (TIMER_TICK_PERIOD);

  //setup callbacks
  spin1_callback_on (DMA_TRANSFER_DONE,dma_transfer_handler,1);
  spin1_callback_on (MCPL_PACKET_RECEIVED,packet_rx,-1);  
  spin1_callback_on (USER_EVENT,data_write,0);
  spin1_callback_on (TIMER_TICK,time_check,1);


  app_init();

  spin1_start (SYNC_WAIT);
  
  app_done ();

}

