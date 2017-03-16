/*
 ============================================================================
 Name        : spinn_IHC_AN_fixed.c
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
#include "spinn_IHC_AN_fixed.h"
#include "spin1_api.h"
#include "math.h"

#define TIMER_TICK_PERIOD  150000//110000
#define TOTAL_TICKS        180
#define PROFILE
//#define PRINT

//=========GLOBAL VARIABLES============//
REAL Fs,dt;
REAL testInput[SEGSIZE];
uint coreID;
uint chipID;
uint test_DMA;
uint seg_index;
uint read_switch;
uint write_switch;
uint processing;
uint index_x;
uint index_y;

float fixer=pow(2,FRACBITS);

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


startupVars generateStartupVars(void)
{
	startupVars out;
	float Gu0,kt0,IHCV,gmax,u0,u1,s1,s0,ga,gk,et,ek,rpc,ekp,gamma,beta,mICaCurr,
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

	//TODO: may need to cast
	out.IHCVnow0=IHCV*fixer;
	out.Ekp=ekp*fixer;
	out.mICaCurr0=mICaCurr*fixer;
	out.CaCurrLSR0=CaCurrLSR*fixer;
	out.CaCurrMSR0=CaCurrMSR*fixer;
	out.CaCurrHSR0=CaCurrHSR*fixer;
	out.ANCleftLSR0=ANCleftLSR*fixer;
	out.ANCleftMSR0=ANCleftMSR*fixer;
	out.ANCleftHSR0=ANCleftHSR*fixer;
	out.ANAvailLSR0=ANAvailLSR*fixer;
	out.ANAvailMSR0=ANAvailMSR*fixer;
	out.ANAvailHSR0=ANAvailHSR*fixer;
	out.ANReproLSR0=ANReproLSR*fixer;
	out.ANReproMSR0=ANReproMSR*fixer;
	out.ANReproHSR0=ANReproHSR*fixer;

	return out;
}

//application initialisation
void app_init(void)
{
	Fs=REAL_CONST(44100.);
	dt=(1.0/Fs);
	seg_index=0;
	read_switch=0;
	write_switch=0;
	/* say hello */
	
	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);
	
	// Allocate buffers somewhere in SDRAM
	
	//output results buffer, 10x fibre outputs
	sdramout_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					 NUMFIBRES * (uint)Fs * sizeof(REAL),
					 coreID,
					 ALLOC_LOCK);	

	sdramin_buffer = (REAL *) sark_xalloc (sv->sdram_heap,
					NUMFIBRES * (uint)Fs *sizeof(REAL),	
					coreID|32,
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
	
		for (uint i=0;i<NUMFIBRES * (uint)Fs;i++)
		{
			sdramout_buffer[i]  = 0;//((uint)1<<31)-1;//0x5a5a5a5a;
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
		io_printf (IO_BUF, "word size=%d\n",sizeof(REAL));
		
	}
	
	//============MODEL INITIALISATION================//
	
	//calculate startup values
	startupValues=generateStartupVars();
	//initialise cilia
	Cilia.tc=(REAL)(0.00012*fixer);
	Cilia.C=(REAL)(0.004*fixer);
	Cilia.u0=(REAL)(5e-9*fixer);
	Cilia.recips0=(REAL)((1/3e-8)*pow(2,33));
	Cilia.u1=(REAL)(1e-9*fixer);
	Cilia.recips1=(REAL)((1/1e-9)*pow(2,33));
	Cilia.Gmax=(REAL)(6e-9*fixer);
	Cilia.Ga=(REAL)(8e-10*fixer);
	Cilia.dtCap=(REAL)((dt/4e-12)*pow(2,33));
	Cilia.Et=(REAL)(0.1*fixer);
	Cilia.Gk=(REAL)(2e-8*fixer);
	Cilia.Ek=(REAL)(-0.08*fixer);
	Cilia.Rpc=(REAL)(0.04*fixer);
	
	//==========Recurring Values=================//
	IHCVnow=startupValues.IHCVnow0;
	mICaCurr=startupValues.mICaCurr0;;
	for(uint i=0;i<NUMLSR;i++)
	{
		CaCurrLSR[i]=startupValues.CaCurrLSR0;
		ANCleftLSR[i]=startupValues.ANCleftLSR0;
		ANAvailLSR[i]=startupValues.ANAvailLSR0;
		ANReproLSR[i]=startupValues.ANReproLSR0;
	}
	for(uint i=0;i<NUMMSR;i++)
	{
		CaCurrMSR[i]=startupValues.CaCurrMSR0;
		ANCleftMSR[i]=startupValues.ANCleftMSR0;
		ANAvailMSR[i]=startupValues.ANAvailMSR0;
		ANReproMSR[i]=startupValues.ANReproMSR0;
	}
	for(uint i=0;i<NUMHSR;i++) 
	{
		CaCurrHSR[i]=startupValues.CaCurrHSR0;
		ANCleftHSR[i]=startupValues.ANCleftHSR0;
		ANAvailHSR[i]=startupValues.ANAvailHSR0;
		ANReproHSR[i]=startupValues.ANReproHSR0;
	}
	
	//=========initialise the pre synapse params========//
	preSyn.recipBeta=2.5e-3*fixer;
	preSyn.gamma=REAL_CONST(130.*fixer);
	preSyn.dtTauM=(dt/1e-4)*fixer;
	preSyn.dtTauCa=(dt/1e-4)*fixer;
	preSyn.Eca=0.066*fixer;
	preSyn.power=REAL_CONST(3.)*fixer;
	//TODO:develop solution to this in fixed point
	preSyn.z=1e38;//reduced from 1e40 to fit within float range (then an additional *100 is added in code)

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

	preSyn.CaTh[0]=1e-14*fixer;
	preSyn.CaTh[1]=1e-14*fixer;
	preSyn.CaTh[2]=1e-13*fixer;
	preSyn.CaTh[3]=1e-13*fixer;
	preSyn.CaTh[4]=1.43e-13*fixer;
	preSyn.CaTh[5]=1.43e-13*fixer;
	preSyn.CaTh[6]=1.43e-13*fixer;
	preSyn.CaTh[7]=1.43e-13*fixer;
	preSyn.CaTh[8]=1.43e-13*fixer;
	preSyn.CaTh[9]=1.43e-13*fixer;

	//=======initialise the synapse params=======//

	Synapse.l=REAL_CONST(2580.*fixer);
	Synapse.y=REAL_CONST(10.*fixer);
	Synapse.x=REAL_CONST(30.*fixer);
	Synapse.r=REAL_CONST(6580.*fixer);
	Synapse.refrac_period=7.5e-4*fixer;

	Synapse.M[0]=REAL_CONST(13.*fixer);
	Synapse.M[1]=REAL_CONST(13.*fixer);
	Synapse.M[2]=REAL_CONST(14.*fixer);
	Synapse.M[3]=REAL_CONST(14.*fixer);
	Synapse.M[4]=REAL_CONST(8.*fixer);
	Synapse.M[5]=REAL_CONST(8.*fixer);
	Synapse.M[6]=REAL_CONST(8.*fixer);
	Synapse.M[7]=REAL_CONST(8.*fixer);
	Synapse.M[8]=REAL_CONST(8.*fixer);
	Synapse.M[9]=REAL_CONST(8.*fixer); 

	  
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
	REAL *dtcm_buffer_in;
	//read from DMA and copy into DTCM
	if(test_DMA == TRUE && ticks<TOTAL_TICKS)
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
	if (ticks >= TOTAL_TICKS) 
	{
		io_printf (IO_BUF, "spinn_exit\n");
		spin1_exit (0); 
	}
}

REAL longMultiply2(REAL val1, REAL val2)
{
	//this function performs rounding on the multiplication
	//of two longs and returns a long output of the product
	REAL product;
	unsigned REAL a,x;//may need to be unsigned
	REAL b,y;
	unsigned REAL LSWMask=(unsigned REAL)OVERFLOWMASK;
	unsigned REAL MSWMask= (unsigned REAL)OVERFLOWMASK<<32;
	//unsigned REAL check= (pow(2,63)-1) - pow(2,63-(63-FRACBITS));
	unsigned REAL check= ((1<<63)-1) - (1<<(63-(63-FRACBITS)));
	//unsigned REAL negcheck=pow(2,63);
	unsigned REAL negcheck=1<<63;
			
	a=val1 & LSWMask;
	b=(REAL)(val1 & MSWMask)>>32;
	x=val2 & LSWMask;
	y=(REAL)(val2 & MSWMask)>>32;

	REAL add=(a*y)+(b*x);
	unsigned REAL add2=(unsigned REAL)a*(unsigned REAL)x;
	REAL top=(b*y);
	//rounding
	add+=(REAL)1<<18;
	add2+=(REAL)1<<50;

	//check top bits of top word for overflow

	//positive number that will overflow
	if((top & check) && !(top & negcheck))
		{
		io_printf (IO_BUF,"Positive Overflow detected!!\n");
			product= MAXSIGNEDL;
		}

	//negative number that will overflow
	else if((top & negcheck) && !(top & check))
		{
		io_printf (IO_BUF,"Negative Overflow detected!!\n");
			product= MINSIGNEDL;
		}

	else
		{
			int topShift=(64-FRACBITS);
			add=add>>(32-topShift);
			add2=add2>>(64-topShift);
			top=top<<topShift;
			product=(top+add+add2);
		}

	//float testconv=(float)product/pow(2,FRACBITS);

return product;
}

REAL lExp(REAL val)
{
	REAL exponential;
	float valconv= (float)val/fixer;
	float fexponent=exp(valconv);

	//check if value is going to overflow when converted to fixed point
	if(fexponent>=MAXEXP)//fexponent==INFINITY)
	{
		exponential=MAXSIGNEDL;//max signed value
	}

	else
	{
		//convert back to fixed point
		exponential=(REAL)(fexponent*fixer);
	}

	return exponential;
}

REAL lApicalExp(REAL ut,REAL u,REAL recips)
{
	REAL exponential;
	REAL temp;
	REAL subt=((REAL)-1*(ut-u));
	//float subtconv=(float)subt/pow(2,FRACBITS);
	
	temp=longMultiply2(subt,recips);
	temp=temp<<FRACBITS-33;//TODO: add overflow check?
	
	exponential=lExp(temp);
	
	return exponential;
}

uint process_chan(REAL *out_buffer,REAL *in_buffer) 
{
	uint segment_offset=NUMFIBRES*SEGSIZE*(seg_index-1);
	REAL ex0,ex1,lone,div,denom,recip,ut,Gu;
	
	for(uint i=0;i<SEGSIZE;i++)
	{
		//=========Apply Scaler==============//
		ut=longMultiply2(in_buffer[i],Cilia.C);
		//utconv[i]=(float)ut[i]/pow(2,FRACBITS);
		

		//========Apical Conductance=========//
		//exponentials
	 	ex0=lApicalExp(ut,Cilia.u0,Cilia.recips0);
		//float exoconv=(float)ex0/pow(2,FRACBITS);

		ex1=lApicalExp(ut,Cilia.u1,Cilia.recips1);
		//float ex1conv=(float)ex1/pow(2,FRACBITS);
		lone=(REAL)1<<FRACBITS;
		denom=longMultiply2(ex0,(ex1+lone));
		//check for maximum overflow
		if (denom==MAXSIGNEDL || denom==MINSIGNEDL)
		{
			div=0;
		}

		else
		{
			denom=denom+lone;
			//float denomconv=(float)denom/pow(2,FRACBITS);
			//TODO:precision is lost here, is there another implementation? Reciprocal approximations
			recip= ((REAL)1<<62)/denom;
			//float recipconv=(float)recip/pow(2,(62-FRACBITS));
			div=longMultiply2(recip<<(FRACBITS-(62-FRACBITS)),Cilia.Gmax);
		}

		//float divconv=(float)div/fixer;
		Gu=Cilia.Ga+div;
		out_buffer[i] = Gu;
		
//		Guconv[i]=(float)Gu[i]/pow(2,FRACBITS);

/*		//========Receptor Potential=========//
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

		//===========mICa=============//
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

		//==============ICa=================//

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

		//==========LSR Fibres==========//
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
		//=========MSR Fibres=========//
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
		//========HSR Fibres=========//
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
			//HSRfibres[j][i]=IHC_RP[i];*/
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
		//increment segment index
		seg_index++;
		
		#ifdef PROFILE
		  int start_count = tc[T2_COUNT];
		#endif
		
		//choose read buffer
		if(!read_switch && !write_switch)
		{
			index_x=process_chan(dtcm_buffer_x,dtcm_buffer_b);
			#ifdef PRINT 
			io_printf (IO_BUF, "buff_b-->buff_x\n");
			#endif
		}
		else if(!read_switch && write_switch)
		{
			index_y=process_chan(dtcm_buffer_y,dtcm_buffer_b);
			#ifdef PRINT
			io_printf (IO_BUF, "buff_b-->buff_y\n");
			#endif
		}
		else if(read_switch && !write_switch)
		{
			index_x=process_chan(dtcm_buffer_x,dtcm_buffer_a);
			#ifdef PRINT
			io_printf (IO_BUF, "buff_a-->buff_x\n");
			#endif	
		}
		else
		{
			index_y=process_chan(dtcm_buffer_y,dtcm_buffer_a);
			#ifdef PRINT
			io_printf (IO_BUF, "buff_a-->buff_y\n");
			#endif
		}		  
			
		#ifdef PROFILE
			int final_count = tc[T2_COUNT];
		#ifdef PRINT 
			io_printf (IO_BUF, "process complete in %d ticks\n",start_count-final_count);
		#endif	
		#endif			
		
		spin1_trigger_user_event(NULL,NULL);
	}
	else if (ttag==DMA_WRITE)
	{
		//flip write buffers
		write_switch=!write_switch;
	//	io_printf(IO_BUF,"[core %d] Write complete\n",coreID);
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
/*
  // report number of DMA transfers
  io_printf (IO_BUF, "[core %d] completed %d DMA transfers\n", coreID, transfers);

  // report number of failed transfers
  io_printf (IO_BUF, "[core %d] failed %d DMA transfers\n", coreID, tfailed);

  if (tfailed)
  {
    io_printf (IO_BUF, "\t%d : %d @ %d\n", tfvald, tfvals, tfaddr);

    io_printf (IO_BUF, "\t%d : %d @ %d\n", dtcm_buffer[tfaddr],
               sdram_buffer[tfaddr], tfaddr);
  }*/

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
  //spin1_callback_on (DMA_TRANSFER_DONE,process_chan,0);
  spin1_callback_on (DMA_TRANSFER_DONE,transfer_handler,0);
  //reads from DMA to DTCM every tick
  spin1_callback_on (TIMER_TICK,data_read,-1);
  spin1_callback_on (USER_EVENT,data_write,0);

  app_init();

  spin1_start (SYNC_WAIT);
  
  app_done ();

}