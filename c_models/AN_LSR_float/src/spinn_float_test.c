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
#include "IHC_AN_softfloat.h"
#include "spin1_api.h"
#include "math.h"

REAL Fs,dt;
REAL testInput[SEGSIZE];
uint coreID;
uint chipID;

IHC_ciliaParams Cilia;

#define TIMER_TICK_PERIOD  1000000 //1 second

uint coreID;
uint chipID;
uint test_DMA;

float *dtcm_buffer;
float *sdramin_buffer;
float *sdramout_buffer;

//application initialisation
void app_init(void)
{
	/* SETUP CONSTANTS */
	Fs=REAL_CONST(44100.);
	dt=(1.0/Fs);

	//initialise the cilia params
	Cilia.C=REAL_CONST(0.004);
	
	/* say hello */
	
	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);
	
	// Allocate buffers somewhere in SDRAM
	
	sdramout_buffer = (float *) sark_xalloc (sv->sdram_heap,
					SEGSIZE * sizeof(float),
					0,
					ALLOC_LOCK);
	
	sdramin_buffer = (float *) sark_xalloc (sv->sdram_heap,
					SEGSIZE * sizeof(float),
					1,
					ALLOC_LOCK);
	
	// and a buffer in DTCM
	
	dtcm_buffer = (float *) sark_alloc (SEGSIZE, sizeof(float));

	if (dtcm_buffer == NULL ||  sdramout_buffer == NULL || sdramin_buffer == NULL)
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
	dtcm_buffer[i]   = i;
	sdramout_buffer[i]  = 0x5a5a5a5a;
	//sdramin_buffer[i]  = 0x3a3a3a3a;
	}
	
	io_printf (IO_BUF, "[core %d] dtcm buffer @ 0x%08x\n", coreID,
			   (uint) dtcm_buffer);
	io_printf (IO_BUF, "[core %d] sdram out buffer @ 0x%08x\n", coreID,
			   (uint) sdramout_buffer);
	io_printf (IO_BUF, "[core %d] sdram in buffer @ 0x%08x\n", coreID,
			   (uint) sdramin_buffer);
	
	}

}

//DMA read
void data_read(uint null_a, uint null_b)
{
	//read from DMA and copy into DTCM
	if(test_DMA == TRUE)
	{
	 spin1_dma_transfer(DMA_READ,sdramin_buffer, dtcm_buffer, DMA_READ,
			   SEGSIZE*sizeof(uint));
	
	 io_printf (IO_BUF, "[core %d] sdram DMA read @ 0x%08x\n", coreID,
			   sdramin_buffer);
	}
}

void process_chan(uint null_a, uint null_b) 
{
	/***********OUTPUTS******************************/
	REAL utconv[SEGSIZE];

	int i;//i corresponds to each sample in the input
	for(i=0;i<SEGSIZE;i++)
	{
		/**********Apply Scaler****************/
		//no casting here as the actual 32bit word saved to 
		//SDRAM is already in single float format
		testInput[i]=dtcm_buffer[i];
		utconv[i]=testInput[i] * Cilia.C;
		/*****write value to SDRAM****/
		//sdramout_buffer[i]  = testInput[i];
		sdramout_buffer[i]  = utconv[i];	

	}

	io_printf (IO_BUF, "[core %d] all samples processed and written to @ 0x%08x\n", coreID,
               sdramout_buffer);
	spin1_exit (0);	

}

void c_main()
{
	// Get core and chip IDs
	
	coreID = spin1_get_core_id ();
	chipID = spin1_get_chip_id ();
	
	//set timer tick
	spin1_set_timer_tick (TIMER_TICK_PERIOD);
	//setup callback
	spin1_callback_on (DMA_TRANSFER_DONE,process_chan,0);
	//reads from DMA every tick
	spin1_callback_on (TIMER_TICK,data_read,1);
	
	app_init();
	
	spin1_start (SYNC_WAIT);

}

