/*
 ============================================================================
 Name        : data_move_test.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Test application to read from SDRAM using DMA into a DTCM buff 
	       then save this buffer into a different SDRAM location. 
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "spin1_api.h"

#define TIMER_TICK_PERIOD  1000000
#define BUFFER_SIZE        50

uint coreID;
uint chipID;
uint test_DMA;

uint *dtcm_buffer;
uint *sdramin_buffer;
uint *sdramout_buffer;

//application initialisation
void app_init(void)
{
  /* say hello */

  io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
  io_printf (IO_BUF, "[core %d] starting simulation\n", coreID);

 // Allocate buffers somewhere in SDRAM

  sdramout_buffer = (uint *) sark_xalloc (sv->sdram_heap,
					BUFFER_SIZE * sizeof(uint),
					0,
					ALLOC_LOCK);
  
  sdramin_buffer = (uint *) sark_xalloc (sv->sdram_heap,
					BUFFER_SIZE * sizeof(uint),
					1,
					ALLOC_LOCK);

  // and a buffer in DTCM

  dtcm_buffer = (uint *) sark_alloc (BUFFER_SIZE, sizeof(uint));

  if (dtcm_buffer == NULL ||  sdramout_buffer == NULL || sdramin_buffer == NULL)
  {
    test_DMA = FALSE;
    io_printf (IO_BUF, "[core %d] error - cannot allocate buffer\n", coreID);
  }
  else
  {
    test_DMA = TRUE;
    // initialize sections of DTCM, system RAM and SDRAM
    for (uint i = 0; i < BUFFER_SIZE; i++)
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
		       BUFFER_SIZE*sizeof(uint));

     io_printf (IO_BUF, "[core %d] sdram DMA read @ 0x%08x\n", coreID,
               sdramin_buffer);
  }
}

//callback on DMA read complete
void data_move(uint null_a, uint null_b)
{
//run a DMA write command to copy contents of the DTCM buffer,
 /*  spin1_dma_transfer(DMA_WRITE, sdram_buffer, dtcm_buffer, DMA_WRITE,
		       BUFFER_SIZE*sizeof(uint));	*/
  for (uint i = 0; i < BUFFER_SIZE; i++)
    {
      sdramout_buffer[i]  = dtcm_buffer[i];
    }

  io_printf (IO_BUF, "[core %d] dtcm to sdram copy from @ 0x%08x to @ 0x%08x\n", coreID,
               (uint) dtcm_buffer,(uint) sdramout_buffer);

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
  spin1_callback_on (DMA_TRANSFER_DONE,data_move,0);
  //reads from DMA every tick
  spin1_callback_on (TIMER_TICK,data_read,1);

  app_init();

  spin1_start (SYNC_WAIT);

}
