/*
 ============================================================================
 Name        : SpiNNakEar_IHCAN.c
 Author      : Robert James
 Version     : 1.0
 Description : Inner Hair Cell + Auditory Nerve model for use in SpiNNakEar system
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdfix.h>
#include "AN_group_node.h"
#include "spin1_api.h"
#include "math.h"
#include "random.h"
#include "stdfix-exp.h"
#include "log.h"
#include <data_specification.h>
#include <simulation.h>
#include <debug.h>

// The parameters to be read from memory
enum params {
    N_IHCS,
    AN_KEY,
    IS_KEY,
    IHC_KEYS
};

uint32_t an_key;
uint32_t is_key;
uint32_t n_ihcs;
static key_mask_table_entry *key_mask_table;
uint32_t coreID,chipID;

//application initialisation
bool app_init(void)
{
	io_printf (IO_BUF, "[core %d] -----------------------\n", coreID);
	io_printf (IO_BUF, "[core %d] an group starting simulation\n", coreID);
	//obtain data spec
	address_t data_address = data_specification_get_data_address();
    io_printf(IO_BUF,"data_address=%d\n",data_address);

    //get parameters
    address_t params = data_specification_get_region(
                                        0,data_address);

    n_ihcs = params[N_IHCS];
    io_printf(IO_BUF,"n_ihcs=%d\n",n_ihcs);

    uint32_t n_ihc_key_bytes =  n_ihcs * sizeof(key_mask_table_entry);
    io_printf(IO_BUF,"n_ihc_key_bytes=%d\n",n_ihc_key_bytes);
    key_mask_table = (key_mask_table_entry *)spin1_malloc(n_ihc_key_bytes);
    if (key_mask_table == NULL){
        io_printf (IO_BUF, "failed to allocate memory\n");
        app_end();
    }

    spin1_memcpy(key_mask_table, &params[IHC_KEYS],n_ihc_key_bytes);

    for (uint i=0;i<n_ihcs;i++){
        io_printf(IO_BUF,"ihc key:%d mask:0x%x\n",key_mask_table[i].key,key_mask_table[i].mask);
    }
    an_key = params[AN_KEY];
    is_key = params[IS_KEY];
    io_printf(IO_BUF,"an_key=%d\n",an_key);
    io_printf(IO_BUF,"is_key=%d\n",is_key);
}

void key_search_and_send(uint32_t spike,uint null){
//    io_printf(IO_BUF,"r\n");
    //search through ihc_keys for the rx key
    uint32_t imin = 0;
    uint32_t imax = n_ihcs;
    uint32_t imid;
    key_mask_table_entry entry;
    while (imin < imax) {
        imid = (imax + imin) >> 1;
        entry = key_mask_table[imid];
        if ((spike & entry.mask) == entry.key){
            uint32_t neuron_id = (entry.n_atoms*imid)+ (spike & ~entry.mask);
//            io_printf(IO_BUF,"id:%d\n",neuron_id);
            while(!spin1_send_mc_packet(an_key|neuron_id,0,NO_PAYLOAD)){
                spin1_delay_us(1);
            }
            break;
        }
        else if (entry.key < spike) {
            // Entry must be in upper part of the table
            imin = imid + 1;
        } else {
            // Entry must be in lower part of the table
            imax = imid;
        }
    }

}

//void user_event_callback(uint null,uint null){
//
//
//}

void spike_rx(uint32_t mc_key, uint null){
//    io_printf(IO_BUF,"k:%d\n",mc_key);
    spin1_schedule_callback(key_search_and_send,mc_key,NULL,1);
//    if (is_key){
//        key_search_and_send(mc_key,NULL);
//        }
}

void app_end()
{
    spin1_exit(0);
    io_printf (IO_BUF, "spinn_exit\n");
}

void mc_payload_rx_callback(uint null_a, uint null_b){
//    app_end();
    spin1_schedule_callback(app_end,NULL,NULL,2);
}

void c_main()
{
  // Get core and chip IDs
  coreID = spin1_get_core_id ();
  chipID = spin1_get_chip_id ();

  app_init();

  //setup callbacks
  spin1_callback_on (MC_PACKET_RECEIVED,spike_rx,-1);
  spin1_callback_on (MCPL_PACKET_RECEIVED,mc_payload_rx_callback,-1);

  spin1_start (SYNC_WAIT);
}
