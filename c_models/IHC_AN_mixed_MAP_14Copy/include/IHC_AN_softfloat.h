/*
 * IHC_AN.h
 *
 *  Created on: 1 Nov 2016
 *      Author: rjames
 */

#ifndef IHC_AN_softfloat_H_
#define IHC_AN_softfloat_H_

#define REAL float
#define REAL_CONST(x) x##f
#define REAL_MAX 1e38
#define random_gen() mars_kiss64_seed(local_seed)//seed_scrambled_rgen()//WELL1024a_seed(seed)//mars_kiss32()//spin1_rand()//WELL1024a_simp()//rand()//
#define RDM_MAX (UINT32_MAX+REAL_CONST(1.0))//(RAND_MAX+REAL_CONST(1.0))//
#define SEED_TYPE mars_kiss64_seed_t//WELL1024a_seed_t
#define MAX_CHIPX 1//255
#define MAX_CHIPY 1//255
#define MAX_COREID 16
//#define SEED_SEL_SIZE (((MAX_CHIPX << 8) | MAX_CHIPY) << 8) | (MAX_COREID << 3)   
#define SEED_SEL_SIZE 1024

#define SEGSIZE 100//200
#define NUMLSR 2//10
#define NUMMSR 0
#define NUMHSR 2//10
#define NUMFIBRES 4//20

#define SAMPLING_FREQUENCY 44100.f
#define RESAMPLE_FACTOR 20//5


#define MAX_SIGNAL_S 1

#define TIMER2_CONF        0x82
#define TIMER2_LOAD        0


#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

typedef struct {
REAL Ga,tc,C,u0,recips0,u1,recips1,Gmax,dtCap,Et,Gk,Ek,Rpc;
} IHC_ciliaParams;

typedef struct {
REAL Eca,recipBeta,gamma,dtTauM,power,z;
REAL TauCa[NUMFIBRES];
REAL GmaxCa[NUMFIBRES];
REAL CaTh[NUMFIBRES];
} IHC_preSynParams;

typedef struct {
//REAL numFibres,refrac_period,TW_delay,spike_Fs,ydt,xdt,ldt,rdt;
REAL numFibres,TW_delay,spike_Fs,ydt,xdt,ldt,rdt;
unsigned accum refrac_period;
REAL M[NUMFIBRES];
} IHC_synParams;

typedef union
{
	unsigned long fract ulf;
	uint32_t ui;
	float sf;
} ui_ulf;


#endif /* IHC_AN_H_ */
