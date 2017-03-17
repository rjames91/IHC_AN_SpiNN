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
#define random_gen() mars_kiss32()//spin1_rand()//WELL1024a_simp()//rand()//
#define RDM_MAX (UINT32_MAX+REAL_CONST(1.0))//(RAND_MAX+REAL_CONST(1.0))//

#define SEGSIZE (uint)100//200
#define NUMLSR 1//2
#define NUMMSR 2
#define NUMHSR 6
#define NUMFIBRES 10

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
REAL Eca,recipBeta,gamma,dtTauM,dtTauCa,power,z;
REAL GmaxCa[10];
REAL CaTh[10];
} IHC_preSynParams;

typedef struct {
REAL numFibres,refrac_period,TW_delay,spike_Fs,y,x,l,r;
REAL M[10];
} IHC_synParams;

#endif /* IHC_AN_H_ */