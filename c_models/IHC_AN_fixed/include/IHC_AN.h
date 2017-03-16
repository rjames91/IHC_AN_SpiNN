/*
 * IHC_AN.h
 *
 *  Created on: 1 Nov 2016
 *      Author: rjames
 */

#ifndef IHC_AN_H_
#define IHC_AN_H_

#define FRACBITS 50//47//50//57
#define MAXSIGNEDL (long)(pow(2,63)-1)
#define MINSIGNEDL (long)-1*(pow(2,63))

#define MASK (short int)pow(2,FRACBITS)

#define OVERFLOWMASK (unsigned int)(pow(2,32)-1)

#define segSize 441//50
#define NUMLSR 2
#define NUMMSR 2
#define NUMHSR 6

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

typedef struct {
long Ga,tc,C,u0,recips0,u1,recips1,Gmax,dtCap,Et,Gk,Ek,Rpc;
//long Ga, ;
} IHC_ciliaParams;

typedef struct {
long Eca,recipBeta,gamma,dtTauM,dtTauCa,power,z;
long GmaxCa[2];
long CaTh[10];
} IHC_preSynParams;

typedef struct {
long numFibres,refrac_period,TW_delay,spike_Fs,y,x,l,r;
long M[10];

} IHC_synParams;





#endif /* IHC_AN_H_ */
