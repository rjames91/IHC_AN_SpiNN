/*
 * ApicalExp.h
 *
 *  Created on: 17 Nov 2016
 *      Author: rjames
 */

#ifndef APICALEXP_H_
#define APICALEXP_H_

#define MAXEXP (long)(pow(2,63)-1)/pow(2,FRACBITS)

typedef struct {
	int result, Qi, Qf;
}fixedPointVars;

int exponentials[100];
fixedPointVars ApicalExpLookup(fixedPointVars ut, int u, int recips);
fixedPointVars ApicalExp(fixedPointVars ut, int u, int recips);
fixedPointVars long2int(long input,int fracbits);
long lpow(long val,long power);
long lExp(long val);
long lApicalExp(long ut,long u,long recips);
long longMultiply(long val1, long val2);
long longMultiply2(long val1, long val2);
long overflowCheck(long temp,int fracbits1, int fracbits2);

#endif /* APICALEXP_H_ */
