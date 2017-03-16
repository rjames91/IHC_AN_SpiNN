#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "IHC_Utils.h"
#include "IHC_AN.h"

fixedPointVars ApicalExp(fixedPointVars ut,int u,int recips)
{
fixedPointVars exponential;
long temp;
float testpre=(float)(ut.result)/pow(2,ut.Qf);

//cast inputs as longs
long lrecips=(long)recips;
long lut=(long)ut.result;
long lu = (long)u;


//TODO:add in constraint to ensure adjustments never lose MSBs
//shift u to match fracbits of ut
int fracbits;
if(ut.Qf>FRACBITS)
{
	//a has larger fractional bits case
	//so remove extra fractional bits
	lu=lu<<(ut.Qf-FRACBITS);
	//ut.result=ut.result>>(ut.Qf-30);
	fracbits=ut.Qf;
}
else if(ut.Qf<FRACBITS)
{
	//b has larger fractional bits case
	//so remove extra fractional bits
	lut=lut<<(FRACBITS-ut.Qf);
	fracbits=FRACBITS;
}
else
{
}
float lutconv=(float)(lut)/pow(2,fracbits);
float luconv=(float)(lu)/pow(2,fracbits);

fixedPointVars subt=long2int(-1*(lut-lu),fracbits);

//temp is in Q30+1.fracbits+1 format

temp = lrecips*(long)(subt.result);
//temp = (long)recips*(long)((-1*ut.result)-u);
//float recipsconv=(float)recips/pow(2,1);
float tempconv= (float)temp/pow(2,subt.Qf+1);
//int input=long2int(temp,1+)
//exp function needs a Q16.15 format
//int input=(int)(temp>>16);
fixedPointVars input=long2int(temp,subt.Qf+1);
float inputconv=(float)input.result/pow(2,input.Qf);
float fexponent=exp(inputconv);

if(fexponent==INFINITY)
{
	exponential.result=pow(2,31)-1;
	exponential.Qf=0;
	exponential.Qi=31;
}

else
{
	//scale exponential based on value of inputconv
	int fracbits=7;
	//exponent is in Q
	long lexponent=fexponent*pow(2,fracbits);
	exponential=long2int(lexponent,fracbits);
}

return exponential;
}

long lpow(long val,long power)
{
	long result;
	float valconv=(float)val/pow(2,FRACBITS);
	float fpow=pow(valconv,power);

	//check if value is going to overflow when converted to fixed point
	if(fpow>=MAXEXP)//fexponent==INFINITY)
	{
		result=MAXSIGNEDL;//max signed value
	}

	else
	{
		//convert back to fixed point
		result=(long)(fpow*pow(2,FRACBITS));
	}
	return result;
}

long lExp(long val)
{
	long exponential;
	float valconv= (float)val/pow(2,FRACBITS);
	float fexponent=exp(valconv);

	//check if value is going to overflow when converted to fixed point
	if(fexponent>=MAXEXP)//fexponent==INFINITY)
	{
		exponential=MAXSIGNEDL;//max signed value
	}

	else
	{
		//convert back to fixed point
		exponential=(long)(fexponent*pow(2,FRACBITS));
	}

	return exponential;
}

long lApicalExp(long ut,long u,long recips)
{
long exponential;
long temp;
long subt=((long)-1*(ut-u));
float subtconv=(float)subt/pow(2,FRACBITS);

//temp is in Q13+30.50+33 format
//temp = longMultiply(subt,recips);
//temp=overflowCheck(temp,FRACBITS,33);
temp=longMultiply2(subt,recips);
temp=temp<<FRACBITS-33;//TODO: add overflow check?

exponential=lExp(temp);
/*
float tempconv= (float)temp/pow(2,FRACBITS);
float fexponent=exp(tempconv);

//check if value is going to overflow when converted to fixed point
if(fexponent>=MAXEXP)//fexponent==INFINITY)
{
	exponential=MAXSIGNEDL;//max signed value
}

else
{
	//convert back to fixed point
	exponential=fexponent*pow(2,FRACBITS);
}
*/
return exponential;
}

long longMultiply(long val1, long val2)
{
	//this function performs rounding on the multiplication
	//of two longs and returns a long output of the product
	long product;
	unsigned long a,x;//may need to be unsigned
	long b,y;
	unsigned long LSWMask=(unsigned long)OVERFLOWMASK;
	unsigned long MSWMask= (unsigned long)OVERFLOWMASK<<32;
	a=val1 & LSWMask;
	b=(long)(val1 & MSWMask)>>32;
	x=val2 & LSWMask;
	y=(long)(val2 & MSWMask)>>32;

	//long test=(unsigned long)a*(long)y;
	unsigned long test=a*x;
	long add=(a*y)+(b*x);
	//rounding
	add+=(long)1<<31;
	//add=add+((add & ((long)1<<(FRACBITS-1)))<<1);
	product=b*y+(add>>32);
	float testconv=(float)product/pow(2,FRACBITS+FRACBITS-64);

return product;
}

long longMultiply2(long val1, long val2)
{
	//this function performs rounding on the multiplication
	//of two longs and returns a long output of the product
	long product;
	unsigned long a,x;//may need to be unsigned
	long b,y;
	unsigned long LSWMask=(unsigned long)OVERFLOWMASK;
	unsigned long MSWMask= (unsigned long)OVERFLOWMASK<<32;
	unsigned long check= (pow(2,63)-1) - pow(2,63-(63-FRACBITS));
	unsigned long negcheck=pow(2,63);

	a=val1 & LSWMask;
	b=(long)(val1 & MSWMask)>>32;
	x=val2 & LSWMask;
	y=(long)(val2 & MSWMask)>>32;

	long add=(a*y)+(b*x);
	unsigned long add2=(unsigned long)a*(unsigned long)x;
	long top=(b*y);
	//rounding
	add+=(long)1<<18;
	add2+=(long)1<<50;

	//check top bits of top word for overflow

	//positive number that will overflow
	if((top & check) && !(top & negcheck))
		{
			printf("Positive Overflow detected!!\n");
			product= MAXSIGNEDL;
		}

	//negative number that will overflow
	else if((top & negcheck) && !(top & check))
		{
			printf("Negative Overflow detected!!\n");
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

	float testconv=(float)product/pow(2,FRACBITS);

return product;
}


fixedPointVars long2int(long input, int fracbits)

{
	fixedPointVars fixedResult;

	unsigned int check=(input>>31) & OVERFLOWMASK;
	int count;

	if(!check || check==pow(2,32)-1)
	{
		//result will be in Q31-fracbits.fracbits format
		fixedResult.result=(int)input;
		fixedResult.Qf=fracbits;
		fixedResult.Qi=31-fixedResult.Qf;
	}
	//otherwise shift right until there is no overflow
	else
	{
		count=0;
		while(check>0 && check<pow(2,32)-1)
		{
			input=input>>1;
			check=(input>>31) & OVERFLOWMASK;
			count++;
		}
		//result has lost count LSR bits
		//this has altered the Q format to Q31-fracbits.fracbits-count-1 format
		//TODO: consider the case when count>fracbits
		fixedResult.result=(int)input;
		fixedResult.Qf=fracbits-count;
		fixedResult.Qi=31-fixedResult.Qf;
	}

	return fixedResult;
}

long overflowCheck(long temp,int fracbits1, int fracbits2 )
{
	//input is assumed to be in Q ?.fracbits1+fracbits2-64 format
	//check that overflow won't occur
	//unsigned long check = pow(2,63)-1-pow(2,64-(fracbits2+fracbits1-64)-1);
	unsigned long check = pow(2,63)-1-pow(2,64-(fracbits1-(fracbits2+fracbits1-64))-1);
	unsigned long negcheck=pow(2,63);
	//positive number that will overflow
	if((temp & check) && !(temp & negcheck))
		{
			printf("Positive Overflow detected!!\n");
			return MAXSIGNEDL;
		}

	//negative number that will overflow
	else if((temp & negcheck) && !(temp & check))
		{
			printf("Negative Overflow detected!!\n");
			return MINSIGNEDL;
		}
	else
		{
			//shift input to give result in Q63-fracbits1.fracbits1 format
			temp=temp<<(fracbits1-(fracbits2+fracbits1-64));
			return temp;
		}
}


fixedPointVars ApicalExpLookup(fixedPointVars ut,int u,int recips)

{
	fixedPointVars output;
/*	int temp;
	//ltemp will be in q25+1.6+30 format
	long ltemp= (long)recips*(long)((-1*ut)-u);
	//convert to a Q16.15 format
	temp=(int)(ltemp>>21);
	//search in lookup table for exponentials
	if(temp<=1146880 && temp>=-1146880)
	{
//		output.result=expTable[temp+35];//index offset
		//Qformat depends on temp
	//	output.Qi=
	//	output.Qf=
	}
	else
	{
		output.result=0xFFFF;
	}
*/
	return output;

}

void generateExpTables(void)

{
	int i;

	for(i=-35;i<=35;i++)
	{
		float exponent=exp((float)i);
		//expTable[i+35]=(int)exp((float)i);
	}
//printf("expTable[0]=%f\n",expTable[0]);
}

