#include <stdint.h>
#include <stdio.h>
#include <arm_neon.h>

/* Name: void Newton_Raph_iter_f(float32_t *d, float32_t *x, int times, int nOP) 
 * Function for testing the Newton-Raphson iteration generating higher accurate 
 * reciprocal result.
*/

void Newton_Raph_iter_f(float32_t *d, float32_t *x, int times, int nOP)
{
	int i,j;

	if(times<=0 || times>100 || nOP!=4) // limit the iteration times and the number of input data
	{
		printf("times = %d , nOP = %d",times, nOP);
		return;
	}
	
	//printf("\nReci_esti: %f  %f  %f %f \n", x[0],x[1],x[2],x[3]);

	for(i=times; i>0; i--)
	{
		
		//printf("round %d  :",times-i+1);
		for(j=0; j<nOP; j++)
		{
			x[j] = x[j]*(2 - d[j]*x[j]);
			//printf(" %f ", x[j]);		
		}
		//printf("\n");
	}	
}


/* Name: void Test_Estimate()
 * A test funtion for calculating reciprocal estimate in 0.5~1.0, and comparing its error 
 * to a normal division result. 
 * Also uses the Newton-Raphson iteration for higher accurate reciprocal estimate result.
*/

void Test_Estimate()
{
	int i = 0;
	float32_t val = 0.5;
	float32x4_t f4_op;
	float32x4_t f4_re;
	float32_t af4_op[4] = {0.5, 0.52, 0.54, 0.56};
	float32_t af4_re[4] = {0,0,0,0};
	float32_t af4_re_it[4] = {0,0,0,0};
	float32_t recip[4] = {0,0,0,0};
	float32_t diff[4] = {0,0,0,0};
	
	printf(" op_val,pc_reciprocal,vrecpeq_estimate,difference, esti_it , diff_it\n");

	do
	{
		for(i=0; i<4; i++)
		{
			af4_op[i] = val; // init val
			recip[i] = 1/af4_op[i]; // calculating the reciprocal with normal division
			
			val += 0.01;//step
		}

		f4_op = vld1q_f32(af4_op);
		f4_re = vrecpeq_f32(f4_op); // floating point vrecpe reciprocal estimate
		vst1q_f32(af4_re, f4_re);

		diff[0] = recip[0]-af4_re[0];
		diff[1] = recip[1]-af4_re[1];
		diff[2] = recip[2]-af4_re[2];
		diff[3] = recip[3]-af4_re[3];
		
		af4_re_it[0] = af4_re[0];
		af4_re_it[1] = af4_re[1];
		af4_re_it[2] = af4_re[2];
		af4_re_it[3] = af4_re[3];

		Newton_Raph_iter_f(af4_op,af4_re_it, 4,4); // using the Newton-Raphson iteration for higher accurate results

		for(i=0; i<4; i++)
		{
			printf("  %f  ,  %f ,  %f  ,  %f ,  %f  ,  %f \n", 
				af4_op[i],recip[i],af4_re[i],diff[i],af4_re_it[i],recip[i]-af4_re_it[i]);	
		}
		
	}while(val < 1.0);
}


/* Name: void UINT_VRECPE_test()
 * Funcion for mixed test of fixed-point unsigned integer type VRECPE and floating point type VRECPE.
*/

void UINT_VRECPE_test()
{
	int i = 0;
	uint32x4_t v4_op;	
	uint32x4_t v4_re;

	uint32_t a4_op[4] = {0x80000000,0xc0000000,0xe0000000,0xf0000000};
	uint32_t a4_re[4] = {0,0,0,0};

	float32x4_t f4_op;
	float32x4_t f4_re;
	float32_t af4_op[4] = {0.5, 0.75, 0.875, 0.9375};
	float32_t af4_re[4] = {0,0,0,0};

	float32_t af4_r[4] = {0,0,0,0};
	float32_t diff[4] = {0,0,0,0};

	/**** uint reciprocal estimate **********/
	v4_op = vld1q_u32(a4_op);
	v4_re = vrecpeq_u32(v4_op);
	vst1q_u32(a4_re, v4_re);

	printf("a4_op(uint32_t)  = 0x%x, 0x%x, 0x%x, 0x%x \n", a4_op[0], a4_op[1],a4_op[2],a4_op[3]);
	printf("Reci_Esti: a4_re = 0x%x, 0x%x, 0x%x, 0x%x \n", a4_re[0], a4_re[1],a4_re[2],a4_re[3]);

	/**** floating point reciprocal estimate **********/
	f4_op = vld1q_f32(af4_op);
	f4_re = vrecpeq_f32(f4_op);
	vst1q_f32(af4_re, f4_re);

	printf("af4_op(float32_t) = %f, %f, %f, %f \n", af4_op[0], af4_op[1],af4_op[2],af4_op[3]);
	printf("Reci_Esti: af4_re = %f, %f, %f, %f \n", af4_re[0], af4_re[1],af4_re[2],af4_re[3]);


	/**reciprocal by normal division**/
	for(i=0; i<4; i++)
	{
		af4_r[i] = 1/af4_op[i];
		diff[i] = af4_r[i] - af4_re[i];
	}
	printf("Reciprocal: af4_r = %f, %f, %f, %f \n", af4_r[0], af4_r[1],af4_r[2],af4_r[3]);
	printf("Difference: diff  = %f, %f, %f, %f \n", diff[0], diff[1],diff[2],diff[3]);


	printf("Test_Estimate() from 0.5 to 1.0: \n");
}

/* Name: int main()
 * calling two test sub-programs.
*/

int main()
{
	UINT_VRECPE_test();
	Test_Estimate();

	return 0;
}

