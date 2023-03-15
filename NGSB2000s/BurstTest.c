/* 
 * File: BurstTest.c     by: Steve Drasco 
 *                                                           
 * STAND ALONE VERSION                                       
 * This program tests the Burst statistic.   
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		LambdaB;
	double		*h;
	int		i, N;


	/* read input */
	N = atof(argv[1]);

	/* allocate memory */
	h = (double *) malloc(N*sizeof(double));
	for(i=0; i < N; i++) h[i] = ((double) i) + 1.0;
	LambdaB = BurstStatistic(N, h);
	printf("\n\t LambdaB = %e\n",LambdaB);
}
