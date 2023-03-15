/* 
 * File: GaussianTest.c     by: Steve Drasco 
 *                                                           
 * STAND ALONE VERSION                                       
 * This program generates data with a Gaussian distribution  
 * with zero mean and users standard deviation and length.   
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		x, y, sigma, Pmin, xmax, ymax, *h;
	double   	pi = 3.14159265358979323846;
	int		i = 0, N;
	FILE 		*fp;
	extern long 	idum;

	/* check inputs */
	if(argc != 5) {
                printf("I need inputs: outfile, length, sigma, Pmin\n");
                return;
        }

	/* open output file */
        fp = fopen(argv[1],"w");
        if(fp == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[1]);
                return;
        }

	/* read other inputs */
	N = atof(argv[2]);
	sigma = atof(argv[3]);
	Pmin = atof(argv[4]);
	printf("\n\tN = %d\n\tsigma = %f\n\tProbability minumum = %f\n\toutput file: %s\n\n",N,sigma,Pmin,argv[1]);

	/* initialize random number generator */
	SeedRand();

	/* compute boundries */
	xmax = sigma * sqrt(  abs(  2*log( Pmin*sigma*sqrt(2*pi) )  )   );
	ymax = 1.0 / ( sigma * sqrt(2*pi) );

	/* main loop */
	while(i < N){

		/* make a new point */
		x = ( 2.0 * xmax * ((double) ran2(&idum)) ) - xmax;
		y = ymax * ((double) ran2(&idum));
	
		/* if point is above curve - keep x */
		if (   y <= exp( - x * x / (2.0 * sigma * sigma) )   /  ( sigma*sqrt(2.0 * pi) )    ) {
			fprintf(fp,"%f\n",x);
			i++;
		}
	}

	/* close the output file */
	fclose(fp);

	/* For testing the library version */
	h = (double *) malloc(N*sizeof(double));
	if (1 == Gaussian(sigma, N, Pmin, h)) {
		printf("\nLibrary version of Gaussian failed. Giving up...\n");
		return;
	} 
	printf("\nResults from Library version:\n");
	for(i = 0; i < N; i++) printf("%f\n",h[i]);
}
