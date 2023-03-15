/*-----------------------------------------------------------*/
/* File: NGSBDataTest.c     by: Steve Drasco 27 July 2000    */
/*                                                           */
/* STAND ALONE VERSION                                       */
/* This program generates detector data with a Non-Gaussian  */
/* stochastic background characterized by user specified     */
/* output file, length (N), sigma, alpha, and xi             */
/*-----------------------------------------------------------*/

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		sigma, alpha, xi, *n, *s, *h;
	int		i, N;
	FILE 		*fp;
	extern long     idum;

	/* check inputs */
	if(argc != 6) {
                printf("I need inputs: outfile, N, sigma, alpha, and xi\n");
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
	alpha = atof(argv[4]);
	xi = atof(argv[5]);

	/* allocate memory for noise and signal */
	n = (double *) malloc(N*sizeof(double));
	s = (double *) malloc(N*sizeof(double));

	/* seed the random number generator */
	SeedRand();

	/* make noise */
	Gaussian(sigma,N,PMIN,n);

	/* make Gaussian background signal */
	Gaussian(alpha,N,PMIN,s);

	/* zero out parts of signal as output is filled */
	for (i = 0; i < N; i++) {
		if(((double) ran2(&idum)) < xi ) {
			fprintf(fp,"%f\t%f\t%f\n",n[i],s[i],n[i]+s[i]);
		} else {
			fprintf(fp,"%f\t%f\t%f\n",n[i],0.0,n[i]);
		}
	}

	/* close the output file */
	fclose(fp);

	/* For testing the library version */
	free(h);
	h = (double *) malloc(N*sizeof(double));
	if (1 == NGSBData(N, sigma, alpha, xi, h)) {
		printf("\nLibrary version of NGSBData failed. Giving up...\n");
		return;
	} 
	printf("\nResults from Library version:\n");
	for(i = 0; i < N; i++) printf("%f\n",h[i]);
}
