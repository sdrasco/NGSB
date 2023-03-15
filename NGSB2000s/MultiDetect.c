/*-----------------------------------------------------------*/
/* File: MultiDetect.c     by: Steve Drasco 27 July 2000     */
/*                                                           */
/* STAND ALONE                                               */
/* This program tries to detect an artificial non-Gaussian   */
/* stochastic background superimposed on artificial data.    */
/* User specifies: AlphaMin, DeltaAlpha, AlphaMax, XiMin,    */
/* DeltaXi, XiMax, and the name of an output file.	     */
/* Mean results of N attempts per point are sent to output.  */
/* 						             */
/*-----------------------------------------------------------*/

#include "NGSB.h"

main(int argc, char *argv[])
{
	double	sigma, *h, LogLambda;
	double	AlphaMin, AlphaMax, DeltaAlpha, XiMin, XiMax, DeltaXi;
	int	i, j, k, N, AlphaSteps, XiSteps;
	FILE 	*fp;

	/* check inputs */
	if(argc != 9) {
                printf("I need inputs: outfile, N, AlphaMin, DeltaAlpha, AlphaMax, XiMin, DeltaXi, and XiMax\n");
                return;
        }

	/* open output file */
        fp = fopen(argv[1],"w");
        if(fp == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[1]);
                return;
        }

	/* read other inputs */
	N = atol(argv[2]);
	AlphaMin = atof(argv[3]);
	DeltaAlpha = atof(argv[4]);
	AlphaMax = atof(argv[5]);
	XiMin = atof(argv[6]);
	DeltaXi = atof(argv[7]);
	XiMax = atof(argv[8]);

	/* allocate memory for artificial data */
	h = (double *) malloc(LENGTH*sizeof(double));

	/* compute steps */
	XiSteps = (XiMax - XiMin) / DeltaXi;
	AlphaSteps = (AlphaMax - AlphaMin) / DeltaAlpha; 

	/* print xi and alpha grids to output file */
	for(i = 0; i < AlphaSteps; i++) {
		for(j = 0; j < XiSteps; j++) {
			fprintf(fp,"%f\t", XiMin + DeltaXi * ((double) j));
			printf("%f\t", XiMin + DeltaXi * ((double) j));
		}
		fprintf(fp,"\n");
		printf("\n");
	}
	fprintf(fp,"\n\n");
	printf("\n\n");
        for(i = 0; i < AlphaSteps; i++){
	                for(j = 0; j < XiSteps; j++){
				fprintf(fp,"%f\t", AlphaMin + DeltaAlpha * ((double) i));
				printf("%f\t", AlphaMin + DeltaAlpha * ((double) i));
	                }
	                fprintf(fp,"\n");
			printf("\n");
        }
        fprintf(fp,"\n\n");
	printf("\n\n");

	/* tell user how many runs */
	printf("\t MultiDetect:  Starting first of %d runs.\n",AlphaSteps*XiSteps*N);
	printf("\t               Runtime ~ %f hours.\n",((double) AlphaSteps*XiSteps*N) * 0.00033333);

	/* main loop for detections */
	for(i = 0; i < AlphaSteps; i++) {

		for(j = 0; j < XiSteps; j++) {

			/* reset average */
			LogLambda = 0.0;

			for(k = 0; k < N; k++) {

           			/* make artificial data */
			        NGSBData(LENGTH, SIGMATRUE, ALPHATRUE, XITRUE, h);

			        /* compute exact likelihood ratio once just to test */
			        LogLambda += ExactLogLikelihood(LENGTH, SIGMATRUE, AlphaMin + DeltaAlpha * ((double) i), XiMin + DeltaXi * ((double) j), h) / ((double) N); 
			}
		
			/* write results to file */
			fprintf(fp,"%f\t",LogLambda);
			printf("%f\t",LogLambda);
		}
		
		/* go to next line in matrix */
		fprintf(fp,"\n");
		printf("\n");
	}
	
	/* close the output file */
	fclose(fp);
}
