/*
 * File: FileMultiDetect.c   by: Steve Drasco 19 Oct 2000    
 *                                                           
 * STAND ALONE                                               
 * Tries to detect an artificial non-Gaussian stochastic background
 * on artificial Gaussian noise data.
 * 
 * User specifies: (input file with) Alpha, Xi pairs, output file, N
 * Mean results of N attempts per point are sent to output.  
 *
 * 19 Oct 2000
 * Now computes first approzimation of likelihood function as well as exact
 * 						             
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	float	*alphaf, *xif;
	double	sigma, snr;
	double	*alpha, *xi, *Sum1, *Sum2, *Sum3, *Sum4;
	double	LogLambda, A1LogLambda, **h;
	int	i, j, k, N, pairs=0,c;
	FILE 	*fpIN, *fpOUT;

	/* check inputs */
	if(argc != 4) {
                printf("I need inputs: infile, outfile, N\n");
                return;
        }

	/* open input file */
        fpIN = fopen(argv[1],"r");
        if(fpIN == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[1]);
                return;
        }

	/* open output file */
        fpOUT = fopen(argv[2],"w");
        if(fpOUT == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[2]);
                return;
        }

	/* print header of output file */
	fprintf(fpOUT,"%% XiTrue      = %f\n",XITRUE);
	fprintf(fpOUT,"%% AlphaTrue   = %f\n",ALPHATRUE);
	fprintf(fpOUT,"%% SigmaTrue   = %f\n",SIGMATRUE);
	fprintf(fpOUT,"%% Data Length = %f\n",LENGTH);
	fprintf(fpOUT,"%% Output format:\n");
	fprintf(fpOUT,"%% XiGuess\tAlphaGuess\tExactLogLambda\tApprox1LogLambda\tSNR\n\n");

	/* read other inputs */
	N = atol(argv[3]);

	/* count alpha,xi pairs */ 
	while( (c = getc(fpIN)) != EOF) if( c == '\n') pairs++;

	/* close and reopen input file */
	fclose(fpIN);
	fpIN = fopen(argv[1],"r");

	/* allocate memory for alpha xi and sums */
	alphaf = (float *)malloc(pairs*sizeof(float));
	xif = (float *)malloc(pairs*sizeof(float));
	alpha = (double *)malloc(pairs*sizeof(double));
	xi = (double *)malloc(pairs*sizeof(double));
	Sum1 = (double *)malloc(pairs*sizeof(double));
	Sum2 = (double *)malloc(pairs*sizeof(double));
	Sum3 = (double *)malloc(pairs*sizeof(double));
	Sum4 = (double *)malloc(pairs*sizeof(double));


	/* fill xi and alpha LEARN HOW TO SCAN IN DOUBLES */
	for(i=0; i < pairs; i++) {
		fscanf(fpIN,"%f\t%f", xif+i,alphaf+i);
		xi[i] = (float) xif[i];
		alpha[i] = (float) alphaf[i];
	}

	/* allocate memory for artificial data */
	h = (double **) malloc(N*sizeof(double *));
	for(i=0;i<N;i++) h[i]=(double *)malloc(LENGTH*sizeof(double));

	/* make artificial data */
	for(i=0;i<N;i++) NGSBData(LENGTH, SIGMATRUE, ALPHATRUE, XITRUE, h[i]);

	/* main loop for detections */
	for(i = 0; i < pairs; i++) {

		/* reset average */
		LogLambda   = 0.0;
		A1LogLambda = 0.0;

		/* compute exact and approximate likelihood ratio */
		for(k = 0; k < N; k++) {
			Approx1(LENGTH, SIGMATRUE, alpha[i], h[k], Sum1+k, Sum2+k, Sum3+k, Sum4+k);
			LogLambda   += ExactLogLikelihood(LENGTH, SIGMATRUE, alpha[i], xi[i], h[k]) / ((double) N);
			A1LogLambda +=  (  xi[i] * Sum1[k]
				     	 - xi[i] * xi[i] * Sum2[k] / 2.0
					 + xi[i] * xi[i] * xi[i] * Sum3[k] / 3.0
					 - xi[i] * xi[i] * xi[i] * xi[i] * Sum4[k] /  4.0
					) / ((double) N);
		}
		
		/* write results to file */
		snr = xi[i] * sqrt( ((double) LENGTH) /2.0 ) * alpha[i] / SIGMATRUE;
		fprintf(fpOUT,"%f\t%f\t%f\t%f\t%f\n", xi[i], alpha[i], LogLambda, A1LogLambda, snr);
		printf("%f\t%f\t%f\t%f\t%f\n", xi[i], alpha[i], LogLambda, A1LogLambda, snr);
	}
	
	/* close the files */
	fclose(fpIN);
	fclose(fpOUT);
}
