/* 
 * File: GeneralMultiDetect.c  (stant alone)
 *
 * This program computes the exact and generalized log of the liklihood function FOR TWO DETECTORS
 * over a grid of (alpha,xi).  The input files must contain the xi values and alpha values.
 *
 * Inputs are: XiInFile, AlphaInFile, OutFile, XiTrue, AlphaTrue, SNR, Sigma1, Sigma2
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	float		*XiF, *AlphaF;
	double		LogLambda, XiTrue, AlphaTrue, Rho, Sigma1, Sigma2;
	double		*h1, *h2, *n1, *n2, *s, *Xi, *Alpha, *theta;
	double		kappa = 100.0;
	int		i, j, k, N, c, AlphaN=0, XiN=0;
	FILE 		*XiInFile, *AlphaInFile, *OutFile;

	/* check inputs */
	if(argc != 9) {
                printf("I need inputs: XiInFile, AlphaInFile, OutFile, XiTrue, AlphaTrue, SNR, Sigma1, Sigma2\n");
                return;
        }

        /* read inputs */
        XiTrue = (double) atof(argv[4]);
        AlphaTrue = (double) atof(argv[5]);
        Rho = (double) atof(argv[6]);
	Sigma1 = (double) atof(argv[7]);
	Sigma2 = (double) atof(argv[8]);

	/* ##### NEEDS TO BE MODIFIED TO REFLECT TWO DETECTORS ? ##### */
	/* compute N and abort if over 1e7 */
	N = (int) ( Rho * Rho * Sigma1 * Sigma1 / (XiTrue * XiTrue * AlphaTrue * AlphaTrue) );
	if(N >= (int) 1e7 || N <= 1) {
		printf("N is unreasonable.  Maximum N is 1e7, you tried: %d\nGiving up...\n",N);
		return;
	}
	printf("This run has N = %d\n",N);

	/* open input files and count entries */
        XiInFile = fopen(argv[1],"r");
        if(XiInFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[1]);
                return;
        }
        AlphaInFile = fopen(argv[2],"r");
        if(AlphaInFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[2]);
                return;
        }
	while( (c = getc(XiInFile)) != EOF) if( c == '\n') XiN++;
	while( (c = getc(AlphaInFile)) != EOF) if( c == '\n') AlphaN++;

        /* close and reopen input files */
        fclose(XiInFile);
	fclose(AlphaInFile);
        XiInFile = fopen(argv[1],"r");
	AlphaInFile = fopen(argv[2],"r");

	/* allocate memory */
	h1 = (double *) malloc(N*sizeof(double));
	h2 = (double *) malloc(N*sizeof(double));
	n1 = (double *) malloc(N*sizeof(double));
	n2 = (double *) malloc(N*sizeof(double));
	s = (double *) malloc(N*sizeof(double));
	theta = (double *) malloc(N*sizeof(double));
        XiF = (float *) malloc(XiN*sizeof(float));
        AlphaF = (float *)malloc(AlphaN*sizeof(double));
        Xi = (double *)malloc(XiN*sizeof(double));
	Alpha = (double *)malloc(AlphaN*sizeof(double));
	printf("We look at %d grid points.\n",XiN*AlphaN);

        /* fill xi and alpha LEARN HOW TO SCAN IN DOUBLES */
        for(i=0; i < XiN; i++) {
                fscanf(XiInFile,"%f", XiF+i);
                Xi[i] = (float) XiF[i];
	}
        for(i=0; i < AlphaN; i++) {
                fscanf(AlphaInFile,"%f", AlphaF+i);
                Alpha[i] = (float) AlphaF[i];
        }

	/* close input files */
	fclose(XiInFile);
        fclose(AlphaInFile);

        /* open output files */
        OutFile = fopen(argv[3],"w");
	if( OutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[3]);
                return;
        }

	/* print headers to output files & screen */
	fprintf(OutFile,"%% Two Detector Log Lambda Results\n");
        fprintf(OutFile,"%% XiTrue      = %e\n",XiTrue);
        fprintf(OutFile,"%% AlphaTrue   = %e\n",AlphaTrue);
        fprintf(OutFile,"%% Sigma1      = %e\n",Sigma1);
        fprintf(OutFile,"%% Sigma2      = %e\n",Sigma2);
        fprintf(OutFile,"%% Data Length = %d\n",N);
	fprintf(OutFile,"%% SNR         = %e\n",Rho);
        fprintf(OutFile,"%% Output format is meshgrid style: (xi,alpha,logLambda)\n");
        printf("%% Two Detector Log Lambda Results\n");
        printf("%% XiTrue      = %e\n",XiTrue);
        printf("%% AlphaTrue   = %e\n",AlphaTrue);
        printf("%% Sigma1      = %e\n",Sigma1);
        printf("%% Sigma2      = %e\n",Sigma2);
        printf("%% Data Length = %d\n",N);
        printf("%% SNR         = %e\n",Rho);
        printf("%% Output format is meshgrid style: (xi,alpha,logLambda)\n");

	/* print xi and alpha grids to output files */
	for(i = 0; i < AlphaN; i++) {
		for(j = 0; j < XiN; j++) fprintf(OutFile,"%e\t", Xi[j]);
		fprintf(OutFile,"\n");
	}
	fprintf(OutFile,"\n\n");
        for(i = 0; i < AlphaN; i++){
	                for(j = 0; j < XiN; j++) fprintf(OutFile,"%e\t", Alpha[i]);
	                fprintf(OutFile,"\n");
        }
        fprintf(OutFile,"\n\n");

	/* seed random number generator */
	SeedRand();

	/* make artificial data */
	PairNGSBData(N, Sigma1, Sigma2, AlphaTrue, XiTrue, h1, h2, n1, n2, s);

	/* main loop for calculating log Lambda */
	for(i = 0; i < AlphaN; i++) {

		for(j = 0; j < XiN; j++) {

		        /* compute exact log Lambda */
		        LogLambda = GeneralLogLikelihood(N, Sigma1, Sigma2, Alpha[i], Xi[j], h1, h2, kappa, theta); 

			/* write results to file */
			fprintf(OutFile,"%f\t",LogLambda);
		}
		
		/* go to next line in matrix */
		fprintf(OutFile,"\n");
	}
	
	/* close the output file */
	fclose(OutFile);

}
