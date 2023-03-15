/* 
 * File: PairLogLogMultiDetect.c  (stant alone)
 *
 * This program computes the exact log of the liklihood function FOR TWO DETECTORS
 * over a grid of (alpha,xi).  The input files must contain the xi values and alpha values.
 *
 * Inputs are: XiInFile, AlphaInFile, OutFile, XiTrue, AlphaTrue, SNR, Sigma1, Sigma2
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	float		*XiF, *AlphaF;
	double		LogLambda, XiTrue, AlphaTrue, Rho, SigmaHat1, SigmaHat2;
	double		SigmaBar1, SigmaBar2, LambdaMax, AlphaBar, XiBar, C;
	double		*h1, *h2, *n1, *n2, *s, *Xi, *Alpha, f=0.5, *XiA, *AlphaA;
	double		**simplex;
	int		i, j, k, N, c, AlphaN=0, XiN=0;
	int		GN = 50;
	FILE 		*XiInFile, *AlphaInFile, *OutFile;


	/* check inputs */
	if(argc != 7) {
                printf("I need inputs: XiInFile, AlphaInFile, OutFile, XiTrue, AlphaTrue, SNR\n");
                return;
        }

        /* read inputs */
        XiTrue = (double) atof(argv[4]); 
        AlphaTrue = (double) atof(argv[5]); 
        Rho = (double) atof(argv[6]);

	/* compute N and abort if over 1e7 */
	N = (int) ( Rho * Rho / (XiTrue * XiTrue * AlphaTrue * AlphaTrue) );
	if(N >= (int) 5e6) {
		printf("Maximum N is 5e6, you tried: %d\nGiving up...\n",N);
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
        XiF = (float *) malloc(XiN*sizeof(float));
        AlphaF = (float *)malloc(AlphaN*sizeof(double));
        Xi = (double *)malloc(XiN*sizeof(double));
	Alpha = (double *)malloc(AlphaN*sizeof(double));
        XiA = (double *) malloc(GN*sizeof(double));
        AlphaA = (double *) malloc(GN*sizeof(double));
	simplex = (double **) malloc(5*sizeof(double *));
	for(i=0; i<5; i++) simplex[i] = (double *)malloc(4*sizeof(double));
	printf("We look at %d grid points.\n",XiN*AlphaN);

        /* fill xi and alpha LEARN HOW TO SCAN IN DOUBLES */
        for(i=0; i < XiN; i++) {
                fscanf(XiInFile,"%e", XiF+i);
                Xi[i] = (float) XiF[i];
	}
        for(i=0; i < AlphaN; i++) {
                fscanf(AlphaInFile,"%e", AlphaF+i);
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
        fprintf(OutFile,"%% XiTrue      = %f\n",XiTrue);
        fprintf(OutFile,"%% AlphaTrue   = %f\n",AlphaTrue);
        fprintf(OutFile,"%% Sigma1      = %f\n",1.0);
        fprintf(OutFile,"%% Sigma2      = %f\n",1.0);
        fprintf(OutFile,"%% Data Length = %d\n",N);
	fprintf(OutFile,"%% SNR         = %f\n",Rho);
        fprintf(OutFile,"%% Output format is meshgrid style: (xi,alpha,logLambda)\n");

        /* seed the random number generator */
	SeedRand();

        /* make artificial signal-FULL data */
	PairNGSBData(N, 1.0, 1.0, AlphaTrue, XiTrue, h1, h2, n1, n2, s); 
	
	/* make artificial signal-less data */
	/*
	Gaussian(1.0, N, PMIN, h1);
	Gaussian(1.0, N, PMIN, h2);
	*/

	/* Compute cross correlation statistic */
	CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C);

	/* make a 'good' guess at parameters */
	simplex[0][0] = 1.0;
        simplex[0][1] = 1.0;
        simplex[0][2] = AlphaTrue; 
        simplex[0][3] = XiTrue; 


	/* Measure ALL PARAMETERS!! */
	printf("ExitStatus: %d\n",
	SimplexMaxLikelihood(&LambdaMax, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, h1, h2, simplex, SigmaHat1, SigmaHat2, C, 7500, 7500, 1e-4, 1e-4));
	
	/* print maximization results to header */
	fprintf(OutFile,"%% MaxLikelihood = %e\n",  LambdaMax);
	fprintf(OutFile,"%% SigmaBar1	  = %e\n",  SigmaBar1);
	fprintf(OutFile,"%% SigmaBar2     = %e\n",  SigmaBar2);
	fprintf(OutFile,"%% XiBar         = %e\n",  XiBar);
	fprintf(OutFile,"%% AlphaBar      = %e\n",  AlphaBar);
	fprintf(OutFile,"%% Lambda@ExPeak = %e\n",  PairExactLogLikelihood(N, 1.0, 1.0, AlphaTrue, XiTrue, SigmaHat1, SigmaHat2, C, h1, h2));
	fflush(NULL);

        /* print xi and alpha grids to output files */
        for(i = 0; i < AlphaN; i++) {
                for(j = 0; j < XiN; j++) fprintf(OutFile,"%f\t", Xi[j]);
                fprintf(OutFile,"\n");
        }
        fprintf(OutFile,"\n\n");
        for(i = 0; i < AlphaN; i++){
                        for(j = 0; j < XiN; j++) fprintf(OutFile,"%f\t", Alpha[i]);
                        fprintf(OutFile,"\n");
        }
        fprintf(OutFile,"\n\n");

        /* main loop for calculating log Lambda on all grid points */
        for(i = 0; i < AlphaN; i++) {

                for(j = 0; j < XiN; j++) {

                        /* compute exact log Lambda */
                        LogLambda = PairExactLogLikelihood(N, SigmaBar1, SigmaBar2, Alpha[i], Xi[j], SigmaHat1, SigmaHat2, C, h1, h2);

                        /* write results to file */
                        fprintf(OutFile,"%f\t",LogLambda);
                }

                /* go to next line in matrix */
                fprintf(OutFile,"\n");
        }


	/* close the output file */
	fclose(OutFile);

}
