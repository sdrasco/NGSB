/* 
 * File: PairLogLogMultiDetect.c  (stant alone)
 *
 * This program computes the exact log of the liklihood function FOR TWO DETECTORS
 * over a grid of (alpha,xi).  The input files must contain the xi values and alpha values.
 *
 * Inputs are: XiInFile, AlphaInFile, OutFile, XiTrue, AlphaTrue, SNR, Sigma1, Sigma2
 *
 * NOTE: we assume unit varriance in detector noise
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	float		*XiF, *AlphaF;
	double		LogLambda, ApproxLogLambda, XiTrue, AlphaTrue, Rho, Sum1, Sum2, Sum3, Sum4;
	double		SigmaHat1, SigmaHat2, C;
	double		*h1, *h2, *n1, *n2, *s, *Xi, *Alpha, f=0.5, *XiA, *AlphaA;
	int		i, j, k, N, c, AlphaN=0, XiN=0;
	int		GN = 50;
	FILE 		*XiInFile, *AlphaInFile, *ExactOutFile, *ApproxOutFile;

	/* check inputs */
	if(argc != 8) {
                printf("I need inputs: XiInFile, AlphaInFile, ExactOutFile, ApproxOutFile, XiTrue, AlphaTrue, SNR\n");
                return;
        }

        /* read inputs */
        XiTrue = (double) atof(argv[5]);
        AlphaTrue = (double) atof(argv[6]);
        Rho = (double) atof(argv[7]);

	/* compute N and abort if over 1e7 */
	N = (int) ( Rho * Rho / (XiTrue * XiTrue * AlphaTrue * AlphaTrue) );
	if(N >= (int) 1e7) {
		printf("Maximum N is 1e7, you tried: %d\nGiving up...\n",N);
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
        ExactOutFile = fopen(argv[3],"w");
	if( ExactOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[3]);
                return;
        }
        ApproxOutFile = fopen(argv[4],"w");
        if( ApproxOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[4]);
                return;
        }


	/* print headers to output files & screen */
	fprintf(ExactOutFile,"%% Two Detector Log Lambda Results\n");
        fprintf(ExactOutFile,"%% XiTrue      = %f\n",XiTrue);
        fprintf(ExactOutFile,"%% AlphaTrue   = %f\n",AlphaTrue);
        fprintf(ExactOutFile,"%% Sigma1      = %f\n",1.0);
        fprintf(ExactOutFile,"%% Sigma2      = %f\n",1.0);
        fprintf(ExactOutFile,"%% Data Length = %d\n",N);
	fprintf(ExactOutFile,"%% SNR         = %f\n",Rho);
        fprintf(ExactOutFile,"%% Output format is meshgrid style: (xi,alpha,logLambda)\n");
        fprintf(ApproxOutFile,"%% Two Detector Approximated Log Lambda Results to order xi^4 \n");
        fprintf(ApproxOutFile,"%% XiTrue      = %f\n",XiTrue);
        fprintf(ApproxOutFile,"%% AlphaTrue   = %f\n",AlphaTrue);
        fprintf(ApproxOutFile,"%% Sigma1      = %f\n",1.0);
        fprintf(ApproxOutFile,"%% Sigma2      = %f\n",1.0);
        fprintf(ApproxOutFile,"%% Data Length = %d\n",N);
        fprintf(ApproxOutFile,"%% SNR         = %f\n",Rho);
        fprintf(ApproxOutFile,"%% Output format is meshgrid style:  (xi,alpha,logLambda) \n");


        /* seed the random number generator */
	SeedRand();

        /* make artificial data */
        PairNGSBData(N, 1.0, 1.0, AlphaTrue, XiTrue, h1, h2, n1, n2, s);

	/* print xi and alpha grids to output files */
	for(i = 0; i < AlphaN; i++) {
		for(j = 0; j < XiN; j++) fprintf(ExactOutFile,"%f\t", Xi[j]);
		fprintf(ExactOutFile,"\n");
	}
	fprintf(ExactOutFile,"\n\n");
        for(i = 0; i < AlphaN; i++){
	                for(j = 0; j < XiN; j++) fprintf(ExactOutFile,"%f\t", Alpha[i]);
	                fprintf(ExactOutFile,"\n");
        }
        fprintf(ExactOutFile,"\n\n");
        for(i = 0; i < AlphaN; i++) {
                for(j = 0; j < XiN; j++) fprintf(ApproxOutFile,"%f\t", Xi[j]);
                fprintf(ApproxOutFile,"\n");
        }
        fprintf(ApproxOutFile,"\n\n");
        for(i = 0; i < AlphaN; i++){
                        for(j = 0; j < XiN; j++) fprintf(ApproxOutFile,"%f\t", Alpha[i]);
                        fprintf(ApproxOutFile,"\n");
        }
        fprintf(ApproxOutFile,"\n\n");


	/* main loop for calculating log Lambda */
	for(i = 0; i < AlphaN; i++) {

		for(j = 0; j < XiN; j++) {

			/* compute cross correlation statistic*/
			CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C);

		        /* compute exact log Lambda */
		        LogLambda = PairExactLogLikelihood(N, SigmaHat1 - Xi[j]*Alpha[i], SigmaHat2 - Xi[j]*Alpha[i], Alpha[i], Xi[j], 
							   SigmaHat1, SigmaHat2, C, h1, h2); 

			/* write results to file */
			fprintf(ExactOutFile,"%f\t",LogLambda);
		}
		
		/* go to next line in matrix */
		fprintf(ExactOutFile,"\n");
	}

	/* main loop for calculating approximated log lambda */
        for(i = 0; i < AlphaN; i++) {

                /* compute Alpha dependent sums */
                Approx2(N, Alpha[i], h1, h2, &Sum1, &Sum2, &Sum3, &Sum4);

                for(j = 0; j < XiN; j++) {

                        /* compute approximate log Lambda */
                        ApproxLogLambda = Xi[j] * Sum1
                                        - Xi[j] * Xi[j] * Sum2 / 2.0
                                        + Xi[j] * Xi[j] * Xi[j] * Sum3 / 3.0
                                        - Xi[j] * Xi[j] * Xi[j] * Xi[j] * Sum4 / 4.0;

                        /* write results to file */
                        fprintf(ApproxOutFile,"%f\t",ApproxLogLambda);
                }

                /* go to next line in matrix */
                fprintf(ApproxOutFile,"\n");
        }

	
	/* close the output files */
	fclose(ExactOutFile);
	fclose(ApproxOutFile);

}
