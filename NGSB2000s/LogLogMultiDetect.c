 /*
  * This program computes the exact and  first approximated log of the liklihood function
  * over a grid of (alpha,xi).  The input files must contain the xi values and alpha values.
  *
  * Inputs are: XiInFile, AlphaInFile, ExactOutFile, ApproxOutFile, XiTrue, AlphaTrue, SNR
  *
  * NOTE: We assume mean squared noise is 1.0
  * 
  */

#include "NGSB.h"

main(int argc, char *argv[])
{
	float	*XiF, *AlphaF;
	double	Sum1, Sum2, Sum3, Sum4, ExactLogLambda, ApproxLogLambda, XiTrue, AlphaTrue, Rho;
	double	 *h, *Xi, *Alpha;
	int	i, j, k, N, c, AlphaN=0, XiN=0;
	FILE 	*XiInFile, *AlphaInFile, *ExactOutFile, *ApproxOutFile;
	

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
	h = (double *)malloc(N*sizeof(double));
        XiF = (float *)malloc(XiN*sizeof(float));
        AlphaF = (float *)malloc(AlphaN*sizeof(double));
        Xi = (double *)malloc(XiN*sizeof(double));
	Alpha = (double *)malloc(AlphaN*sizeof(double));

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

	/* print headers to output files */
	fprintf(ExactOutFile,"%% Exact Log Lambda Results\n");
        fprintf(ExactOutFile,"%% XiTrue      = %f\n",XiTrue);
        fprintf(ExactOutFile,"%% AlphaTrue   = %f\n",AlphaTrue);
        fprintf(ExactOutFile,"%% SigmaTrue   = %f\n",1.0);
        fprintf(ExactOutFile,"%% Data Length = %d\n",N);
	fprintf(ExactOutFile,"%% SNR         = %f\n",Rho);
        fprintf(ExactOutFile,"%% Output format is meshgrid style: (xi,alpha,logLambda)\n");
	fprintf(ApproxOutFile,"%% Approx Log Lambda Results\n");
        fprintf(ApproxOutFile,"%% XiTrue      = %f\n",XiTrue);
        fprintf(ApproxOutFile,"%% AlphaTrue   = %f\n",AlphaTrue);
        fprintf(ApproxOutFile,"%% SigmaTrue   = %f\n",1.0);
        fprintf(ApproxOutFile,"%% Data Length = %d\n",N);
	fprintf(ApproxOutFile,"%% SNR         = %f\n",Rho);
        fprintf(ApproxOutFile,"%% Output format is meshgrid style: (xi,alpha,logLambda)\n");


	/* print xi and alpha grids to output files */
	for(i = 0; i < AlphaN; i++) {
		for(j = 0; j < XiN; j++) {
			fprintf(ExactOutFile,"%f\t", Xi[j]);
			fprintf(ApproxOutFile,"%f\t", Xi[j]);
		}
		fprintf(ExactOutFile,"\n");
		fprintf(ApproxOutFile,"\n");
	}
	fprintf(ExactOutFile,"\n\n");
	fprintf(ApproxOutFile,"\n\n");
        for(i = 0; i < AlphaN; i++){
	                for(j = 0; j < XiN; j++){
				fprintf(ExactOutFile,"%f\t", Alpha[i]);
				fprintf(ApproxOutFile,"%f\t", Alpha[i]);
	                }
	                fprintf(ExactOutFile,"\n");
			fprintf(ApproxOutFile,"\n");
        }
        fprintf(ExactOutFile,"\n\n");
	fprintf(ApproxOutFile,"\n\n");

	/* make artificial data */
	NGSBData(N, 1.0, AlphaTrue, XiTrue, h);


	/* main loop for calculating exact log Lambda */
	for(i = 0; i < AlphaN; i++) {

		for(j = 0; j < XiN; j++) {

		        /* compute exact log Lambda */
		        ExactLogLambda = ExactLogLikelihood(N, 1.0, Alpha[i], Xi[j], h); 

			/* write results to file */
			fprintf(ExactOutFile,"%f\t",ExactLogLambda);
		}
		
		/* go to next line in matrix */
		fprintf(ExactOutFile,"\n");
	}
	
	/* close the output file */
	fclose(ExactOutFile);

	/* main loop for calculating approximate log Lambda */
        for(i = 0; i < AlphaN; i++) {

		/* compute Alpha dependent sums */
		Approx1(N, 1.0, Alpha[i], h, &Sum1, &Sum2, &Sum3, &Sum4);

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

	/* close the output file */
	fclose(ApproxOutFile);


}
