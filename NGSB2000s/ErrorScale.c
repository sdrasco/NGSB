/* 
 * File: ErrorScale.c  (stant alone)
 *
 * This program computes error in the xi-estimator for a variety
 * of N values.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		XiTrue, AlphaTrue, RhoTrue, SigmaHat1, SigmaHat2;
	double		SigmaBar1, SigmaBar2, LambdaMax, AlphaBar;
	double		ExpectedXi, DeltaXi, DeltaXiOverXi;
	double		*h1, *h2, *n1, *n2, *s, *XiBar;
	double		**simplex, **sum;
	double		f = 0.5;
	int		i, j, k, z, MinN, MaxN, DeltaN, Nsteps;
	int             *N;
	int		AlphaN = 0, XiN = 0, GN = 50, AveScale = 10, SmallAveScale = 5;
	char		*OutFileName;
	FILE 		*OutFile;
        time_t          tp;
        struct  tm      *TimeStruct;


	/* check inputs */
	if(argc != 6) {
                printf("I need inputs: MinN, DeltaN, MaxN, XiTrue, RhoTrue\n");
                return;
        }

        /* read inputs */
	MinN = (int) atoi(argv[1]);
	DeltaN = (int) atoi(argv[2]);
        MaxN = (int) atoi(argv[3]);
        XiTrue = (double) atof(argv[4]); 
        RhoTrue = (double) atof(argv[5]); 

        /* compute N info */
        Nsteps = ((int) ((MaxN - MinN) / DeltaN)) + 1;
        MaxN = MinN + Nsteps*DeltaN;

	/* abort if N is too big */
	if(MaxN >= (int) 5e6) {
		printf("Maximum N is 5e6, you tried: %d\nGiving up...\n",MaxN);
		return;
	}

	/* allocate memory */
	h1 = (double *) malloc(MaxN*sizeof(double));
	h2 = (double *) malloc(MaxN*sizeof(double));
	n1 = (double *) malloc(MaxN*sizeof(double));
	n2 = (double *) malloc(MaxN*sizeof(double));
	OutFileName = (char *) malloc(100*sizeof(char));
	s = (double *) malloc(MaxN*sizeof(double));
        XiBar = (double *)malloc(AveScale*sizeof(double));
	simplex = (double **) malloc(5*sizeof(double *));
	for(i=0; i<5; i++) simplex[i] = (double *)malloc(4*sizeof(double));
        sum = (double **) malloc(9*sizeof(double *));
        for(i=0; i<9; i++) sum[i] = (double *)malloc(9*sizeof(double));
        N = (int *) malloc(Nsteps*sizeof(int));
        for(i=0; i<Nsteps; i++) N[i] = MinN + i*DeltaN;


        /* open output files */
	sprintf(OutFileName,"../results/ErrorScale-");
        time(&tp);
        TimeStruct = localtime(&tp);
        strftime(OutFileName+22,7,"%d%m%y",TimeStruct);
        sprintf(OutFileName+28,"-%d-%d-%d-%d",MinN,MaxN,DeltaN,time(NULL));
        OutFile = fopen(OutFileName,"w");
	if( OutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",OutFileName);
                return;
        }

	/* print headers to output files & screen */
        fprintf(OutFile,"%% XiTrue      = %f\n",XiTrue);
        fprintf(OutFile,"%% AlphaTrue   = %f\n",AlphaTrue);
        fprintf(OutFile,"%% Sigma1      = %f\n",1.0);
        fprintf(OutFile,"%% Sigma2      = %f\n",1.0);
	fprintf(OutFile,"%% RhoTrue     = %f\n",RhoTrue);
        fprintf(OutFile,"%% Output format: [DeltaXi / XiTrue , N; ...]\n");

        /* seed the random number generator */
	SeedRand();

	/* main loop */
        for(i = 0; i < Nsteps; i++) {

		/* small average to smooth results */
		DeltaXiOverXi = 0.0;
		for(z = 0; z < SmallAveScale; z++){

			/* loop for averaging to get means and fluctuations*/
			AveScale = 10;
			for(j = 0; j < AveScale; j++){
	
			        /* make artificial data */
				PairNGSBData(N[i], 1.0, 1.0, AlphaTrue, XiTrue, h1, h2, n1, n2, s); 
					
				/* zero all sums */
				sum[2][0] = 0.0;
				sum[0][2] = 0.0;
				sum[1][1] = 0.0;
		
				/* compute sums */
				for(k=0; k<N[i]; k++){
				sum[2][0] += h1[k]*h1[k];
				sum[0][2] += h2[k]*h2[k];
				sum[1][1] += h1[k]*h2[k];
				}
		
				/* compute AlphaTrue[i] */
				AlphaTrue = RhoTrue / (   XiTrue * sqrt( ((double) N[i]) )   );
		
				/* compute statistic (no need to store value) */
				simplex[0][0] = 1.0;
				simplex[0][1] = 1.0;
				simplex[0][2] = AlphaTrue;
				simplex[0][3] = XiTrue;
				SimplexMaxLikelihood(&LambdaMax, &SigmaBar1, &SigmaBar2, &AlphaBar, XiBar+j, N[i], h1, h2,
						simplex, sum[2][0]/((double) N[i]), sum[0][2]/((double) N[i]), sum[1][1],
						7500, 7500, 1e-4, 1e-4);
			}

			/* compute DeltaXi/Xi */
			ExpectedXi = DeltaXi = 0.0;
			for(j = 0; j < AveScale; j++) ExpectedXi += XiBar[j] / ((double) AveScale);
			for(j = 0; j < AveScale; j++) DeltaXi += (XiBar[j] - ExpectedXi)*(XiBar[j] - ExpectedXi);
			DeltaXi = sqrt(DeltaXi / ((double) AveScale));
			DeltaXiOverXi += DeltaXi/(ExpectedXi * ((double)  SmallAveScale));
		}

		/* output results */
		fprintf(OutFile,"%e\t%d\n", DeltaXiOverXi, N[i]);
		fflush(NULL);
	
        }

	/* close the output file */
	fclose(OutFile);

}
