/* 
 * File: timing.c  (stand alone)     By: Steve Drasco
 *
 * This program finds out how computation time scales with N for 
 * rho_detectable vs. N plots. 
 *
 * NOTE: time recorded in seconds.
 *
 * 	 Scale RunTime by X to estimate time cost for a false 
 *       dismissal versus false alarm curve made from X trials.
 *
 *       Scale result by Y to estimate time cost for rho_detectable
 *       computation based on Y iterations.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		Rho, f = 0.5, xi = 1, alpha = 0.01;
	double		stat, SigmaBar1, SigmaBar2, AlphaBar, XiBar, RunTime;
	double		*h1, *h2, *n1, *n2, *s;
	double		**simplex, **sum;
	int		i, j, k, MinN, DeltaN, MaxN, Nsteps;
	int		*N;
	char            *CFileName, *AFileName, *MFileName;
	FILE		*COutFile, *AOutFile, *MOutFile;
	time_t          tp;
	struct  tm      *TimeStruct;
	

	/* check inputs */
	if(argc != 4) {
                printf("I need inputs: MinN DeltaN MaxN \n");
                return;
        }

        /* read inputs */
	MinN = (int) atoi(argv[1]);
	DeltaN = (int) atoi(argv[2]);
        MaxN = (int) atoi(argv[3]);

	/* compute N info */
	Nsteps = ((int) ((MaxN - MinN) / DeltaN)) + 1;
	MaxN = MinN + Nsteps*DeltaN;

        /* allocate memory */
        h1 = (double *) malloc(MaxN*sizeof(double));
        h2 = (double *) malloc(MaxN*sizeof(double));
        n1 = (double *) malloc(MaxN*sizeof(double));
        n2 = (double *) malloc(MaxN*sizeof(double));
        s = (double *) malloc(MaxN*sizeof(double));
        MFileName = (char *) malloc(100*sizeof(char));
	AFileName = (char *) malloc(100*sizeof(char));
	CFileName = (char *) malloc(100*sizeof(char));
        simplex = (double **) malloc(5*sizeof(double *));
        for(i=0; i<5; i++) simplex[i] = (double *)malloc(4*sizeof(double));
        sum = (double **) malloc(9*sizeof(double *));
        for(i=0; i<9; i++) sum[i] = (double *)malloc(9*sizeof(double));
	N = (int *) malloc(Nsteps*sizeof(int));
	for(i=0; i<Nsteps; i++) N[i] = MinN + i*DeltaN;

	/* make output file name strings */
	sprintf(MFileName,"../results/M-Timing-");
	sprintf(AFileName,"../results/A-Timing-");
	sprintf(CFileName,"../results/C-Timing-");
        time(&tp);
        TimeStruct = localtime(&tp);
	strftime(MFileName+20,7,"%d%m%y",TimeStruct);
	strftime(AFileName+20,7,"%d%m%y",TimeStruct);
	strftime(CFileName+20,7,"%d%m%y",TimeStruct);
	sprintf(MFileName+26,"-%d-%d-%d-%d",MinN,MaxN,DeltaN,time(NULL));
	sprintf(AFileName+26,"-%d-%d-%d-%d",MinN,MaxN,DeltaN,time(NULL));
	sprintf(CFileName+26,"-%d-%d-%d-%d",MinN,MaxN,DeltaN,time(NULL));

        /* open output files */
	MOutFile = fopen(MFileName,"w");
	AOutFile = fopen(AFileName,"w");
	COutFile = fopen(CFileName,"w");
	if( MOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",MFileName);
                return;
        }
        if( AOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",AFileName);
                return;
        }
	if( COutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",CFileName);
                return;
        }

	/* print headers to output files */
	
        fprintf(MOutFile,"%% Maximum Likelihood results\n");
        fprintf(MOutFile,"%% Sigma1      = 1.0\n");
        fprintf(MOutFile,"%% Sigma2      = 1.0\n");
        fprintf(MOutFile,"%% xi          = 0.01\n");
        fprintf(MOutFile,"%% alpha       = 0.01\n");
        fprintf(MOutFile,"%% Output format is style: N RunTime\n");
        fprintf(AOutFile,"%% approximate statistic results\n");
        fprintf(AOutFile,"%% Sigma1      = 1.0\n");
        fprintf(AOutFile,"%% Sigma2      = 1.0\n");
        fprintf(AOutFile,"%% xi          = 0.01\n");
        fprintf(AOutFile,"%% alpha       = 0.01\n");
        fprintf(AOutFile,"%% Output format is style: N RunTime\n");
	fprintf(COutFile,"%% cross correlation results\n");
        fprintf(COutFile,"%% Sigma1      = 1.0\n");
        fprintf(COutFile,"%% Sigma2      = 1.0\n");
        fprintf(COutFile,"%% xi          = 0.01\n");
        fprintf(COutFile,"%% alpha       = 0.01\n");
	fflush(NULL);

	/* seed the random number generator */
	SeedRand();

        /* time the full statistic */
        for(i = 0; i < Nsteps; i++) {

                /* start time (we average over 10) */
		RunTime = ((double) time(NULL));
                for(j = 0; j < 10; j++) {

                        /* make artificial signal-less data for detector 1 */
                        Gaussian(1.0, N[i], PMIN, h1);

                        /* make artificial signal-less data for detector 2 */
                        Gaussian(1.0, N[i], PMIN, h2);

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

                        /* compute statistic (no need to store value) */
	                simplex[0][0] = 1.0;
	                simplex[0][1] = 1.0;
			simplex[0][2] = alpha;
			simplex[0][3] = xi;
			SimplexMaxLikelihood(&stat, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N[i], h1, h2,
                                             simplex, sum[2][0]/((double) N[i]), sum[0][2]/((double) N[i]), sum[1][1],
                                             7500, 7500, 1e-4, 1e-4);

                        /* make artificial signal-full data */
			PairNGSBData(N[i], 1.0, 1.0, alpha, xi, h1, h2, n1, n2, s);

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

			/* compute statistic (no need to store value) */
                        simplex[0][0] = 1.0;
                        simplex[0][1] = 1.0;
                        simplex[0][2] = alpha;
                        simplex[0][3] = xi;
                        SimplexMaxLikelihood(&stat, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N[i], h1, h2,
                                             simplex, sum[2][0]/((double) N[i]), sum[0][2]/((double) N[i]), sum[1][1],
                                             7500, 7500, 1e-4, 1e-4);

                }

                /* compute average run time */
                RunTime = (((double) time(NULL)) - RunTime) / 10;

                /* output results */
                fprintf(MOutFile,"%d\t%e\n", N[i], RunTime);
                fflush(NULL);
        }

        /* time the approximate statistic*/
        for(i = 0; i < Nsteps; i++) {

                /* start time (we average over 10) */
                RunTime = ((double) time(NULL));
                for(j = 0; j < 10; j++) {

                        /* make artificial signal-less data for detector 1 */
                        Gaussian(1.0, N[i], PMIN, h1);

                        /* make artificial signal-less data for detector 2 */
                        Gaussian(1.0, N[i], PMIN, h2);

	                /* zero all sums */
	                sum[2][0] = 0.0;
	                sum[0][2] = 0.0;
	                sum[1][1] = 0.0;
	                sum[4][0] = 0.0;
	                sum[0][4] = 0.0;
	                sum[3][1] = 0.0;
		        sum[1][3] = 0.0;
        	        sum[2][2] = 0.0;
       	        	sum[6][0] = 0.0;
                	sum[0][6] = 0.0;
                	sum[5][1] = 0.0;
                	sum[1][5] = 0.0;
                	sum[4][2] = 0.0;
                	sum[2][4] = 0.0;
                	sum[3][3] = 0.0;
                	sum[8][0] = 0.0;
                	sum[0][8] = 0.0;
                	sum[7][1] = 0.0;
                	sum[1][7] = 0.0;
                	sum[6][2] = 0.0;
                	sum[2][6] = 0.0;
                	sum[5][3] = 0.0;
               		sum[3][5] = 0.0;
                	sum[4][4] = 0.0;

                	/* compute sums */
                	for(k=0; k<N[i]; k++){
                        	sum[2][0] += Power(h1[k],2)*Power(h2[k],0);
                        	sum[0][2] += Power(h1[k],0)*Power(h2[k],2);
                        	sum[1][1] += Power(h1[k],1)*Power(h2[k],1);
                        	sum[4][0] += Power(h1[k],4)*Power(h2[k],0);
                        	sum[0][4] += Power(h1[k],0)*Power(h2[k],4);
                        	sum[3][1] += Power(h1[k],3)*Power(h2[k],1);
                       		sum[1][3] += Power(h1[k],1)*Power(h2[k],3);
                       		sum[2][2] += Power(h1[k],2)*Power(h2[k],2);
                        	sum[6][0] += Power(h1[k],6)*Power(h2[k],0);
                        	sum[0][6] += Power(h1[k],0)*Power(h2[k],6);
                        	sum[5][1] += Power(h1[k],5)*Power(h2[k],1);
                        	sum[1][5] += Power(h1[k],1)*Power(h2[k],5);
                        	sum[4][2] += Power(h1[k],4)*Power(h2[k],2);
                        	sum[2][4] += Power(h1[k],2)*Power(h2[k],4);
                        	sum[3][3] += Power(h1[k],3)*Power(h2[k],3);
                        	sum[8][0] += Power(h1[k],8)*Power(h2[k],0);
                        	sum[0][8] += Power(h1[k],0)*Power(h2[k],8);
                        	sum[7][1] += Power(h1[k],7)*Power(h2[k],1);
                        	sum[1][7] += Power(h1[k],1)*Power(h2[k],7);
                        	sum[6][2] += Power(h1[k],6)*Power(h2[k],2);
                        	sum[2][6] += Power(h1[k],2)*Power(h2[k],6);
                        	sum[5][3] += Power(h1[k],5)*Power(h2[k],3);
                        	sum[3][5] += Power(h1[k],3)*Power(h2[k],5);
                        	sum[4][4] += Power(h1[k],4)*Power(h2[k],4);
                	}



                        /* compute statistic (no need to store value) */
	                simplex[0][0] = 1.0;
	                simplex[0][1] = 1.0;
	                simplex[0][2] = alpha;
	                simplex[0][3] = xi;
	                SimplexMaxLikelihood2(&stat, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N[i], simplex, sum,
                                               7500, 7500, 1e-4, 1e-4);

                        /* make artificial signal-full data */
                        PairNGSBData(N[i], 1.0, 1.0, alpha, xi, h1, h2, n1, n2, s);

                        /* zero all sums */
                        sum[2][0] = 0.0;
                        sum[0][2] = 0.0;
                        sum[1][1] = 0.0;
                        sum[4][0] = 0.0;
                        sum[0][4] = 0.0;
                        sum[3][1] = 0.0;
                        sum[1][3] = 0.0;
                        sum[2][2] = 0.0;
                        sum[6][0] = 0.0;
                        sum[0][6] = 0.0;
                        sum[5][1] = 0.0;
                        sum[1][5] = 0.0;
                        sum[4][2] = 0.0;
                        sum[2][4] = 0.0;
                        sum[3][3] = 0.0;
                        sum[8][0] = 0.0;
                        sum[0][8] = 0.0;
                        sum[7][1] = 0.0;
                        sum[1][7] = 0.0;
                        sum[6][2] = 0.0;
                        sum[2][6] = 0.0;
                        sum[5][3] = 0.0;
                        sum[3][5] = 0.0;
                        sum[4][4] = 0.0;

                        /* compute sums */
                        for(k=0; k<N[i]; k++){
                                sum[2][0] += Power(h1[k],2)*Power(h2[k],0);
                                sum[0][2] += Power(h1[k],0)*Power(h2[k],2);
                                sum[1][1] += Power(h1[k],1)*Power(h2[k],1);
                                sum[4][0] += Power(h1[k],4)*Power(h2[k],0);
                                sum[0][4] += Power(h1[k],0)*Power(h2[k],4);
                                sum[3][1] += Power(h1[k],3)*Power(h2[k],1);
                                sum[1][3] += Power(h1[k],1)*Power(h2[k],3);
                                sum[2][2] += Power(h1[k],2)*Power(h2[k],2);
                                sum[6][0] += Power(h1[k],6)*Power(h2[k],0);
                                sum[0][6] += Power(h1[k],0)*Power(h2[k],6);
                                sum[5][1] += Power(h1[k],5)*Power(h2[k],1);
                                sum[1][5] += Power(h1[k],1)*Power(h2[k],5);
                                sum[4][2] += Power(h1[k],4)*Power(h2[k],2);
                                sum[2][4] += Power(h1[k],2)*Power(h2[k],4);
                                sum[3][3] += Power(h1[k],3)*Power(h2[k],3);
                                sum[8][0] += Power(h1[k],8)*Power(h2[k],0);
                                sum[0][8] += Power(h1[k],0)*Power(h2[k],8);
                                sum[7][1] += Power(h1[k],7)*Power(h2[k],1);
                                sum[1][7] += Power(h1[k],1)*Power(h2[k],7);
                                sum[6][2] += Power(h1[k],6)*Power(h2[k],2);
                                sum[2][6] += Power(h1[k],2)*Power(h2[k],6);
                                sum[5][3] += Power(h1[k],5)*Power(h2[k],3);
                                sum[3][5] += Power(h1[k],3)*Power(h2[k],5);
                                sum[4][4] += Power(h1[k],4)*Power(h2[k],4);
                        }

                        /* compute statistic (no need to store value) */
                        simplex[0][0] = 1.0;
                        simplex[0][1] = 1.0;
                        simplex[0][2] = alpha;
                        simplex[0][3] = xi;
                        SimplexMaxLikelihood2(&stat, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N[i], simplex, sum,
                                               7500, 7500, 1e-4, 1e-4);

                }

                /* compute average run time */
		RunTime = (((double) time(NULL)) - RunTime) / 10;

                /* output results */
                fprintf(AOutFile,"%d\t%e\n", N[i], RunTime);
                fflush(NULL);
        }
	
	/* time the cross-correlation statistic*/
	for(i = 0; i < Nsteps; i++) {

		/* start time (we average over 10) */
		RunTime = ((double) time(NULL));
		for(j = 0; j < 10; j++) {

			/* make artificial signal-less data for detector 1 */
	                Gaussian(1.0, N[i], PMIN, h1);
	
	                /* make artificial signal-less data for detector 2 */
	                Gaussian(1.0, N[i], PMIN, h2);
	
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
	
	                /* statistic (no need to store value) */
	                if(sum[1][1] > 0)
				stat = sum[1][1] / sqrt(sum[2][0] * sum[0][2]);
	                else
				stat = 0.0;
	
			/* make artificial signal-full data */
			PairNGSBData(N[i], 1.0, 1.0, alpha, xi, h1, h2, n1, n2, s);
	
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
	
                	/* compute statistic (no need to store value) */
                	if(sum[1][1] > 0)
                       		stat = sum[1][1] / sqrt(sum[2][0] * sum[0][2]);
                	else
                        	stat = 0.0;
		}

		/* compute average run time */
		RunTime = (((double) time(NULL)) - RunTime) / 10;
	
		/* output results */	
                fprintf(COutFile,"%d\t%e\n", N[i], RunTime);
		fflush(NULL);
        }

	/* close the output files */
	fclose(MOutFile);
        fclose(AOutFile);
	fclose(COutFile);

}
