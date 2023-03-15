/* 
 * File: BurstDistribution.c  (stand alone)     By: Steve Drasco
 *
 * Computes the --> ML <--  statistic many times.  All results are kept so as to be 
 * turned into a distribution.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		*h1, *h2;
	double          **simplex;
	double          SigmaHat1, SigmaHat2,  C, stat;
	double          SigmaBar1, SigmaBar2, AlphaBar, XiBar;
	int		i, N, NTrials;
	char            *FileName;
	FILE		*OutFile;
	time_t          tp;
	struct  tm      *TimeStruct;

	/* check inputs */
	if(argc != 3) {
		printf("I need inputs: N, NTrials\n");
                return;
        }

        /* read inputs */
	N = (int) atol(argv[1]);
	NTrials = (int) atol(argv[2]);

        /* allocate memory */
        h1 = (double *) malloc(N*sizeof(double));
	h2 = (double *) malloc(N*sizeof(double));
	FileName = (char *) malloc(100*sizeof(char));
	simplex = (double **) malloc(5*sizeof(double *));
        for(i=0; i<5; i++) simplex[i] = (double *)malloc(4*sizeof(double));

	/* make output file name strings */
	sprintf(FileName,"../results/Distribution-");
        time(&tp);        
        TimeStruct = localtime(&tp);
	strftime(FileName+24,7,"%d%m%y",TimeStruct);
	sprintf(FileName+30,"-%d-%d-%d",N,NTrials,time(NULL));

        /* open output files */
	OutFile = fopen(FileName,"w");
	if( OutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",FileName);
                return;
        }

	/* print headers to output files */
	fprintf(OutFile,"%% ML statistic values for:\n");
        fprintf(OutFile,"%% Sigma  = 1.0\n");
	fprintf(OutFile,"%% xi = alpha = 0");
	fprintf(OutFile,"%% N = %d\n",N);
	fprintf(OutFile,"%% trials      = %d\n",NTrials);

	/* seed the random number generator */
	SeedRand();


	/* main loop  */
	for(i = 0; i < NTrials; i++) {

		/* make signal-less data */
                Gaussian(1.0,N,PMIN,h1);
		Gaussian(1.0,N,PMIN,h2);

		/* compute statistic */
                simplex[0][0] = 1.0;
                simplex[0][1] = 1.0;
                simplex[0][2] = 1.0;
                simplex[0][3] = 0.5;
                CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C);
                SimplexMaxLikelihood(&stat, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, h1, h2,
                                     simplex, SigmaHat1, SigmaHat2, C, 7500, 7500, 1e-4, 1e-4);

		/* output results */
		fprintf(OutFile,"%e\n", stat );
		fflush(NULL);

	} 

	/* close the output file */
	fclose(OutFile);

}
