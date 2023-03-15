/* 
 * File: StandardPerformance.c  (stand alone)     By: Steve Drasco
 *
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	float   	*ThreshF;
	double		Rho, SigmaHat1, SigmaHat2, C, xi, alpha;
	double		*Thresh, **h1, **h2, *n1, *n2, *s;
	int		i, j, k, N, c, ThreshN=0, FalseAlarm, FalseDismisal, Trials;
	FILE 		*OutFile, *InFile;
	

	/* check inputs */
	if(argc != 7) {
                printf("I need inputs: InFile, OutFile, xi, alpha, SNR, Trials\n");
                return;
        }

        /* read inputs */
	xi = (double) atof(argv[3]);
	alpha = (double) atof(argv[4]);
        Rho = (double) atof(argv[5]);
	Trials = (int) atol(argv[6]);

	/* compute N */
	N = (int) ( Rho * Rho / (xi * xi * alpha * alpha) );
	if(N >= (int) 1e7 || N <= 1) {
		printf("N is unreasonable.  Maximum N is 1e7, you tried: %d\nGiving up...\n",N);
		return;
	}
	printf("starting a run with N = %d\n",N);

        /* open input file and count entries */
        InFile = fopen(argv[1],"r");
        if(InFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[1]);
                return;
        }
	while( (c = getc(InFile)) != EOF) if( c == '\n') ThreshN++;
        if(N*Trials >= (int) 1.7e7) {
                printf("Maximum N*Trials is 1e7, you tried: %d\nGiving up...\n",N*Trials);
                return;
        }

        /* close and reopen input file */
        fclose(InFile);
	InFile = fopen(argv[1],"r");

        /* allocate memory */
        h1 = (double **) malloc(Trials*sizeof(double));
        h2 = (double **) malloc(Trials*sizeof(double *));
	for(i = 0; i < Trials; i++) {
		h1[i] = (double *) malloc(N*sizeof(double));
		h2[i] = (double *) malloc(N*sizeof(double));
	}
        ThreshF = (float *) malloc(ThreshN*sizeof(float));
        Thresh = (double *) malloc(ThreshN*sizeof(double));
        n1 = (double *) malloc(N*sizeof(double));
        n2 = (double *) malloc(N*sizeof(double));
        s = (double *) malloc(N*sizeof(double));


        /* fill thresh LEARN HOW TO SCAN IN DOUBLES */
        for(i=0; i < ThreshN; i++) {
                fscanf(InFile,"%e", ThreshF+i);
                Thresh[i] = (double) ThreshF[i];
        }

	/* close input */
	fclose(InFile);

        /* open output files */
        OutFile = fopen(argv[2],"w");
	if( OutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",argv[3]);
                return;
        }

	/* print headers to output files & screen */
	fprintf(OutFile,"%% Cross Correlation Results\n");
        fprintf(OutFile,"%% Sigma1      = 1.0\n");
        fprintf(OutFile,"%% Sigma2      = 1.0\n");
	fprintf(OutFile,"%% xi          = %e\n",xi);
	fprintf(OutFile,"%% alpha       = %e\n",alpha);
	fprintf(OutFile,"%% SNR         = %e\n",Rho);
	fprintf(OutFile,"%% data length = %d\n",N);
	fprintf(OutFile,"%% trials      = %d\n",Trials);
        fprintf(OutFile,"%% Output format is style: [false alarm probability ... false dismisal probability]\n");

	/* seed the random number generator */
	SeedRand();

	/* make signal-less data */
	for(i = 0; i < Trials; i++) {
	                        
			/* make artificial data for detector 1 */
                        Gaussian(1.0, N, PMIN, h1[i]);

                        /* make artificial data for detector 2 */
                        Gaussian(1.0, N, PMIN, h2[i]);

	}


	/* false alarm loop */
	for(i = 0; i < ThreshN; i++) {

		FalseAlarm = 0;
		
		/* compare thresh to statistic for each of the data sets */
		for(j = 0; j < Trials; j++) {
		        if(CrossCorr(N, h1[j], h2[j], &SigmaHat1, &SigmaHat2, &C) > Thresh[i]) FalseAlarm++;
		}

		/* write results to file */
		fprintf(OutFile,"%e\n", ((double) FalseAlarm) / ((double) Trials));

	}

        /* make signal data */
        for(i = 0; i < Trials; i++) {

		/* make artificial data */
		PairNGSBData(N, 1.0, 1.0, alpha, xi, h1[i], h2[i], n1, n2, s);
                                
        }


        /* false dismisal */
        for(i = 0; i < ThreshN; i++) {

                FalseDismisal = 0;

                for(j = 0; j < Trials; j++) {
                        /* compare to thresh exact log Lambda */
                        if(CrossCorr(N, h1[j], h2[j], &SigmaHat1, &SigmaHat2, &C) < Thresh[i]) FalseDismisal++;
                }

                /* write results to file */
                fprintf(OutFile,"%e\n", ((double) FalseDismisal) / ((double) Trials));

        }


	/* close the output file */
	fclose(OutFile);

}
