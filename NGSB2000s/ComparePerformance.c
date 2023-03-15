/* 
 * File: ComparePerformance.c  (stand alone)     By: Steve Drasco
 *
 * This compares performance of cross correlation method to performance of taylored 
 * maxium likelihood method.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		Rho, SigmaHat1, SigmaHat2, f = 0.5, xi, alpha;
	double		LambdaMax, SigmaBar1, SigmaBar2, AlphaBar, XiBar, C;
	double		*h1, *h2, *n1, *n2, *s;
	double		**simplex;
	int		i, j, N, Trials, status = 1;
	char            *GFileName, *CFileName, *MFileName;
	FILE		*GOutFile, *COutFile, *MOutFile;
	time_t          tp;
	struct  tm      *TimeStruct;
	

	/* check inputs */
	if(argc != 5) {
                printf("I need inputs: xi, alpha, N, Trials\n");
                return;
        }

        /* read inputs */
	xi = (double) atof(argv[1]);
	alpha = (double) atof(argv[2]);
        N = (int) atoi(argv[3]);
	Trials = (int) atol(argv[4]);
	
        /* compute rho */
        Rho = sqrt(((double) N)) * xi * alpha; 

        /* allocate memory */
        h1 = (double *) malloc(N*sizeof(double));
        h2 = (double *) malloc(N*sizeof(double));
        n1 = (double *) malloc(N*sizeof(double));
        n2 = (double *) malloc(N*sizeof(double));
        s = (double *) malloc(N*sizeof(double));
	GFileName = (char *) malloc(100*sizeof(char));
	CFileName = (char *) malloc(100*sizeof(char));
        MFileName = (char *) malloc(100*sizeof(char));
        simplex = (double **) malloc(5*sizeof(double *));
        for(i=0; i<5; i++) simplex[i] = (double *)malloc(4*sizeof(double));

	/* make output file name strings */
	sprintf(GFileName,"../results/G-");
	sprintf(CFileName,"../results/C-");
	sprintf(MFileName,"../results/M-");
        time(&tp);
        TimeStruct = localtime(&tp);
	strftime(GFileName+13,7,"%d%m%y",TimeStruct);
	strftime(CFileName+13,7,"%d%m%y",TimeStruct);
	strftime(MFileName+13,7,"%d%m%y",TimeStruct);
	sprintf(GFileName+19,"-%f-%f-%f-%d-%d-%d",xi,alpha,Rho,N,Trials,time(NULL));
	sprintf(CFileName+19,"-%f-%f-%f-%d-%d-%d",xi,alpha,Rho,N,Trials,time(NULL));
	sprintf(MFileName+19,"-%f-%f-%f-%d-%d-%d",xi,alpha,Rho,N,Trials,time(NULL));

        /* open output files */
	GOutFile = fopen(GFileName,"w");
	COutFile = fopen(CFileName,"w");
	MOutFile = fopen(MFileName,"w");
	if( GOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",GFileName);
                return;
        }
        if( COutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",CFileName);
                return;
        }
        if( MOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",MFileName);
                return;
        }

	/* print headers to output files */
	
	fprintf(GOutFile,"%% Gaussian Optimal Statistic Results\n");
        fprintf(GOutFile,"%% Sigma1      = 1.0\n");
        fprintf(GOutFile,"%% Sigma2      = 1.0\n");
	fprintf(GOutFile,"%% xi          = %e\n",xi);
	fprintf(GOutFile,"%% alpha       = %e\n",alpha);
	fprintf(GOutFile,"%% SNR         = %e\n",Rho);
	fprintf(GOutFile,"%% data length = %d\n",N);
	fprintf(GOutFile,"%% trials      = %d\n",Trials);
        fprintf(GOutFile,"%% Output format is style: [NoSignalResponse_1 SignalResponse_1 NoSignalResponse_2 SignalResponse_2 ...]\n");
        fprintf(COutFile,"%% Cross Correlation Results\n");
        fprintf(COutFile,"%% Sigma1      = 1.0\n");
        fprintf(COutFile,"%% Sigma2      = 1.0\n");
        fprintf(COutFile,"%% xi          = %e\n",xi);
        fprintf(COutFile,"%% alpha       = %e\n",alpha);
        fprintf(COutFile,"%% SNR         = %e\n",Rho);
        fprintf(COutFile,"%% data length = %d\n",N);
        fprintf(COutFile,"%% trials      = %d\n",Trials);
        fprintf(COutFile,"%% Output format is style: [NoSignalResponse_1 SignalResponse_1 NoSignalResponse_2 SignalResponse_2 ...]\n");
        fprintf(MOutFile,"%% Non-Gaussian Optimal Statistic (with Gaussian noise varriance estimator) Results\n");
        fprintf(MOutFile,"%% Sigma1      = 1.0\n");
        fprintf(MOutFile,"%% Sigma2      = 1.0\n");
        fprintf(MOutFile,"%% xi          = %e\n",xi);
        fprintf(MOutFile,"%% alpha       = %e\n",alpha);
        fprintf(MOutFile,"%% SNR         = %e\n",Rho);
        fprintf(MOutFile,"%% data length = %d\n",N);
        fprintf(MOutFile,"%% trials      = %d\n",Trials);
        fprintf(MOutFile,"%% Output format is style: [NoSignalResponse_1 SignalResponse_1 NoSignalResponse_2 SignalResponse_2 ...]\n");
	fflush(NULL);

	/* seed the random number generator */
	SeedRand();
	
	/* main loop over trials */
	for(i = 0; i < Trials; i++) {

		/* make artificial signal-less data for detector 1 */
                Gaussian(1.0, N, PMIN, h1);

                /* make artificial signal-less data for detector 2 */
                Gaussian(1.0, N, PMIN, h2);

		/* compute detection statistics */
		fprintf(GOutFile,"%e\n", CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C));
		if(C > 0.0)
			fprintf(COutFile,"%e\n", C);
		else
			fprintf(COutFile,"%e\n", 0.0);
		fflush(NULL);
	        simplex[0][0] = 1.0;
	        simplex[0][1] = 1.0;
	        simplex[0][2] = alpha;
	        simplex[0][3] = xi;
		while(status) status = SimplexMaxLikelihood(&LambdaMax, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, h1, h2, 
							    simplex, SigmaHat1, SigmaHat2, C, 7500, 7500, 1e-4, 1e-4);
		status = 1;
                fprintf(MOutFile,"%e\t%% %e\t%e\t%e\t%e \n",  LambdaMax,SigmaBar1,SigmaBar2,AlphaBar,XiBar);
		fflush(NULL);

		/* make artificial signal-full data */
		PairNGSBData(N, 1.0, 1.0, alpha, xi, h1, h2, n1, n2, s);

                /* compute detection statistics */
                fprintf(GOutFile,"%e\n", CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C));
                if(C > 0.0)
                        fprintf(COutFile,"%e\n", C);
                else
                        fprintf(COutFile,"%e\n", 0.0);
		fflush(NULL);
                simplex[0][0] = 1.0;
                simplex[0][1] = 1.0;
                simplex[0][2] = alpha;
                simplex[0][3] = xi;
                while(status) status = SimplexMaxLikelihood(&LambdaMax, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, h1, h2,
                                                            simplex, SigmaHat1, SigmaHat2, C, 7500, 7500, 1e-4, 1e-4);
		status = 1;
                fprintf(MOutFile,"%e\t%% %e\t%e\t%e\t%e \n",  LambdaMax,SigmaBar1,SigmaBar2,AlphaBar,XiBar);
		fflush(NULL);

        }

	/* close the output files */
	fclose(GOutFile);
        fclose(COutFile);
        fclose(MOutFile);

}
