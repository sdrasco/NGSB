/* 
 * File: ComparePerformance2.c  (stand alone)     By: Steve Drasco
 *
 * This compares performance of cross correlation method to performance of taylored 
 * maxium likelihood method ... and now the approximate version of the latter.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		Rho, f = 0.5, xi, alpha;
	double		LambdaMax, SigmaBar1, SigmaBar2, AlphaBar, XiBar;
	double		*h1, *h2, *n1, *n2, *s;
	double		**simplex, **sum;
	int		i, j, N, Trials, status = 1;
	char            *MFileName, *AFileName;
	FILE		*MOutFile, *AOutFile;
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
        MFileName = (char *) malloc(100*sizeof(char));
	AFileName = (char *) malloc(100*sizeof(char));
        simplex = (double **) malloc(5*sizeof(double *));
        for(i=0; i<5; i++) simplex[i] = (double *)malloc(4*sizeof(double));
        sum = (double **) malloc(9*sizeof(double *));
        for(i=0; i<9; i++) sum[i] = (double *)malloc(9*sizeof(double));

	/* make output file name strings */
	sprintf(MFileName,"../results/M-");
	sprintf(AFileName,"../results/A-");
        time(&tp);
        TimeStruct = localtime(&tp);
	strftime(MFileName+13,7,"%d%m%y",TimeStruct);
	strftime(AFileName+13,7,"%d%m%y",TimeStruct);
	sprintf(MFileName+19,"-%f-%f-%f-%d-%d-%d",xi,alpha,Rho,N,Trials,time(NULL));
	sprintf(AFileName+19,"-%f-%f-%f-%d-%d-%d",xi,alpha,Rho,N,Trials,time(NULL));

        /* open output files */
	MOutFile = fopen(MFileName,"w");
	AOutFile = fopen(AFileName,"w");
        if( MOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",MFileName);
                return;
        }
	if( AOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",AFileName);
                return;
        }

	/* print headers to output files */
	
        fprintf(MOutFile,"%% Non-Gaussian Optimal Statistic (with Gaussian noise varriance estimator) Results\n");
        fprintf(MOutFile,"%% Sigma1      = 1.0\n");
        fprintf(MOutFile,"%% Sigma2      = 1.0\n");
        fprintf(MOutFile,"%% xi          = %e\n",xi);
        fprintf(MOutFile,"%% alpha       = %e\n",alpha);
        fprintf(MOutFile,"%% SNR         = %e\n",Rho);
        fprintf(MOutFile,"%% data length = %d\n",N);
        fprintf(MOutFile,"%% trials      = %d\n",Trials);
        fprintf(MOutFile,"%% Output format is style: [NoSignalResponse_1 SignalResponse_1 NoSignalResponse_2 SignalResponse_2 ...]\n");
        fprintf(AOutFile,"%% APPROXIMATE! Non-Gaussian Optimal Statistic (with Gaussian noise varriance estimator) Results\n");
        fprintf(AOutFile,"%% Sigma1      = 1.0\n");
        fprintf(AOutFile,"%% Sigma2      = 1.0\n");
        fprintf(AOutFile,"%% xi          = %e\n",xi);
        fprintf(AOutFile,"%% alpha       = %e\n",alpha);
        fprintf(AOutFile,"%% SNR         = %e\n",Rho);
        fprintf(AOutFile,"%% data length = %d\n",N);
        fprintf(AOutFile,"%% trials      = %d\n",Trials);
        fprintf(AOutFile,"%% Output format is style: [NoSignalResponse_1 SignalResponse_1 NoSignalResponse_2 SignalResponse_2 ...]\n");
	fflush(NULL);


	/* seed the random number generator */
	SeedRand();
	
	/* main loop over trials */
	for(i = 0; i < Trials; i++) {

		/* make artificial signal-less data for detector 1 */
                Gaussian(1.0, N, PMIN, h1);

                /* make artificial signal-less data for detector 2 */
                Gaussian(1.0, N, PMIN, h2);

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
	        for(j=0; j<N; j++){
	                sum[2][0] += Power(h1[j],2)*Power(h2[j],0);
	                sum[0][2] += Power(h1[j],0)*Power(h2[j],2);
	                sum[1][1] += Power(h1[j],1)*Power(h2[j],1);
	                sum[4][0] += Power(h1[j],4)*Power(h2[j],0);
	                sum[0][4] += Power(h1[j],0)*Power(h2[j],4);
	                sum[3][1] += Power(h1[j],3)*Power(h2[j],1);
	                sum[1][3] += Power(h1[j],1)*Power(h2[j],3);
	                sum[2][2] += Power(h1[j],2)*Power(h2[j],2);
	                sum[6][0] += Power(h1[j],6)*Power(h2[j],0);
	                sum[0][6] += Power(h1[j],0)*Power(h2[j],6);
	                sum[5][1] += Power(h1[j],5)*Power(h2[j],1);
	                sum[1][5] += Power(h1[j],1)*Power(h2[j],5);
	                sum[4][2] += Power(h1[j],4)*Power(h2[j],2);
	                sum[2][4] += Power(h1[j],2)*Power(h2[j],4);
	                sum[3][3] += Power(h1[j],3)*Power(h2[j],3);
	                sum[8][0] += Power(h1[j],8)*Power(h2[j],0);
	                sum[0][8] += Power(h1[j],0)*Power(h2[j],8);
	                sum[7][1] += Power(h1[j],7)*Power(h2[j],1);
	                sum[1][7] += Power(h1[j],1)*Power(h2[j],7);
	                sum[6][2] += Power(h1[j],6)*Power(h2[j],2);
	                sum[2][6] += Power(h1[j],2)*Power(h2[j],6);
	               	sum[5][3] += Power(h1[j],5)*Power(h2[j],3);
        	        sum[3][5] += Power(h1[j],3)*Power(h2[j],5);
        	        sum[4][4] += Power(h1[j],4)*Power(h2[j],4);
	        }

                /* full maximum likelihood statistic */
                simplex[0][0] = 1.0;
                simplex[0][1] = 1.0;
                simplex[0][2] = alpha;
                simplex[0][3] = xi;
		status = 1;
                status = SimplexMaxLikelihood(&LambdaMax, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, h1, h2,
                                                             simplex, sum[2][0]/((double) N), sum[0][2]/((double) N), sum[1][1],
                                                             7500, 7500, 1e-4, 1e-4);
                fprintf(MOutFile,"%e\t%% %e\t%e\t%e\t%e \n",  LambdaMax,SigmaBar1,SigmaBar2,AlphaBar,XiBar);
		fflush(NULL);

                /* approximate maximum likelihood statistic! */
                simplex[0][0] = 1.0;
                simplex[0][1] = 1.0;
                simplex[0][2] = alpha;
                simplex[0][3] = xi;
                status = SimplexMaxLikelihood2(&LambdaMax, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, simplex, sum,
                                               7500, 7500, 1e-4, 1e-4);
                if(status == 0)
                        fprintf(AOutFile,"%e\t%% %e\t%e\t%e\t%e \n",  LambdaMax,SigmaBar1,SigmaBar2,AlphaBar,XiBar);
                else
                        fprintf(AOutFile,"%e\t%% %e\t%e\t%e\t%e %% WARNING !! DID NOT CONVERGE\n",  LambdaMax,SigmaBar1,SigmaBar2,AlphaBar,XiBar);
		fflush(NULL);
		

		/* make artificial signal-full data */
		PairNGSBData(N, 1.0, 1.0, alpha, xi, h1, h2, n1, n2, s);

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
                for(j=0; j<N; j++){
                        sum[2][0] += Power(h1[j],2)*Power(h2[j],0);
                        sum[0][2] += Power(h1[j],0)*Power(h2[j],2);
                        sum[1][1] += Power(h1[j],1)*Power(h2[j],1);
                        sum[4][0] += Power(h1[j],4)*Power(h2[j],0);
                        sum[0][4] += Power(h1[j],0)*Power(h2[j],4);
                        sum[3][1] += Power(h1[j],3)*Power(h2[j],1);
                        sum[1][3] += Power(h1[j],1)*Power(h2[j],3);
                        sum[2][2] += Power(h1[j],2)*Power(h2[j],2);
                        sum[6][0] += Power(h1[j],6)*Power(h2[j],0);
                        sum[0][6] += Power(h1[j],0)*Power(h2[j],6);
                        sum[5][1] += Power(h1[j],5)*Power(h2[j],1);
                        sum[1][5] += Power(h1[j],1)*Power(h2[j],5);
                        sum[4][2] += Power(h1[j],4)*Power(h2[j],2);
                        sum[2][4] += Power(h1[j],2)*Power(h2[j],4);
                        sum[3][3] += Power(h1[j],3)*Power(h2[j],3);
                        sum[8][0] += Power(h1[j],8)*Power(h2[j],0);
                        sum[0][8] += Power(h1[j],0)*Power(h2[j],8);
                        sum[7][1] += Power(h1[j],7)*Power(h2[j],1);
                        sum[1][7] += Power(h1[j],1)*Power(h2[j],7);
                        sum[6][2] += Power(h1[j],6)*Power(h2[j],2);
                        sum[2][6] += Power(h1[j],2)*Power(h2[j],6);
                        sum[5][3] += Power(h1[j],5)*Power(h2[j],3);
                        sum[3][5] += Power(h1[j],3)*Power(h2[j],5);
                        sum[4][4] += Power(h1[j],4)*Power(h2[j],4);
                }

                /* full maximum likelihood statistic */
                simplex[0][0] = 1.0;
                simplex[0][1] = 1.0;
                simplex[0][2] = alpha;
                simplex[0][3] = xi;
		status = 1;
                status = SimplexMaxLikelihood(&LambdaMax, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, h1, h2,
                                                             simplex, sum[2][0]/((double) N), sum[0][2]/((double) N), sum[1][1],
                                                             7500, 7500, 1e-4, 1e-4) ;
                fprintf(MOutFile,"%e\t%% %e\t%e\t%e\t%e \n",  LambdaMax,SigmaBar1,SigmaBar2,AlphaBar,XiBar);
		fflush(NULL);

                /* approximate maximum likelihood statistic! */
                simplex[0][0] = 1.0;
                simplex[0][1] = 1.0;
                simplex[0][2] = alpha;
                simplex[0][3] = xi;
                status = SimplexMaxLikelihood2(&LambdaMax, &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, simplex, sum,
                                               7500, 7500, 1e-4, 1e-4);
                if(status == 0)
                        fprintf(AOutFile,"%e\t%% %e\t%e\t%e\t%e \n",  LambdaMax,SigmaBar1,SigmaBar2,AlphaBar,XiBar);
                else
                        fprintf(AOutFile,"%e\t%% %e\t%e\t%e\t%e  %% WARNING !! DID NOT CONVERGE\n",  LambdaMax,SigmaBar1,SigmaBar2,AlphaBar,XiBar);
		fflush(NULL);

        }

	/* close the output files */
        fclose(MOutFile);
	fclose(AOutFile);

}
