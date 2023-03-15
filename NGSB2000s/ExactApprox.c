/* 
 * File: ExactApprox.c  (stand alone)     By: Steve Drasco
 *
 * Compares the exact and approximate non-Gaussian ML statistics.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		xi, alpha, sigma1=1.0, sigma2=1.0, rho;
	double		AlphaMin, AlphaMax, ExactStat, ApproxStat;
	double		*h1, *h2, *n1, *n2, *s, *AlphaA;
	double		**sum;
	int		i, N, NAlpha;
	char            *FileName;
	FILE		*OutFile;
	time_t          tp;
	struct  tm      *TimeStruct;
	

	/* check inputs */
	if(argc != 7) {
                printf("I need inputs: xi, alpha, N, AlphaMin, AlphaMax, NAlpha\n");
                return;
        }

        /* read inputs */
	xi = (double) atof(argv[1]);
	alpha = (double) atof(argv[2]);
        N = (int) atoi(argv[3]);
	AlphaMin = (double) atof(argv[4]);
	AlphaMax = (double) atof(argv[5]);
	NAlpha = (int) atoi(argv[6]);

        /* compute rho */
        rho = sqrt(((double) N)) * xi * alpha; 

        /* allocate memory */
        h1 = (double *) malloc(N*sizeof(double));
        h2 = (double *) malloc(N*sizeof(double));
        n1 = (double *) malloc(N*sizeof(double));
        n2 = (double *) malloc(N*sizeof(double));
        s = (double *) malloc(N*sizeof(double));
	FileName = (char *) malloc(100*sizeof(char));
        sum = (double **) malloc(9*sizeof(double *));
        for(i=0; i<9; i++) sum[i] = (double *)calloc(9,sizeof(double));  

	/* make output file name string */
	sprintf(FileName,"../results/X-");
        time(&tp);
        TimeStruct = localtime(&tp);
	strftime(FileName+13,7,"%d%m%y",TimeStruct);
	sprintf(FileName+19,"-%f-%f-%f-%d-%f-%f-%d",xi,alpha,rho,N,AlphaMin,AlphaMax,NAlpha);

        /* open output file */
	OutFile = fopen(FileName,"w");
	if( OutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",FileName);
                return;
        }

	/* print headers to output files */
	
	fprintf(OutFile,"%% Comparison of exact and approximate non-Gaussian ML statistics\n");
        fprintf(OutFile,"%% Sigma1      = 1.0\n");
        fprintf(OutFile,"%% Sigma2      = 1.0\n");
	fprintf(OutFile,"%% xi          = %e\n",xi);
	fprintf(OutFile,"%% alpha       = %e\n",alpha);
	fprintf(OutFile,"%% SNR         = %e\n",rho);
	fprintf(OutFile,"%% data length = %d\n",N);
	fprintf(OutFile,"%% AlphaMin    = %e\n",AlphaMin);
	fprintf(OutFile,"%% AlphaMax    = %e\n",AlphaMax);
	fprintf(OutFile,"%% NAlpha      = %d\n",NAlpha);
        fprintf(OutFile,"%% Output format is style: [alpha; exact; approx]\n");

	/* seed the random number generator */
	SeedRand();

	/* make artificial signal-full data */
	PairNGSBData(N, 1.0, 1.0, alpha, xi, h1, h2, n1, n2, s);

	/* compute sums */
	for(i=0; i<N; i++){
		sum[2][0] += Power(h1[i],2)*Power(h2[i],0);
                sum[0][2] += Power(h1[i],0)*Power(h2[i],2);
                sum[1][1] += Power(h1[i],1)*Power(h2[i],1);
                sum[4][0] += Power(h1[i],4)*Power(h2[i],0);
                sum[0][4] += Power(h1[i],0)*Power(h2[i],4);
                sum[3][1] += Power(h1[i],3)*Power(h2[i],1);
                sum[1][3] += Power(h1[i],1)*Power(h2[i],3);
                sum[2][2] += Power(h1[i],2)*Power(h2[i],2);
                sum[6][0] += Power(h1[i],6)*Power(h2[i],0);
                sum[0][6] += Power(h1[i],0)*Power(h2[i],6);
                sum[5][1] += Power(h1[i],5)*Power(h2[i],1);
                sum[1][5] += Power(h1[i],1)*Power(h2[i],5);
                sum[4][2] += Power(h1[i],4)*Power(h2[i],2);
                sum[2][4] += Power(h1[i],2)*Power(h2[i],4);
                sum[3][3] += Power(h1[i],3)*Power(h2[i],3);
                sum[8][0] += Power(h1[i],8)*Power(h2[i],0);
                sum[0][8] += Power(h1[i],0)*Power(h2[i],8);
                sum[7][1] += Power(h1[i],7)*Power(h2[i],1);
                sum[1][7] += Power(h1[i],1)*Power(h2[i],7);
                sum[6][2] += Power(h1[i],6)*Power(h2[i],2);
                sum[2][6] += Power(h1[i],2)*Power(h2[i],6);
                sum[5][3] += Power(h1[i],5)*Power(h2[i],3);
                sum[3][5] += Power(h1[i],3)*Power(h2[i],5);
                sum[4][4] += Power(h1[i],4)*Power(h2[i],4);
	}

	/* main loop over alpha array */
	for(i = 0; i < NAlpha; i++) {

		/* compute alpha */
		alpha = pow(10.0,log10(AlphaMin) + (((double) i) * (log10(AlphaMax) - log10(AlphaMin)) / (NAlpha-1.0) ) );

		/* compute exact statistic */
		ExactStat = PairExactLogLikelihood(N, sigma1, sigma2, alpha, xi, 
						   sum[2][0] / ((double) N), sum[0][2] / ((double) N), 
						   sum[1][1], h1, h2);

		/* compute approximate statistic */
		ApproxStat = PairApproxLogLikelihood(N, sigma1, sigma2, alpha, xi, sum);

		/* print to output file */
		fprintf(OutFile,"%e\t%e\t%e\n",alpha,ExactStat,ApproxStat);

        }

	/* close the output file */
	fclose(OutFile);

	/* exit */
	printf("ExactApprox done! We used:\n\n");
        printf("Sigma1      = 1.0\n");
        printf("Sigma2      = 1.0\n");
        printf("xi          = %e\n",xi);
        printf("alpha       = %e\n",(double) atof(argv[2]));
        printf("SNR         = %e\n",rho);
        printf("data length = %d\n",N);
        printf("AlphaMin    = %e\n",AlphaMin);
        printf("AlphaMin    = %e\n",AlphaMin);
        printf("NAlpha      = %d\n\n",NAlpha);

}
