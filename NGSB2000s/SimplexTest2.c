/* 
 * File: SimplexTest2.c     by: Steve Drasco 
 *                                                           
 * STAND ALONE VERSION                                       
 * Just a forum to test the simplex maximium likelihood (3) routine.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		SigmaBar1, SigmaBar2, XiBar, AlphaBar, Lambda, HatLambda;
	double		Sigma1, Sigma2, xi, alpha, *n1, *n2, *s, *h1, *h2;
	double		**simplex, **sum;
	int		i, j, N;

	/* check inputs */
	if(argc != 4) {
                printf("I need inputs: xi, alpha, N\n");
                return;
        } 

        /* read inputs */
	Sigma1 = Sigma2 = 1.0;
        xi = (double) atof(argv[1]);
        alpha = (double) atof(argv[2]);
        N = (int) atoi(argv[3]);
	printf("Signal to noise ratio: %f\n",sqrt(((double) N)) * xi * alpha);

	/* allocate memory */
        h1 = (double *) malloc(N*sizeof(double));
        h2 = (double *) malloc(N*sizeof(double));
        n1 = (double *) malloc(N*sizeof(double));
        n2 = (double *) malloc(N*sizeof(double));
        s = (double *) malloc(N*sizeof(double));
	simplex = (double **)malloc(4*sizeof(double *));
	for(i=0; i<4; i++) simplex[i] = (double *)malloc(3*sizeof(double));
        sum = (double **) malloc(9*sizeof(double *));
        for(i=0; i<9; i++) sum[i] = (double *)malloc(9*sizeof(double));

	/* make artificial signal-less data for detector 1 
	Gaussian(Sigma1, N, PMIN, h1); */

	/* make artificial signal-less data for detector 2 
	Gaussian(Sigma2, N, PMIN, h2); */

	/* make artificial signal-full data */
	PairNGSBData(N, Sigma1, Sigma2, alpha, xi, h1, h2, n1, n2, s);

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

        /* fill simplex */
        simplex[0][0] = 1.2345; /* Sigma1 */
        simplex[0][1] = 678.90; /* Sigma2 */
        simplex[0][2] = 0.1234; /* xi */

	/* alpha test  THIS WORKS: peak = 10.002
	printf("GetAlpha test: alpha = %f\n",GetAlpha(-15,200,0.1,3)); */

	/* compute statistic */
	/* printf("returned: %d\n(0 means converged)\n",SimplexMaxLikelihood3(&SigmaBar1, &SigmaBar2, &XiBar, N, simplex, sum, 7500, 7500, 1e-4, 1e-4)); 
	LambdaMax = PairApproxLogLikelihood2(N, SigmaBar1, SigmaBar2, &AlphaBar, XiBar, sum); 
	LambdaMax = PairApproxLogLikelihood2(N, 1.0, 1.0, &AlphaBar, xi, sum);
	printf("Lambda = %e at Sigma1 = %e\tSigma2 = %e\txi = %e\talpha = %e\n",LambdaMax,SigmaBar1,SigmaBar2,XiBar,AlphaBar); */

	HatLambda = PairApproxLogLikelihood2(N, 1.0, 1.0, &AlphaBar, xi, sum);
	Lambda = PairApproxLogLikelihood(N, 1.0, 1.0, alpha, xi, sum);
	printf("Lambda = %e\t HatLambda = %e \n",Lambda, HatLambda);
	
}
