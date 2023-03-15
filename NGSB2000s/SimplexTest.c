/* 
 * File: SimplexTest.c     by: Steve Drasco 
 *                                                           
 * STAND ALONE VERSION                                       
 * Just a forum to test the simplex maximium likelihood routine.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		h1, h2, MaxLogLambda, SigmaBar1, SigmaBar2, XiBar, AlphaBar;
	double		**simplex;
	int		i;

	/* check inputs */
	if(argc != 5) {
                printf("I need inputs: Sigma1, Sigma2, xi, and alpha\n");
                return;
        }

	/* allocate memory */
	simplex = (double **)malloc(5*sizeof(double *));
	for(i=0; i<5; i++) simplex[i] = (double *)malloc(4*sizeof(double));

	/* fill simplex */
	simplex[0][0] = atof(argv[1]);
        simplex[0][1] = atof(argv[2]);
        simplex[0][2] = atof(argv[3]);
        simplex[0][3] = atof(argv[4]);

	while(SimplexMaxLikelihood(&MaxLogLambda, &SigmaBar1, &SigmaBar2, &XiBar, &AlphaBar, 1.0, &h1, &h2, simplex, 1.0, 1.0, 1.0, 1000, 1000, 1e-4, 1e-4));

	printf("\tMaxLogLambda = %e at: Sigma1 = %e\tSigma2 = %e\txi = %e\talpha = %e\n",MaxLogLambda,SigmaBar1,SigmaBar2,XiBar,AlphaBar);
}
