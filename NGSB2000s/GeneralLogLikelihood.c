/*
 * File: GeneralLogLikelihood.c  by: Steve Drasco 
 *                                                           
 * This program estimates the log of the likelihood   
 * function for a non-gaussian stochastic background         
 * in Gaussian detector noise.
 *
 * This is an estimate because we compute the I_j / C_j 
 * by clipping and then doing lowest order (2 step) numerical 
 * integration.
 *
 */
 

#include "NGSB.h"

double GeneralLogLikelihood(int N, double sigma1, double sigma2, double alpha, double xi, double *h1, double *h2, double kappa, double *theta)
{
	double	A, B, C1, C2, D;
	double	pi = 3.14159265358979323846;
	double	LogLambda, S1 = 0.0, S2 = 0.0;
        int     j;

        /* check data pointer */
        if(h1 == NULL || h2 == NULL || theta == NULL) {
                printf("data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* makes constants to reduce flops */
	A = xi * sqrt(kappa/pi) * (1.0 - (2.0*alpha*kappa));
	B = alpha*kappa * ( (1.0/sigma1) + (1.0/sigma2)  );
	C1 = sqrt(2*kappa*alpha) / sigma1;
	C2 = sqrt(2*kappa*alpha) / sigma2;
	D = 1 - xi;

	/* main loop */
	for(j = 0; j < N; j++) {
		theta[j] = (C1*h1[j]) + (C2*h2[j]);
		S1 += theta[j];
		S2 += log( (A * (1 + exp(-2*theta[j]) + 2*exp(-B+theta[j]) + 2*exp(-0.5*(B + theta[j]))  + 2*exp(-0.5*B - 1.5*theta[j])  ) )   +   (D * exp(-B -theta[j])) );
	}
	LogLambda += N*B;

	/* exit */
	return LogLambda;
}
