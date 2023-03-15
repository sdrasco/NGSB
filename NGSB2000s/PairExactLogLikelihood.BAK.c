/*
 * File: ExactLogLikelihood.c  by: Steve Drasco 
 *                                                           
 * This program computes the *Exact* log of the likelihood  
 * function for a non-gaussian stochastic background       
 * characterized by user specified sigma1, sigma2, alpha, and xi.  
 * User also gives detector data & length.               
 * 
 * Remember this program is specifically for a pair of detectors!
 * 
 * NOTE:
 * 	T2 = N * detector 1 output varriance =  SumOverj( h1[j] h1[j] )
 * 	T3 = N * detector 2 output varriance =  SumOverj( h2[j] h2[j] )
 * 	T4 = cross correlation  = SumOverj( h1[j] h2[j] )
 */

#include "NGSB.h"

double PairExactLogLikelihood(int N, double sigma1, double sigma2, double alpha, double xi, double T2, double T3, double T4, double *h1, double *h2)
{
	double	A, B, C, D, E, T1 = 0.0, LogLambda;
        int     i;

        /* check data pointer */
        if(h1 == NULL || h2 == NULL) {
                printf("Data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* compute constants A, B, & C for faster product */
	A = xi / (  sqrt(  (alpha/sigma1) + (alpha/sigma2) + 1.0  )    );
	B = 1.0 / ( 2.0 * sigma1 * sigma1 * ( (1.0/sigma1) + (1.0/sigma2) + (1.0/alpha) )  );
	C = 1.0 / ( 2.0 * sigma2 * sigma2 * ( (1.0/sigma1) + (1.0/sigma2) + (1.0/alpha) )  );
	D = 1.0 / ( sigma1 + sigma2 + (sigma1*sigma2/alpha)  );
	E = 1.0 - xi;
	
	/* main loop for product */
	for(i = 0; i < N; i++) {
		T1 += log(  A + (E * exp(- B*h1[i]*h1[i] - C*h2[i]*h2[i] - D*h1[i]*h2[i]))  );
	}
	T2 *= B;
	T3 *= C;
	T4 *= D;

	/* exit */
	LogLambda = T1 + T2 + T3 + T4;
	return LogLambda;
}
