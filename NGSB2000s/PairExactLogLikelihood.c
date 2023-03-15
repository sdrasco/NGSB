/*
 * File: ExactLogLikelihood.c  by: Steve Drasco 
 *                                                           
 * This program computes the log of the two-detector likelihood  
 * function for a non-Gaussian stochastic background        
 * (characterized by user specified sigma1, sigma2, alpha, and xi)
 * on Gaussian Noise.  
 *      
 * NOTE:
 * 	name of this file is misleading.  uses 'maximum likelihood' approximated denominator integral (the no signal integral)
 * 	includes ALL noise terms (i.e. treats noise parameters as unknown)
 * 	SigmaHat1 = detector 1 output varriance =  SumOverj( h1[j] h1[j] ) / N
 * 	SigmaHat2 = detector 2 output varriance =  SumOverj( h2[j] h2[j] ) / N
 * 	c = cross correlation  = SumOverj( h1[j] h2[j] )
 */

#include "NGSB.h"

double PairExactLogLikelihood(int N, double Sigma1, double Sigma2, double alpha, double xi, double SigmaHat1, double SigmaHat2, double c, double *h1, double *h2)
{
	double	A, B = 0.0, C, D, E, F, G;
        int     j;

        /* check data pointer */
        if(h1 == NULL || h2 == NULL) {
                printf("Data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* compute constants to speed up the sum */
	A  = SigmaHat1/(Sigma1*Sigma1)  + SigmaHat2/(Sigma2*Sigma2) + 2.0*c/( ((double) N)*Sigma1*Sigma2 );
	A /= 1.0/Sigma1 + 1.0/Sigma2 + 1.0/alpha;
	A += 2.0 - SigmaHat1/Sigma1  - SigmaHat2/Sigma2;
	A *= 0.5*((double) N);
	C = xi / sqrt( Sigma1*Sigma2/(SigmaHat1*SigmaHat2) + Sigma1*alpha/(SigmaHat1*SigmaHat2) + Sigma2*alpha/(SigmaHat1*SigmaHat2) );
	D = (1.0 - xi) * sqrt( SigmaHat1*SigmaHat2/(Sigma1*Sigma2) );
	E = - 0.5 / (Sigma1 * Sigma1 * (1.0/Sigma1 + 1.0/Sigma2 + 1.0/alpha));
	F = - 0.5 / (Sigma2 * Sigma2 * (1.0/Sigma1 + 1.0/Sigma2 + 1.0/alpha));
	G = - 1.0 / (Sigma1 * Sigma2 * (1.0/Sigma1 + 1.0/Sigma2 + 1.0/alpha));

	/* main loop for sum */
	for(j = 0; j < N; j++) {
		B += log( C + D * exp(E*h1[j]*h1[j] + F*h2[j]*h2[j] + G*h1[j]*h2[j])  );
	}

	/* exit */
	return A+B;
}
