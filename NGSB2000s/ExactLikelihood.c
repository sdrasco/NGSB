/* 
 * File: ExactLikelihood.c     by: Steve Drasco 
 *                                                           
 * This program computes the *Exact* likelihood function     
 * (lambda) for a non-gaussian stochastic background         
 * characterized by user specified sigma, alpha, and xi.     
 * User also gives detector data & length.                   
 *
 */

#include "NGSB.h"

float ExactLikelihood(int N, float sigma, float alpha, float xi, float *h)
{
	float	A, B, C, lambda;
        int     i;

        /* check data pointer */
        if(h == NULL) {
                printf("Data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* compute constants A, B, & C for faster product */ 
	A = 1.0 - xi;
	B = xi / sqrt(   1 + (  alpha / sigma  )   );
	C = 1 /  (   2 * ( sigma + ((sigma*sigma) / alpha)  )   );

	/* initialize lambda */
	lambda = 1.0;

	/* main loop for product */
	for(i = 0; i < N; i++) lambda *= A + ( B * exp(C * h[i] * h[i]) );

	/* exit */
	return lambda;
}
