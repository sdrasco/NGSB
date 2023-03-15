/*-----------------------------------------------------------*/
/* File: ExactLogLikelihood.c  by: Steve Drasco 27 July 2000 */
/*                                                           */
/* This program computes the *Exact* log of the likelihood   */
/* function for a non-gaussian stochastic background         */
/* characterized by user specified sigma, alpha, and xi.     */
/* User also gives detector data & length.                   */
/*-----------------------------------------------------------*/

#include "NGSB.h"

double ExactLogLikelihood(int N, double sigma, double alpha, double xi, double *h)
{
	double	A, B, C, LogLambda;
        int     i;

        /* check data pointer */
        if(h == NULL) {
                printf("Data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* compute constants A, B, & C for faster product */ 
	A = 1.0 - xi;
	B = xi / sqrt(   1.0 + (  (alpha*alpha) / (sigma*sigma)  )   );
	C = 1.0 /  (   2.0 * ( (sigma*sigma) + ((sigma*sigma*sigma*sigma) / (alpha*alpha))  )   );

	/* initialize lambda */
	LogLambda = 0.0;

	/* main loop for product */
	for(i = 0; i < N; i++) LogLambda += C * h[i] * h[i]  +  log(  B + ( A * exp(-C * h[i] * h[i]) )  );

	/* exit */
	return LogLambda;
}
