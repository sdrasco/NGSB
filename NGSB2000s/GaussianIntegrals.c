/*
 * File: GaussianIntegrals.c  by: Steve Drasco  
 *
 * Given data, length, and alpha this program returns the I_j integrals
 * approximated by the specified parameter kappa.
 *                                                           
 *
 */
 

#include "NGSB.h"

double GeneralLogLikelihood(int N, double alpha, double *h, double *I, double kappa)
{
	double	pi = 3.14159265358979323846;
        int     j;

        /* check data pointer */
        if(h == NULL || I == NULL) {
                printf("Data or I pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* main loop for product */
	for(j = 0; j < N; j++)  I[j] = 

	/* exit */
	return LogLambda;
}
