/*-----------------------------------------------------------*/
/* File: NGSBData.c     by: Steve Drasco 27 July 2000        */
/*                                                           */
/* This program generates detector data with a Non-Gaussian  */
/* stochastic background characterized by user specified     */
/* length (N), sigma, alpha, and xi                          */
/*-----------------------------------------------------------*/

#include "NGSB.h"

int NGSBData(int N, double sigma, double alpha, double xi, double *h)
{
        double  	*n, *s;
        int     	i;
	extern long     idum;

        /* check output pointer */
        if(h == NULL) {
                printf("Output pointer was NULL.  Giving up ...\n");
                return 1;
        }

        /* allocate memory for noise and signal */
        n = (double *) malloc(N*sizeof(double));
        s = (double *) malloc(N*sizeof(double));

        /* make noise */
        Gaussian(sigma,N,PMIN,n);

	/* make Gaussian background signal */
        Gaussian(alpha,N,PMIN,s);

        /* zero out parts of signal as output is filled */
        for (i = 0; i < N; i++) {
	
                if( ((double) ran2(&idum)) < xi ) {
                        h[i] = n[i]+s[i];
                } else {
                        h[i] = n[i];
                }
        }

	/* check out free */	
	free(n);
	free(s);

	/* normal exit */
	return 0;
}
