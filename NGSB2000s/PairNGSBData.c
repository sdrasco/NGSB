/*
 *  File: PairNGSBData.c     by: Steve Drasco     
 *                                                           
 * This program generates detector data with a Non-Gaussian  
 * stochastic background characterized by user specified    
 * length (N), sigma, alpha, and xi                          
 *
 * built to conserve memory...
 * 
 */

#include "NGSB.h"

int PairNGSBData(int N, double sigma1, double sigma2, double alpha, double xi, double *h1, double *h2, double  *n1, double *n2, double *s)
{
        int     	i;
	extern long     idum;

        /* check output pointer */
        if(h1 == NULL || h2 == NULL) {
                printf("Output pointer was NULL.  Giving up ...\n");
                return 1;
        }

        /* make noise for detector 1 */
	Gaussian(sigma1,N,PMIN,n1);

	/* make noise for detector 2 */
        Gaussian(sigma2,N,PMIN,n2);

	/* make Gaussian background signal */
        Gaussian(alpha,N,PMIN,s);

        /* zero out parts of signal as output is filled */
        for (i = 0; i < N; i++) {
	
                if( ((double) ran2(&idum)) < xi ) {
			h1[i] = n1[i]+s[i]; 
			h2[i] = n2[i]+s[i]; 
                } else {
                        h1[i] = n1[i]; 
			h2[i] = n2[i]; 
                }
        }

	/* normal exit */
	return 0;
}
