/*
 * File: CrossCorr.c  by: Steve Drasco 
 *                                                           
 * This program computes the standard cross correlation statistic
 * for detection stochastic background signals.
 *
 */
 

#include "NGSB.h"

double CrossCorr(int N, double *h1, double *h2, double *SigmaHat1, double *SigmaHat2, double *C)
{
	double	stat;
        int     j;

        /* check data pointer */
        if(h1 == NULL || h2 == NULL) {
                printf("data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* initialize terms */
	*SigmaHat1 = *SigmaHat2 = *C = 0.0;

	/* compute statistic */
	for(j = 0; j < N; j++) {
		*SigmaHat1 += h1[j]*h1[j];
		*SigmaHat2 += h2[j]*h2[j];
		*C         += h1[j]*h2[j];
	}
	stat = *C / sqrt(*SigmaHat1*(*SigmaHat2));

	/* estimate detector varriances */
	*SigmaHat1 /= ((double) N);
	*SigmaHat2 /= ((double) N);

	/* exit */
	if(*C > 0) return stat;
	else return 0.0;
}
