/*
 * File: BurstStatistic.c  by: Steve Drasco 
 *                                                           
 * This program computes a naive statistic for detecting
 * spike like bursts in one detector:
 *
 * [output] = max_k  | h_1^k |
 * 
 */
 

#include "NGSB.h"

double BurstStatistic(int N, double *h)
{
	double	stat = 0.0;
        int     j;

        /* check data pointer */
        if(h == NULL ) {
                printf("data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* compute statistic */
	for(j = 0; j < N; j++) if(stat < fabs(h[j])) stat = fabs(h[j]); 

	/* normal exit */
	return stat;
}
