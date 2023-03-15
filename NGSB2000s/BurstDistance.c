/*
 *  File: BurstDistance.c               By: Steve Drasco
 *
 *  Computes the distance (in false dismissal probability) between 
 *  a given alpha and a given and the critical alpha.
 *
 *
 *  NOTE: ran2 random number generator must be seeded by parent 
 *        function.
 *
 *  NOTE: uses the naive burst-detecting statistic.
 *                                                           
 */

#include "NGSB.h"

void BurstDistance(int N, double Sigma1, double Sigma2, double xi, 
                       double alpha, int NTrials, double *Lambda0, 
		       double *Lambda1, double *h, double *n,double *s, 
                       double ThreshDelta, double PFAc, 
		       double PFDc, double *distance)
{
	double		MeanThresh;
	double		DeltaPFA = 1.0, MinThresh = 10.0, MaxThresh = -10.0;
	int		i, j;
	int		counter = 0;
	

	/* main loop over trials */
	for(i = 0; i < NTrials; i++) {

		/* make artificial signal-less data */
                Gaussian(1.0, N, PMIN, h);

		/* compute detection statistics */
		Lambda0[i] = BurstStatistic(N, h);
		
		/* make artificial signal-full data 
		PairNGSBData(N, 1.0, 1.0, alpha, xi, h, h, n, n, s);
		*/

		/* ################################################################# */
		/* BEGIN: clip from PairNGSBData */

		        /* make noise */
		        Gaussian(1.0,N,PMIN,n);

		        /* make Gaussian background signal */
		        Gaussian(alpha,N,PMIN,s);

		        /* zero out parts of signal as output is filled */
		        for (j = 0; j < N; j++) {

		                if( ((double) ran2(&idum)) < xi ) {
		                        h[j] = n[j]+s[j];
        		        } else {
                		        h[j] = n[j];
               		 	}
        		}

		/* END: clip from PairNGSBData */
		/* ################################################################# */

                /* compute detection statistics */
		Lambda1[i] = BurstStatistic(N, h);

		/* update minimum and maximum threshold (greater than 1) */
		if(Lambda0[i] > MaxThresh) MaxThresh = Lambda0[i];
		if(Lambda1[i] > MaxThresh) MaxThresh = Lambda1[i];
		if(Lambda0[i] < MinThresh && Lambda0[i] != 1.0) MinThresh = Lambda0[i];
		if(Lambda1[i] < MinThresh && Lambda1[i] != 1.0) MinThresh = Lambda1[i];

        }

	/* find critical threshold PFA(CritThresh) = PFAc */  
	while(fabs(DeltaPFA) > ThreshDelta && counter < 1000000) {  /* the 1000000 is a quick fix to what is not likely a problem */

		/* reset DeltaPFA */
		DeltaPFA = PFAc;

		/* find middle threshold */
		MeanThresh = 0.5 * (MaxThresh + MinThresh);

		/* compute deltaPFA  = PFAc - PFA for middle threshold */
		for(i = 0; i < NTrials; i++) if(Lambda0[i] > MeanThresh) DeltaPFA -= 1.0 / ((double) NTrials);

		/* update bracket */
		if(DeltaPFA > 0 ) MaxThresh = MeanThresh;
		else MinThresh = MeanThresh;

		/* increment counter */
		counter++;

	}

	/* print an error if we did not get below threshold */
	if(counter >= 1000000) {
		printf("COUNTER OVERFLOW ERROR!!\n");
		printf("Did not find a threshold with |DeltaPFA| < %f.\n", ThreshDelta);
		printf("Best threshold had DeltaPFA = %f", DeltaPFA);
	}
	

	/* compute distance = PFDc - PFD */
	*distance = PFDc;
	for(i = 0; i < NTrials; i++) if(Lambda1[i] < MeanThresh) *distance -= 1.0 / ((double) NTrials);

}
