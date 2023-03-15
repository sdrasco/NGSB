/*
 *  File: GaussianAlarmDismissal.c               By: Steve Drasco
 *
 *  Computes the distance (in false dismissal probability) between 
 *  a given alpha and a given and the critical alpha.
 *
 *
 *  NOTE: ran2 random number generator must be seeded by parent 
 *        function.
 *
 *  NOTE: uses the maximum likelihood detection statistic for a 
 *  	  Gaussian signal.
 *                                                           
 */

#include "NGSB.h"

void GaussianDistance(int N, double Sigma1, double Sigma2, double xi, 
                       double alpha, int NTrials, double *Lambda0, 
		       double *Lambda1, double *h1, double *h2, double *n1, 
		       double *n2, double *s, double ThreshDelta, double PFAc, 
		       double PFDc, double *distance)
{
	double		SigmaHat1, SigmaHat2,  C, stat, MeanThresh;
	double		DeltaPFA = 1.0, MinThresh = 10.0, MaxThresh = -10.0;
	int		i;
	int		counter = 0;
	

	/* main loop over trials */
	for(i = 0; i < NTrials; i++) {

		/* make artificial signal-less data for detector 1 */
                Gaussian(1.0, N, PMIN, h1);

                /* make artificial signal-less data for detector 2 */
                Gaussian(1.0, N, PMIN, h2);

		/* compute detection statistics */
		stat = CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C);
		Lambda0[i] = 1.0/(1.0 - (stat*stat));
		
		/* make artificial signal-full data */
		PairNGSBData(N, 1.0, 1.0, alpha, xi, h1, h2, n1, n2, s);

                /* compute detection statistics */
                stat = CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C);
		Lambda1[i] = 1.0/(1.0 - (stat*stat));

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
