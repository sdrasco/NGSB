/*
 *  File: NonGaussianAlarmDismissal.c               By: Steve Drasco
 *
 *  Computes the distance (in false dismissal probability) between 
 *  a given alpha and a given and the critical alpha.
 *
 *
 *  NOTE: ran2 random number generator must be seeded by parent 
 *        function.
 *
 *  NOTE: uses the maximum likelihood detection statistic for a 
 *  	  non-Gaussian signal.
 *                                                           
 */

#include "NGSB.h"

void NonGaussianDistance(int N, double Sigma1, double Sigma2, double xi, 
			double alpha, int NTrials, double *Lambda0, 
			double *Lambda1, double *h1, double *h2, double *n1, 
			double *n2, double *s, double ThreshDelta, double PFAc, 
			double PFDc, double *distance)
{
	double		SigmaHat1, SigmaHat2,  C, stat, MeanThresh;
	double		DeltaPFA = 1.0, MinThresh = 10.0, MaxThresh = -10.0;
	double		SigmaBar1, SigmaBar2, AlphaBar, XiBar;
	double          **simplex;
	int		i;
	int		counter = 0;
	

	/* allocate memory for the simplex (will leak memory - but hardly any) */
        simplex = (double **) malloc(5*sizeof(double *));
        for(i=0; i<5; i++) simplex[i] = (double *)malloc(4*sizeof(double));

	/* main loop over trials */
	for(i = 0; i < NTrials; i++) {

		/* make artificial signal-less data for detector 1 */
                Gaussian(1.0, N, PMIN, h1);

                /* make artificial signal-less data for detector 2 */
                Gaussian(1.0, N, PMIN, h2);

		/* compute detection statistics */
                simplex[0][0] = 1.0;
                simplex[0][1] = 1.0;
                simplex[0][2] = alpha;
                simplex[0][3] = xi;
		CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C);
                SimplexMaxLikelihood(&Lambda0[i], &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, h1, h2,
                                     simplex, SigmaHat1, SigmaHat2, C, 7500, 7500, 1e-4, 1e-4);
		
		/* make artificial signal-full data */
		PairNGSBData(N, 1.0, 1.0, alpha, xi, h1, h2, n1, n2, s);

                /* compute detection statistics */
		simplex[0][0] = 1.0;
                simplex[0][1] = 1.0;
                simplex[0][2] = alpha;
                simplex[0][3] = xi;
		CrossCorr(N, h1, h2, &SigmaHat1, &SigmaHat2, &C);
                SimplexMaxLikelihood(&Lambda1[i], &SigmaBar1, &SigmaBar2, &AlphaBar, &XiBar, N, h1, h2,
                                     simplex, SigmaHat1, SigmaHat2, C, 7500, 7500, 1e-4, 1e-4);

		/* update minimum and maximum threshold (greater than 1) */
		if(Lambda0[i] > MaxThresh) MaxThresh = Lambda0[i];
		if(Lambda1[i] > MaxThresh) MaxThresh = Lambda1[i];
		if(Lambda0[i] < MinThresh && Lambda0[i] != 1.0) MinThresh = Lambda0[i];
		if(Lambda1[i] < MinThresh && Lambda1[i] != 1.0) MinThresh = Lambda1[i];

        }

	/* find critical threshold PFA(CritThresh) = PFAc */
	while(fabs(DeltaPFA) > ThreshDelta && counter < 1000000) { /* note 1000000 is a quick fix to what is not likely a problem */

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
