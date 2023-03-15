/*
 * File: GaussianStatistic.c  by: Steve Drasco 
 *                                                           
 * This program computes the log of the two-detector likelihood  
 * function for a Gaussian stochastic background on Gaussian Noise.
 *      
 * NOTE:
 * 	name of this file is misleading.  uses 'maximum likelihood' approximated denominator integral (the no signal integral)
 * 	includes ALL noise terms (i.e. treats noise parameters as unknown)
 * 	SigmaHat1 = detector 1 output varriance =  SumOverj( h1[j] h1[j] )
 * 	SigmaHat2 = detector 2 output varriance =  SumOverj( h2[j] h2[j] )
 * 	c = cross correlation  = SumOverj( h1[j] h2[j] )
 */

#include "NGSB.h"

double GaussianStatistic(int N, double SigmaHat1, double SigmaHat2, double c)
{
	double	stat, alpha, epsilon;

	/* estimate signal varriance from cross correlation */
	alpha = c / ((double) N);

	/* compute statistic */
	epsilon = alpha*alpha / (SigmaHat1*SigmaHat2);
	stat = 1.0/(1.0 - epsilon);

	/* exit */
	return stat;
}
