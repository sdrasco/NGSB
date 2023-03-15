/*
 * File: SimplexMaxLikelihood.c  by: Steve Drasco 
 *                                                           
 * This program estimates the maximum of the likelihood function 
 * for detection of non-Gaussian stochastic background signals.  It 
 * treats variance of noise (1,2), variance of signal, and Gaussianity
 * as unknowns.  Hence there is a 4-D maximization.  The maximization 
 * is done with a simplex (Nelder and Mead). We speed things up by using 
 * good guesses.  
 *
 * NOTE: simplex is a 5x4 matrix with starting guess in first row with
 * 	 format (sigma1, sigma2, xi, alpha).
 *
 * NOTE: I'm basicly ripping the simplex bit off from the routine fminsearch 
 * 	 in Matlab.
 *
 */
 

#include "NGSB.h"

int SimplexMaxLikelihood(double *MaxLogLambda, double *SigmaBar1, double *SigmaBar2, double *AlphaBar, 
			double *XiBar, int N, double *h1, double *h2, double **simplex,
			double SigmaHat1, double SigmaHat2, double C, int MaxIterate, 
			int MaxCall, double TolX,double TolF)
{
	double	*NegLogLambda, *AveragePoint, *ReflectionPoint, *ExpansionPoint, *ContractionPoint;
	double	UsualDelta, ZeroTermDelta, rho, chi, psi, gamma, x, f;
	double	AverageValue, ReflectionValue, ExpansionValue, ContractionValue;
        int     i, j, call, shrink, iterate;

	/* check data pointer */
        if(h1 == NULL || h2 == NULL || simplex == NULL) {
                printf("MaxLikelihood reports: An input pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* initialize parameters */
	rho = 1.0;
	chi = 2.0;
	psi = 0.5;
	gamma = 0.5;
	shrink = 0;
	UsualDelta = 5e-2;
	ZeroTermDelta = 2.5e-4;

	/* allocate memory */
	NegLogLambda = (double *)malloc(5*sizeof(double));
	AveragePoint = (double *)malloc(4*sizeof(double));
	ReflectionPoint = (double *)malloc(4*sizeof(double));
	ExpansionPoint = (double *)malloc(4*sizeof(double));
	ContractionPoint = (double *)malloc(4*sizeof(double));

	/* give birth to the simplex */
	for(i=0; i<4; i++) {
		for(j=0; j<4; j++){ 
			if(i == j) { 
				if(simplex[0][j] == 0.0) {
					simplex[i+1][j] = ZeroTermDelta;
				} else { 
					simplex[i+1][j] = (1.0 + UsualDelta) * simplex[0][j];	
				}
			} else {
			simplex[i+1][j] =  simplex[0][j];
			}
		}
	}

	/* evaluate likelihood on simplex */	
	for(i=0; i<5; i++) {
		BoundaryCheck(simplex[i]);
		NegLogLambda[i] = - PairExactLogLikelihood(N, simplex[i][0], simplex[i][1], simplex[i][2], simplex[i][3], SigmaHat1, 
						      	   SigmaHat2, C, h1, h2);
		/* ################# uncomment on this to test on parabaloid ################
		NegLogLambda[i] =    (simplex[i][0]-0.5)*(simplex[i][0]-0.5) 
                                 + (simplex[i][1]-0.5)*(simplex[i][1]-0.5) 
                                 + (simplex[i][2]-0.5)*(simplex[i][2]-0.5) 
                                 + (simplex[i][3]-0.5)*(simplex[i][3]-0.5);
		*/
	}
	call = 5;

	/* sort the simplex */
	SimplexSort(simplex,NegLogLambda);
	iterate = 1;

	/* main algorithm (self contained) */
	while(1) {

		/* check stop conditions */
		if(call > MaxCall || iterate > MaxIterate) return 1;
		x = f = 0.0;
		for(i=1; i<5; i++) for(j=0; j<4; j++) {
			if(fabs(simplex[0][j] - simplex[i][j]) >= x) x = fabs(simplex[0][j] - simplex[i][j]);
		}
                for(i=1; i<5; i++) {
                        if(  fabs(NegLogLambda[0] - NegLogLambda[i])  >= f) f = fabs(NegLogLambda[0] - NegLogLambda[i]);
                }
		if( x < TolX && f < TolF) return 0;

		/* compute reflection point */
		for(i=0; i<4; i++) AveragePoint[i] = 0.0;
		for(i=0; i<4; i++) for(j=0; j<4; j++) AveragePoint[i] += simplex[j][i] / 4.0; 
		for(i=0; i<4; i++) ReflectionPoint[i] = AveragePoint[i] + rho*(AveragePoint[i] - simplex[4][i]);

		/* evaluate likelihood at reflection point */
		BoundaryCheck(ReflectionPoint); 
		ReflectionValue = - PairExactLogLikelihood(N, ReflectionPoint[0], ReflectionPoint[1], ReflectionPoint[2], ReflectionPoint[3], 
							 SigmaHat1, SigmaHat2, C, h1, h2);
		/* ################# uncomment on this to test on parabaloid ################
		ReflectionValue = (ReflectionPoint[0]-0.5)*(ReflectionPoint[0]-0.5)
				 + (ReflectionPoint[1]-0.5)*(ReflectionPoint[1]-0.5) 
				 + (ReflectionPoint[2]-0.5)*(ReflectionPoint[2]-0.5) 
				 + (ReflectionPoint[3]-0.5)*(ReflectionPoint[3]-0.5);
		*/
		call++;

		/* if reflection is better than best, check expanded reflection, otherwise compare to worst point and reflect or contract */
		if(ReflectionValue < NegLogLambda[0]) { 
	
			/* calculate the expansion point */
			for(i=0; i<4; i++) ExpansionPoint[i] = AveragePoint[i] + rho*chi*(AveragePoint[i] - simplex[4][i]);

			/* evaluate likelihood at expansion point */
			BoundaryCheck(ExpansionPoint); 
			ExpansionValue = - PairExactLogLikelihood(N, ExpansionPoint[0], ExpansionPoint[1], ExpansionPoint[2], 
								ExpansionPoint[3],SigmaHat1, SigmaHat2, C, h1, h2);
			/* ################# uncomment on this to test on parabaloid ################
			ExpansionValue=  (ExpansionPoint[0]-0.5)*(ExpansionPoint[0]-0.5)
                                 	+ (ExpansionPoint[1]-0.5)*(ExpansionPoint[1]-0.5)
                                 	+ (ExpansionPoint[2]-0.5)*(ExpansionPoint[2]-0.5)
                                	+ (ExpansionPoint[3]-0.5)*(ExpansionPoint[3]-0.5);
			*/
			call++;

			/* keep whichever point is better (expansion or reflection) */
			if(ExpansionValue < ReflectionValue) {
				for(i=0; i<4; i++) simplex[4][i] = ExpansionPoint[i];
				NegLogLambda[4] = ExpansionValue;
			} else {
				for(i=0; i<4; i++) simplex[4][i] = ReflectionPoint[i];
				NegLogLambda[4] = ReflectionValue;
			}

		} else {
			
			/* reflection point was not better than best.  compare to second worst point and then reflect or contract*/
			if(ReflectionValue < NegLogLambda[3]) {
                                for(i=0; i<4; i++) simplex[4][i] = ReflectionPoint[i];
                                NegLogLambda[4] = ReflectionValue;
			} else {
						
				/* reflection is worse than second worst point, if better(worse) than worst contract or shrink inside(outside) */
				if(ReflectionValue < NegLogLambda[4]) {

					/* calculate outside contraction point */
					for(i=0; i<4; i++) ContractionPoint[i] = AveragePoint[i] + rho*psi*(AveragePoint[i] - simplex[4][i]);
	
					/* evaluate likeilhood at contraction point */
					BoundaryCheck(ContractionPoint);
					ContractionValue = - PairExactLogLikelihood(N, ContractionPoint[0], ContractionPoint[1], ContractionPoint[2],
										  ContractionPoint[3],SigmaHat1, SigmaHat2, C, h1, h2);
					/* ################# uncomment on this to test on parabaloid ################
					ContractionValue =  (ContractionPoint[0]-0.5)*(ContractionPoint[0]-0.5)
                                        		   + (ContractionPoint[1]-0.5)*(ContractionPoint[1]-0.5)
                                        		   + (ContractionPoint[2]-0.5)*(ContractionPoint[2]-0.5)
                                        		   + (ContractionPoint[3]-0.5)*(ContractionPoint[3]-0.5);
					*/
					call++;

					/* contract or shrink */
					if(ContractionValue <= ReflectionValue) {
						for(i=0; i<4; i++) simplex[4][i] = ContractionPoint[i];
						NegLogLambda[4] = ContractionValue;
					} else shrink = 1;
				} else {

					/* calculate inside contraction point */
                                        for(i=0; i<4; i++) ContractionPoint[i] = AveragePoint[i] - psi*(AveragePoint[i] - simplex[4][i]);

                                        /* evaluate likeilhood at contraction point */
					BoundaryCheck(ContractionPoint);
                                        ContractionValue = - PairExactLogLikelihood(N, ContractionPoint[0], ContractionPoint[1], ContractionPoint[2],
                                                                                    ContractionPoint[3],SigmaHat1, SigmaHat2, C, h1, h2);
					/* ################# uncomment on this to test on parabaloid ################
                                        ContractionValue =  (ContractionPoint[0]-0.5)*(ContractionPoint[0]-0.5)
                                                           + (ContractionPoint[1]-0.5)*(ContractionPoint[1]-0.5)
                                                           + (ContractionPoint[2]-0.5)*(ContractionPoint[2]-0.5)
                                                           + (ContractionPoint[3]-0.5)*(ContractionPoint[3]-0.5);
					*/
                                        call++;

                                        /* contract or shrink */
                                        if(ContractionValue < NegLogLambda[4]) {
                                                for(i=0; i<4; i++) simplex[4][i] = ContractionPoint[i];
                                                NegLogLambda[4] = ContractionValue;
                                        } else shrink = 1;
			
				}

			}

			/* shrink if necessary */
			if(shrink) {

				for(i=1; i<5; i++) {
					for(j=0; j<4; j++) simplex[i][j] = simplex[0][j] + gamma*(simplex[i][j] - simplex[0][j]);
					BoundaryCheck(simplex[i]); 
					NegLogLambda[i] = - PairExactLogLikelihood(N, simplex[i][0], simplex[i][1], simplex[i][2], simplex[i][3], SigmaHat1,
                                                      			      	   SigmaHat2, C, h1, h2);
					/* ################# uncomment on this to test on parabaloid ################
					NegLogLambda[i] =  (simplex[i][0]-0.5)*(simplex[i][0]-0.5)
                                                       + (simplex[i][1]-0.5)*(simplex[i][1]-0.5)
                                                       + (simplex[i][2]-0.5)*(simplex[i][2]-0.5)
                                                       + (simplex[i][3]-0.5)*(simplex[i][3]-0.5);
					*/
				}
				call += 4;
				shrink = 0;
			}

			
		}

		/* sort the simplex */
		SimplexSort(simplex,NegLogLambda);

		/* increment iteration counter */
		iterate++;

		/* update best value and point */
		*MaxLogLambda = - NegLogLambda[0];
		*SigmaBar1 = simplex[0][0];
		*SigmaBar2 = simplex[0][1];
		*AlphaBar = simplex[0][2];
		*XiBar = simplex[0][3];

	} /* end of self contained main algorithm */

}
