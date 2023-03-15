/*
 * File: SimplexMaxLikelihood3.c  by: Steve Drasco 
 *                                                           
 *
 * see SimplexMaxLikelihood.c ... the version here uses the approximate statistic - where AlphaBar is found by other means.
 *
 * NOTE:
 *      You *MUST* call MaxLogLambda = PairApproxLogLikelihood2(N, SigmaBar1, SigmaBar2, &AlphaBar, XiBar, sum) after 
 *      using this to get the appropriate estimate of AlphaBar! 
 *
 */
 

#include "NGSB.h"

int SimplexMaxLikelihood3(double *SigmaBar1, double *SigmaBar2, double *XiBar, 
			int N, double **simplex, double **sum, int MaxIterate,
			int MaxCall, double TolX,double TolF)
{
	double	*NegLogLambda, *AveragePoint, *ReflectionPoint, *ExpansionPoint, *ContractionPoint;
	double	UsualDelta, ZeroTermDelta, rho, chi, psi, gamma, x, f, DummyAlphaBar;
	double	AverageValue, ReflectionValue, ExpansionValue, ContractionValue;
        int     i, j, call, shrink, iterate;

	/* initialize parameters */
	rho = 1.0;
	chi = 2.0;
	psi = 0.5;
	gamma = 0.5;
	shrink = 0;
	UsualDelta = 5e-2;
	ZeroTermDelta = 2.5e-4;

	/* allocate memory */
	NegLogLambda = (double *)malloc(4*sizeof(double));
	AveragePoint = (double *)malloc(3*sizeof(double));
	ReflectionPoint = (double *)malloc(3*sizeof(double));
	ExpansionPoint = (double *)malloc(3*sizeof(double));
	ContractionPoint = (double *)malloc(3*sizeof(double));

	/* give birth to the simplex */
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++){ 
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
	for(i=0; i<4; i++) {  
		BoundaryCheck2(simplex[i]);
		NegLogLambda[i] = - PairApproxLogLikelihood2(N, simplex[i][0], simplex[i][1], &DummyAlphaBar, simplex[i][2], sum); 
		/* ################# uncomment on this to test on parabaloid ################ 
		NegLogLambda[i] =    (simplex[i][0]-0.5)*(simplex[i][0]-0.5) 
                                 + (simplex[i][1]-0.5)*(simplex[i][1]-0.5) 
                                 + (simplex[i][2]-0.5)*(simplex[i][2]-0.5); */
	}
	call = 4; 

	/* sort the simplex */
	SimplexSort2(simplex,NegLogLambda);
	iterate = 1;

	/* main algorithm (self contained) */
	while(1) {

		/* check stop conditions */
		if(call > MaxCall || iterate > MaxIterate) return 1;
		x = f = 0.0;
		for(i=1; i<4; i++) for(j=0; j<3; j++) {
			if(fabs(simplex[0][j] - simplex[i][j]) >= x) x = fabs(simplex[0][j] - simplex[i][j]);
		}
                for(i=1; i<4; i++) {
                        if(  fabs(NegLogLambda[0] - NegLogLambda[i])  >= f) f = fabs(NegLogLambda[0] - NegLogLambda[i]);
                }
		if( x < TolX && f < TolF) return 0;

		/* compute reflection point */
		for(i=0; i<3; i++) AveragePoint[i] = 0.0;
		for(i=0; i<3; i++) for(j=0; j<3; j++) AveragePoint[i] += simplex[j][i] / 3.0; 
		for(i=0; i<3; i++) ReflectionPoint[i] = AveragePoint[i] + rho*(AveragePoint[i] - simplex[3][i]);

		/* evaluate likelihood at reflection point */
		BoundaryCheck2(ReflectionPoint); 
		ReflectionValue = - PairApproxLogLikelihood2(N, ReflectionPoint[0], ReflectionPoint[1], &DummyAlphaBar, ReflectionPoint[2], sum);
		/* ################# uncomment on this to test on parabaloid ################ 
		ReflectionValue = (ReflectionPoint[0]-0.5)*(ReflectionPoint[0]-0.5)
				 + (ReflectionPoint[1]-0.5)*(ReflectionPoint[1]-0.5) 
				 + (ReflectionPoint[2]-0.5)*(ReflectionPoint[2]-0.5) ; */
		call++;

		/* if reflection is better than best, check expanded reflection, otherwise compare to worst point and reflect or contract */
		if(ReflectionValue < NegLogLambda[0]) { 
	
			/* calculate the expansion point */
			for(i=0; i<3; i++) ExpansionPoint[i] = AveragePoint[i] + rho*chi*(AveragePoint[i] - simplex[3][i]);

			/* evaluate likelihood at expansion point */
			BoundaryCheck2(ExpansionPoint); 
			ExpansionValue = - PairApproxLogLikelihood2(N, ExpansionPoint[0], ExpansionPoint[1], &DummyAlphaBar, ExpansionPoint[2], sum);
			/* ################# uncomment on this to test on parabaloid ################ 
			ExpansionValue=  (ExpansionPoint[0]-0.5)*(ExpansionPoint[0]-0.5)
                                 	+ (ExpansionPoint[1]-0.5)*(ExpansionPoint[1]-0.5)
                                 	+ (ExpansionPoint[2]-0.5)*(ExpansionPoint[2]-0.5); */
			call++;

			/* keep whichever point is better (expansion or reflection) */
			if(ExpansionValue < ReflectionValue) {
				for(i=0; i<3; i++) simplex[3][i] = ExpansionPoint[i];
				NegLogLambda[3] = ExpansionValue;
			} else {
				for(i=0; i<3; i++) simplex[3][i] = ReflectionPoint[i];
				NegLogLambda[3] = ReflectionValue;
			}

		} else {
			
			/* reflection point was not better than best.  compare to second worst point and then reflect or contract*/
			if(ReflectionValue < NegLogLambda[2]) {
                                for(i=0; i<3; i++) simplex[3][i] = ReflectionPoint[i];
                                NegLogLambda[3] = ReflectionValue;
			} else {
						
				/* reflection is worse than second worst point, if better(worse) than worst contract or shrink inside(outside) */
				if(ReflectionValue < NegLogLambda[3]) {

					/* calculate outside contraction point */
					for(i=0; i<3; i++) ContractionPoint[i] = AveragePoint[i] + rho*psi*(AveragePoint[i] - simplex[3][i]);
	
					/* evaluate likeilhood at contraction point */
					BoundaryCheck2(ContractionPoint);
					ContractionValue = - PairApproxLogLikelihood2(N, ContractionPoint[0], ContractionPoint[1], &DummyAlphaBar,
										  ContractionPoint[2], sum);
					/* ################# uncomment on this to test on parabaloid ################
					ContractionValue =  (ContractionPoint[0]-0.5)*(ContractionPoint[0]-0.5)
                                        		   + (ContractionPoint[1]-0.5)*(ContractionPoint[1]-0.5)
                                        		   + (ContractionPoint[2]-0.5)*(ContractionPoint[2]-0.5); */
					call++;

					/* contract or shrink */
					if(ContractionValue <= ReflectionValue) {
						for(i=0; i<3; i++) simplex[3][i] = ContractionPoint[i];
						NegLogLambda[3] = ContractionValue;
					} else shrink = 1;
				} else {

					/* calculate inside contraction point */
                                        for(i=0; i<3; i++) ContractionPoint[i] = AveragePoint[i] - psi*(AveragePoint[i] - simplex[3][i]);

                                        /* evaluate likeilhood at contraction point */
					BoundaryCheck2(ContractionPoint);
                                        ContractionValue = - PairApproxLogLikelihood2(N, ContractionPoint[0], ContractionPoint[1], &DummyAlphaBar,
                                                                                    ContractionPoint[2], sum);
					/* ################# uncomment on this to test on parabaloid ################ 
                                        ContractionValue =  (ContractionPoint[0]-0.5)*(ContractionPoint[0]-0.5)
                                                           + (ContractionPoint[1]-0.5)*(ContractionPoint[1]-0.5)
                                                           + (ContractionPoint[2]-0.5)*(ContractionPoint[2]-0.5); */
                                        call++;

                                        /* contract or shrink */
                                        if(ContractionValue < NegLogLambda[3]) {
                                                for(i=0; i<3; i++) simplex[3][i] = ContractionPoint[i];
                                                NegLogLambda[3] = ContractionValue;
                                        } else shrink = 1;
			
				}

			}

			/* shrink if necessary */
			if(shrink) {

				for(i=1; i<4; i++) {
					for(j=0; j<3; j++) simplex[i][j] = simplex[0][j] + gamma*(simplex[i][j] - simplex[0][j]);
					BoundaryCheck2(simplex[i]); 
					NegLogLambda[i] = - PairApproxLogLikelihood2(N, simplex[i][0], simplex[i][1], &DummyAlphaBar, simplex[i][2], sum);
					/* ################# uncomment on this to test on parabaloid ################ 
					NegLogLambda[i] =  (simplex[i][0]-0.5)*(simplex[i][0]-0.5)
                                                       + (simplex[i][1]-0.5)*(simplex[i][1]-0.5)
                                                       + (simplex[i][2]-0.5)*(simplex[i][2]-0.5); */
				}
				call += 3;
				shrink = 0;
			}

			
		}

		/* sort the simplex */
		SimplexSort2(simplex,NegLogLambda);

		/* increment iteration counter */
		iterate++;

		/* update best value and point */
		*SigmaBar1 = simplex[0][0];
		*SigmaBar2 = simplex[0][1];
		*XiBar = simplex[0][2];

	} /* end of self contained main algorithm */

}
