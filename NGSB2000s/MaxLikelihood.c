/*
 * File: MaxLikelihood.c  by: Steve Drasco 
 *                                                           
 * This program estimates the maximum of the likelihood function 
 * for detection stochastic background signals.
 *
 */
 

#include "NGSB.h"

int MaxLikelihood(double *LambdaMax, double *XiBar, double *AlphaBar, int N, double SigmaHat1, double SigmaHat2, double C, double *h1, 
		  double *h2, double *xi, double *alpha, int GN, double f)
{
	double	XiMin, XiMax, AlphaMin, AlphaMax, DeltaXi, DeltaAlpha, Lambda, GridN;
        int     i, j;

        /* check data pointer */
        if(h1 == NULL || h2 == NULL || xi == NULL || alpha == NULL) {
                printf("MaxLikelihood reports: An input pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* make double version of Grid N  and random boundary for course search */
	GridN = (double) GN;
	XiMin = - (Krand() + 10.0);
	AlphaMin = - (Krand() + 10.0);

	/* generate coarse grid (GNxGN) of xi and alpha values */
	/* NOTE: coarse grid is a fixed window */
	for(j = 0; j < GN; j++) {
		xi[j] = pow(10.0,XiMin - (((double) j) * XiMin / (GridN-1.0) ) );
		alpha[j] = pow(10.0,AlphaMin - (((double) j) * (AlphaMin-1.0) / (GridN-1.0) ) );
	}

	/* initialize LambdaMax */
	*LambdaMax = -1e100;

	/* maximize over grid */
	for(i = 0; i < GN; i++) {
		for(j = 0; j < GN; j++) {
			Lambda = PairExactLogLikelihood(N, SigmaHat1 - xi[i]*alpha[j], SigmaHat2 - xi[i]*alpha[j], 
							alpha[j], xi[i], SigmaHat1, SigmaHat2, C, h1, h2);
			if( Lambda > *LambdaMax) {
				*AlphaBar = alpha[j];
				*XiBar = xi[i];
				*LambdaMax = Lambda;
			}
		}
	}

        /* compute boundaries of dynamic fine grid window */
        XiMax = log10(*XiBar)*(1.1);
        XiMin = log10(*XiBar)*(0.9);
        DeltaXi = XiMax - XiMin;
        AlphaMin = log10(*AlphaBar)*(1.1);
        AlphaMax = log10(*AlphaBar)*(0.9);
        DeltaAlpha = AlphaMax - AlphaMin;

        /* generate fine grid of xi and alpha values */
        for(j = 0; j < GN; j++) {
                xi[j] = pow(   10.0,  XiMin + (((double) j) * DeltaXi / (GridN-1.0) )   );
                alpha[j] = pow(   10.0,  AlphaMin + (((double) j) * DeltaAlpha / (GridN-1.0) )   );
        }


        /* maximize over grid */
        for(i = 0; i < GN; i++) {
                for(j = 0; j < GN; j++) {
                        Lambda = PairExactLogLikelihood(N, SigmaHat1 - xi[i]*alpha[j], SigmaHat2 - xi[i]*alpha[j], 
                                                        alpha[j], xi[i], SigmaHat1, SigmaHat2, C, h1, h2);
                        if( Lambda > *LambdaMax) {
                                *AlphaBar = alpha[j];
                                *XiBar = xi[i];
                                *LambdaMax = Lambda;
                        }
                }
        }



	/* compute boundaries of dynamic ultra fine grid window */
        XiMax = *XiBar*(1.0+f);
        XiMin = *XiBar*(1.0-f);
	DeltaXi = XiMax - XiMin;
        AlphaMin = *AlphaBar*(1.0+f);
        AlphaMax = *AlphaBar*(1.0-f);
	DeltaAlpha = AlphaMax - AlphaMin;
       	
	/* generate ultra fine grid of xi and alpha values */
        for(j = 0; j < GN; j++) {
                xi[j] = XiMin + (((double) j) * DeltaXi / (GridN-1.0) );
                alpha[j] = AlphaMin + (((double) j) * DeltaAlpha / (GridN-1.0) );
        }

        /* maximize over grid */
        for(i = 0; i < GN; i++) {
                for(j = 0; j < GN; j++) {
                        Lambda = PairExactLogLikelihood(N, SigmaHat1 - xi[i]*alpha[j], SigmaHat2 - xi[i]*alpha[j], 
                                                        alpha[j], xi[i], SigmaHat1, SigmaHat2, C, h1, h2);
                        if( Lambda > *LambdaMax) {
                                *AlphaBar = alpha[j];
                                *XiBar = xi[i];
                                *LambdaMax = Lambda;
                        }
                }
        }

	/* normal exit */	
	return 0;
}
