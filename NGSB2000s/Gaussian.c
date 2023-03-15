/*-----------------------------------------------------------*/
/* File: Gaussian.c     by: Steve Drasco                     */
/*                                                           */
/* This program generates data with a Gaussian distribution  */
/* with zero mean and users standard deviation and length.   */
/*-----------------------------------------------------------*/

#include "NGSB.h"

int Gaussian(double sigma, int N, double Pmin, double *out)
{
        double   	x, y, xmax, ymax;
        double  	pi = 3.14159265358979323846;
        int     	i = 0;
	extern long     idum;

        /* check output pointer */
        if(out == NULL) {
                printf("Output pointer was NULL.  Giving up ...\n");
                return 1;
        }

        /* compute boundries */
        xmax = sqrt(     fabs(   2.0 * sigma * log( Pmin * sqrt(2*pi*sigma) )   )       );
        ymax = 1.0 / sqrt(2.0*pi*sigma) ;

        /* main loop */
        while(i < N){

                /* make a new point */
                x = ( 2.0 * xmax * ((double) ran2(&idum)) ) - xmax;
                y = ymax * ((double) ran2(&idum)); 

                /* if point is below curve - keep x */
                if (   y <= exp( - x * x / (2.0 * sigma) )   /   sqrt(2.0 * pi * sigma)     ) {
                        out[i]=x;
                        i++;
                }
        }

	/* normal exit */
	return 0;
}
