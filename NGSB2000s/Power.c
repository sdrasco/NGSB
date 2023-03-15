/*
 * File: Power.c  by: Steve Drasco 
 *
 * Simple program.  Input x (double) and a (INTEGER!).  Output x^a.
 *                                                           
 */

#include "NGSB.h"

double Power(double x, int a)
{
	double	result=1.0;
        int     i;

	/* main loop */
	for(i=0; i<a; i++) result *= x;

	/* exit */
	return result;
}
