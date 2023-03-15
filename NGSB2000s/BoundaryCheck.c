/*
 * File: BoundaryCheck.c  by: Steve Drasco 
 *
 * Use this to make sure simplex maximum likelihood search
 * stays within it's boundries when searching so as not generate 
 * NaN's etc.
 *
 */
 

#include "NGSB.h"

int BoundaryCheck(double *v)
{
	if(v[0] < 0.0) v[0] = 0.00025*exp(v[0]);
        if(v[1] < 0.0) v[1] = 0.00025*exp(v[1]);
        if(v[2] < 0.0) v[2] = 0.00025*exp(v[2]);
        if(v[3] < 0.0) v[3] = 0.00025*exp(v[3]);
        if(v[3] > 1.0) v[3] = 1.0 - 0.00025*exp(1.0-v[3]);
	return 0;
}
