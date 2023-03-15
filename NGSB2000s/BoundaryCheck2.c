/*
 * File: BoundaryCheck2.c  by: Steve Drasco 
 *
 * Use this to make sure simplex maximum likelihood search
 * stays within it's boundries when searching so as not generate 
 * NaN's etc.
 *
 * NOTE: 
 * 	this one is for SimplexMaxLikelihood3
 *
 */
 

#include "NGSB.h"

int BoundaryCheck2(double *v)
{
	if(v[0] < 0.0) v[0] = 0.00025*exp(v[0]);
        if(v[1] < 0.0) v[1] = 0.00025*exp(v[1]);
        if(v[2] < 0.0) v[2] = 0.00025*exp(v[2]);
        if(v[2] > 1.0) v[2] = 1.0 - 0.00025*exp(1.0-v[2]);
	return 0;
}
