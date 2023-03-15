/*
 * File: GetAlpha.c  by: Steve Drasco 
 *
 * Given 4 coefficients c4-c3 defining a 4th order polynomial, this finds the extremum nearest zero.
 *
 *  c4 alpha^4 + c3 alpha^3 + c2 alpha^2 + c1 alpha + c0  = (4th order polynomial)
 *
 * of course we don't need c0.
 *
 */
 

#include "NGSB.h"

double GetAlpha(double c4, double c3, double c2, double c1)
{
	double	a1, a2, a3, Q, R, D, S, T, alpha, theta;
	double	x1, x2, x3, small, medium, large;

	/* compute "schaum-racine" coefficients */
	a1 = 0.75 * c3 / c4;
	a2 = 0.5 * c2 / c4;
	a3 = c1 / (4*c4);

	/* compute discriminant */
	Q = (3.0*a2 - a1*a1) / 9.0;
	R = (9.0*a1*a2 - 27.0*a3 - 2.0*a1*a1*a1) / 54.0;
	D = Q*Q*Q + R*R;

	/* could be one two or three unique real extrema */
	if(D > 0.0) {
		
		/* one real extrema */
		S = ((R+sqrt(D))/fabs(R+sqrt(D)))  * pow(fabs(R+sqrt(D)),1.0/3.0);  /* being careful of nans */
		T = ((R+sqrt(D))/fabs(R+sqrt(D)))  * pow(fabs(R-sqrt(D)),1.0/3.0);  /* being careful of nans */
		alpha = S + T - a1/3.0;
		if(alpha < 0.0) alpha = 0.0;
		
	} else if(D == 0.0) {

		/* three real extrema, at least two are equal */
		S = (R / fabs(R)) * pow(fabs(R),1.0/3.0); /* being careful of nans */
                T = (R / fabs(R)) * pow(fabs(R),1.0/3.0); /* being careful of nans */
		x1 = S + T - a1/3.0;
		x2 = -0.5*(S + T) - a1/3.0;

		/* sort them */
		if(x1 <= x2) {
			small = x1;
			large = x2;
		} else {
			small = x2;
			large = x1;
		}

                /* choose smallest positive one */
                if(small >= 0.0) {
                        alpha = small;
                } else if(large >= 0.0) {
                        alpha = large;
                } else {
                        alpha = 0.0;
                }

	} else {
		/* three unique real extrema */
		theta = atan(R / sqrt(-D)) / 3.0;
		x1 = - cos(theta) - a1/3.0;
		x2 = (D - R*R)*(cos(theta) + sqrt(3.0)*sin(theta)) - a1/3.0;
		x3 = (D - R*R)*(cos(theta) - sqrt(3.0)*sin(theta)) - a1/3.0;

		/* sort them */
		if(x1 <=  x2 && x1 <= x3){
			small = x1;
			if(x2 < x3) {
				medium = x2;
				large = x3;
			} else {
				medium = x3;
				large = x2;
			}
		} else if(x2 <=  x1 && x2 <= x3) {
                        small = x2;
                        if(x1 < x3) {
                                medium = x1;
                                large = x3;
                        } else {
                                medium = x3;
                                large = x1;
                        }
		} else if(x3 <=  x2 && x3 <= x1) {
                        small = x3;
                        if(x2 < x1) {
                                medium = x2;
                                large = x1;
                        } else {
                                medium = x1;
                                large = x2;
                        }
		}

		/* choose smallest positive one */
		if(small >= 0.0) {
			alpha = small;
		} else if(medium >= 0.0) {
			alpha = medium;
		} else if(large >= 0.0) {
			alpha = large;
		} else {
			alpha = 0.0;
		}

	} 

	/* normal exit */	
	return alpha;
}
