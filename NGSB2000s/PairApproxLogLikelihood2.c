/*
 * File: PairApproxLogLikelihood2.c  by: Steve Drasco 
 *                                                           
 * This program computes the log of the two-detector likelihood  
 * function of user specified sigma1, sigma2, and xi to fourth order in alpha 
 * for a non-Gaussian stochastic background on Gaussian Noise.  
 * 
 * We taylor expand to fourth order in alpha and use the extremum at the 
 * smallest positive value of alpha, or we set alpha = 0.
 *      
 * NOTE:
 *      sums is a 8x8 matrix with components:  sums[a][b] = SumOverj( h1[j]^a h2[j]^b ) 
 *	(Need only components with a+b = 2,4,6,8  that is need only 24 sums)
 *
 * NOTE:
 *	The alpha we use is returned as AlphaBar.
 */

#include "NGSB.h"


double PairApproxLogLikelihood2(int N, double Sigma1, double Sigma2, double *AlphaBar, double xi, double **sum)
{
        double  stat, c0, c1, c2, c3, c4, SigmaHat1, SigmaHat2;

	/* define SigmaHat1,2 */
	SigmaHat1 = sum[2][0] / ((double) N);
	SigmaHat2 = sum[0][2] / ((double) N);

	/* compute coefficients  (from CForm in Mathematica)*/
c0 = N - sum[0][2]/(2.*Sigma2) - sum[2][0]/(2.*Sigma1) + 
   (N*log((SigmaHat1*SigmaHat2)/(Sigma1*Sigma2)))/2.; 

c1 = -(N*xi)/(2.*Sigma1) - (N*xi)/(2.*Sigma2) + 
   (sum[0][2]*xi)/(2.*Power(Sigma2,2)) + 
   (sum[1][1]*xi)/(Sigma1*Sigma2) + 
   (sum[2][0]*xi)/(2.*Power(Sigma1,2));

c2 = -(N*(-3 + xi)*xi)/(8.*Power(Sigma1,2)) - 
   (N*(-3 + xi)*xi)/(8.*Power(Sigma2,2)) - 
   (N*(-3 + xi)*xi)/(4.*Sigma1*Sigma2) + 
   (sum[0][2]*(-3 + xi)*xi)/(4.*Power(Sigma2,3)) + 
   (sum[0][2]*(-3 + xi)*xi)/(4.*Sigma1*Power(Sigma2,2)) + 
   (sum[1][1]*(-3 + xi)*xi)/(2.*Sigma1*Power(Sigma2,2)) + 
   (sum[1][1]*(-3 + xi)*xi)/(2.*Power(Sigma1,2)*Sigma2) + 
   (sum[2][0]*(-3 + xi)*xi)/(4.*Power(Sigma1,3)) + 
   (sum[2][0]*(-3 + xi)*xi)/(4.*Power(Sigma1,2)*Sigma2) - 
   (sum[0][4]*(-1 + xi)*xi)/(8.*Power(Sigma2,4)) - 
   (sum[1][3]*(-1 + xi)*xi)/(2.*Sigma1*Power(Sigma2,3)) - 
   (3*sum[2][2]*(-1 + xi)*xi)/
    (4.*Power(Sigma1,2)*Power(Sigma2,2)) - 
   (sum[3][1]*(-1 + xi)*xi)/(2.*Power(Sigma1,3)*Sigma2) - 
   (sum[4][0]*(-1 + xi)*xi)/(8.*Power(Sigma1,4));

c3 = (15*sum[0][2]*xi)/(16.*Power(Sigma1,2)*Power(Sigma2,2)) + 
   (15*sum[1][1]*xi)/(8.*Power(Sigma1,3)*Sigma2) + 
   (15*sum[2][0]*xi)/(16.*Power(Sigma1,4)) - 
   (15*sum[2][2]*xi)/(8.*Power(Sigma1,2)*Power(Sigma2,3)) - 
   (5*sum[3][1]*xi)/(4.*Power(Sigma1,3)*Power(Sigma2,2)) - 
   (5*sum[4][0]*xi)/(16.*Power(Sigma1,4)*Sigma2) - 
   (5*sum[2][4]*(1 - 2*xi)*(-1 + xi)*xi)/
    (16.*Power(Sigma1,2)*Power(Sigma2,4)) - 
   (5*sum[3][3]*(1 - 2*xi)*(-1 + xi)*xi)/
    (11.99999999999999*Power(Sigma1,3)*Power(Sigma2,3))\
    - (9*sum[0][2]*Power(xi,2))/
    (16.*Power(Sigma1,2)*Power(Sigma2,2)) - 
   (9*sum[1][1]*Power(xi,2))/(8.*Power(Sigma1,3)*Sigma2) - 
   (9*sum[2][0]*Power(xi,2))/(16.*Power(Sigma1,4)) + 
   (21*sum[2][2]*Power(xi,2))/
    (8.*Power(Sigma1,2)*Power(Sigma2,3)) + 
   (7*sum[3][1]*Power(xi,2))/
    (4.*Power(Sigma1,3)*Power(Sigma2,2)) + 
   (7*sum[4][0]*Power(xi,2))/(16.*Power(Sigma1,4)*Sigma2) + 
   (sum[0][2]*Power(xi,3))/
    (8.*Power(Sigma1,2)*Power(Sigma2,2)) + 
   (sum[1][1]*Power(xi,3))/(4.*Power(Sigma1,3)*Sigma2) + 
   (sum[2][0]*Power(xi,3))/(8.*Power(Sigma1,4)) - 
   (3*sum[2][2]*Power(xi,3))/
    (4.*Power(Sigma1,2)*Power(Sigma2,3)) - 
   (sum[3][1]*Power(xi,3))/
    (2.*Power(Sigma1,3)*Power(Sigma2,2)) - 
   (sum[4][0]*Power(xi,3))/(8.*Power(Sigma1,4)*Sigma2) - 
   (sum[0][4]*(-1 + xi)*xi*(-5 + 2*xi))/
    (16.*Sigma1*Power(Sigma2,4)) - 
   (sum[1][3]*(-1 + xi)*xi*(-5 + 2*xi))/
    (4.*Sigma1*Power(Sigma2,4)) - 
   (sum[1][3]*(-1 + xi)*xi*(-5 + 2*xi))/
    (4.*Power(Sigma1,2)*Power(Sigma2,3)) - 
   (sum[4][0]*(-1 + xi)*xi*(-5 + 2*xi))/
    (16.*Power(Sigma2,5)) + 
   (N*xi*(-15 + 9*xi - 2*Power(xi,2)))/
    (47.99999999999999*Power(Sigma2,3)) - 
   (N*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (47.99999999999999*Power(Sigma1,3)) - 
   (N*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (16.*Sigma1*Power(Sigma2,2)) - 
   (N*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (16.*Power(Sigma1,2)*Sigma2) + 
   (sum[0][2]*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (16.*Power(Sigma2,4)) + 
   (sum[0][2]*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (8.*Sigma1*Power(Sigma2,3)) + 
   (sum[1][1]*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (8.*Sigma1*Power(Sigma2,3)) + 
   (sum[1][1]*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (4.*Power(Sigma1,2)*Power(Sigma2,2)) + 
   (sum[2][0]*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (16.*Power(Sigma1,2)*Power(Sigma2,2)) + 
   (sum[2][0]*xi*(15 - 9*xi + 2*Power(xi,2)))/
    (8.*Power(Sigma1,3)*Sigma2) - 
   (3*sum[2][2]*xi*(5 - 7*xi + 2*Power(xi,2)))/
    (8.*Power(Sigma1,3)*Power(Sigma2,2)) - 
   (sum[3][1]*xi*(5 - 7*xi + 2*Power(xi,2)))/
    (4.*Power(Sigma1,4)*Sigma2) - 
   (sum[4][0]*xi*(5 - 7*xi + 2*Power(xi,2)))/
    (16.*Power(Sigma1,5)) + 
   (sum[1][5]*xi*(1 - 3*xi + 2*Power(xi,2)))/
    (8.*Sigma1*Power(Sigma2,5)) + 
   (5*sum[4][2]*xi*(1 - 3*xi + 2*Power(xi,2)))/
    (16.*Power(Sigma1,4)*Power(Sigma2,2)) + 
   (sum[5][1]*xi*(1 - 3*xi + 2*Power(xi,2)))/
    (8.*Power(Sigma1,5)*Sigma2) + 
   (sum[6][0]*xi*(1 - 3*xi + 2*Power(xi,2)))/
    (47.99999999999999*Power(Sigma1,6)) + 
   (sum[6][0]*xi*(1 - 3*xi + 2*Power(xi,2)))/
    (47.99999999999999*Power(Sigma2,6));

c4 = (35*N*xi)/(128.*Power(Sigma1,4)) - 
   (105*sum[0][2]*xi)/
    (32.*Power(Sigma1,2)*Power(Sigma2,3)) - 
   (35*sum[0][2]*xi)/(32.*Power(Sigma1,3)*Power(Sigma2,2)) - 
   (105*sum[1][1]*xi)/
    (16.*Power(Sigma1,3)*Power(Sigma2,2)) - 
   (35*sum[1][1]*xi)/(16.*Power(Sigma1,4)*Sigma2) - 
   (35*sum[2][0]*xi)/(32.*Power(Sigma1,5)) - 
   (105*sum[2][0]*xi)/(32.*Power(Sigma1,4)*Sigma2) + 
   (105*sum[2][2]*xi)/
    (32.*Power(Sigma1,2)*Power(Sigma2,4)) + 
   (105*sum[2][2]*xi)/
    (16.*Power(Sigma1,3)*Power(Sigma2,3)) + 
   (35*sum[3][1]*xi)/(16.*Power(Sigma1,3)*Power(Sigma2,3)) + 
   (35*sum[3][1]*xi)/(8.*Power(Sigma1,4)*Power(Sigma2,2)) + 
   (35*sum[4][0]*xi)/(64.*Power(Sigma1,4)*Power(Sigma2,2)) + 
   (35*sum[4][0]*xi)/(32.*Power(Sigma1,5)*Sigma2) - 
   (35*sum[0][4]*(-1 + xi)*xi)/
    (64.*Power(Sigma1,2)*Power(Sigma2,4)) - 
   (35*sum[1][3]*(-1 + xi)*xi)/
    (16.*Power(Sigma1,3)*Power(Sigma2,3)) - 
   (105*sum[2][2]*(-1 + xi)*xi)/
    (32.*Power(Sigma1,4)*Power(Sigma2,2)) + 
   (35*sum[2][4]*(-1 + xi)*xi)/
    (32.*Power(Sigma1,2)*Power(Sigma2,5)) - 
   (35*sum[3][1]*(-1 + xi)*xi)/
    (16.*Power(Sigma1,5)*Sigma2) + 
   (35*sum[3][3]*(-1 + xi)*xi)/
    (23.99999999999999*Power(Sigma1,3)*Power(Sigma2,4))\
    - (35*sum[4][0]*(-1 + xi)*xi)/(64.*Power(Sigma1,6)) + 
   (35*sum[4][2]*(-1 + xi)*xi)/
    (32.*Power(Sigma1,4)*Power(Sigma2,3)) + 
   (7*sum[5][1]*(-1 + xi)*xi)/
    (16.*Power(Sigma1,5)*Power(Sigma2,2)) + 
   (7*sum[6][0]*(-1 + xi)*xi)/(96.*Power(Sigma1,6)*Sigma2) - 
   (29*N*Power(xi,2))/(128.*Power(Sigma1,4)) + 
   (87*sum[0][2]*Power(xi,2))/
    (32.*Power(Sigma1,2)*Power(Sigma2,3)) + 
   (29*sum[0][2]*Power(xi,2))/
    (32.*Power(Sigma1,3)*Power(Sigma2,2)) + 
   (87*sum[1][1]*Power(xi,2))/
    (16.*Power(Sigma1,3)*Power(Sigma2,2)) + 
   (29*sum[1][1]*Power(xi,2))/(16.*Power(Sigma1,4)*Sigma2) + 
   (29*sum[2][0]*Power(xi,2))/(32.*Power(Sigma1,5)) + 
   (87*sum[2][0]*Power(xi,2))/(32.*Power(Sigma1,4)*Sigma2) - 
   (183*sum[2][2]*Power(xi,2))/
    (32.*Power(Sigma1,2)*Power(Sigma2,4)) - 
   (183*sum[2][2]*Power(xi,2))/
    (16.*Power(Sigma1,3)*Power(Sigma2,3)) - 
   (61*sum[3][1]*Power(xi,2))/
    (16.*Power(Sigma1,3)*Power(Sigma2,3)) - 
   (61*sum[3][1]*Power(xi,2))/
    (8.*Power(Sigma1,4)*Power(Sigma2,2)) - 
   (61*sum[4][0]*Power(xi,2))/
    (64.*Power(Sigma1,4)*Power(Sigma2,2)) - 
   (61*sum[4][0]*Power(xi,2))/(32.*Power(Sigma1,5)*Sigma2) + 
   (13*sum[0][4]*(-1 + xi)*Power(xi,2))/
    (32.*Power(Sigma1,2)*Power(Sigma2,4)) + 
   (13*sum[1][3]*(-1 + xi)*Power(xi,2))/
    (8.*Power(Sigma1,3)*Power(Sigma2,3)) + 
   (39*sum[2][2]*(-1 + xi)*Power(xi,2))/
    (16.*Power(Sigma1,4)*Power(Sigma2,2)) - 
   (45*sum[2][4]*(-1 + xi)*Power(xi,2))/
    (16.*Power(Sigma1,2)*Power(Sigma2,5)) + 
   (13*sum[3][1]*(-1 + xi)*Power(xi,2))/
    (8.*Power(Sigma1,5)*Sigma2) - 
   (15*sum[3][3]*(-1 + xi)*Power(xi,2))/
    (4.*Power(Sigma1,3)*Power(Sigma2,4)) + 
   (13*sum[4][0]*(-1 + xi)*Power(xi,2))/
    (32.*Power(Sigma1,6)) - 
   (45*sum[4][2]*(-1 + xi)*Power(xi,2))/
    (16.*Power(Sigma1,4)*Power(Sigma2,3)) - 
   (9*sum[5][1]*(-1 + xi)*Power(xi,2))/
    (8.*Power(Sigma1,5)*Power(Sigma2,2)) - 
   (3*sum[6][0]*(-1 + xi)*Power(xi,2))/
    (16.*Power(Sigma1,6)*Sigma2) + 
   (3*N*Power(xi,3))/(32.*Power(Sigma1,4)) - 
   (9*sum[0][2]*Power(xi,3))/
    (8.*Power(Sigma1,2)*Power(Sigma2,3)) - 
   (3*sum[0][2]*Power(xi,3))/
    (8.*Power(Sigma1,3)*Power(Sigma2,2)) - 
   (9*sum[1][1]*Power(xi,3))/
    (4.*Power(Sigma1,3)*Power(Sigma2,2)) - 
   (3*sum[1][1]*Power(xi,3))/(4.*Power(Sigma1,4)*Sigma2) - 
   (3*sum[2][0]*Power(xi,3))/(8.*Power(Sigma1,5)) - 
   (9*sum[2][0]*Power(xi,3))/(8.*Power(Sigma1,4)*Sigma2) + 
   (3*sum[2][2]*Power(xi,3))/
    (Power(Sigma1,2)*Power(Sigma2,4)) + 
   (6*sum[2][2]*Power(xi,3))/
    (Power(Sigma1,3)*Power(Sigma2,3)) + 
   (2*sum[3][1]*Power(xi,3))/
    (Power(Sigma1,3)*Power(Sigma2,3)) + 
   (4*sum[3][1]*Power(xi,3))/
    (Power(Sigma1,4)*Power(Sigma2,2)) + 
   (sum[4][0]*Power(xi,3))/
    (2.*Power(Sigma1,4)*Power(Sigma2,2)) + 
   (sum[4][0]*Power(xi,3))/(Power(Sigma1,5)*Sigma2) - 
   (3*sum[0][4]*(-1 + xi)*Power(xi,3))/
    (32.*Power(Sigma1,2)*Power(Sigma2,4)) - 
   (3*sum[1][3]*(-1 + xi)*Power(xi,3))/
    (8.*Power(Sigma1,3)*Power(Sigma2,3)) - 
   (9*sum[2][2]*(-1 + xi)*Power(xi,3))/
    (16.*Power(Sigma1,4)*Power(Sigma2,2)) + 
   (15*sum[2][4]*(-1 + xi)*Power(xi,3))/
    (16.*Power(Sigma1,2)*Power(Sigma2,5)) - 
   (3*sum[3][1]*(-1 + xi)*Power(xi,3))/
    (8.*Power(Sigma1,5)*Sigma2) + 
   (5*sum[3][3]*(-1 + xi)*Power(xi,3))/
    (4.*Power(Sigma1,3)*Power(Sigma2,4)) - 
   (3*sum[4][0]*(-1 + xi)*Power(xi,3))/
    (32.*Power(Sigma1,6)) + 
   (15*sum[4][2]*(-1 + xi)*Power(xi,3))/
    (16.*Power(Sigma1,4)*Power(Sigma2,3)) + 
   (3*sum[5][1]*(-1 + xi)*Power(xi,3))/
    (8.*Power(Sigma1,5)*Power(Sigma2,2)) + 
   (sum[6][0]*(-1 + xi)*Power(xi,3))/
    (16.*Power(Sigma1,6)*Sigma2) - 
   (N*Power(xi,4))/(64.*Power(Sigma1,4)) + 
   (3*sum[0][2]*Power(xi,4))/
    (16.*Power(Sigma1,2)*Power(Sigma2,3)) + 
   (sum[0][2]*Power(xi,4))/
    (16.*Power(Sigma1,3)*Power(Sigma2,2)) + 
   (3*sum[1][1]*Power(xi,4))/
    (8.*Power(Sigma1,3)*Power(Sigma2,2)) + 
   (sum[1][1]*Power(xi,4))/(8.*Power(Sigma1,4)*Sigma2) + 
   (sum[2][0]*Power(xi,4))/(16.*Power(Sigma1,5)) + 
   (3*sum[2][0]*Power(xi,4))/(16.*Power(Sigma1,4)*Sigma2) - 
   (9*sum[2][2]*Power(xi,4))/
    (16.*Power(Sigma1,2)*Power(Sigma2,4)) - 
   (9*sum[2][2]*Power(xi,4))/
    (8.*Power(Sigma1,3)*Power(Sigma2,3)) - 
   (3*sum[3][1]*Power(xi,4))/
    (8.*Power(Sigma1,3)*Power(Sigma2,3)) - 
   (3*sum[3][1]*Power(xi,4))/
    (4.*Power(Sigma1,4)*Power(Sigma2,2)) - 
   (3*sum[4][0]*Power(xi,4))/
    (32.*Power(Sigma1,4)*Power(Sigma2,2)) - 
   (3*sum[4][0]*Power(xi,4))/(16.*Power(Sigma1,5)*Sigma2) - 
   (sum[0][4]*(-1 + xi)*xi*(35 - 26*xi + 6*Power(xi,2)))/
    (64.*Power(Sigma2,6)) - 
   (sum[0][4]*(-1 + xi)*xi*(35 - 26*xi + 6*Power(xi,2)))/
    (32.*Sigma1*Power(Sigma2,5)) - 
   (sum[1][3]*(-1 + xi)*xi*(35 - 26*xi + 6*Power(xi,2)))/
    (16.*Sigma1*Power(Sigma2,5)) - 
   (sum[1][3]*(-1 + xi)*xi*(35 - 26*xi + 6*Power(xi,2)))/
    (8.*Power(Sigma1,2)*Power(Sigma2,4)) + 
   (sum[0][6]*(-1 + xi)*xi*(7 - 18*xi + 6*Power(xi,2)))/
    (96.*Power(Sigma2,7)) + 
   (sum[0][6]*(-1 + xi)*xi*(7 - 18*xi + 6*Power(xi,2)))/
    (96.*Sigma1*Power(Sigma2,6)) + 
   (sum[1][5]*(-1 + xi)*xi*(7 - 18*xi + 6*Power(xi,2)))/
    (16.*Sigma1*Power(Sigma2,6)) + 
   (sum[1][5]*(-1 + xi)*xi*(7 - 18*xi + 6*Power(xi,2)))/
    (16.*Power(Sigma1,2)*Power(Sigma2,5)) + 
   (5*sum[2][4]*(-1 + xi)*xi*(7 - 18*xi + 6*Power(xi,2)))/
    (32.*Power(Sigma1,3)*Power(Sigma2,4)) + 
   (5*sum[3][3]*(-1 + xi)*xi*(7 - 18*xi + 6*Power(xi,2)))/
    (23.99999999999999*Power(Sigma1,4)*Power(Sigma2,3))\
    - (7*sum[2][6]*(-1 + xi)*xi*(1 - 6*xi + 6*Power(xi,2)))/
    (96.*Power(Sigma1,2)*Power(Sigma2,6)) - 
   (7*sum[3][5]*(-1 + xi)*xi*(1 - 6*xi + 6*Power(xi,2)))/
    (47.99999999999999*Power(Sigma1,3)*Power(Sigma2,5))\
    - (35*sum[4][4]*(-1 + xi)*xi*(1 - 6*xi + 6*Power(xi,2)))/
    (191.9999999999999*Power(Sigma1,4)*Power(Sigma2,4))\
    - (7*sum[5][3]*(-1 + xi)*xi*(1 - 6*xi + 6*Power(xi,2)))/
    (47.99999999999999*Power(Sigma1,5)*Power(Sigma2,3))\
    - (N*xi*(-35 + 29*xi - 12*Power(xi,2) + 
        2*Power(xi,3)))/(128.*Power(Sigma2,4)) - 
   (N*xi*(-35 + 29*xi - 12*Power(xi,2) + 2*Power(xi,3)))/
    (32.*Sigma1*Power(Sigma2,3)) - 
   (3*N*xi*(-35 + 29*xi - 12*Power(xi,2) + 
        2*Power(xi,3)))/
    (64.*Power(Sigma1,2)*Power(Sigma2,2)) - 
   (N*xi*(-35 + 29*xi - 12*Power(xi,2) + 2*Power(xi,3)))/
    (32.*Power(Sigma1,3)*Sigma2) + 
   (sum[0][2]*xi*(-35 + 29*xi - 12*Power(xi,2) + 
        2*Power(xi,3)))/(32.*Power(Sigma2,5)) + 
   (3*sum[0][2]*xi*(-35 + 29*xi - 12*Power(xi,2) + 
        2*Power(xi,3)))/(32.*Sigma1*Power(Sigma2,4)) + 
   (sum[1][1]*xi*(-35 + 29*xi - 12*Power(xi,2) + 
        2*Power(xi,3)))/(16.*Sigma1*Power(Sigma2,4)) + 
   (3*sum[1][1]*xi*(-35 + 29*xi - 12*Power(xi,2) + 
        2*Power(xi,3)))/
    (16.*Power(Sigma1,2)*Power(Sigma2,3)) + 
   (sum[2][0]*xi*(-35 + 29*xi - 12*Power(xi,2) + 
        2*Power(xi,3)))/
    (32.*Power(Sigma1,2)*Power(Sigma2,3)) + 
   (3*sum[2][0]*xi*(-35 + 29*xi - 12*Power(xi,2) + 
        2*Power(xi,3)))/
    (32.*Power(Sigma1,3)*Power(Sigma2,2)) + 
   (5*sum[4][2]*xi*(-7 + 25*xi - 24*Power(xi,2) + 
        6*Power(xi,3)))/
    (32.*Power(Sigma1,5)*Power(Sigma2,2)) + 
   (sum[5][1]*xi*(-7 + 25*xi - 24*Power(xi,2) + 
        6*Power(xi,3)))/(16.*Power(Sigma1,6)*Sigma2) + 
   (sum[6][0]*xi*(-7 + 25*xi - 24*Power(xi,2) + 
        6*Power(xi,3)))/(96.*Power(Sigma1,7)) - 
   (sum[0][8]*xi*(-1 + 7*xi - 12*Power(xi,2) + 
        6*Power(xi,3)))/
    (383.9999999999999*Power(Sigma2,8)) - 
   (sum[1][7]*xi*(-1 + 7*xi - 12*Power(xi,2) + 
        6*Power(xi,3)))/
    (47.99999999999999*Sigma1*Power(Sigma2,7)) - 
   (7*sum[6][2]*xi*(-1 + 7*xi - 12*Power(xi,2) + 
        6*Power(xi,3)))/
    (96.*Power(Sigma1,6)*Power(Sigma2,2)) - 
   (sum[7][1]*xi*(-1 + 7*xi - 12*Power(xi,2) + 
        6*Power(xi,3)))/
    (47.99999999999999*Power(Sigma1,7)*Sigma2) - 
   (sum[8][0]*xi*(-1 + 7*xi - 12*Power(xi,2) + 
        6*Power(xi,3)))/
    (383.9999999999999*Power(Sigma1,8));

	/*&&&&&&&&&
	if(c1-c1 || c2-c2 || c3-c3 || c4-c4 || c0-c0) printf("%e\t%e\t%e\t%e\t%d\t%e\t%e\t%e\n",c4,c3,c2,c1,N,xi,Sigma1,Sigma2); */

	/* compute AlphaBar */
	*AlphaBar = GetAlpha(c4,c3,c2,c1);

	/* return statistic */
	return  c0 
		+ c1 * (*AlphaBar) 
		+ c2 * (*AlphaBar) * (*AlphaBar)
		+ c3 * (*AlphaBar) * (*AlphaBar) * (*AlphaBar)
		+ c4 * (*AlphaBar) * (*AlphaBar) * (*AlphaBar) * (*AlphaBar);

}
