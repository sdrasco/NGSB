/*
 * File: ApproxLogLikelihood.c  by: Steve Drasco 
 *                                                           
 * This program computes the log of the two-detector likelihood  
 * function for a non-Gaussian stochastic background        
 * (characterized by user specified sigma1, sigma2, alpha, and xi)
 * on Gaussian Noise.  
 *      
 * NOTE:
 *      sums is a 8x8 matrix with components:  sums[a][b] = SumOverj( h1[j]^a h2[j]^b ) 
 *	(Need only components with a+b = 2,4,6,8  that is need only 24 sums)
 */

#include "NGSB.h"

double PairApproxLogLikelihood(int N, double Sigma1, double Sigma2, double alpha, double xi, double **sum)
{
        double  stat;
	
/* here's the monster! */
stat=(192*Power(Sigma1,7)*Power(Sigma2,7)*
      (2*N*Sigma1*Sigma2 - Sigma1*sum[0][2] - Sigma2*sum[2][0])\
      + alpha*(192*Power(Sigma1,6)*Power(Sigma2,6)*
         (-(Sigma1*(N*Sigma1*Sigma2 + 
                N*Power(Sigma2,2) - Sigma1*sum[0][2] - 
                2*Sigma2*sum[1][1])) + Power(Sigma2,2)*sum[2][0])
          + 48*alpha*Power(Sigma1,4)*Power(Sigma2,4)*
         (Sigma1*(3*N*Power(Sigma1,3)*Power(Sigma2,2) + 
              6*N*Power(Sigma1,2)*Power(Sigma2,3) + 
              3*N*Sigma1*Power(Sigma2,4) - 
              6*Power(Sigma1,3)*Sigma2*sum[0][2] - 
              6*Power(Sigma1,2)*Power(Sigma2,2)*sum[0][2] + 
              Power(Sigma1,3)*sum[0][4] - 
              12*Power(Sigma1,2)*Power(Sigma2,2)*sum[1][1] - 
              12*Sigma1*Power(Sigma2,3)*sum[1][1] + 
              4*Power(Sigma1,2)*Sigma2*sum[1][3] - 
              6*Sigma1*Power(Sigma2,3)*sum[2][0] - 
              6*Power(Sigma2,4)*sum[2][0] + 
              6*Sigma1*Power(Sigma2,2)*sum[2][2] + 
              4*Power(Sigma2,3)*sum[3][1]) + 
           Power(Sigma2,4)*sum[4][0]) + 
        8*Power(alpha,2)*Power(Sigma1,2)*Power(Sigma2,2)*
         (-(Sigma1*Sigma2*
              (15*N*Power(Sigma1,5)*Power(Sigma2,2) + 
                45*N*Power(Sigma1,4)*Power(Sigma2,3) + 
                45*N*Power(Sigma1,3)*Power(Sigma2,4) + 
                15*N*Power(Sigma1,2)*Power(Sigma2,5) - 
                45*Power(Sigma1,5)*Sigma2*sum[0][2] - 
                90*Power(Sigma1,4)*Power(Sigma2,2)*
                 sum[0][2] - 
                45*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[0][2] + 
                15*Power(Sigma1,4)*Sigma2*sum[0][4] - 
                90*Power(Sigma1,4)*Power(Sigma2,2)*
                 sum[1][1] - 
                180*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[1][1] - 
                90*Power(Sigma1,2)*Power(Sigma2,4)*
                 sum[1][1] + 
                60*Power(Sigma1,4)*Sigma2*sum[1][3] + 
                60*Power(Sigma1,3)*Power(Sigma2,2)*
                 sum[1][3] - 6*Power(Sigma1,4)*sum[1][5] - 
                45*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[0][2] - 
                90*Power(Sigma1,2)*Power(Sigma2,4)*
                 sum[0][2] - 
                45*Sigma1*Power(Sigma2,5)*sum[0][2] + 
                90*Power(Sigma1,3)*Power(Sigma2,2)*
                 sum[2][2] + 
                90*Power(Sigma1,2)*Power(Sigma2,3)*
                 sum[2][2] - 
                15*Power(Sigma1,3)*Sigma2*sum[2][4] + 
                60*Power(Sigma1,2)*Power(Sigma2,3)*
                 sum[3][1] + 
                60*Sigma1*Power(Sigma2,4)*sum[3][1] - 
                20*Power(Sigma1,2)*Power(Sigma2,2)*
                 sum[3][3] + 15*Power(Sigma1,5)*sum[4][0] + 
                15*Sigma1*Power(Sigma2,4)*sum[4][0] + 
                15*Power(Sigma2,5)*sum[4][0] - 
                15*Sigma1*Power(Sigma2,3)*sum[4][2] - 
                6*Power(Sigma2,4)*sum[5][1])) + 
           (Power(Sigma1,6) + Power(Sigma2,6))*sum[6][0]) + 
        Power(alpha,3)*
         (Sigma1*(105*N*Power(Sigma1,7)*
               Power(Sigma2,4) + 
              420*N*Power(Sigma1,6)*Power(Sigma2,5) + 
              630*N*Power(Sigma1,5)*Power(Sigma2,6) + 
              420*N*Power(Sigma1,4)*Power(Sigma2,7) + 
              105*N*Power(Sigma1,3)*Power(Sigma2,8) - 
              420*Power(Sigma1,7)*Power(Sigma2,3)*
               sum[0][2] - 
              1260*Power(Sigma1,6)*Power(Sigma2,4)*
               sum[0][2] - 
              1260*Power(Sigma1,5)*Power(Sigma2,5)*
               sum[0][2] - 
              420*Power(Sigma1,4)*Power(Sigma2,6)*
               sum[0][2] + 
              210*Power(Sigma1,7)*Power(Sigma2,2)*
               sum[0][4] + 
              420*Power(Sigma1,6)*Power(Sigma2,3)*
               sum[0][4] + 
              210*Power(Sigma1,5)*Power(Sigma2,4)*
               sum[0][4] - 28*Power(Sigma1,7)*Sigma2*sum[0][6] - 
              28*Power(Sigma1,6)*Power(Sigma2,2)*sum[0][6] + 
              Power(Sigma1,7)*sum[0][8] - 
              840*Power(Sigma1,6)*Power(Sigma2,4)*
               sum[1][1] - 
              2520*Power(Sigma1,5)*Power(Sigma2,5)*
               sum[1][1] - 
              2520*Power(Sigma1,4)*Power(Sigma2,6)*
               sum[1][1] - 
              840*Power(Sigma1,3)*Power(Sigma2,7)*
               sum[1][1] + 
              840*Power(Sigma1,6)*Power(Sigma2,3)*
               sum[1][3] + 
              1680*Power(Sigma1,5)*Power(Sigma2,4)*
               sum[1][3] + 
              840*Power(Sigma1,4)*Power(Sigma2,5)*
               sum[1][3] - 
              168*Power(Sigma1,6)*Power(Sigma2,2)*
               sum[1][5] - 
              168*Power(Sigma1,5)*Power(Sigma2,3)*
               sum[1][5] + 8*Power(Sigma1,6)*Sigma2*sum[1][7] - 
              420*Power(Sigma1,5)*Power(Sigma2,5)*
               sum[0][2] - 
              1260*Power(Sigma1,4)*Power(Sigma2,6)*
               sum[0][2] - 
              1260*Power(Sigma1,3)*Power(Sigma2,7)*
               sum[0][2] - 
              420*Power(Sigma1,2)*Power(Sigma2,8)*
               sum[0][2] + 
              1260*Power(Sigma1,5)*Power(Sigma2,4)*
               sum[2][2] + 
              2520*Power(Sigma1,4)*Power(Sigma2,5)*
               sum[2][2] + 
              1260*Power(Sigma1,3)*Power(Sigma2,6)*
               sum[2][2] - 
              420*Power(Sigma1,5)*Power(Sigma2,3)*
               sum[2][4] - 
              420*Power(Sigma1,4)*Power(Sigma2,4)*
               sum[2][4] + 
              28*Power(Sigma1,5)*Power(Sigma2,2)*sum[2][6] + 
              840*Power(Sigma1,4)*Power(Sigma2,5)*
               sum[3][1] + 
              1680*Power(Sigma1,3)*Power(Sigma2,6)*
               sum[3][1] + 
              840*Power(Sigma1,2)*Power(Sigma2,7)*
               sum[3][1] - 
              560*Power(Sigma1,4)*Power(Sigma2,4)*
               sum[3][3] - 
              560*Power(Sigma1,3)*Power(Sigma2,5)*
               sum[3][3] + 
              56*Power(Sigma1,4)*Power(Sigma2,3)*sum[3][5] + 
              210*Power(Sigma1,3)*Power(Sigma2,6)*
               sum[4][0] + 
              420*Power(Sigma1,2)*Power(Sigma2,7)*
               sum[4][0] + 
              210*Sigma1*Power(Sigma2,8)*sum[4][0] - 
              420*Power(Sigma1,3)*Power(Sigma2,5)*
               sum[4][2] - 
              420*Power(Sigma1,2)*Power(Sigma2,6)*
               sum[4][2] + 
              70*Power(Sigma1,3)*Power(Sigma2,4)*sum[4][4] - 
              168*Power(Sigma1,2)*Power(Sigma2,6)*
               sum[5][1] - 
              168*Sigma1*Power(Sigma2,7)*sum[5][1] + 
              56*Power(Sigma1,2)*Power(Sigma2,5)*sum[5][3] - 
              28*Sigma1*Power(Sigma2,7)*sum[6][0] - 
              28*Power(Sigma2,8)*sum[6][0] + 
              28*Sigma1*Power(Sigma2,6)*sum[6][2] + 
              8*Power(Sigma2,7)*sum[7][1]) + 
           Power(Sigma2,8)*sum[8][0]))*xi - 
     Power(alpha,2)*(48*Power(Sigma1,4)*Power(Sigma2,4)*
         (Sigma1*(N*Power(Sigma1,3)*Power(Sigma2,2) + 
              2*N*Power(Sigma1,2)*Power(Sigma2,3) + 
              N*Sigma1*Power(Sigma2,4) - 
              2*Power(Sigma1,3)*Sigma2*sum[0][2] - 
              2*Power(Sigma1,2)*Power(Sigma2,2)*sum[0][2] + 
              Power(Sigma1,3)*sum[0][4] - 
              4*Power(Sigma1,2)*Power(Sigma2,2)*sum[1][1] - 
              4*Sigma1*Power(Sigma2,3)*sum[1][1] + 
              4*Power(Sigma1,2)*Sigma2*sum[1][3] - 
              2*Sigma1*Power(Sigma2,3)*sum[0][2] - 
              2*Power(Sigma2,4)*sum[0][2] + 
              6*Sigma1*Power(Sigma2,2)*sum[2][2] + 
              4*Power(Sigma2,3)*sum[3][1]) + 
           Power(Sigma2,4)*sum[4][0]) + 
        24*alpha*Power(Sigma1,2)*Power(Sigma2,2)*
         (-(Sigma1*Sigma2*
              (3*N*Power(Sigma1,5)*Power(Sigma2,2) + 
                9*N*Power(Sigma1,4)*Power(Sigma2,3) + 
                9*N*Power(Sigma1,3)*Power(Sigma2,4) + 
                3*N*Power(Sigma1,2)*Power(Sigma2,5) - 
                9*Power(Sigma1,5)*Sigma2*sum[0][2] - 
                18*Power(Sigma1,4)*Power(Sigma2,2)*
                 sum[0][2] - 
                9*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[0][2] + 
                7*Power(Sigma1,4)*Sigma2*sum[0][4] - 
                18*Power(Sigma1,4)*Power(Sigma2,2)*
                 sum[1][1] - 
                36*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[1][1] - 
                18*Power(Sigma1,2)*Power(Sigma2,4)*
                 sum[1][1] + 
                28*Power(Sigma1,4)*Sigma2*sum[1][3] + 
                28*Power(Sigma1,3)*Power(Sigma2,2)*
                 sum[1][3] - 6*Power(Sigma1,4)*sum[1][5] - 
                9*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[0][2] - 
                18*Power(Sigma1,2)*Power(Sigma2,4)*
                 sum[0][2] - 
                9*Sigma1*Power(Sigma2,5)*sum[0][2] + 
                42*Power(Sigma1,3)*Power(Sigma2,2)*
                 sum[2][2] + 
                42*Power(Sigma1,2)*Power(Sigma2,3)*
                 sum[2][2] - 
                15*Power(Sigma1,3)*Sigma2*sum[2][4] + 
                28*Power(Sigma1,2)*Power(Sigma2,3)*
                 sum[3][1] + 
                28*Sigma1*Power(Sigma2,4)*sum[3][1] - 
                20*Power(Sigma1,2)*Power(Sigma2,2)*
                 sum[3][3] + 7*Power(Sigma1,5)*sum[4][0] + 
                7*Sigma1*Power(Sigma2,4)*sum[4][0] + 
                7*Power(Sigma2,5)*sum[4][0] - 
                15*Sigma1*Power(Sigma2,3)*sum[4][2] - 
                6*Power(Sigma2,4)*sum[5][1])) + 
           (Power(Sigma1,6) + Power(Sigma2,6))*sum[6][0]) + 
        Power(alpha,2)*
         (Sigma1*(87*N*Power(Sigma1,7)*Power(Sigma2,4) + 
              348*N*Power(Sigma1,6)*Power(Sigma2,5) + 
              522*N*Power(Sigma1,5)*Power(Sigma2,6) + 
              348*N*Power(Sigma1,4)*Power(Sigma2,7) + 
              87*N*Power(Sigma1,3)*Power(Sigma2,8) - 
              348*Power(Sigma1,7)*Power(Sigma2,3)*
               sum[0][2] - 
              1044*Power(Sigma1,6)*Power(Sigma2,4)*
               sum[0][2] - 
              1044*Power(Sigma1,5)*Power(Sigma2,5)*
               sum[0][2] - 
              348*Power(Sigma1,4)*Power(Sigma2,6)*
               sum[0][2] + 
              366*Power(Sigma1,7)*Power(Sigma2,2)*
               sum[0][4] + 
              732*Power(Sigma1,6)*Power(Sigma2,3)*
               sum[0][4] + 
              366*Power(Sigma1,5)*Power(Sigma2,4)*
               sum[0][4] - 
              100*Power(Sigma1,7)*Sigma2*sum[0][6] - 
              100*Power(Sigma1,6)*Power(Sigma2,2)*
               sum[0][6] + 7*Power(Sigma1,7)*sum[0][8] - 
              696*Power(Sigma1,6)*Power(Sigma2,4)*
               sum[1][1] - 
              2088*Power(Sigma1,5)*Power(Sigma2,5)*
               sum[1][1] - 
              2088*Power(Sigma1,4)*Power(Sigma2,6)*
               sum[1][1] - 
              696*Power(Sigma1,3)*Power(Sigma2,7)*
               sum[1][1] + 
              1464*Power(Sigma1,6)*Power(Sigma2,3)*
               sum[1][3] + 
              2928*Power(Sigma1,5)*Power(Sigma2,4)*
               sum[1][3] + 
              1464*Power(Sigma1,4)*Power(Sigma2,5)*
               sum[1][3] - 
              600*Power(Sigma1,6)*Power(Sigma2,2)*
               sum[1][5] - 
              600*Power(Sigma1,5)*Power(Sigma2,3)*
               sum[1][5] + 56*Power(Sigma1,6)*Sigma2*sum[1][7] - 
              348*Power(Sigma1,5)*Power(Sigma2,5)*
               sum[0][2] - 
              1044*Power(Sigma1,4)*Power(Sigma2,6)*
               sum[0][2] - 
              1044*Power(Sigma1,3)*Power(Sigma2,7)*
               sum[0][2] - 
              348*Power(Sigma1,2)*Power(Sigma2,8)*
               sum[0][2] + 
              2196*Power(Sigma1,5)*Power(Sigma2,4)*
               sum[2][2] + 
              4392*Power(Sigma1,4)*Power(Sigma2,5)*
               sum[2][2] + 
              2196*Power(Sigma1,3)*Power(Sigma2,6)*
               sum[2][2] - 
              1500*Power(Sigma1,5)*Power(Sigma2,3)*
               sum[2][4] - 
              1500*Power(Sigma1,4)*Power(Sigma2,4)*
               sum[2][4] + 
              196*Power(Sigma1,5)*Power(Sigma2,2)*
               sum[2][6] + 
              1464*Power(Sigma1,4)*Power(Sigma2,5)*
               sum[3][1] + 
              2928*Power(Sigma1,3)*Power(Sigma2,6)*
               sum[3][1] + 
              1464*Power(Sigma1,2)*Power(Sigma2,7)*
               sum[3][1] - 
              2000*Power(Sigma1,4)*Power(Sigma2,4)*
               sum[3][3] - 
              2000*Power(Sigma1,3)*Power(Sigma2,5)*
               sum[3][3] + 
              392*Power(Sigma1,4)*Power(Sigma2,3)*
               sum[3][5] + 
              366*Power(Sigma1,3)*Power(Sigma2,6)*
               sum[4][0] + 
              732*Power(Sigma1,2)*Power(Sigma2,7)*
               sum[4][0] + 
              366*Sigma1*Power(Sigma2,8)*sum[4][0] - 
              1500*Power(Sigma1,3)*Power(Sigma2,5)*
               sum[4][2] - 
              1500*Power(Sigma1,2)*Power(Sigma2,6)*
               sum[4][2] + 
              490*Power(Sigma1,3)*Power(Sigma2,4)*
               sum[4][4] - 
              600*Power(Sigma1,2)*Power(Sigma2,6)*
               sum[5][1] - 
              600*Sigma1*Power(Sigma2,7)*sum[5][1] + 
              392*Power(Sigma1,2)*Power(Sigma2,5)*
               sum[5][3] - 
              100*Sigma1*Power(Sigma2,7)*sum[6][0] - 
              100*Power(Sigma2,8)*sum[6][0] + 
              196*Sigma1*Power(Sigma2,6)*sum[6][2] + 
              56*Power(Sigma2,7)*sum[7][1]) + 
           7*Power(Sigma2,8)*sum[8][0]))*Power(xi,2) + 
     4*Power(alpha,3)*(4*Power(Sigma1,2)*Power(Sigma2,2)*
         (-(Sigma1*Sigma2*
              (N*Power(Sigma1,5)*Power(Sigma2,2) + 
                3*N*Power(Sigma1,4)*Power(Sigma2,3) + 
                3*N*Power(Sigma1,3)*Power(Sigma2,4) + 
                N*Power(Sigma1,2)*Power(Sigma2,5) - 
                3*Power(Sigma1,5)*Sigma2*sum[0][2] - 
                6*Power(Sigma1,4)*Power(Sigma2,2)*
                 sum[0][2] - 
                3*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[0][2] + 
                3*Power(Sigma1,4)*Sigma2*sum[0][4] - 
                6*Power(Sigma1,4)*Power(Sigma2,2)*
                 sum[1][1] - 
                12*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[1][1] - 
                6*Power(Sigma1,2)*Power(Sigma2,4)*
                 sum[1][1] + 
                12*Power(Sigma1,4)*Sigma2*sum[1][3] + 
                12*Power(Sigma1,3)*Power(Sigma2,2)*
                 sum[1][3] - 6*Power(Sigma1,4)*sum[1][5] - 
                3*Power(Sigma1,3)*Power(Sigma2,3)*
                 sum[0][2] - 
                6*Power(Sigma1,2)*Power(Sigma2,4)*
                 sum[0][2] - 
                3*Sigma1*Power(Sigma2,5)*sum[0][2] + 
                18*Power(Sigma1,3)*Power(Sigma2,2)*
                 sum[2][2] + 
                18*Power(Sigma1,2)*Power(Sigma2,3)*
                 sum[2][2] - 
                15*Power(Sigma1,3)*Sigma2*sum[2][4] + 
                12*Power(Sigma1,2)*Power(Sigma2,3)*
                 sum[3][1] + 
                12*Sigma1*Power(Sigma2,4)*sum[3][1] - 
                20*Power(Sigma1,2)*Power(Sigma2,2)*
                 sum[3][3] + 3*Power(Sigma1,5)*sum[4][0] + 
                3*Sigma1*Power(Sigma2,4)*sum[4][0] + 
                3*Power(Sigma2,5)*sum[4][0] - 
                15*Sigma1*Power(Sigma2,3)*sum[4][2] - 
                6*Power(Sigma2,4)*sum[5][1])) + 
           (Power(Sigma1,6) + Power(Sigma2,6))*sum[6][0]) + 
        3*alpha*(Sigma1*
            (3*N*Power(Sigma1,7)*Power(Sigma2,4) + 
              12*N*Power(Sigma1,6)*Power(Sigma2,5) + 
              18*N*Power(Sigma1,5)*Power(Sigma2,6) + 
              12*N*Power(Sigma1,4)*Power(Sigma2,7) + 
              3*N*Power(Sigma1,3)*Power(Sigma2,8) - 
              12*Power(Sigma1,7)*Power(Sigma2,3)*sum[0][2] - 
              36*Power(Sigma1,6)*Power(Sigma2,4)*sum[0][2] - 
              36*Power(Sigma1,5)*Power(Sigma2,5)*sum[0][2] - 
              12*Power(Sigma1,4)*Power(Sigma2,6)*sum[0][2] + 
              16*Power(Sigma1,7)*Power(Sigma2,2)*sum[0][4] + 
              32*Power(Sigma1,6)*Power(Sigma2,3)*sum[0][4] + 
              16*Power(Sigma1,5)*Power(Sigma2,4)*sum[0][4] - 
              8*Power(Sigma1,7)*Sigma2*sum[0][6] - 
              8*Power(Sigma1,6)*Power(Sigma2,2)*sum[0][6] + 
              Power(Sigma1,7)*sum[0][8] - 
              24*Power(Sigma1,6)*Power(Sigma2,4)*sum[1][1] - 
              72*Power(Sigma1,5)*Power(Sigma2,5)*sum[1][1] - 
              72*Power(Sigma1,4)*Power(Sigma2,6)*sum[1][1] - 
              24*Power(Sigma1,3)*Power(Sigma2,7)*sum[1][1] + 
              64*Power(Sigma1,6)*Power(Sigma2,3)*sum[1][3] + 
              128*Power(Sigma1,5)*Power(Sigma2,4)*
               sum[1][3] + 
              64*Power(Sigma1,4)*Power(Sigma2,5)*sum[1][3] - 
              48*Power(Sigma1,6)*Power(Sigma2,2)*sum[1][5] - 
              48*Power(Sigma1,5)*Power(Sigma2,3)*sum[1][5] + 
              8*Power(Sigma1,6)*Sigma2*sum[1][7] - 
              12*Power(Sigma1,5)*Power(Sigma2,5)*sum[0][2] - 
              36*Power(Sigma1,4)*Power(Sigma2,6)*sum[0][2] - 
              36*Power(Sigma1,3)*Power(Sigma2,7)*sum[0][2] - 
              12*Power(Sigma1,2)*Power(Sigma2,8)*sum[0][2] + 
              96*Power(Sigma1,5)*Power(Sigma2,4)*sum[2][2] + 
              192*Power(Sigma1,4)*Power(Sigma2,5)*
               sum[2][2] + 
              96*Power(Sigma1,3)*Power(Sigma2,6)*sum[2][2] - 
              120*Power(Sigma1,5)*Power(Sigma2,3)*
               sum[2][4] - 
              120*Power(Sigma1,4)*Power(Sigma2,4)*
               sum[2][4] + 
              28*Power(Sigma1,5)*Power(Sigma2,2)*sum[2][6] + 
              64*Power(Sigma1,4)*Power(Sigma2,5)*sum[3][1] + 
              128*Power(Sigma1,3)*Power(Sigma2,6)*
               sum[3][1] + 
              64*Power(Sigma1,2)*Power(Sigma2,7)*sum[3][1] - 
              160*Power(Sigma1,4)*Power(Sigma2,4)*
               sum[3][3] - 
              160*Power(Sigma1,3)*Power(Sigma2,5)*
               sum[3][3] + 
              56*Power(Sigma1,4)*Power(Sigma2,3)*sum[3][5] + 
              16*Power(Sigma1,3)*Power(Sigma2,6)*sum[4][0] + 
              32*Power(Sigma1,2)*Power(Sigma2,7)*sum[4][0] + 
              16*Sigma1*Power(Sigma2,8)*sum[4][0] - 
              120*Power(Sigma1,3)*Power(Sigma2,5)*
               sum[4][2] - 
              120*Power(Sigma1,2)*Power(Sigma2,6)*
               sum[4][2] + 
              70*Power(Sigma1,3)*Power(Sigma2,4)*sum[4][4] - 
              48*Power(Sigma1,2)*Power(Sigma2,6)*sum[5][1] - 
              48*Sigma1*Power(Sigma2,7)*sum[5][1] + 
              56*Power(Sigma1,2)*Power(Sigma2,5)*sum[5][3] - 
              8*Sigma1*Power(Sigma2,7)*sum[6][0] - 
              8*Power(Sigma2,8)*sum[6][0] + 
              28*Sigma1*Power(Sigma2,6)*sum[6][2] + 
              8*Power(Sigma2,7)*sum[7][1]) + 
           Power(Sigma2,8)*sum[8][0]))*Power(xi,3) - 
     6*Power(alpha,4)*(N*Power(Sigma1,4)*Power(Sigma2,4)*
         Power(Sigma1 + Sigma2,4) - 
        Sigma1*(4*Power(Sigma1,7)*Power(Sigma2,3)*
            sum[0][2] + 12*Power(Sigma1,6)*Power(Sigma2,4)*
            sum[0][2] + 12*Power(Sigma1,5)*Power(Sigma2,5)*
            sum[0][2] + 4*Power(Sigma1,4)*Power(Sigma2,6)*
            sum[0][2] - 6*Power(Sigma1,7)*Power(Sigma2,2)*
            sum[0][4] - 12*Power(Sigma1,6)*Power(Sigma2,3)*
            sum[0][4] - 6*Power(Sigma1,5)*Power(Sigma2,4)*
            sum[0][4] + 4*Power(Sigma1,7)*Sigma2*sum[0][6] + 
           4*Power(Sigma1,6)*Power(Sigma2,2)*sum[0][6] - 
           Power(Sigma1,7)*sum[0][8] + 
           8*Power(Sigma1,6)*Power(Sigma2,4)*sum[1][1] + 
           24*Power(Sigma1,5)*Power(Sigma2,5)*sum[1][1] + 
           24*Power(Sigma1,4)*Power(Sigma2,6)*sum[1][1] + 
           8*Power(Sigma1,3)*Power(Sigma2,7)*sum[1][1] - 
           24*Power(Sigma1,6)*Power(Sigma2,3)*sum[1][3] - 
           48*Power(Sigma1,5)*Power(Sigma2,4)*sum[1][3] - 
           24*Power(Sigma1,4)*Power(Sigma2,5)*sum[1][3] + 
           24*Power(Sigma1,6)*Power(Sigma2,2)*sum[1][5] + 
           24*Power(Sigma1,5)*Power(Sigma2,3)*sum[1][5] - 
           8*Power(Sigma1,6)*Sigma2*sum[1][7] + 
           4*Power(Sigma1,5)*Power(Sigma2,5)*sum[0][2] + 
           12*Power(Sigma1,4)*Power(Sigma2,6)*sum[0][2] + 
           12*Power(Sigma1,3)*Power(Sigma2,7)*sum[0][2] + 
           4*Power(Sigma1,2)*Power(Sigma2,8)*sum[0][2] - 
           36*Power(Sigma1,5)*Power(Sigma2,4)*sum[2][2] - 
           72*Power(Sigma1,4)*Power(Sigma2,5)*sum[2][2] - 
           36*Power(Sigma1,3)*Power(Sigma2,6)*sum[2][2] + 
           60*Power(Sigma1,5)*Power(Sigma2,3)*sum[2][4] + 
           60*Power(Sigma1,4)*Power(Sigma2,4)*sum[2][4] - 
           28*Power(Sigma1,5)*Power(Sigma2,2)*sum[2][6] - 
           24*Power(Sigma1,4)*Power(Sigma2,5)*sum[3][1] - 
           48*Power(Sigma1,3)*Power(Sigma2,6)*sum[3][1] - 
           24*Power(Sigma1,2)*Power(Sigma2,7)*sum[3][1] + 
           80*Power(Sigma1,4)*Power(Sigma2,4)*sum[3][3] + 
           80*Power(Sigma1,3)*Power(Sigma2,5)*sum[3][3] - 
           56*Power(Sigma1,4)*Power(Sigma2,3)*sum[3][5] - 
           6*Power(Sigma1,3)*Power(Sigma2,6)*sum[4][0] - 
           12*Power(Sigma1,2)*Power(Sigma2,7)*sum[4][0] - 
           6*Sigma1*Power(Sigma2,8)*sum[4][0] + 
           60*Power(Sigma1,3)*Power(Sigma2,5)*sum[4][2] + 
           60*Power(Sigma1,2)*Power(Sigma2,6)*sum[4][2] - 
           70*Power(Sigma1,3)*Power(Sigma2,4)*sum[4][4] + 
           24*Power(Sigma1,2)*Power(Sigma2,6)*sum[5][1] + 
           24*Sigma1*Power(Sigma2,7)*sum[5][1] - 
           56*Power(Sigma1,2)*Power(Sigma2,5)*sum[5][3] + 
           4*Sigma1*Power(Sigma2,7)*sum[6][0] + 
           4*Power(Sigma2,8)*sum[6][0] - 
           28*Sigma1*Power(Sigma2,6)*sum[6][2] - 
           8*Power(Sigma2,7)*sum[7][1]) + 
        Power(Sigma2,8)*sum[8][0])*Power(xi,4) + 
     192*N*Power(Sigma1,8)*Power(Sigma2,8)*
      log((  (sum[2][0]/N)*(sum[0][2]/N)    )/(Sigma1*Sigma2)))/
   (384*Power(Sigma1,8)*Power(Sigma2,8));

	/* exit */
	return stat;
}
