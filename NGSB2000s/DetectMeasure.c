/*
 * File: DetectMeasure.c     by: Steve Drasco 
 *                                                           
 * STAND ALONE                                               
 * This program tries to detect an artificial non-Gaussian   
 * stochastic background superimposed on artificial data.   
 * User specifies output file, true values for sigma,        
 * alpha, and xi.                                            
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double	sigma, alpha, xi, *h, lambda, LambdaThresh, Dsigma, Dalpha, Dxi;
	int	i, detect=0;

	/* check inputs */
	if(argc != 4) {
                printf("I need inputs: sigma, alpha, and xi\n");
                return;
        }

	/* read other inputs */
	sigma = atof(argv[1]);
	alpha = atof(argv[2]);
	xi = atof(argv[3]);

	/* allocate memory for artificial data */
	h = (double *) malloc(LENGTH*sizeof(double));

	/* make artificial data */
	NGSBData(LENGTH, SIGMATRUE, ALPHATRUE, XITRUE, h); 

	/* compute exact likelihood ratio once just to test */
	lambda = ExactLogLikelihood(LENGTH, sigma, alpha, xi, h);	

	/* compute threshold */
	LambdaThresh = PTHRESH / (1.0 - PTHRESH);

	/* check for detection */
	if(lambda > LambdaThresh) {
		detect=1;
		Dsigma = sigma;
		Dalpha = alpha;
		Dxi = xi;
	}

        /* write results to screen */
        printf("\t--------------------------------- Detect & Measure Results -----------------------------------\n\n");
        if(detect) {
                printf("\t\tDETECTION !!\n");
		printf("\t\tLambda threshold:\t%f\n",LambdaThresh);
		printf("\t\tLambda approximately:\t%f\n",lambda);
                printf("\t\tProbability threshold:\t%f\n",PTHRESH);
                printf("\t\tDetection probability approximately:\t%f\n",lambda/(1.0 + lambda));
                printf("\t\tGiven:\tsigma = %f\talpha = %f\txi = %f\n",sigma,alpha,xi);
                printf("\t\tMeasured:\tsigma = %f\talpha = %f\txi = %f\n\n",Dsigma,Dalpha,Dxi);
        } else {
                printf("\t\tNO DETECTION\n");
		printf("\t\tGiven:\tsigma = %f\talpha = %f\txi = %f\n",sigma,alpha,xi);
                printf("\t\tLambda threshold:\t%f\n",LambdaThresh);
		printf("\t\tLambda approximately:\t%f\n",lambda);
		printf("\t\tProbability threshold:\t%f\n",PTHRESH);
                printf("\t\tDetection probability approximately:\t%f\n\n",lambda/(1.0 + lambda));
        }
        printf("\n\t----------------------------------------------------------------------------------------------\n");

}
