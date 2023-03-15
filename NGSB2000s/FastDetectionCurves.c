/* 
 * File: FastDetectionCurves.c  (stand alone)     By: Steve Drasco
 *
 * Generates detection curve (alpha detectable vs. N) data. 
 * It is "fast" because it only treats the cross correlation statistic, 
 * and the simple burst-detecting statistic.
 *
 */

#include "NGSB.h"

main(int argc, char *argv[])
{
	double		*h1, *h2, *n1, *n2, *s, *Lambda0, *Lambda1;
	double		alpha0, alpha1, MeanAlpha, xi, distance, val0, val1;
	double		PFAc = 0.1, PFDc = 0.1, Sigma1 = 1.0, Sigma2 = 1.0, ThreshDistance = 0.01;
	int		i, j, *N, NTrials, iterate, MinIndex;
	int		NThresh = 2000, DetectionPoints = 1, MaxIterate = 10, MaxN = 20000;
	char            *GFileName, *BFileName;
	FILE		*GOutFile, *BOutFile;
	time_t          tp;
	struct  tm      *TimeStruct;

	/* check inputs */
	if(argc != 3) {
		printf("I need inputs: xi,  NTrials\n");
                return;
        }

        /* read inputs */
	xi = (double) atof(argv[1]);
	NTrials = (int) atol(argv[2]);

        /* allocate memory */
        h1 = (double *) malloc(MaxN*sizeof(double));
        h2 = (double *) malloc(MaxN*sizeof(double));
        n1 = (double *) malloc(MaxN*sizeof(double));
        n2 = (double *) malloc(MaxN*sizeof(double));
	s = (double *) malloc(MaxN*sizeof(double));
	Lambda0 = (double *) malloc(NTrials*sizeof(double));
	Lambda1 = (double *) malloc(NTrials*sizeof(double));
	GFileName = (char *) malloc(100*sizeof(char));
	BFileName = (char *) malloc(100*sizeof(char));
	N = (int *) malloc(DetectionPoints*sizeof(int));

	/* make output file name strings */
	sprintf(GFileName,"../results/G-DetectionCurve-");
	sprintf(BFileName,"../results/B-DetectionCurve-");
        time(&tp);
        TimeStruct = localtime(&tp);
	strftime(GFileName+28,7,"%d%m%y",TimeStruct);
	strftime(BFileName+28,7,"%d%m%y",TimeStruct);
	sprintf(GFileName+34,"-%f-%d-%d",xi,NTrials,time(NULL));
	sprintf(BFileName+34,"-%f-%d-%d",xi,NTrials,time(NULL));

        /* open output files */
	GOutFile = fopen(GFileName,"w");
	if( GOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",GFileName);
                return;
        }
        BOutFile = fopen(BFileName,"w");
        if( BOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",BFileName);
                return;
        }

	/* print headers to output files */
	fprintf(GOutFile,"%% Gaussian Detection Curve Results\n");
        fprintf(GOutFile,"%% Sigma1      = 1.0\n");
        fprintf(GOutFile,"%% Sigma2      = 1.0\n");
	fprintf(GOutFile,"%% xi          = %e\n",xi);
	fprintf(GOutFile,"%% max data length = %d\n",MaxN);
	fprintf(GOutFile,"%% (PFDc,PFAc) = (%f,%f)\n",PFDc,PFAc);
	fprintf(GOutFile,"%% trials      = %d\n",NTrials);
        fprintf(GOutFile,"%% Output format is style: [AlphaDetectable;  N; DeltaAlpha]\n");
        fprintf(BOutFile,"%% Burst statistic Detection Curve Results\n");
        fprintf(BOutFile,"%% Sigma1      = 1.0\n");
        fprintf(BOutFile,"%% Sigma2      = 1.0\n");
        fprintf(BOutFile,"%% xi          = %e\n",xi);
        fprintf(BOutFile,"%% max data length = %d\n",MaxN);
        fprintf(BOutFile,"%% (PFDc,PFAc) = (%f,%f)\n",PFDc,PFAc);
        fprintf(BOutFile,"%% trials      = %d\n",NTrials);
        fprintf(BOutFile,"%% Output format is style: [AlphaDetectable;  N; DeltaAlpha]\n");
	fflush(NULL);

	/* seed the random number generator */
	SeedRand();

	/* make N array NOTE: CHANGE DetectionPoints to 10 when using this */
	/*
        N[0] = 4000;
        N[1] = 4667;
        N[2] = 5333;
        N[3] = 6000;
        N[4] = 6667;
        N[5] = 7333;
        N[6] = 8000;
        N[7] = 8667;
        N[8] = 9333;
        N[9] = 10000; 
	*/

	/* use only N = 10000 NOTE: CHANGE DetectionPoints to 1 when using this */ 
	N[0] = 10000;

	/* give an initial alpha bracket  (between SNR of 1.0 and 10)*/
	alpha0 = 1.0 / ( xi * sqrt( ((double) N[0]) ) );
	alpha1 = 10.0 / ( xi * sqrt( ((double) N[0]) ) );


	/* loop for Gaussian detection points  */
	for(i = 0; i < DetectionPoints; i++) {

		if(i != 0) {
			alpha0 = 1.0 / ( xi * sqrt( ((double) N[i]) ) );
			alpha1 = MeanAlpha;
		} 

                val1 = val0 = 1.0;
                while( val0*val1 >= 0.0 ) {
                        GaussianDistance(N[i], Sigma1, Sigma2, xi, alpha0, NTrials, Lambda0, Lambda1, h1, h2, n1,
					 n2, s, ThreshDistance/100.0, PFAc, PFDc, &val0);
                        GaussianDistance(N[i], Sigma1, Sigma2, xi, alpha1, NTrials, Lambda0, Lambda1, h1, h2, n1,
					 n2, s, ThreshDistance/100.0, PFAc, PFDc, &val1);

                        if(val0 <= 0.0 && val1 <= 0.0){
				alpha0 = alpha1;
				alpha1 = 10*alpha1;
			} else if(val0 > 0.0 && val1 > 0.0) {
				alpha1 = alpha0;
				alpha0 = alpha0 / 10.0;
			}

                }

		distance = 2.0*ThreshDistance;
		iterate = 0;

		while(fabs(distance) > ThreshDistance   &&   iterate < MaxIterate){

			MeanAlpha = 0.5 * (alpha0 + alpha1);

			GaussianDistance(N[i], Sigma1, Sigma2, xi, MeanAlpha, NTrials, Lambda0, Lambda1, h1, h2, n1,
					n2, s, ThreshDistance/100.0, PFAc, PFDc, &distance);

			if(distance > 0.0) alpha1 = MeanAlpha;
			else alpha0 = MeanAlpha;

                	iterate++;
	
		}

		if(fabs(distance) > ThreshDistance){
			printf("WARNING! Did not find a detectable alpha for Gaussian statistic using N = %d.\n", N[i]);
			printf("Closest value was alpha = %f for distance = %f \n", MeanAlpha, distance);
			printf("ThreshDistance = %f. \n", ThreshDistance);
			printf("We will use estimate: alpha = %f\n",0.5 * (alpha0 + alpha1));
			printf("moving on...\n");
			fprintf(GOutFile,"%f\t%d\t%e %% WARNING! distance = %f\n", 0.5 * (alpha0 + alpha1), N[i], 0.5*fabs(alpha1-alpha0), distance);
			fflush(NULL);
		} else {
			fprintf(GOutFile,"%f\t%d\t%e\n", 0.5 * (alpha0 + alpha1), N[i], 0.5*fabs(alpha1-alpha0) );
			fflush(NULL);
		}

	} 

        /* loop for Burst detection points  */
        for(i = 0; i < DetectionPoints; i++) {

                if(i != 0) {
                        alpha0 = 1.0 / ( xi * sqrt( ((double) N[i]) ) );
                        alpha1 = MeanAlpha;
                }

                val1 = val0 = 1.0;
                while( val0*val1 >= 0.0 ) {
                        BurstDistance(N[i], Sigma1, Sigma2, xi, alpha0, NTrials, Lambda0, Lambda1, h1, n1,
                                         s, ThreshDistance/100.0, PFAc, PFDc, &val0);
                        BurstDistance(N[i], Sigma1, Sigma2, xi, alpha1, NTrials, Lambda0, Lambda1, h1, n1,
                                         s, ThreshDistance/100.0, PFAc, PFDc, &val1);

                        if(val0 <= 0.0 && val1 <= 0.0){
                                alpha0 = alpha1;
                                alpha1 = 10*alpha1;
                        } else if(val0 > 0.0 && val1 > 0.0) {
                                alpha1 = alpha0;
                                alpha0 = alpha0 / 10.0;
                        }
                }

                distance = 2.0*ThreshDistance;
                iterate = 0;

                while(fabs(distance) > ThreshDistance   &&   iterate < MaxIterate){

                        MeanAlpha = 0.5 * (alpha0 + alpha1);

                        BurstDistance(N[i], Sigma1, Sigma2, xi, MeanAlpha, NTrials, Lambda0, Lambda1, h1, n1,
                                        s, ThreshDistance/100.0, PFAc, PFDc, &distance);

                        if(distance > 0.0) alpha1 = MeanAlpha;
                        else alpha0 = MeanAlpha;

                        iterate++;

                }

                if(fabs(distance) > ThreshDistance){
                        printf("WARNING! Did not find a detectable alpha for Gaussian statistic using N = %d.\n", N[i]);
                        printf("Closest value was alpha = %f for distance = %f \n", MeanAlpha, distance);
                        printf("ThreshDistance = %f. \n", ThreshDistance);
                        printf("We will use estimate: alpha = %f\n",0.5 * (alpha0 + alpha1));
                        printf("moving on...\n");
                        fprintf(BOutFile,"%f\t%d\t%e %% WARNING! distance = %f\n", 0.5 * (alpha0 + alpha1), N[i], 0.5*fabs(alpha1-alpha0), distance);
                        fflush(NULL);
                } else {
                        fprintf(BOutFile,"%f\t%d\t%e\n", 0.5 * (alpha0 + alpha1), N[i], 0.5*fabs(alpha1-alpha0) );
                        fflush(NULL);
                }

        }

	/* close the output files */
	fclose(GOutFile);
	fclose(BOutFile);

}
