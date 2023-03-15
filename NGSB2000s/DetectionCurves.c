/* 
 * File: DetectionCurves.c  (stand alone)     By: Steve Drasco
 *
 * Generates detection curve (alpha detectable vs. N) data.  
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
	char            *GFileName, *NGFileName;
	FILE		*GOutFile, *NGOutFile;
	time_t          tp;
	struct  tm      *TimeStruct;

	/* ### DEBUG ### 
        char            *GDFileName, *NGDFileName;
        FILE            *GDOutFile, *NGDOutFile;
	*/

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
	NGFileName = (char *) malloc(100*sizeof(char));

	/* ### DEBUG ### 
        GDFileName = (char *) malloc(100*sizeof(char));
        NGDFileName = (char *) malloc(100*sizeof(char));
	*/

	N = (int *) malloc(DetectionPoints*sizeof(int));

	/* make output file name strings */
	sprintf(GFileName,"../results/G-DetectionCurve-");
	sprintf(NGFileName,"../results/NG-DetectionCurve-");
        time(&tp);
        TimeStruct = localtime(&tp);
	strftime(GFileName+28,7,"%d%m%y",TimeStruct);
	strftime(NGFileName+29,7,"%d%m%y",TimeStruct);
	sprintf(GFileName+34,"-%f-%d-%d",xi,NTrials,time(NULL));
	sprintf(NGFileName+35,"-%f-%d-%d",xi,NTrials,time(NULL));

	/* ### DEBUG ### 
        sprintf(GDFileName,"../results/G-DetectionCurve-");
        sprintf(NGDFileName,"../results/NG-DetectionCurve-");
        strftime(GDFileName+28,7,"%d%m%y",TimeStruct);
        strftime(NGDFileName+29,7,"%d%m%y",TimeStruct);
        sprintf(GDFileName+34,"-%f-%d-%d-DEBUG",xi,NTrials,time(NULL));
        sprintf(NGDFileName+35,"-%f-%d-%d-DEBUG",xi,NTrials,time(NULL));
	*/

        /* open output files */
	GOutFile = fopen(GFileName,"w");
	if( GOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",GFileName);
                return;
        }
        NGOutFile = fopen(NGFileName,"w");
        if( NGOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",NGFileName);
                return;
        }

	/* ### DEBUG ### 
        GDOutFile = fopen(GDFileName,"w");
        if( GDOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",GDFileName);
                return;
        }
        NGDOutFile = fopen(NGDFileName,"w");
        if( NGDOutFile == NULL) {
                printf("Cannot open file: %s\nGiving up...\n",NGDFileName);
                return;
        }
	*/

	/* print headers to output files */
	fprintf(GOutFile,"%% Gaussian Detection Curve Results\n");
        fprintf(GOutFile,"%% Sigma1      = 1.0\n");
        fprintf(GOutFile,"%% Sigma2      = 1.0\n");
	fprintf(GOutFile,"%% xi          = %e\n",xi);
	fprintf(GOutFile,"%% max data length = %d\n",MaxN);
	fprintf(GOutFile,"%% (PFDc,PFAc) = (%f,%f)\n",PFDc,PFAc);
	fprintf(GOutFile,"%% trials      = %d\n",NTrials);
        fprintf(GOutFile,"%% Output format is style: [AlphaDetectable;  N; DeltaAlpha]\n");
        fprintf(NGOutFile,"%% Non-Gaussian Detection Curve Results\n");
        fprintf(NGOutFile,"%% Sigma1      = 1.0\n");
        fprintf(NGOutFile,"%% Sigma2      = 1.0\n");
        fprintf(NGOutFile,"%% xi          = %e\n",xi);
        fprintf(NGOutFile,"%% max data length = %d\n",MaxN);
        fprintf(NGOutFile,"%% (PFDc,PFAc) = (%f,%f)\n",PFDc,PFAc);
        fprintf(NGOutFile,"%% trials      = %d\n",NTrials);
        fprintf(NGOutFile,"%% Output format is style: [AlphaDetectable;  N; DeltaAlpha]\n");
	fflush(NULL);

	/* seed the random number generator */
	SeedRand();

	/* make N array 
        N[0] = 4000;
        N[1] = 4667;
        N[2] = 5333;
        N[3] = 6000;
        N[4] = 6667;
        N[5] = 7333;
        N[6] = 8000;
        N[7] = 8667;
        N[8] = 9333;
        N[9] = 10000; */

	/* make N array top 6 only 
        N[0] = 6667;
        N[1] = 7333;
        N[2] = 8000;
        N[3] = 8667;
        N[4] = 9333;
        N[5] = 10000; */

	/* make N array top 1 only */
	N[0] = 10000;
	

	/* give an initial alpha bracket  (between SNR of 1.0 and 10)*/
	alpha0 = 1.0 / ( xi * sqrt( ((double) N[0]) ) );
	alpha1 = 10.0 / ( xi * sqrt( ((double) N[0]) ) );

	/* ########### Gaussian detection points ############# */
	/* loop for Gaussian detection points  
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

	} */
	/* ########### END Gaussian detection points ############# */
	
        /* initial alpha bracket  for non-Gaussian detection (between SNR of 1.0 and 10) */
        alpha0 = 1.0 / ( xi * sqrt( ((double) N[0]) ) );
        alpha1 = 10.0 / ( xi * sqrt( ((double) N[0]) ) );

        /* loop for non-Gaussian detection points */
        for(i = 0; i < DetectionPoints; i++) {

                /* give new alpha bracket */
                if(i != 0) {
                        alpha0 = 1.0 / ( xi * sqrt( ((double) N[i]) ) );
                        alpha1 = MeanAlpha;
                }

                /* ### DEBUG ### 
                fprintf(NGDOutFile,"checking initial bracket\n");
		fflush(NULL);
		*/

		/* confirm the bracket is good */
                val1 = val0 = 1.0;
                while( val0*val1 >= 0.0 ) {
                        NonGaussianDistance(N[i], Sigma1, Sigma2, xi, alpha0, NTrials, Lambda0, Lambda1, h1, h2, n1,
				 	    n2, s, ThreshDistance/100.0, PFAc, PFDc, &val0);
                        NonGaussianDistance(N[i], Sigma1, Sigma2, xi, alpha1, NTrials, Lambda0, Lambda1, h1, h2, n1,
                                            n2, s, ThreshDistance/100.0, PFAc, PFDc, &val1);

                        /* ### DEBUG ### 
                        fprintf(NGDOutFile,"alpha0 = %e\tval0 = %e\talpha1 = %e\tval1 = %e\n",alpha0,val0,alpha1,val1);
			fflush(NULL);
			*/

                        if(val0 <= 0.0 && val1 <= 0.0){
                                alpha0 = alpha1;
                                alpha1 = 10*alpha1;
                        } else if(val0 > 0.0 && val1 > 0.0) {
                                alpha1 = alpha0;
                                alpha0 = alpha0 / 10.0;
                        }
                }

                /* reset distance and iteration counter */
                distance = 2.0*ThreshDistance;
                iterate = 0;

                /* ### DEBUG ### 
                fprintf(NGDOutFile,"bracket good. start root finding\n");
		fflush(NULL);
		*/

                while(fabs(distance) > ThreshDistance   &&   iterate < MaxIterate){

                        /* find middle of alpha bracket */
                        MeanAlpha = 0.5 * (alpha0 + alpha1);

                        /* compute new distance */
                        NonGaussianDistance(N[i], Sigma1, Sigma2, xi, MeanAlpha, NTrials, Lambda0, Lambda1, h1, h2, n1,
                                            n2, s, ThreshDistance/100.0, PFAc, PFDc, &distance);

                        /* ### DEBUG ### 
                        fprintf(NGDOutFile,"alpha0 = %e\talpha1 = %e\tMeanAlpha = %e\tdistance = %e\n",alpha0,alpha1,MeanAlpha,distance);
			fflush(NULL);
			*/

                        /* update alpha bracket */
                        if(distance > 0.0) alpha1 = MeanAlpha;
                        else alpha0 = MeanAlpha;

                        /* increment iteration counter */
                        iterate++;

                }

                /* if there's a new detection point print it to output file */
                if(fabs(distance) > ThreshDistance){
                        printf("WARNING! Did not find a detectable alpha for non-Gaussian statistic using N = %d.\n", N[i]);
                        printf("Colsest value was alpha = %f for distance = %f \n", MeanAlpha, distance);
                        printf("ThreshDistance = %f.\n", ThreshDistance);
			printf("We will use estimate: alpha = %f\n",0.5 * (alpha0 + alpha1));
                        fprintf(NGOutFile,"%f\t%d\t%e %% WARNING! distance = %f\n", 0.5 * (alpha0 + alpha1), N[i], 0.5*fabs(alpha1-alpha0), distance);
			fflush(NULL);
                } else {
                        fprintf(NGOutFile,"%f\t%d\t%e\n", 0.5 * (alpha0 + alpha1), N[i], 0.5*fabs(alpha1-alpha0) );
			fflush(NULL);
                }

        }


	/* close the output files */
	fclose(GOutFile);
	fclose(NGOutFile);

	/* ### DEBUG ### 
	fclose(GDOutFile);
        fclose(NGDOutFile);
	*/

}
