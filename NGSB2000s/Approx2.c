/* File: Approx1LogLikelihood.c by: Steve Drasco 
*                                                            
* Computes sums necessary for Xi expanded log of the likelihood function
* (for a pair of detectors)
*
* NOTE: we assume unit varriance detector noise.
* 
*/

#include "NGSB.h"

double Approx2(int N, double alpha, double *h1, double *h2, double *Sum1, double *Sum2, double *Sum3, double *Sum4)
{
	double	s1=0.0, s2=0.0, s3=0.0, s4=0.0, A, B, TwoB, ThreeB, FourB;
        int     i;

        /* check data pointer */
        if(h1 == NULL || h2==NULL) {
                printf("Data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* constants */
	A = 1.0 / sqrt( 1.0 + 2.0*alpha );
	B = 1.0 /  (   4.0 + (2.0/alpha)  );
        TwoB = 2.0*B;
        ThreeB = 3.0*B;
        FourB = 4.0*B;

	/* do the length N sums */
        for(i=0; i < N; i++) {
                s1 += exp(     B * (h1[i] + h2[i]) * (h1[i] + h2[i]) );
                s2 += exp(  TwoB * (h1[i] + h2[i]) * (h1[i] + h2[i]) );
                s3 += exp(ThreeB * (h1[i] + h2[i]) * (h1[i] + h2[i]) );
                s4 += exp( FourB * (h1[i] + h2[i]) * (h1[i] + h2[i]) );
        }

	/* compute sums over powers of v_j */
	*Sum1 = A * s1 
              - ((double) N);
	*Sum2 = A * A * s2
              - 2.0 * A * s1
              + ((double) N);
	*Sum3 = A * A * A * s3 
              - 3.0 * A * A *s2
              + 3.0 * A * s1
              - ((double) N);
	*Sum4 = A * A * A * A * s4
              - 4.0 * A * A * A * s3
              + 6.0 * A * A * s2
              - 4.0 * A * s1
              + ((double) N);

	/* exit */
	return 0;
}
