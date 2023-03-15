/* File: Approx1LogLikelihood.c by: Steve Drasco 27 Nov 2000 
*                                                            
* Computes sums necessary for Xi expanded log of the likelihood function
* 
*/

#include "NGSB.h"

int Approx1(int N, double sigma, double alpha, double *h, double *Sum1, double *Sum2, double *Sum3, double *Sum4)
{
	double	s1=0.0, s2=0.0, s3=0.0, s4=0.0, A, B, TwoB, ThreeB, FourB;
        int     i;

        /* check data pointer */
        if(h == NULL) {
                printf("Data pointer was NULL.  Giving up ...\n");
                return 1;
        }

	/* constants */
	A = 1 / sqrt( 1 + (alpha/sigma) );
	B = 1 /  (   2.0 * ( sigma + ((sigma*sigma) / alpha)  )   );
	TwoB = 2*B;
	ThreeB = 3*B;
	FourB = 4*B;

	/* do the length N sums */
	for(i=0; i < N; i++) {
		s1 += exp(B * h[i] * h[i]);
		s2 += exp(TwoB * h[i] * h[i]);
		s3 += exp(ThreeB * h[i] * h[i]);
		s4 += exp(FourB * h[i] * h[i]);
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
