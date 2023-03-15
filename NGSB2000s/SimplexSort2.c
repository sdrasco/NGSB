/*
 * File: SimplexSort.c  by: Steve Drasco 
 *
 * Sorts a 4-vector (LogLambda) and then sorts the rows of a 
 * 3 x 4 matrix (simplex) in the same manner.  Probably there is a 
 * much better way, but probably doesn't matter.
 *
 */
 

#include "NGSB.h"

int SimplexSort2(double **simplex, double *LogLambda)
{
	double	SimplexTemp[4][3], LogLambdaTemp[4];
	int	order[2][4]; /* change 2 to 1? */
        int     i, j;

	/* store temporary values */
	for(i = 0; i < 4; i++) {
		order[0][i] = order[1][i] = i;
		LogLambdaTemp[i] = LogLambda[i];
		for(j = 0; j < 3; j++) SimplexTemp[i][j] = simplex[i][j];
	}

	/* sort the five-vector */
	for(i = 0; i < 3; i++) {
		for(j = i; j < 4; j++) {
			if(LogLambda[j] < LogLambda[i]) {
				LogLambda[i] = LogLambdaTemp[j];
				LogLambda[j] = LogLambdaTemp[i];
				LogLambdaTemp[i] = LogLambda[i];
				LogLambdaTemp[j] = LogLambda[j];
				order[0][i] = order[1][j];
				order[0][j] = order[1][i];
				order[1][i] = order[0][i];
				order[1][j] = order[0][j];
			}
		}
	}
	
	/* sort the rows of the matrix */
	for(i = 0; i < 4; i++) for(j = 0; j < 3; j++) {
			simplex[i][j] = SimplexTemp[order[0][i]][j];
	}

	/* normal exit */	
	return 0;
}
