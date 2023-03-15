/*
 * File: SimplexSort.c  by: Steve Drasco 
 *
 * Sorts a five-vector (LogLambda) and then sorts the rows of a 
 * 4 x 5 matrix (simplex) in the same manner.  Probably there is a 
 * much better way, but probably doesn't matter.
 *
 */
 

#include "NGSB.h"

int SimplexSort(double **simplex, double *LogLambda)
{
	double	SimplexTemp[5][4], LogLambdaTemp[5];
	int	order[2][5];
        int     i, j;

	/* store temporary values */
	for(i = 0; i < 5; i++) {
		order[0][i] = order[1][i] = i;
		LogLambdaTemp[i] = LogLambda[i];
		for(j = 0; j < 4; j++) SimplexTemp[i][j] = simplex[i][j];
	}

	/* sort the five-vector */
	for(i = 0; i < 4; i++) {
		for(j = i; j < 5; j++) {
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
	for(i = 0; i < 5; i++) for(j = 0; j < 4; j++) {
			simplex[i][j] = SimplexTemp[order[0][i]][j];
	}

	/* normal exit */	
	return 0;
}
