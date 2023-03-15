/* 
 *
 * File: SeedRand.c        by: Steve Drasco      
 *                                                           
 * Seeds the external varriable for ran2.
 *
 */

#include "NGSB.h"

int SeedRand()
{
	extern long idum;

	idum = (long) time(NULL);	
	if(idum >= 0) idum *= -1.0;

}
