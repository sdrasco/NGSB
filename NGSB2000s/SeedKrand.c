/* 
 *
 * File: SeedKrand.c        by: Steve Drasco      
 *                                                           
 * Seeds the external varriable for Krand.
 *
 */

#include "NGSB.h"

int SeedKrand()
{
	extern unsigned int idum;

	idum = (unsigned int) time(NULL);	

}
