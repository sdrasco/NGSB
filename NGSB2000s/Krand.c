/* 
 *
 * File: Krand.c        by: Steve Drasco      
 *                                                           
 * "an even quicker random number generator" - NR            
 *
 */

#include "NGSB.h"

double Krand()
{
	extern unsigned int idum;
	unsigned int itemp;
        static unsigned int jflone=0x3F800000;
        static unsigned int jflmsk=0x007fffff;

	/* make new random integer */
        idum = 1664525 * idum + 1013904223;

	/* turn it into a double between 0 and 1 */
        itemp = jflone | (jflmsk & idum );
        return (double) (*(float *)&itemp)-1.0;
}
