/*  -- translated by f2c (version 20200916).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/*     Last change:  S    14 Apr 2007   10:24 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: combinatorial                                      CCC */
/* CC                                                                  CCC */
/* CC   Calculation of binomial coefficients (combinatorial number)    CCC */
/* CC   based on the relation:                                         CCC */
/* CC                                                                  CCC */
/* CC                  (n+1; p+1) = (n; p) + (n; p+1)                  CCC */
/* CC                                                                  CCC */
/* CC   The result is stored in the array Cnp such that:               CCC */
/* CC                                                                  CCC */
/* CC                   Cnp[n*(n+1)/2 + p] = (n; p)                    CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int combinatorial_(integer *nu, doublereal *cnp)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer n, np;

    i__1 = *nu;
    for (n = 0; n <= i__1; ++n) {
	cnp[n * (n + 1) / 2] = 1.;
	cnp[n * (n + 1) / 2 + n] = 1.;
	i__2 = n - 1;
	for (np = 1; np <= i__2; ++np) {
	    cnp[n * (n + 1) / 2 + np] = cnp[n * (n - 1) / 2 + np - 1] + cnp[n 
		    * (n - 1) / 2 + np];
	}
    }
    return 0;
} /* combinatorial_ */

