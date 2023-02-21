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

/*     Last change:  S    14 Apr 2007   10:38 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: power                                              CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates all powers of a number x up to the   CCC */
/* CC   nth power, and stores the result in the array Xpn such that    CCC */
/* CC   Xpn[i] = x**i                                                  CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int power_(integer *n, doublereal *x, doublereal *xpn)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    xpn[0] = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xpn[i__] = *x * xpn[i__ - 1];
    }
    return 0;
} /* power_ */

