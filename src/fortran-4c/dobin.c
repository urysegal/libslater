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
/* CC   Function: Dobin                                                CCC */
/* CC                                                                  CCC */
/* CC   Dobin(n,m) = (n)!! / (n-2m)!!                                  CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal dobin_(integer *n, integer *m)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__;

    ret_val = 1.;
    i__1 = *m - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	ret_val *= (doublereal) (*n - (i__ << 1));
    }
    return ret_val;
} /* dobin_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: TDobin                                             CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates Dobin(n,m) = (n)!!/(n-2m)!!          CCC */
/* CC   and stores the result in the array Dob[i], i=0...m-1           CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int tdobin_(integer *n, integer *m, doublereal *dob)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    dob[0] = 1.;
    i__1 = *m - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	dob[i__ + 1] = dob[i__] * (doublereal) (*n - (i__ << 1));
    }
    return 0;
} /* tdobin_ */

