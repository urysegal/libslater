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

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Function: DX                                                   CCC */
/* CC                                                                  CCC */
/* CC   Evaluation of the kth derivative of the function (x)^{n} at    CCC */
/* CC   the point x.                                                   CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - n : power of the function                                 CCC */
/* CC      - x : point at which we evaluate the function               CCC */
/* CC      - k : order of the derivative                               CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal dx_(integer *n, doublereal *x, integer *k)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer j;
    static doublereal terme;

    if (*k > *n) {
	ret_val = 0.;
	return ret_val;
    }
    if (*k == 0) {
	ret_val = pow_di(x, n);
    } else {
	terme = 1.;
	i__1 = *k;
	for (j = 1; j <= i__1; ++j) {
	    terme = (doublereal) (*n - j - 1) * terme;
/* L5: */
	}
	i__1 = *n - *k;
	ret_val = terme * pow_di(x, &i__1);
    }
    return ret_val;
} /* dx_ */

