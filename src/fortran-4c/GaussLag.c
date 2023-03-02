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

/*     Last change:  H     3 Mar 2009    1:28 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: GaussLag                                           CCC */
/* CC                                                                  CCC */
/* CC   Purpose : This routine returns arrays X(1:n) and W(1:n)        CCC */
/* CC             containing the abscissas and weights of the          CCC */
/* CC             n-point Gauss-Laguerre quadrature formula.           CCC */
/* CC             The smallest abscissa is returned in X(1),           CCC */
/* CC             the largest in X(n).                                 CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - n : Order of the Legendre polynomial                      CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC      - X(n) : Zeros of the Laguerre polynomial                   CCC */
/* CC      - W(n) : Corresponding weighting coefficients               CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int gausslag_(integer *n, doublereal *x, doublereal *w)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal z__, p1, p2, p3, z1, ai, pp;

    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == 1) {
	    z__ = 3. / (*n * 2.4 + 1.);
	} else if (i__ == 2) {
	    z__ += 15. / (*n * 2.5 + 1.);
	} else {
	    ai = (doublereal) (i__ - 2);
	    z__ += (ai * 2.55 + 1.) / (ai * 1.9) * (z__ - x[i__ - 2]);
	}
L10:
	p1 = 1.;
	p2 = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    p3 = p2;
	    p2 = p1;
	    p1 = (((doublereal) ((j << 1) - 1) - z__) * p2 - (doublereal) (j 
		    - 1) * p3) / (doublereal) j;
	}
	pp = (doublereal) (*n) * (p1 - p2) / z__;
	z1 = z__;
	z__ = z1 - p1 / pp;
	if ((d__1 = z__ - z1, abs(d__1)) > z__ * 3e-13) {
	    goto L10;
	}
/* L1: */
	x[i__] = z__;
	w[i__] = -1. / (pp * *n * p2);
    }
    return 0;
} /* gausslag_ */

