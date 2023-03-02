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
/* CC   Subroutine: GaussLeg                                           CCC */
/* CC                                                                  CCC */
/* CC   Purpose : Compute the zeros of Legendre polynomial Pn(x)       CCC */
/* CC             in the interval [-1,1], and the corresponding        CCC */
/* CC             weighting coefficients for Gauss-Legendre            CCC */
/* CC             integration                                          CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - n : Order of the Legendre polynomial                      CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC      - X(n) : Zeros of the Legendre polynomial                   CCC */
/* CC      - W(n) : Corresponding weighting coefficients               CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int gaussleg_(integer *n, doublereal *x, doublereal *w)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal p, q, z__, f0, f1;
    static integer n0;
    static doublereal z0, fd, gd, pd, pf;
    static integer nr;
    static doublereal wp;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*  Programmer : Hassan Safouhi                  sept  25 2002 */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* .... pi */
    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
/* .... Tiny   => the smallest number */
/* .... Huge   => the largest number */
/* .... Eps    => Epsilon used in general */
/* .... EpsV   => limit to consider v as zero */
/* .... EpsLu  => Epsilon used in LU subroutine */
/* .... EpsSS  =>  limit used for the pre-determined accuracy for the sum using infinite series */
/* .... EpsDB  =>  limit used for the pre-determined accuracy for \bar{D} */
/* .... EpsEP  =>  limit used for the pre-determined accuracy for epsilon algorithm */
/* .... EpsUL  =>  limit used for the pre-determined accuracy for Levin's u */
/* .... EpsHD  =>  limit used for the pre-determined accuracy for H\bar{D} */
/* .... EpsER  =>  limit used for the pre-determined accuracy for epsilon RECURSIVE */
/* .... EpsUR  =>  limit used for the pre-determined accuracy for Levin's u RECURSIVE */
/* .... EpsSD  =>  limit used for the pre-determined accuracy for S\bar{D} */
/* .... EpsWD  =>  limit used for the pre-determined accuracy for W */
/* .... EpsSC  =>  limit used for the pre-determined accuracy for SD Sin Cos */
/* .....FOR THE CALCULATION IN bESSEL SIN-COS, USE THE FOLLOWING */
/* ......... data EpsSC, EpsSD, EpsWD / 1.0d-13, 1.0d-14, 1.0d-14/ */
/* .... Cutoff */
/*      data Cutoff/1.0d-12/ */
/* .... RADEG => Convert degrees to radians */
/*      data RADEG /1.74532925199432958d-02/ */
/* .... nombre maximal d'iteration dans le cycle SCF */
/*      data MAX_ITER / 2000/ */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
    }
    n0 = (*n + 1) / 2;
    i__1 = n0;
    for (nr = 1; nr <= i__1; ++nr) {
	z__ = cos(pi * (doublereal) ((nr << 2) - 1) / (doublereal) (*n << 2));
L10:
	z0 = z__;
	p = 1.;
	i__2 = nr - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    p *= z__ - x[i__];
/* L15: */
	}
	f0 = 1.;
	if (nr == n0 && *n != *n / 2 << 1) {
	    z__ = 0.;
	}
	f1 = z__;
	i__2 = *n;
	for (k = 2; k <= i__2; ++k) {
	    pf = (2. - 1. / k) * z__ * f1 - (1. - 1. / k) * f0;
	    pd = k * (f1 - z__ * pf) / (1. - z__ * z__);
	    f0 = f1;
	    f1 = pf;
/* L20: */
	}
	if (z__ == 0.) {
	    goto L40;
	}
	fd = pf / p;
	q = 0.;
	i__2 = nr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wp = 1.;
	    i__3 = nr;
	    for (j = 1; j <= i__3; ++j) {
		if (j != i__) {
		    wp *= z__ - x[j];
		}
/* L30: */
	    }
	    q += wp;
/* L35: */
	}
	gd = (pd - q * fd) / p;
	z__ -= fd / gd;
	if ((d__1 = z__ - z0, abs(d__1)) > abs(z__) * 1e-15) {
	    goto L10;
	}
L40:
	x[nr] = -z__;
	x[*n + 1 - nr] = z__;
	w[nr] = 2. / ((1. - z__ * z__) * pd * pd);
	w[*n + 1 - nr] = w[nr];
/* L45: */
    }
    return 0;
} /* gaussleg_ */

