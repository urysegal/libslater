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

/* Table of constant values */

static integer c__200 = 200;
static integer c__15 = 15;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: SPHJ                                               CCC */
/* CC                                                                  CCC */
/* CC   Purpose: Compute spherical Bessel functions j_{n}(x) and their CCC */
/* CC            derivatives                                           CCC */
/* CC   Input :  x --- Argument of jn(x)                               CCC */
/* CC            n --- Order of jn(x)  ( n = 0,1,... )                 CCC */
/* CC   Output:  SJ(n) --- j_n(x)                                      CCC */
/* CC            DJ(n) --- j_n'(x)                                     CCC */
/* CC            nm --- Highest order computed                         CCC */
/* CC   Routines called:                                               CCC */
/* CC            msta1 and msta2 for computing the starting            CCC */
/* CC            point for backward recurrence                         CCC */
/* CC                                                                  CCC */
/* CC   Note :                                                         CCC */
/* CC            If the arguument x is greater then the order n, then  CCC */
/* CC            we may proceed with the recursion in the forward      CCC */
/* CC            direction. Otherwise the recusion relations are       CCC */
/* CC            unstable in the forward, so instead we proceed in the CCC */
/* CC            bacwards direction.                                   CCC */
/* CC                                                                  CCC */
/* CC   Acknowledgement :                                              CCC */
/* CC            The majority of the code in this file is copywrited   CCC */
/* CC            and comes from the accompaniement to the book         CCC */
/* CC            "Computation of Special Functions" by Shanjie Zhang   CCC */
/* CC            and Jianming Jin.  The original code can be found on  CCC */
/* CC            the World Wide Web at :                               CCC */
/* CC            http://jin.ece.uiuc.edu/routines/routines.html        CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int sphj_(integer *n, doublereal *x, integer *nm, doublereal 
	*sj, doublereal *dj)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal f;
    static integer i__, k, m;
    static doublereal f0, f1, sa, sb, cs;
    extern integer msta1_(doublereal *, integer *), msta2_(doublereal *, 
	    integer *, integer *);

    *nm = *n;
/* .... If the argument is very small, we will give approximate values */
/* .... in order to avoid overflow */
    if (abs(*x) < 1e-100) {
	i__1 = *n;
	for (k = 0; k <= i__1; ++k) {
	    sj[k] = 0.;
	    dj[k] = 0.;
	}
	sj[0] = 1.;
	dj[1] = .3333333333333333;
	return 0;
    }
/* .... Order zero spherical Bessel function: j_{0}(x) = sin(x)/x */
    sj[0] = sin(*x) / *x;
/* .... This test checks if the argument is larger than the order. */
/* .... If so, then we proceed with the recusion in the forward direction */
/* .... in order to save CPU time, since in such a case the recursion is */
/* .... not unstable. */
    if (*x > (doublereal) (*n) && *n > 0) {
	sj[1] = (sj[0] - cos(*x)) / *x;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    sj[i__] = (doublereal) ((i__ << 1) - 1) * sj[i__ - 1] / *x - sj[
		    i__ - 2];
	}
	return 0;
    } else {
	if (*n == 1) {
	    sj[1] = (sj[0] - cos(*x)) / *x;
	    return 0;
	}
	if (*n >= 2) {
	    sj[1] = (sj[0] - cos(*x)) / *x;
	    sa = sj[0];
	    sb = sj[1];
	    m = msta1_(x, &c__200);
	    if (m < *n) {
		*nm = m;
	    } else {
		m = msta2_(x, n, &c__15);
	    }
	    f0 = 0.;
	    f1 = 1e-100;
	    for (k = m; k >= 0; --k) {
		f = (k * 2. + 3.) * f1 / *x - f0;
		if (k <= *nm) {
		    sj[k] = f;
		}
		f0 = f1;
/* L15: */
		f1 = f;
	    }
	    if (abs(sa) > abs(sb)) {
		cs = sa / f;
	    }
	    if (abs(sa) <= abs(sb)) {
		cs = sb / f0;
	    }
	    i__1 = *nm;
	    for (k = 0; k <= i__1; ++k) {
/* L20: */
		sj[k] = cs * sj[k];
	    }
	}
	dj[0] = (cos(*x) - sin(*x) / *x) / *x;
    }
    i__1 = *nm;
    for (k = 1; k <= i__1; ++k) {
/* L25: */
	dj[k] = sj[k - 1] - (k + 1.) * sj[k] / *x;
    }
    return 0;
} /* sphj_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC     Purpose: Determine the starting point for backward           CCC */
/* CC                recurrence such that the magnitude of             CCC */
/* CC              Jn(x) at that point is about 10^(-mp)               CCC */
/* CC     Input :  x     --- Argument of Jn(x)                         CCC */
/* CC              mp    --- Value of magnitude                        CCC */
/* CC     Output:  msta1 --- Starting point                            CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
integer msta1_(doublereal *x, integer *mp)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static doublereal f, a0, f0, f1;
    static integer n0, n1, nn, it;
    extern doublereal envj_(integer *, doublereal *);

    a0 = abs(*x);
    n0 = (integer) (a0 * 1.1f) + 1;
    f0 = envj_(&n0, &a0) - *mp;
    n1 = n0 + 5;
    f1 = envj_(&n1, &a0) - *mp;
    for (it = 1; it <= 20; ++it) {
	nn = (integer) (n1 - (n1 - n0) / (1. - f0 / f1));
	f = envj_(&nn, &a0) - *mp;
	if ((i__1 = nn - n1, abs(i__1)) < 1) {
	    goto L20;
	}
	n0 = n1;
	f0 = f1;
	n1 = nn;
/* L10: */
	f1 = f;
    }
L20:
    ret_val = nn;
    return ret_val;
} /* msta1_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC     Purpose: Determine the starting point for backward           CCC */
/* CC              recurrence such that all Jn(x) has MP               CCC */
/* CC              significant digits                                  CCC */
/* CC     Input :  x  --- Argument of Jn(x)                            CCC */
/* CC              n  --- Order of Jn(x)                               CCC */
/* CC              mp --- Significant digit                            CCC */
/* CC     Output:  msta2 --- Starting point                            CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
integer msta2_(doublereal *x, integer *n, integer *mp)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static doublereal f, a0, f0, f1;
    static integer n0, n1, nn, it;
    static doublereal obj, ejn, hmp;
    extern doublereal envj_(integer *, doublereal *);

    a0 = abs(*x);
    hmp = *mp * .5;
    ejn = envj_(n, &a0);
    if (ejn <= hmp) {
	obj = (doublereal) (*mp);
	n0 = (integer) (a0 * 1.1f) + 1;
    } else {
	obj = hmp + ejn;
	n0 = *n;
    }
    f0 = envj_(&n0, &a0) - obj;
    n1 = n0 + 5;
    f1 = envj_(&n1, &a0) - obj;
    for (it = 1; it <= 20; ++it) {
	nn = (integer) (n1 - (n1 - n0) / (1. - f0 / f1));
	f = envj_(&nn, &a0) - obj;
	if ((i__1 = nn - n1, abs(i__1)) < 1) {
	    goto L20;
	}
	n0 = n1;
	f0 = f1;
	n1 = nn;
/* L10: */
	f1 = f;
    }
L20:
    ret_val = nn + 10;
    return ret_val;
} /* msta2_ */

doublereal envj_(integer *n, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double d_lg10(doublereal *);

    d__1 = *n * 6.28;
    d__2 = *x * 1.36 / *n;
    ret_val = d_lg10(&d__1) * .5 - *n * d_lg10(&d__2);
    return ret_val;
} /* envj_ */

