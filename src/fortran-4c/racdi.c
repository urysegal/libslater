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

/*     Last change:  H     4 Jul 2010    2:08 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC   Programmer: Hassan Safouhi                                     CCC */
/* CC   Subroutine: racdi                                              CCC */
/* CC                                                                  CCC */
/* CC   This subroutine returns the position of the root of the        CCC */
/* CC   spherical Bessel function of order n, between a and b using    CCC */
/* CC   the bisection method.                                          CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int racdi_(integer *n, doublereal *a, doublereal *b, 
	doublereal *xn, integer *kmax, integer *ier)
{
    /* Initialized data */

    static doublereal eps = 1e-15;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k;
    static doublereal an, bn;
    static integer nm;
    static doublereal dja[1001], sja[1001], djx[1001], sjx[1001];
    extern /* Subroutine */ int sphj_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*  Programmer : Hassan Safouhi                  sept  25 2002 */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* .... pi */
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
    *ier = 0;
    an = *a;
    bn = *b;
    i__1 = *kmax;
    for (k = 1; k <= i__1; ++k) {
	*xn = (an + bn) / 2.;
	sphj_(n, &an, &nm, sja, dja);
	sphj_(n, xn, &nm, sjx, djx);
	if (sja[*n] * sjx[*n] > 0.) {
	    an = *xn;
	    goto L50;
	} else if (sja[*n] * sjx[*n] < eps) {
	    bn = *xn;
	    goto L50;
	} else if (sja[*n] * sjx[*n] == 0.) {
	    if ((d__1 = sjx[*n], abs(d__1)) < 0.) {
		return 0;
	    }
	} else {
	    *xn = an;
	    return 0;
	}
L50:
	if ((d__1 = sjx[*n], abs(d__1)) < eps) {
	    return 0;
	}
/* L5: */
    }
    *ier = 1;
    return 0;
} /* racdi_ */

