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

/*     Last change:  S    25 Apr 2007    2:58 am */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: harmonique                                         CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the spherical harmonics at the point CCC */
/* CC   cos(tetha) for all combinations of l, m with l <= lmax, and    CCC */
/* CC   the result is stored in the array Ylm such that :              CCC */
/* CC   -  Ylm(l*(l+1) + m) = Y(l,m)        (m = -l,.,-1,0,1,.,l)      CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int harmonique_(doublereal *theta, integer *l, doublereal *
	ylm)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal fact, dmun, fact1, fact2, dbfact, ctheta, stheta, 
	    sthetai;

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
    ctheta = cos(*theta);
    stheta = sin(*theta);
    ylm[0] = 1. / sqrt(pi * 4.);
    ylm[2] = sqrt(3. / (pi * 4.)) * ctheta;
    ylm[6] = ctheta * sqrt(15.) / 2. * ylm[2] - sqrt(5.) / 2. * ylm[0];
    dbfact = .5;
    sthetai = stheta;
    dmun = -1.;
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ylm[i__ * (i__ + 1) - i__] = sqrt(((i__ << 1) + 1) * dbfact / (pi * 
		4.)) * sthetai;
	ylm[i__ * (i__ + 1) + i__] = dmun * ylm[i__ * (i__ + 1) - i__];
	dbfact = dbfact * (doublereal) ((i__ << 1) + 1) / (doublereal) ((i__ 
		<< 1) + 2);
	sthetai *= stheta;
	dmun *= -1.;
/* L5: */
    }
    i__1 = *l;
    for (i__ = 2; i__ <= i__1; ++i__) {
	fact = ctheta * sqrt((doublereal) (((i__ << 1) + 1) * ((i__ << 1) - 1)
		) / (doublereal) ((i__ << 1) - 1));
/*         fact = ctheta * dsqrt( dble(2*i+1) ) */
	ylm[i__ * (i__ + 1) + i__ - 1] = fact * ylm[i__ * (i__ - 1) + i__ - 1]
		;
	ylm[i__ * (i__ + 1) - i__ + 1] = fact * ylm[i__ * (i__ - 1) - i__ + 1]
		;
/* L10: */
    }
    i__1 = *l;
    for (i__ = 3; i__ <= i__1; ++i__) {
	i__2 = i__ - 2;
	for (j = -i__ + 2; j <= i__2; ++j) {
	    fact1 = sqrt((doublereal) (((i__ << 1) + 1) * ((i__ << 1) - 1)) / 
		    (doublereal) ((i__ + j) * (i__ - j)));
	    fact2 = sqrt((doublereal) (((i__ << 1) + 1) * (i__ + j - 1) * (
		    i__ - j - 1)) / (doublereal) (((i__ << 1) - 3) * (i__ + j)
		     * (i__ - j)));
	    ylm[i__ * (i__ + 1) + j] = ctheta * fact1 * ylm[i__ * (i__ - 1) + 
		    j] - fact2 * ylm[(i__ - 2) * (i__ - 1) + j];
/* L20: */
	}
/* L15: */
    }
    return 0;
} /* harmonique_ */

