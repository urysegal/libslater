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

/*     Last change:  S     6 Apr 2007    2:46 am */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: Cartco                                             CCC */
/* CC                                                                  CCC */
/* CC   This subroutine calculates the spherical polar coordinates of  CCC */
/* CC   \vec{AB} from the cartesian coordinates of \vec{A} and \vec{B}.CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - Xa, Ya, Za : cartesian coordinates of \vec{A}             CCC */
/* CC      - Xb, Yb, Zb : cartesian coordinates of \vec{B}             CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC      - AB, Phiab, Thetab : radial, azimuthal, and polar          CCC */
/* CC                            coordinates of \vec{AB}               CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int cartco_(doublereal *xa, doublereal *ya, doublereal *za, 
	doublereal *xb, doublereal *yb, doublereal *zb, doublereal *ab, 
	doublereal *thetab, doublereal *phiab)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal eps = 1e-15;

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal), cos(doublereal), atan(
	    doublereal);

    /* Local variables */
    static doublereal abx, aby, abz, sgnx;

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
    abx = *xb - *xa;
    aby = *yb - *ya;
    abz = *zb - *za;
    *ab = sqrt(abx * abx + aby * aby + abz * abz);
    if (*ab < eps) {
	*ab = 0.;
    }
    *thetab = 0.;
    *phiab = 0.;
    if (*ab != 0.) {
	*thetab = acos(abz / *ab);
	if ((d__1 = 1. - cos(*thetab), abs(d__1)) < eps) {
	    *thetab = 0.;
	}
	if ((d__1 = cos(*thetab) + 1., abs(d__1)) < eps) {
	    *thetab = pi;
	}
	if (abs(abx) < eps) {
	    abx = 0.;
	}
	if (abx == 0.) {
	    if (abs(aby) < eps) {
		aby = 0.;
	    }
	    if (aby != 0.) {
		*phiab = abs(aby) / aby * pi / 2.;
	    }
	} else {
	    sgnx = abs(abx) / abx;
	    *phiab = (d__1 = sgnx - 1., abs(d__1)) * pi / 2. + atan(aby / abx)
		    ;
	}
    }
    return 0;
} /* cartco_ */

