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

/*     Last change:  H     3 Mar 2009    2:01 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: LU                                                 CCC */
/* CC                                                                  CCC */
/* CC   Resolution of the linear system A.x = b                        CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - A[][] : system matrix                                     CCC */
/* CC      - b[]   : vector containing the constants                   CCC */
/* CC      - x[]   : vector containing the unknowns                    CCC */
/* CC      - m     : dimension of the system                           CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int lu_(doublereal *a, doublereal *b, doublereal *x, integer 
	*m, integer *ier, doublereal *det)
{
    extern /* Subroutine */ int triang_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *), resolv_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/* .... MAX_REP_SS : Maximal number of iterations for the subroutine using the infinite series */
/* .... MAX_REP_UL :   Maximal number of iterations for the subroutine using the Levin's u transform */
/* .... MAX_REP_EP : Maximal number of iterations for the subroutine using the epsilon algorithm of Wynn */
/* .... MAX_REP_UR :   Maximal number of iterations for the subroutine using the Levin's u recurrence */
/* .... MAX_REP_ER : Maximal number of iterations for the subroutine using the epsilon recurrence */
/* .... MAX_REP_HD :  Maximal number of iterations for the subroutine using H\bar{D} */
/* .... MAX_REP_DB :   Maximal number of iterations for the subroutine using \bar{D} */
/* .... MAX_REP_SD :  Maximal number of iterations for the subroutine using S\bar{D} */
/* .... MAX_REP_WD : Maximal number of iterations for the subroutine using H\bar{D} and the W-algorithm of Sid
i */
/* .... MAX_REP_SC : Maximal number of iterations for the subroutine using Sin-Cos-S\bar{D} */
/* .... MAX_REP_TC : Maximal number of iterations for the subroutine using Sin-Cos-S\bar{D} */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC  Here we define the maxima of all parameters used in the code.   CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* .... MAX_N : maximal value of the principal quantum number n */
/* .... MAX_L : maximal value of the quantum number l */
/* .... For the Four-center Two-electron Coulomb integral */
/* .... MAX_bsj : Maximum order of spherical Bessel functions in zero_bsj (lambda = 30) */
/* .... MAX_Fact : Size of arrays Fact, Dfact */
/* .... MAX_Cnp : Size of arrays Cnp */
/* .... MAX_LEG : Maximum order of the Gauss-Legendre quadrature */
/* .... MAX_RAC : Maximum order of the Gauss-Legendre quadrature */
/* .... MAX_LAG : Maximum order of the Gauss-Laguerre quadrature */
/* .... For the Four-center Two-electron Coulomb integral */
/* .... MAX_DEV : Maximum number of terms in infinite series */
/* .... MAX_DIV : Maximal number of subdivisions of the finite intervals when using Gauss-Legendre */
/* .... MAX_DB : Maximum order for Dbar */
/* .... MAX_HD : Maximum order for HDbar equal to MAX_DB for LU to operate correctly */
/* .... MAX_UL : maximal order of Levin u */
/* .... MAX_EP : maximal order of epsilon algorithm */
/* .... MAX_UR : maximal order of Levin u using recurrence formulae */
/* .... MAX_ER : maximal order of Epsilon using recurrence formulae */
/* .... MAX_SD : Maximum order for SDbar */
/* .... MAX_WD  : Maximum order for GREP-W */
    *ier = 0;
    *det = 1.;
    triang_(a, b, m, ier, det);
    resolv_(a, b, x, m, ier);
    return 0;
} /* lu_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int triang_(doublereal *a, doublereal *b, integer *m, 
	integer *ier, doublereal *det)
{
    /* Initialized data */

    static doublereal tiny = 1e-150;
    static doublereal eps = 1e-15;

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal am, gik, piv;

/* .... MAX_REP_SS : Maximal number of iterations for the subroutine using the infinite series */
/* .... MAX_REP_UL :   Maximal number of iterations for the subroutine using the Levin's u transform */
/* .... MAX_REP_EP : Maximal number of iterations for the subroutine using the epsilon algorithm of Wynn */
/* .... MAX_REP_UR :   Maximal number of iterations for the subroutine using the Levin's u recurrence */
/* .... MAX_REP_ER : Maximal number of iterations for the subroutine using the epsilon recurrence */
/* .... MAX_REP_HD :  Maximal number of iterations for the subroutine using H\bar{D} */
/* .... MAX_REP_DB :   Maximal number of iterations for the subroutine using \bar{D} */
/* .... MAX_REP_SD :  Maximal number of iterations for the subroutine using S\bar{D} */
/* .... MAX_REP_WD : Maximal number of iterations for the subroutine using H\bar{D} and the W-algorithm of Sid
i */
/* .... MAX_REP_SC : Maximal number of iterations for the subroutine using Sin-Cos-S\bar{D} */
/* .... MAX_REP_TC : Maximal number of iterations for the subroutine using Sin-Cos-S\bar{D} */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC  Here we define the maxima of all parameters used in the code.   CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* .... MAX_N : maximal value of the principal quantum number n */
/* .... MAX_L : maximal value of the quantum number l */
/* .... For the Four-center Two-electron Coulomb integral */
/* .... MAX_bsj : Maximum order of spherical Bessel functions in zero_bsj (lambda = 30) */
/* .... MAX_Fact : Size of arrays Fact, Dfact */
/* .... MAX_Cnp : Size of arrays Cnp */
/* .... MAX_LEG : Maximum order of the Gauss-Legendre quadrature */
/* .... MAX_RAC : Maximum order of the Gauss-Legendre quadrature */
/* .... MAX_LAG : Maximum order of the Gauss-Laguerre quadrature */
/* .... For the Four-center Two-electron Coulomb integral */
/* .... MAX_DEV : Maximum number of terms in infinite series */
/* .... MAX_DIV : Maximal number of subdivisions of the finite intervals when using Gauss-Legendre */
/* .... MAX_DB : Maximum order for Dbar */
/* .... MAX_HD : Maximum order for HDbar equal to MAX_DB for LU to operate correctly */
/* .... MAX_UL : maximal order of Levin u */
/* .... MAX_EP : maximal order of epsilon algorithm */
/* .... MAX_UR : maximal order of Levin u using recurrence formulae */
/* .... MAX_ER : maximal order of Epsilon using recurrence formulae */
/* .... MAX_SD : Maximum order for SDbar */
/* .... MAX_WD  : Maximum order for GREP-W */
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
    *det = 1.;
    i__1 = *m - 1;
    for (k = 0; k <= i__1; ++k) {
	am = (d__1 = a[k + k * 51], abs(d__1));
	if (am < eps) {
	    *ier = 1;
	}
	piv = a[k + k * 51];
	*det *= piv;
	i__2 = *m;
	for (i__ = k + 1; i__ <= i__2; ++i__) {
	    gik = a[i__ + k * 51] / piv;
	    i__3 = *m;
	    for (j = k + 1; j <= i__3; ++j) {
		a[i__ + j * 51] -= gik * a[k + j * 51];
		if ((d__1 = a[i__ + j * 51], abs(d__1)) < tiny) {
		    a[i__ + j * 51] = tiny;
		}
/* L35: */
	    }
	    b[i__] -= gik * b[k];
/* L30: */
	}
/* L10: */
    }
    *det *= a[*m + *m * 51];
    return 0;
} /* triang_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int resolv_(doublereal *a, doublereal *b, doublereal *x, 
	integer *m, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static doublereal s;

/* .... MAX_REP_SS : Maximal number of iterations for the subroutine using the infinite series */
/* .... MAX_REP_UL :   Maximal number of iterations for the subroutine using the Levin's u transform */
/* .... MAX_REP_EP : Maximal number of iterations for the subroutine using the epsilon algorithm of Wynn */
/* .... MAX_REP_UR :   Maximal number of iterations for the subroutine using the Levin's u recurrence */
/* .... MAX_REP_ER : Maximal number of iterations for the subroutine using the epsilon recurrence */
/* .... MAX_REP_HD :  Maximal number of iterations for the subroutine using H\bar{D} */
/* .... MAX_REP_DB :   Maximal number of iterations for the subroutine using \bar{D} */
/* .... MAX_REP_SD :  Maximal number of iterations for the subroutine using S\bar{D} */
/* .... MAX_REP_WD : Maximal number of iterations for the subroutine using H\bar{D} and the W-algorithm of Sid
i */
/* .... MAX_REP_SC : Maximal number of iterations for the subroutine using Sin-Cos-S\bar{D} */
/* .... MAX_REP_TC : Maximal number of iterations for the subroutine using Sin-Cos-S\bar{D} */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC  Here we define the maxima of all parameters used in the code.   CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* .... MAX_N : maximal value of the principal quantum number n */
/* .... MAX_L : maximal value of the quantum number l */
/* .... For the Four-center Two-electron Coulomb integral */
/* .... MAX_bsj : Maximum order of spherical Bessel functions in zero_bsj (lambda = 30) */
/* .... MAX_Fact : Size of arrays Fact, Dfact */
/* .... MAX_Cnp : Size of arrays Cnp */
/* .... MAX_LEG : Maximum order of the Gauss-Legendre quadrature */
/* .... MAX_RAC : Maximum order of the Gauss-Legendre quadrature */
/* .... MAX_LAG : Maximum order of the Gauss-Laguerre quadrature */
/* .... For the Four-center Two-electron Coulomb integral */
/* .... MAX_DEV : Maximum number of terms in infinite series */
/* .... MAX_DIV : Maximal number of subdivisions of the finite intervals when using Gauss-Legendre */
/* .... MAX_DB : Maximum order for Dbar */
/* .... MAX_HD : Maximum order for HDbar equal to MAX_DB for LU to operate correctly */
/* .... MAX_UL : maximal order of Levin u */
/* .... MAX_EP : maximal order of epsilon algorithm */
/* .... MAX_UR : maximal order of Levin u using recurrence formulae */
/* .... MAX_ER : maximal order of Epsilon using recurrence formulae */
/* .... MAX_SD : Maximum order for SDbar */
/* .... MAX_WD  : Maximum order for GREP-W */
    *ier = 0;
    x[*m] = b[*m] / a[*m + *m * 51];
    for (i__ = *m - 1; i__ >= 0; --i__) {
	s = b[i__];
	i__1 = *m;
	for (j = i__ + 1; j <= i__1; ++j) {
	    s -= a[i__ + j * 51] * x[j];
/* L10: */
	}
	x[i__] = s / a[i__ + i__ * 51];
/* L5: */
    }
    return 0;
} /* resolv_ */

