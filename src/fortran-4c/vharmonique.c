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

static doublereal c_b2 = -1.;

/*     Last change:  H    26 Jul 2009    0:47 am */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: vharmonique1                                       CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the Ylmv(thetav)                     CCC */
/* CC   appearing in the outer integral over s for the                 CCC */
/* CC   Fourier transform approach of three center:                    CCC */
/* CC   - nuclear attraction integrals                                 CCC */
/* CC   - Coulomb integrals                                            CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - lmax : the max value of lambda appearing in the integral  CCC */
/* CC      - nrac : order of the Gauss-Legendre quadrature             CCC */
/* CC      - Thetav[] : array of thetav for each quadrature point      CCC */
/* CC                                                                  CCC */
/* CC   ouput:                                                         CCC */
/* CC      - Ylmv[] : array of Ylmv(thetav)                            CCC */
/* CC                Ylmv(nrac*(l*(l+1)+m) + i-1) == Ylm(Thetav[i])    CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int vharmonique1_(integer *lmax, integer *nrac, doublereal *
	thetav, doublereal *ylmv)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal eps = 1e-15;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sqrt(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    extern /* Subroutine */ int legendre_(doublereal *, integer *, doublereal 
	    *);
    static integer i__, l, m;
    static doublereal dcosthetav, plm[6904], poch;
    static integer nylm1, nylm2;

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
    /* Parameter adjustments */
    --thetav;

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
    i__1 = *nrac;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* .... Associated Legendre polynomials */
	dcosthetav = cos(thetav[i__]);
	if (abs(dcosthetav) < eps) {
	    dcosthetav = 0.;
	}
	legendre_(&dcosthetav, lmax, plm);
	i__2 = *lmax;
	for (l = 0; l <= i__2; ++l) {
	    poch = 1.;
	    i__3 = l;
	    for (m = 0; m <= i__3; ++m) {
		nylm1 = *nrac * (l * (l + 1) + m) + i__ - 1;
		nylm2 = *nrac * (l * (l + 1) - m) + i__ - 1;
		ylmv[nylm1] = pow_di(&c_b2, &m) * plm[l * (l + 1) / 2 + m] * 
			sqrt((doublereal) ((l << 1) + 1) / (pi * 4.)) / poch;
		ylmv[nylm2] = pow_di(&c_b2, &m) * ylmv[nylm1];
		poch *= sqrt((doublereal) ((l + m + 1) * (l - m)));
	    }
	}
    }
    return 0;
} /* vharmonique1_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the Ylmv(thetav)                     CCC */
/* CC   appearing in the outer integral over s for the                 CCC */
/* CC   Fourier transform approach of four center Coulomb integral     CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - lmax : the max value of lambda appearing in the integral  CCC */
/* CC      - nrac : order of the Gauss-Legendre quadrature             CCC */
/* CC      - Thetav[] : array of thetav for each quadrature point      CCC */
/* CC                                                                  CCC */
/* CC   ouput:                                                         CCC */
/* CC      - Ylmv[] : array of Ylmv(thetav)                            CCC */
/* CC          Ylmv(nrac*nrac*(l*(l+1)+m) + nst-1) == Ylm(Thetav[nst]) CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int vharmonique2_(integer *lmax, integer *nrac, doublereal *
	thetav, doublereal *ylmv)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal eps = 1e-15;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double cos(doublereal), sqrt(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    extern /* Subroutine */ int legendre_(doublereal *, integer *, doublereal 
	    *);
    static integer l, m, ns, nt;
    static doublereal dcosthetav, plm[6904];
    static integer nst;
    static doublereal poch;
    static integer nylm1, nylm2;

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
    /* Parameter adjustments */
    --thetav;

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
    i__1 = *nrac;
    for (ns = 1; ns <= i__1; ++ns) {
	i__2 = *nrac;
	for (nt = 1; nt <= i__2; ++nt) {
	    nst = *nrac * (ns - 1) + nt;
/* .... Associated Legendre polynomials */
	    dcosthetav = cos(thetav[nst]);
	    if (abs(dcosthetav) < eps) {
		dcosthetav = 0.;
	    }
	    legendre_(&dcosthetav, lmax, plm);
	    i__3 = *lmax;
	    for (l = 0; l <= i__3; ++l) {
		poch = 1.f;
		i__4 = l;
		for (m = 0; m <= i__4; ++m) {
/*                  nylm1 = nrac*(l*(l+1) + m) + nst - 1 */
/*                  nylm2 = nrac*(l*(l+1) - m) + nst - 1 */
		    nylm1 = *nrac * *nrac * (l * (l + 1) + m) + nst - 1;
		    nylm2 = *nrac * *nrac * (l * (l + 1) - m) + nst - 1;
		    ylmv[nylm1] = pow_di(&c_b2, &m) * plm[l * (l + 1) / 2 + m]
			     * sqrt((doublereal) ((l << 1) + 1) / (pi * 4.)) /
			     poch;
		    ylmv[nylm2] = pow_di(&c_b2, &m) * ylmv[nylm1];
		    poch *= sqrt((doublereal) ((l + m + 1) * (l - m)));
		}
	    }
	}
    }
    return 0;
} /* vharmonique2_ */

