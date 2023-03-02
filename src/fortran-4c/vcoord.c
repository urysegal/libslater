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

/*     Last change:  S    22 Apr 2007    4:25 am */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: vcoord1                                            CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the coordinates of the vector V      CCC */
/* CC   appearing in the outer integral over s for the                 CCC */
/* CC   Fourier transform approach of three center:                    CCC */
/* CC   - nuclear attraction integrals                                 CCC */
/* CC   - Coulomb integrals                                            CCC */
/* CC                                                                  CCC */
/* CC        \vec{V} = (1-s) \vec{ab} - \vec{ac}                       CCC */
/* CC   Input :                                                        CCC */
/* CC      - ab, thetab, phiab : the coordinates of \vec{ab}           CCC */
/* CC      - ac, thetac, phiac : the coordinates of \vec{ac}           CCC */
/* CC      - nrac : order of the Gauss-Legendre quadrature             CCC */
/* CC      - Xbar : the array containing the Gauss-Legendre abscissas  CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC      - V[], Thetav, Phiv[] : array of coordinates                CCC */
/* CC        of \vec{v} for each point of the quadrature.              CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int vcoord1_(doublereal *ab, doublereal *thetab, doublereal *
	phiab, doublereal *ac, doublereal *thetac, doublereal *phiac, integer 
	*nrac, doublereal *xbar, doublereal *v, doublereal *thetav, 
	doublereal *phiv)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal eps = 1e-15;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal), acos(
	    doublereal), atan(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s, vx, vy, vz, sgnx, vnorm;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*  Programmer : Hassan Safouhi                  sept  25 2002 */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* .... pi */
    /* Parameter adjustments */
    --phiv;
    --thetav;
    --v;
    --xbar;

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
	s = (xbar[i__] + 1.) * .5;
/* .... evaluation of the spherical coordinates of */
/* .... \vec{v} = (1-s) \vec{ab} - \vec{ac} */
	vx = (1. - s) * *ab * sin(*thetab) * cos(*phiab) - *ac * sin(*thetac) 
		* cos(*phiac);
	if (abs(vx) < eps) {
	    vx = 0.;
	}
	vy = (1. - s) * *ab * sin(*thetab) * sin(*phiab) - *ac * sin(*thetac) 
		* sin(*phiac);
	if (abs(vy) < eps) {
	    vy = 0.;
	}
	vz = (1. - s) * *ab * cos(*thetab) - *ac * cos(*thetac);
	if (abs(vz) < eps) {
	    vz = 0.;
	}
	vnorm = sqrt(vx * vx + vy * vy + vz * vz);
	v[i__] = vnorm;
	if (abs(vnorm) < eps) {
	    v[i__] = 0.;
	}
	thetav[i__] = 0.;
	phiv[i__] = 0.;
	if (vnorm != 0.) {
	    thetav[i__] = acos(vz / vnorm);
	    if ((d__1 = 1 - cos(thetav[i__]), abs(d__1)) < eps) {
		thetav[i__] = 0.;
	    }
	    if ((d__1 = cos(thetav[i__]) + 1, abs(d__1)) < eps) {
		thetav[i__] = pi;
	    }
	    if (vx == 0.) {
		if (vy == 0.) {
		    phiv[i__] = 0.;
		} else {
		    phiv[i__] = abs(vy) / vy * pi / 2.;
		}
	    } else {
		sgnx = abs(vx) / vx;
		phiv[i__] = (d__1 = sgnx - 1., abs(d__1)) * pi / 2. + atan(vy 
			/ vx);
	    }
	}
    }
    return 0;
} /* vcoord1_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the coordinates of the vector V      CCC */
/* CC   appearing in the outer integral on s & t for the               CCC */
/* CC   Fourier transform approach of the three center exchange        CCC */
/* CC   integral                                                       CCC */
/* CC                                                                  CCC */
/* CC   \vec{v} = (1-s) \vec{ab} - t \vec{ac}                          CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - ab, thetab, phiab : the coordinates of \vec{ab}           CCC */
/* CC      - ac, thetac, phiac : the coordinates of \vec{ac}           CCC */
/* CC      - nrac : order of the Gauss-Legendre quadrature             CCC */
/* CC      - Xbar : the array containing the Gauss-Legendre abscissas  CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC      - V[], Thetav, Phiv[] : array of coordinates                CCC */
/* CC        of \vec{v} for each point of the quadrature.              CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int vcoord2_(doublereal *ab, doublereal *thetab, doublereal *
	phiab, doublereal *ac, doublereal *thetac, doublereal *phiac, integer 
	*nrac, doublereal *xbar, doublereal *v, doublereal *thetav, 
	doublereal *phiv)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal eps = 1e-15;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal), acos(
	    doublereal), atan(doublereal);

    /* Local variables */
    static doublereal s, t;
    static integer ns, nt;
    static doublereal vx, vy, vz;
    static integer nst;
    static doublereal sgnx, vnorm;

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
    --phiv;
    --thetav;
    --v;
    --xbar;

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
	    if (ns <= *nrac) {
		s = (xbar[ns] + 1.) / 20.;
	    } else if (ns <= *nrac + *nrac) {
		s = (xbar[ns] * 8. + 10.) / 20.;
	    } else {
		s = (xbar[ns] + 19.) / 20.;
	    }
	    if (nt <= *nrac) {
		t = (xbar[nt] + 1.) / 20.;
	    } else if (nt <= *nrac + *nrac) {
		t = (xbar[nt] * 8. + 10.) / 20.;
	    } else {
		t = (xbar[nt] + 19.) / 20.;
	    }
/* .... evaluation of the spherical coordinates of */
/* .... \vec{v} = (1-s) \vec{ab} - t \vec{ac} */
	    vx = (1. - s) * *ab * sin(*thetab) * cos(*phiab) - t * *ac * sin(*
		    thetac) * cos(*phiac);
/*         if(dabs(vx) .lt. Eps) vx = 0.0d0 */
	    vy = (1. - s) * *ab * sin(*thetab) * sin(*phiab) - t * *ac * sin(*
		    thetac) * sin(*phiac);
/*         if(dabs(vy) .lt. Eps) vy = 0.0d0 */
	    vz = (1. - s) * *ab * cos(*thetab) - t * *ac * cos(*thetac);
/*         if(dabs(vz) .lt. Eps) vz = 0.0d0 */
	    vnorm = sqrt(vx * vx + vy * vy + vz * vz);
	    nst = *nrac * (ns - 1) + nt;
	    v[nst] = vnorm;
/*         if(dabs(vnorm) .lt. Eps) V(i) = 0.0d0 */
	    thetav[nst] = 0.;
	    phiv[nst] = 0.;
	    if (vnorm != 0.) {
		thetav[nst] = acos(vz / vnorm);
		if ((d__1 = 1 - cos(thetav[nst]), abs(d__1)) < eps) {
		    thetav[nst] = 0.;
		}
		if ((d__1 = cos(thetav[nst]) + 1, abs(d__1)) < eps) {
		    thetav[nst] = pi;
		}
		if (vx == 0.) {
		    if (vy == 0.) {
			phiv[nst] = 0.;
		    } else {
			phiv[nst] = abs(vy) / vy * pi / 2.;
		    }
		} else {
		    sgnx = abs(vx) / vx;
		    phiv[nst] = (d__1 = sgnx - 1., abs(d__1)) * pi / 2. + 
			    atan(vy / vx);
		}
	    }
/* L5: */
	}
    }
    return 0;
} /* vcoord2_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the coordinates of the vector V      CCC */
/* CC   appearing in the outer integral on s & t for the               CCC */
/* CC   Fourier transform approach of the four center Coulomb          CCC */
/* CC   integral.                                                      CCC */
/* CC                                                                  CCC */
/* CC   \vec{v} = (1-s) \vec{ab} + (1-t) \vec{cd} - \vec{ad}           CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - ac, thetac, phiac : the coordinates of \vec{ac}           CCC */
/* CC      - cd, thetcd, phicd : the coordinates of \vec{cd}           CCC */
/* CC      - ad, thetad, phiad : the coordinates of \vec{ad}           CCC */
/* CC      - nrac : order of the Gauss-Legendre quadrature             CCC */
/* CC      - Xbar : the array containing the Gauss-Legendre abscissas  CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC      - V[], Thetav, Phiv[] : array of coordinates                CCC */
/* CC        of \vec{v} for each point of the quadrature.              CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int vcoord3_(doublereal *ab, doublereal *thetab, doublereal *
	phiab, doublereal *cd, doublereal *thetcd, doublereal *phicd, 
	doublereal *ad, doublereal *thetad, doublereal *phiad, integer *nrac, 
	doublereal *xbar, doublereal *v, doublereal *thetav, doublereal *phiv)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal eps = 1e-15;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal), acos(
	    doublereal), atan(doublereal);

    /* Local variables */
    static doublereal s, t;
    static integer ns, nt;
    static doublereal vx, vy, vz;
    static integer nst;
    static doublereal sgnx, vnorm;

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
    --phiv;
    --thetav;
    --v;
    --xbar;

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
	    s = (xbar[ns] + 1.) * .5;
	    t = (xbar[nt] + 1.) * .5;
/* .... evaluation of the spherical coordinates of */
/* .... \vec{v} = (1-s) \vec{ab} + (1-t) \vec{cd} - \vec{ad} */
	    vx = (1. - s) * *ab * sin(*thetab) * cos(*phiab) + (1. - t) * *cd 
		    * sin(*thetcd) * cos(*phicd) - *ad * sin(*thetad) * cos(*
		    phiad);
/*         if(dabs(vx) .lt. Eps) vx = 0.0d0 */
	    vy = (1. - s) * *ab * sin(*thetab) * sin(*phiab) + (1. - t) * *cd 
		    * sin(*thetcd) * sin(*phicd) - *ad * sin(*thetad) * sin(*
		    phiad);
/*         if(dabs(vy) .lt. Eps) vy = 0.0d0 */
	    vz = (1. - s) * *ab * cos(*thetab) + (1. - t) * *cd * cos(*thetcd)
		     - *ad * cos(*thetad);
/*         if(dabs(vz) .lt. Eps) vz = 0.0d0 */
	    vnorm = sqrt(vx * vx + vy * vy + vz * vz);
	    nst = *nrac * (ns - 1) + nt;
	    v[nst] = vnorm;
/*         if(dabs(vnorm) .lt. Eps) V(i) = 0.0d0 */
	    thetav[nst] = 0.;
	    phiv[nst] = 0.;
	    if (vnorm != 0.) {
		thetav[nst] = acos(vz / vnorm);
		if ((d__1 = 1 - cos(thetav[nst]), abs(d__1)) < eps) {
		    thetav[nst] = 0.;
		}
		if ((d__1 = cos(thetav[nst]) + 1, abs(d__1)) < eps) {
		    thetav[nst] = pi;
		}
		if (vx == 0.) {
		    if (vy == 0.) {
			phiv[nst] = 0.;
		    } else {
			phiv[nst] = abs(vy) / vy * pi / 2.;
		    }
		} else {
		    sgnx = abs(vx) / vx;
		    phiv[nst] = (d__1 = sgnx - 1., abs(d__1)) * pi / 2. + 
			    atan(vy / vx);
		}
	    }
/* L5: */
	}
    }
    return 0;
} /* vcoord3_ */

