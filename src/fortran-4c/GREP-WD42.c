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

/*     Last change:  S    22 Apr 2007    1:29 am */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: GREPW4                                             CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the semi-infinite integrals          CCC */
/* CC   occurring in the analytic expression of four-center Coulomb    CCC */
/* CC   integral over B functions, using the W-algorithm of Sidi and   CCC */
/* CC   the second order differential equation of H\bar{D}.  The       CCC */
/* CC   algorithm is recursive. The accuracy is determined according   CCC */
/* CC   to EpsW.                                                       CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - nk, nx, nu, ng, lambda : parameters of the integral       CCC */
/* CC      - r    : distance between the centers.                      CCC */
/* CC      - v    : norm of the vector v                               CCC */
/* CC      - a, b : parameters for gamma.                              CCC */
/* CC      - X    : zeros of the spherical Beseel function             CCC */
/* CC      - nrac, Xbar[], Xh[] : order, roots and weights of G-Leg    CCC */
/* CC      - nlag, Rlag, Wlag   : order, roots and weights of G-Lag    CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC     - valW  : value of the integral                              CCC */
/* CC     - maxns : order at which the desired precision is attained   CCC */
/* CC     - t     : calculation time                                   CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int grepwd42_(integer *nx, integer *lambda, doublereal *v, 
	integer *nu12, integer *ng12, doublereal *ab, doublereal *a12, 
	doublereal *b12, integer *nu34, integer *ng34, doublereal *cd, 
	doublereal *a34, doublereal *b34, integer *nrac, doublereal *xbar, 
	doublereal *xh, integer *nlag, doublereal *rlag, doublereal *wlag, 
	doublereal *x, doublereal *valw, integer *maxns, doublereal *t)
{
    /* Initialized data */

    static doublereal epswd = 1e-15;
    static doublereal tiny = 1e-150;
    static doublereal epsv = 1e-15;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), pow_di(doublereal *, integer *), exp(doublereal);

    /* Local variables */
    extern /* Subroutine */ int cpu_time__(doublereal *);
    static doublereal d__[102];
    extern doublereal g_(doublereal *, doublereal *, doublereal *);
    static integer i__, j, k;
    static doublereal q[10201]	/* was [101][101] */, s[10201]	/* was [101][
	    101] */, constants[101], ui, xi;
    static integer ns;
    static doublereal qx[10201]	/* was [101][101] */, sx[10201]	/* was [101][
	    101] */;
    extern doublereal fa4_(integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    static integer iii;
    static doublereal z12i, z34i, uij, xij, xx21, err1, err2, err3, calf;
    extern doublereal hatk_(integer *, doublereal *);
    static integer nint;
    static doublereal size, psix, xpnx, xxss;
    extern doublereal hfa4w_(integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    static doublereal psix1, psix2, b_inf__, b_sup__, zpng12, zpng34, start, 
	    psixx;
    extern doublereal hfa4wx_(integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    static doublereal psix1x, psix2x, finish, weight;

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
    --wlag;
    --rlag;
    --xh;
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
    cpu_time__(&start);
    for (iii = 1; iii <= 1; ++iii) {
/* .... v close to zero */
/* .... calculation using Gauss-Laguerre */
/* .... Taylor development of sin(vx)/vx of order 1 was used */
	if (*v < epsv) {
	    if (*lambda == 0) {
		*valw = 0.;
		weight = sqrt(*a12 + *b12) + sqrt(*a34 + *b34);
		i__1 = *nlag;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    xi = rlag[i__] / weight;
		    z12i = g_(a12, b12, &xi);
		    z34i = g_(a34, b34, &xi);
		    xpnx = pow_di(&xi, nx);
		    zpng12 = pow_di(&z12i, ng12);
		    zpng34 = pow_di(&z34i, ng34);
		    d__1 = *ab * z12i;
		    d__2 = *cd * z34i;
		    calf = hatk_(nu12, &d__1) / zpng12 * hatk_(nu34, &d__2) / 
			    zpng34 * xpnx;
		    ui = wlag[i__] * exp(rlag[i__]) / weight * calf;
		    *valw += ui;
		}
		goto L400;
	    } else {
		*valw = 0.;
		goto L400;
	    }
	}
/* .... Evaluation of the first two finite integrals using Gauss_Legendre of order nrac=nrac3 */
/* .....subdivision of the interval for Gauss-Legendre */
/* Computing MIN */
/* Computing MAX */
	i__2 = (integer) (1 / *v);
	i__1 = max(i__2,1);
	nint = min(i__1,100);
	size = (x[1] - x[0]) / (nint * *v);
	constants[0] = 0.;
	b_inf__ = 0.;
	i__1 = nint;
	for (j = 1; j <= i__1; ++j) {
	    b_sup__ = b_inf__ + size;
	    uij = 0.;
	    i__2 = *nrac;
	    for (k = 1; k <= i__2; ++k) {
		xij = ((b_sup__ - b_inf__) * xbar[k] + b_sup__ + b_inf__) * 
			.5;
		uij += size * .5 * xh[k] * fa4_(nx, lambda, v, nu12, ng12, ab,
			 a12, b12, nu34, ng34, cd, a34, b34, &xij);
	    }
	    constants[0] += uij;
	    b_inf__ = b_sup__;
	}
	size = (x[2] - x[1]) / (nint * *v);
	b_inf__ = x[1] / *v;
	constants[1] = constants[0];
	i__1 = nint;
	for (j = 1; j <= i__1; ++j) {
	    b_sup__ = b_inf__ + size;
	    uij = 0.;
	    i__2 = *nrac;
	    for (k = 1; k <= i__2; ++k) {
		xij = ((b_sup__ - b_inf__) * xbar[k] + b_sup__ + b_inf__) * 
			.5;
		uij += size * .5 * xh[k] * fa4_(nx, lambda, v, nu12, ng12, ab,
			 a12, b12, nu34, ng34, cd, a34, b34, &xij);
	    }
	    constants[1] += uij;
	    b_inf__ = b_sup__;
	}
/* .....Computing the recurrence relations using GREP */
/* .....calculate the first three elements of the table */
	d__1 = x[1] / *v;
	psix1 = hfa4w_(nx, lambda, v, nu12, ng12, ab, a12, b12, nu34, ng34, 
		cd, a34, b34, &d__1);
	q[0] = constants[0] / psix1;
	s[0] = 1. / psix1;
	d__1 = x[1] / *v;
	psix1x = hfa4wx_(nx, lambda, v, nu12, ng12, ab, a12, b12, nu34, ng34, 
		cd, a34, b34, &d__1);
	qx[0] = constants[0] / psix1x;
	sx[0] = 1. / psix1x;
/* .....This test is to avoid getting -NaN */
	if ((d__1 = q[0] / qx[0] - s[0] / sx[0], abs(d__1)) < tiny || (d__2 = 
		qx[0] / q[0] - sx[0] / s[0], abs(d__2)) < tiny) {
	    d__[0] = qx[0] / sx[0];
	} else {
	    d__[0] = q[0] / s[0];
	}
	d__1 = x[2] / *v;
	psix2 = hfa4w_(nx, lambda, v, nu12, ng12, ab, a12, b12, nu34, ng34, 
		cd, a34, b34, &d__1);
	q[101] = constants[1] / psix2;
	s[101] = 1. / psix2;
	d__1 = x[1] / *v;
	psix2x = hfa4wx_(nx, lambda, v, nu12, ng12, ab, a12, b12, nu34, ng34, 
		cd, a34, b34, &d__1);
	qx[101] = constants[1] / psix2x;
	sx[101] = 1. / psix2x;
	xx21 = *v / x[2] - *v / x[1];
	q[1] = (q[101] - q[0]) / xx21;
	s[1] = (s[101] - s[0]) / xx21;
	qx[1] = (qx[101] - qx[0]) / xx21;
	sx[1] = (sx[101] - sx[0]) / xx21;
/* .....This test is to avoid getting -NaN */
	if ((d__1 = q[1] / qx[1] - s[1] / sx[1], abs(d__1)) < tiny || (d__2 = 
		qx[1] / q[1] - sx[1] / s[1], abs(d__2)) < tiny) {
	    d__[1] = qx[1] / sx[1];
	} else {
	    d__[1] = q[1] / s[1];
	}
/* .....Loop through each row until the desired accuracy is reached */
	err2 = 1e100;
	err3 = 1e100;
	ns = 1;
	while(err3 > epswd && ns <= 99) {
/* ....    Evaluation of the first two finite integrals using Gauss_Legendre of order nrac=nrac3 */
	    size = (x[ns + 2] - x[ns + 1]) / (nint * *v);
	    b_inf__ = x[ns + 1] / *v;
	    ui = constants[ns];
	    i__1 = nint;
	    for (j = 1; j <= i__1; ++j) {
		b_sup__ = b_inf__ + size;
		uij = 0.;
		i__2 = *nrac;
		for (k = 1; k <= i__2; ++k) {
		    xij = ((b_sup__ - b_inf__) * xbar[k] + b_sup__ + b_inf__) 
			    * .5;
		    uij += size * .5 * xh[k] * fa4_(nx, lambda, v, nu12, ng12,
			     ab, a12, b12, nu34, ng34, cd, a34, b34, &xij);
		}
		ui += uij;
		b_inf__ = b_sup__;
	    }
	    constants[ns + 1] = ui;
	    d__1 = x[ns + 2] / *v;
	    psix = hfa4w_(nx, lambda, v, nu12, ng12, ab, a12, b12, nu34, ng34,
		     cd, a34, b34, &d__1);
	    d__1 = x[ns + 2] / *v;
	    psixx = hfa4wx_(nx, lambda, v, nu12, ng12, ab, a12, b12, nu34, 
		    ng34, cd, a34, b34, &d__1);
	    q[(ns + 1) * 101] = constants[ns + 1] / psix;
	    qx[(ns + 1) * 101] = constants[ns + 1] / psixx;
	    s[(ns + 1) * 101] = 1. / psix;
	    sx[(ns + 1) * 101] = 1. / psixx;
	    for (i__ = ns; i__ >= 0; --i__) {
		xxss = *v / x[ns + 2] - *v / x[i__ + 1];
		q[ns - i__ + i__ * 101 + 1] = (q[ns - i__ - 1 + (i__ + 1) * 
			101 + 1] - q[ns - i__ - 1 + i__ * 101 + 1]) / xxss;
		qx[ns - i__ + i__ * 101 + 1] = (qx[ns - i__ - 1 + (i__ + 1) * 
			101 + 1] - qx[ns - i__ - 1 + i__ * 101 + 1]) / xxss;
		s[ns - i__ + i__ * 101 + 1] = (s[ns - i__ - 1 + (i__ + 1) * 
			101 + 1] - s[ns - i__ - 1 + i__ * 101 + 1]) / xxss;
		sx[ns - i__ + i__ * 101 + 1] = (sx[ns - i__ - 1 + (i__ + 1) * 
			101 + 1] - sx[ns - i__ - 1 + i__ * 101 + 1]) / xxss;
	    }
/* .....This test is to avoid getting -NaN */
	    if ((d__1 = q[ns + 1] / qx[ns + 1] - s[ns + 1] / sx[ns + 1], abs(
		    d__1)) < tiny || (d__2 = qx[ns + 1] / q[ns + 1] - sx[ns + 
		    1] / s[ns + 1], abs(d__2)) < tiny) {
		d__[ns + 1] = qx[ns + 1] / sx[ns + 1];
	    } else {
		d__[ns + 1] = q[ns + 1] / s[ns + 1];
	    }
	    err1 = err2;
	    err2 = err3;
	    err3 = (d__1 = d__[ns + 1] - d__[ns], abs(d__1));
/*         Err3 = dabs((D(ns)-D(ns-1))/D(ns)) */
/* .....Test: if the error increases twice then stop */
	    if (err2 / err1 > 1. && err3 / err2 > 1.) {
		*maxns = ns - 3;
		goto L40;
	    }
	    *maxns = ns;
	    ++ns;
	}
L40:
	*valw = d__[*maxns + 1];
L400:
	;
    }
    cpu_time__(&finish);
    *t = finish - start;
    return 0;
} /* grepwd42_ */

