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

/*     Last change:  H     2 Mar 2009    1:35 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: sdbar4_lrzv                                        CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the semi-infinie integral over x     CCC */
/* CC   found within the integral over s of eq(49) using S\bar{D},     CCC */
/* CC   for the four-center Coulomb integral                           CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - nx     : nx = (l1+l2+l3+l4) - (l'1+l'2+l'3+l'4)           CCC */
/* CC      - lambda : lambda == order of the spherical Bessel function CCC */
/* CC      - v      : |v| = |(1-s)*ab - (1-t)*dc - ad|                 CCC */
/* CC                                                                  CCC */
/* CC      - nu12 : nu12 == n1 + n2 + l1 + l2 - l - j12                CCC */
/* CC      - ng12 : ng12 == 2(n1+l1+n2+l2)-(l'1+l'2+l)+1               CCC */
/* CC      - ab   : position of the function B translated wrt A        CCC */
/* CC      - a12  : a12 = (1-s)*zeta1*zeta1 + s * zeta2*zeta2          CCC */
/* CC      - b12  : b12 = s*(1-s)                                      CCC */
/* CC                                                                  CCC */
/* CC      - nu34 : nu34 == n3 + n4 + l3 + l4 - l' - j34               CCC */
/* CC      - ng34 : ng34 == 2(n3+l3+n4+l4)-(l'3+l'4+l')+1              CCC */
/* CC      - cd   : position of the function C translated wrt D        CCC */
/* CC      - a34  : a34 = (1-t)*zeta4*zeta4 + t * zeta3*zeta3          CCC */
/* CC      - b34  : b34 = t*(1-t)                                      CCC */
/* CC                                                                  CCC */
/* CC      - nsd  : order of S\bar{D}.                                 CCC */
/* CC      - jc   : the j appearing in equation (50)                   CCC */
/* CC      - nrac : order of the Gauss_Legendre quadrature             CCC */
/* CC      - xbar[], xh[]: roots and weights of the G-L quadrature     CCC */
/* CC      - Cnp  : array containing the binomial coefficients         CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC      - sd : value of the S\bar{D}  approximation                 CCC */
/* CC      - t  : calculation time.                                    CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int sdbar42_lrzv__(integer *nx, integer *lambda, doublereal *
	v, integer *nu12, integer *ng12, doublereal *ab, doublereal *a12, 
	doublereal *b12, integer *nu34, integer *ng34, doublereal *cd, 
	doublereal *a34, doublereal *b34, integer *jc, integer *nrac, 
	doublereal *xbar, doublereal *xh, integer *nlag, doublereal *rlag, 
	doublereal *wlag, doublereal *cnp, doublereal *sd, integer *nconv, 
	doublereal *t)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal epssd = 1e-15;
    static doublereal tiny = 1e-150;
    static doublereal epsv = 1e-15;

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), pow_di(doublereal *, integer *), exp(doublereal),
	     sin(doublereal);

    /* Local variables */
    extern /* Subroutine */ int cpu_time__(doublereal *);
    extern doublereal g_(doublereal *, doublereal *, doublereal *);
    static integer i__, j, k, n;
    static doublereal v1[101], ui, xi;
    static integer iii;
    static doublereal z12i, z34i, uij, sd_0__, calf, dijc, fact, dnjc, sd_n__;
    extern doublereal hatk_(integer *, doublereal *);
    static doublereal vden[101], xijc, gxij, xden, vplb, xijk;
    static integer nint;
    static doublereal size, vnum[101], xnum;
    extern doublereal calf4_(integer *, integer *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *);
    static doublereal xpnx, b_inf__, dxijc, b_sup__, vdenp[101], zpng12, 
	    terme, zpng34, xdenp, start, vnump[101], xnump, finish, weight;

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
	*nconv = 1;
/* .... v close to zero */
/* .... calculation using Gauss-Laguerre */
/* .... Taylor development of sin(vx)/vx of order 1 was used */
	if (*v < epsv) {
	    if (*lambda == 0) {
		*sd = 0.;
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
		    *sd += ui;
		}
		goto L400;
	    } else {
		*sd = 0.;
		goto L400;
	    }
	}
/* .... calculation by Gauss_Legendre */
	i__1 = *lambda + 1;
	vplb = pow_di(v, &i__1);
/* .... initialization nsd = 0 */
	terme = 0.;
	b_sup__ = 0.;
/* Computing MAX */
/* Computing MIN */
	i__2 = (integer) (1 / *v);
	i__1 = min(i__2,100);
	nint = max(i__1,1);
/*      nint = max(min(1/v**2, MAX_DIV), 1) */
/*      nint = max(min(1/v**3, MAX_DIV), 1) */
	size = pi / (nint * *v);
	i__1 = *jc + 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    v1[i__] = terme;
/* .... calculation of F[xi] */
	    ui = 0.;
	    i__2 = nint;
	    for (j = 1; j <= i__2; ++j) {
		b_inf__ = b_sup__;
		b_sup__ = b_inf__ + size;
		uij = 0.;
		i__3 = *nrac;
		for (k = 1; k <= i__3; ++k) {
		    xijk = ((b_sup__ - b_inf__) * xbar[k] + b_sup__ + b_inf__)
			     * .5;
		    uij += size * .5 * xh[k] * sin(xijk * *v) * calf4_(nx, 
			    lambda, &xijk, cnp, nu12, ng12, ab, a12, b12, 
			    nu34, ng34, cd, a34, b34);
		}
		ui += uij;
/* L10: */
	    }
	    v1[i__] += ui / vplb;
	    terme = v1[i__];
/* L5: */
	}
/* .... U0 */
	i__ = 0;
	xijc = (doublereal) (i__ + 1 + *jc) * pi / *v;
	dxijc = xijc * xijc;
	gxij = calf4_(nx, lambda, &xijc, cnp, nu12, ng12, ab, a12, b12, nu34, 
		ng34, cd, a34, b34);
	vnum[i__] = v1[i__ + *jc] / (dxijc * gxij);
	vden[i__] = 1. / (dxijc * gxij);
	vnump[i__] = v1[i__ + *jc] / dxijc;
	vdenp[i__] = 1. / dxijc;
/* .... U1 */
	i__ = 1;
	xijc = (doublereal) (i__ + 1 + *jc) * pi / *v;
	dxijc = xijc * xijc;
	gxij = calf4_(nx, lambda, &xijc, cnp, nu12, ng12, ab, a12, b12, nu34, 
		ng34, cd, a34, b34);
	vnum[i__] = v1[i__ + *jc] / (dxijc * gxij);
	vden[i__] = 1. / (dxijc * gxij);
	vnump[i__] = v1[i__ + *jc] / dxijc;
	vdenp[i__] = 1. / dxijc;
/* .... SD_0 */
	xnum = vnum[0] + vnum[1];
	xden = vden[0] + vden[1];
	xnump = vnump[0] + vnump[1];
	xdenp = vdenp[0] + vdenp[1];
	if ((d__1 = xnump / xnum - xdenp / xden, abs(d__1)) < tiny || (d__2 = 
		xnum / xnump - xden / xdenp, abs(d__2)) < tiny) {
	    sd_n__ = xnump / xdenp;
	    goto L25;
	} else {
	    sd_0__ = xnum / xden;
	}
/* .... Un */
	for (n = 1; n <= 99; ++n) {
	    dnjc = (doublereal) (n + 2 + *jc);
	    dijc = pow_di(&dnjc, &n);
	    xijc = dnjc * pi / *v;
	    dxijc = xijc * xijc;
	    gxij = calf4_(nx, lambda, &xijc, cnp, nu12, ng12, ab, a12, b12, 
		    nu34, ng34, cd, a34, b34);
/* .... calculation of F[x_{n+1+jc}] */
	    v1[n + 1 + *jc] = v1[n + *jc];
	    ui = 0.;
	    i__1 = nint;
	    for (j = 1; j <= i__1; ++j) {
		b_inf__ = b_sup__;
		b_sup__ = b_inf__ + size;
		uij = 0.;
		i__2 = *nrac;
		for (k = 1; k <= i__2; ++k) {
		    xijk = ((b_sup__ - b_inf__) * xbar[k] + b_sup__ + b_inf__)
			     * .5;
		    uij += size * .5 * xh[k] * sin(xijk * *v) * calf4_(nx, 
			    lambda, &xijk, cnp, nu12, ng12, ab, a12, b12, 
			    nu34, ng34, cd, a34, b34);
		}
		ui += uij;
/* L20: */
	    }
	    v1[n + 1 + *jc] += ui / vplb;
	    vnum[n + 1] = dijc * v1[n + 1 + *jc] / (dxijc * gxij);
	    vden[n + 1] = dijc / (dxijc * gxij);
	    vnump[n + 1] = dijc * v1[n + 1 + *jc] / dxijc;
	    vdenp[n + 1] = dijc / dxijc;
/* .... SD_n */
	    xnum = vnum[n + 1];
	    xden = vden[n + 1];
	    xnump = vnump[n + 1];
	    xdenp = vdenp[n + 1];
	    i__1 = n;
	    for (i__ = 0; i__ <= i__1; ++i__) {
		fact = (doublereal) ((i__ + 1 + *jc) * (n + 1)) / (doublereal)
			 (n + 1 - i__);
		vnum[i__] = fact * vnum[i__];
		vden[i__] = fact * vden[i__];
		vnump[i__] = fact * vnump[i__];
		vdenp[i__] = fact * vdenp[i__];
		xnum += vnum[i__];
		xden += vden[i__];
		xnump += vnump[i__];
		xdenp += vdenp[i__];
	    }
	    if ((d__1 = xnump / xnum - xdenp / xden, abs(d__1)) < tiny || (
		    d__2 = xnum / xnump - xden / xdenp, abs(d__2)) < tiny) {
		sd_n__ = xnump / xdenp;
		goto L25;
	    } else {
		sd_n__ = xnum / xden;
	    }
	    *nconv = n;
	    if ((d__1 = sd_0__ - sd_n__, abs(d__1)) < epssd) {
		goto L25;
	    }
	    sd_0__ = sd_n__;
/* L15: */
	}
	*nconv = n - 1;
L25:
	*sd = sd_n__;
L400:
	;
    }
    cpu_time__(&finish);
    *t = finish - start;
    return 0;
} /* sdbar42_lrzv__ */

