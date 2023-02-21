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
/* CC   Function calf4                                                 CCC */
/* CC                                                                  CCC */
/* CC   This function evaluates the integrand of the semi-infinite     CCC */
/* CC   integral contained in the expression of the four-center        CCC */
/* CC   coulomb integral over B functions, as needed for the S\bar{D}  CCC */
/* CC   approach.                                                      CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - nx     : nx = (l1+l2+l3+l4) - (l'1+l'2+l'3+l'4)           CCC */
/* CC      - lambda : order of the spherical bessel function           CCC */
/* CC      - x      : point at which we evaluate the integrand         CCC */
/* CC      - Cnp    : binomial coefficients                            CCC */
/* CC      - nu12   : nu12 = n1 + n2 + l1 + l2 - l - j12               CCC */
/* CC      - ng12   : ng12 = 2(n1+l1+n2+l2)-(l'1+l'2+l)+1              CCC */
/* CC      - ab     : position of B wrt A                              CCC */
/* CC      - a12    : a12 = (1-s)*zeta1*zeta1 + s * zeta2*zeta2        CCC */
/* CC      - b12    : b12 = s*(1-s)                                    CCC */
/* CC      - nu34   : nu34 = n3 + n4 + l3 + l4 - l' - j34              CCC */
/* CC      - ng34   : ng34 = 2(n3+l3+n4+l4)-(l'3+l'4+l')+1             CCC */
/* CC      - cd     : position of C wrt D                              CCC */
/* CC      - a34    : a34 = (1-t)*zeta4*zeta4 + t * zeta3*zeta3        CCC */
/* CC      - b34    : b34 = t*(1-t)                                    CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal calf4_(integer *nx, integer *lambda, doublereal *x, doublereal *
	cnp, integer *nu12, integer *ng12, doublereal *ab, doublereal *a12, 
	doublereal *b12, integer *nu34, integer *ng34, doublereal *cd, 
	doublereal *a34, doublereal *b34)
{

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    extern doublereal g_(doublereal *, doublereal *, doublereal *);
    static doublereal z12, z34;
    extern /* Subroutine */ int bessel_red__(integer *, doublereal *, 
	    doublereal *);
    extern doublereal hatk_(integer *, doublereal *);
    static doublereal xpnx;
    extern doublereal sinf4_(integer *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal besk12[241], besk34[241], zpng12, zpng34;

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
    z12 = g_(a12, b12, x);
    z34 = g_(a34, b34, x);
    if (*lambda == 0) {
	i__1 = *nx - 1;
	xpnx = pow_di(x, &i__1);
	zpng12 = pow_di(&z12, ng12);
	zpng34 = pow_di(&z34, ng34);
	d__1 = *ab * z12;
	d__2 = *cd * z34;
	ret_val = hatk_(nu12, &d__1) / zpng12 * hatk_(nu34, &d__2) / zpng34 * 
		xpnx;
    } else {
	i__1 = *nu12 + *lambda;
	d__1 = z12 * *ab;
	bessel_red__(&i__1, &d__1, besk12);
	i__1 = *nu34 + *lambda;
	d__1 = z34 * *cd;
	bessel_red__(&i__1, &d__1, besk34);
	ret_val = sinf4_(nx, lambda, x, cnp, besk12, besk34, nu12, ng12, ab, &
		z12, b12, nu34, ng34, cd, &z34, b34);
    }
    return ret_val;
} /* calf4_ */

