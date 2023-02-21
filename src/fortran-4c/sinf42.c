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

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Function: sinf4                                                CCC */
/* CC                                                                  CCC */
/* CC   This function evaluates the derivative:                        CCC */
/* CC                   [d/(x dx)]^nj (x^{nj-1} g(x))                  CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - nx     : nx = (l1+l2+l3+l4) - (l'1+l'2+l'3+l'4)           CCC */
/* CC      - lambda : lambda == order of the spherical Bessel function CCC */
/* CC      - x      : point at which we evaluate the integrand         CCC */
/* CC      - Cnp    : array containing the binomial coefficients       CCC */
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
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal sinf4_(integer *nx, integer *nj, doublereal *x, doublereal *cnp, 
	doublereal *besk12, doublereal *besk34, integer *nu12, integer *ng12, 
	doublereal *ab, doublereal *z12, doublereal *b12, integer *nu34, 
	integer *ng34, doublereal *cd, doublereal *z34, doublereal *b34)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal dob[117], xbg[117], dob12[117], dob34[117], xbg34[117], 
	    cste, binoj, binol, binoi, binok, dmupi, termi, termj, termk, 
	    dmupk;
    extern /* Subroutine */ int power_(integer *, doublereal *, doublereal *);
    static doublereal cste_j__, cste_l__, cste_i__, cste_k__;
    extern /* Subroutine */ int tdobin_(integer *, integer *, doublereal *);

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
/* .... initialization of sinf4 */
    ret_val = 0.;
    i__1 = *nx - *nj - 1;
    cste = pow_di(x, &i__1) / (pow_di(z12, ng12) * pow_di(z34, ng34));
    i__1 = *nx + *nj - 1;
    tdobin_(&i__1, nj, dob);
    d__1 = *b34 * *x * *x / (*z34 * *z34);
    power_(nj, &d__1, xbg34);
    d__1 = *b12 * *z34 * *z34 / (*b34 * *z12 * *z12);
    power_(nj, &d__1, xbg);
    if ((*nu12 << 1) + 1 != *ng12) {
	i__1 = (*nu12 << 1) + 1 - *ng12;
	tdobin_(&i__1, nj, dob12);
    }
    if ((*nu34 << 1) + 1 != *ng34) {
	i__1 = (*nu34 << 1) + 1 - *ng34;
	tdobin_(&i__1, nj, dob34);
    }
    i__1 = *nj;
    for (l = 0; l <= i__1; ++l) {
	binol = cnp[*nj * (*nj + 1) / 2 + l];
	cste_l__ = binol * xbg34[l] * dob[*nj - l];
	termj = 0.;
	i__2 = l;
	for (j = 0; j <= i__2; ++j) {
	    binoj = cnp[l * (l + 1) / 2 + j];
	    cste_j__ = binoj * xbg[j];
	    if ((*nu12 << 1) + 1 == *ng12) {
		termi = pow_di(&c_b2, &j) * besk12[*nu12 + j];
	    } else {
		termi = 0.;
		dmupi = 1.;
		i__3 = j;
		for (i__ = 0; i__ <= i__3; ++i__) {
		    binoi = cnp[j * (j + 1) / 2 + i__];
		    cste_i__ = dmupi * binoi * dob12[j - i__] * besk12[*nu12 
			    + i__];
		    dmupi = -dmupi;
		    termi += cste_i__;
		}
	    }
	    if ((*nu34 << 1) + 1 == *ng34) {
		i__3 = l - j;
		termk = pow_di(&c_b2, &i__3) * besk34[*nu34 + l - j];
	    } else {
		termk = 0.;
		dmupk = 1.;
		i__3 = l - j;
		for (k = 0; k <= i__3; ++k) {
		    binok = cnp[(l - j) * (l - j + 1) / 2 + k];
		    cste_k__ = dmupk * binok * dob34[l - j - k] * besk34[*
			    nu34 + k];
		    dmupk = -dmupk;
		    termk += cste_k__;
		}
	    }
	    termj += cste_j__ * termi * termk;
	}
	ret_val += cste_l__ * termj;
    }
    ret_val = cste * ret_val;
    return ret_val;
} /* sinf4_ */

