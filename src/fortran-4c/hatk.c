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
/* CC   Function: Hatk                                                 CCC */
/* CC                                                                  CCC */
/* CC   This function evaluates the reduced Bessel function (RBF) of   CCC */
/* CC   half-integral order nu=n+1/2, n \in Z.  The function is        CCC */
/* CC   evaluated at the point x.                                      CCC */
/* CC                                                                  CCC */
/* CC   The RBFs are calculated using a three-term recurrence formula, CCC */
/* CC   stable in the forward direction.                               CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal hatk_(integer *n, doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal terme_0__, terme_1__, terme_i__;

    terme_0__ = exp(-(*x));
    if (*n == 0) {
	terme_i__ = terme_0__;
    }
    terme_1__ = exp(-(*x)) * (*x + 1.);
    if (*n == 1) {
	terme_i__ = terme_1__;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	terme_i__ = (doublereal) ((i__ << 1) - 1) * terme_1__ + *x * *x * 
		terme_0__;
	terme_0__ = terme_1__;
	terme_1__ = terme_i__;
    }
    ret_val = terme_i__;
/* L10: */
    return ret_val;
} /* hatk_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Function: tildk                                                CCC */
/* CC                                                                  CCC */
/* CC   This is a renormalized RBF introduced by Weniger and Steinborn CCC */
/* CC   in order to avoid overflow.  It is defined as :                CCC */
/* CC                                                                  CCC */
/* CC        \tilde{k}_{n-1/2}(x) = (2^n n!)^{-1} \hat{k}_{n-1/2}(x)   CCC */
/* CC                                                                  CCC */
/* CC   See: E. Joachim Weniger and E. Otto Steinborn, Numerical       CCC */
/* CC        properties of the convolution theorems of B functions,    CCC */
/* CC        Phys. Rev. A, Vol 28, No 4, 1983. p 2028,                 CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal tildk_(integer *nmax, doublereal *z__)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal tildk_0__, tildk_1__, tildk_i__;

    tildk_0__ = exp(-(*z__)) / 2.;
    if (*nmax == 0) {
	tildk_i__ = tildk_0__;
    }
    tildk_1__ = (*z__ + 1.) * exp(-(*z__)) / 8.;
    if (*nmax == 1) {
	tildk_i__ = tildk_1__;
    }
    i__1 = *nmax;
    for (i__ = 2; i__ <= i__1; ++i__) {
	tildk_i__ = ((doublereal) ((i__ << 1) - 1) * tildk_1__ + *z__ * *z__ /
		 (doublereal) (i__ << 1) * tildk_0__) / (doublereal) ((i__ + 1) << 1);
    }
    ret_val = tildk_i__;
    return ret_val;
} /* tildk_ */

