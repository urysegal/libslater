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
/* CC   Function: Dtk                                                  CCC */
/* CC                                                                  CCC */
/* CC   Evaluation of the kth derivative of the function:              CCC */
/* CC          \frac{\hat{k}_{\nu}(\gamma(x).r)}{\gamma(x)^{2\nu}}     CCC */
/* CC   at the point x.                                                CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - nu   : order of the reduced Bessel function               CCC */
/* CC      - r    : position of B2 in intnuc and B3/B4 in bielectronic CCC */
/* CC      - a, b : parameters of gamma                                CCC */
/* CC      - x    : point at which we evaluate the derivative          CCC */
/* CC      - k    : order of the derivative                            CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal dtk_(integer *nu, doublereal *r__, doublereal *a, doublereal *b, 
	doublereal *x, integer *k)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    extern doublereal g_(doublereal *, doublereal *, doublereal *);
    static doublereal z__;
    extern doublereal tk_(integer *, doublereal *, doublereal *);

    z__ = g_(a, b, x);
    if (*k == 0) {
	ret_val = tk_(nu, r__, &z__);
    } else if (*k == 1) {
	i__1 = *nu + 1;
	ret_val = -(*b) * *x * tk_(&i__1, r__, &z__);
    } else if (*k == 2) {
	i__1 = *nu + 1;
/* Computing 2nd power */
	d__1 = *b * *x;
	i__2 = *nu + 2;
	ret_val = -(*b) * tk_(&i__1, r__, &z__) + d__1 * d__1 * tk_(&i__2, 
		r__, &z__);
    } else if (*k == 3) {
/* Computing 2nd power */
	d__1 = *b;
	i__1 = *nu + 2;
/* Computing 3rd power */
	d__2 = *b * *x;
	i__2 = *nu + 3;
	ret_val = d__1 * d__1 * 3. * *x * tk_(&i__1, r__, &z__) - d__2 * (
		d__2 * d__2) * tk_(&i__2, r__, &z__);
    } else if (*k == 4) {
/* Computing 2nd power */
	d__1 = *b;
	i__1 = *nu + 2;
/* Computing 3rd power */
	d__2 = *b;
/* Computing 2nd power */
	d__3 = *x;
	i__2 = *nu + 3;
/* Computing 4th power */
	d__4 = *b * *x, d__4 *= d__4;
	i__3 = *nu + 4;
	ret_val = d__1 * d__1 * 3. * tk_(&i__1, r__, &z__) - d__2 * (d__2 * 
		d__2) * 6. * (d__3 * d__3) * tk_(&i__2, r__, &z__) + d__4 * 
		d__4 * tk_(&i__3, r__, &z__);
    } else if (*k == 5) {
/* Computing 3rd power */
	d__1 = *b;
	i__1 = *nu + 3;
/* Computing 4th power */
	d__2 = *b, d__2 *= d__2;
/* Computing 3rd power */
	d__3 = *x;
	i__2 = *nu + 4;
/* Computing 5th power */
	d__4 = *b * *x, d__5 = d__4, d__4 *= d__4;
	i__3 = *nu + 5;
	ret_val = d__1 * (d__1 * d__1) * -15. * *x * tk_(&i__1, r__, &z__) + 
		d__2 * d__2 * 10. * (d__3 * (d__3 * d__3)) * tk_(&i__2, r__, &
		z__) - d__5 * (d__4 * d__4) * tk_(&i__3, r__, &z__);
    }
    return ret_val;
} /* dtk_ */

/* .... tk(n,r,z) : This function occurs in the analytic expression of */
/* .... the integrand of the semi-infinite integral with the spherical Bessel function */
doublereal tk_(integer *nu, doublereal *r__, doublereal *z__)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    extern doublereal hatk_(integer *, doublereal *);

    d__1 = *r__ * *z__;
    i__1 = (*nu << 1) + 1;
    ret_val = hatk_(nu, &d__1) / pow_di(z__, &i__1);
    return ret_val;
} /* tk_ */

