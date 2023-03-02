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
/* CC   Function: DG                                                   CCC */
/* CC                                                                  CCC */
/* CC   Evaluation of the kth dedrivative of the function \gamma(x)^n  CCC */
/* CC   at the point x.                                                CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - n    : power of the function G  (\gamma)                  CCC */
/* CC      - a, b : parameters of \gamma                               CCC */
/* CC      - x    : point at which we evaluate the derivative          CCC */
/* CC      - k    : order of the derivative                            CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal dg_(doublereal *a, doublereal *b, doublereal *x, integer *n, 
	integer *k)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    extern doublereal g_(doublereal *, doublereal *, doublereal *);

    if (*k == 0) {
	d__1 = g_(a, b, x);
	ret_val = pow_di(&d__1, n);
    } else if (*k == 1) {
	d__1 = g_(a, b, x);
	i__1 = *n - 2;
	ret_val = (doublereal) (*n) * *b * *x * pow_di(&d__1, &i__1);
    } else if (*k == 2) {
	d__1 = g_(a, b, x);
	i__1 = *n - 2;
/* Computing 2nd power */
	d__2 = *b * *x;
	d__3 = g_(a, b, x);
	i__2 = *n - 4;
	ret_val = (doublereal) (*n) * *b * pow_di(&d__1, &i__1) + (doublereal)
		 (*n * (*n - 2)) * (d__2 * d__2) * pow_di(&d__3, &i__2);
    } else if (*k == 3) {
/* Computing 2nd power */
	d__1 = *b;
	d__2 = g_(a, b, x);
	i__1 = *n - 4;
/* Computing 3rd power */
	d__3 = *b * *x;
	d__4 = g_(a, b, x);
	i__2 = *n - 6;
	ret_val = (doublereal) (*n * 3 * (*n - 2)) * (d__1 * d__1) * *x * 
		pow_di(&d__2, &i__1) + (doublereal) (*n * (*n - 2) * (*n - 4))
		 * (d__3 * (d__3 * d__3)) * pow_di(&d__4, &i__2);
    } else if (*k == 4) {
/* Computing 2nd power */
	d__1 = *b;
	d__2 = g_(a, b, x);
	i__1 = *n - 4;
/* Computing 3rd power */
	d__3 = *b;
/* Computing 2nd power */
	d__4 = *x;
	d__5 = g_(a, b, x);
	i__2 = *n - 6;
/* Computing 4th power */
	d__6 = *b * *x, d__6 *= d__6;
	d__7 = g_(a, b, x);
	i__3 = *n - 8;
	ret_val = (doublereal) (*n * 3 * (*n - 2)) * (d__1 * d__1) * pow_di(&
		d__2, &i__1) + (doublereal) (*n * 6 * (*n - 2) * (*n - 4)) * (
		d__3 * (d__3 * d__3)) * (d__4 * d__4) * pow_di(&d__5, &i__2) 
		+ (doublereal) (*n * (*n - 2) * (*n - 4) * (*n - 6)) * (d__6 *
		 d__6) * pow_di(&d__7, &i__3);
    } else if (*k == 5) {
/* Computing 3rd power */
	d__1 = *b;
	d__2 = g_(a, b, x);
	i__1 = *n - 6;
/* Computing 4th power */
	d__3 = *b, d__3 *= d__3;
/* Computing 3rd power */
	d__4 = *x;
	d__5 = g_(a, b, x);
	i__2 = *n - 8;
/* Computing 5th power */
	d__6 = *b * *x, d__7 = d__6, d__6 *= d__6;
	d__8 = g_(a, b, x);
	i__3 = *n - 10;
	ret_val = (doublereal) (*n * 15 * (*n - 2) * (*n - 4)) * (d__1 * (
		d__1 * d__1)) * *x * pow_di(&d__2, &i__1) + (doublereal) (*n *
		 10 * (*n - 2) * (*n - 4) * (*n - 6)) * (d__3 * d__3) * (d__4 
		* (d__4 * d__4)) * pow_di(&d__5, &i__2) + (doublereal) (*n * (
		*n - 2) * (*n - 4) * (*n - 6) * (*n - 8)) * (d__7 * (d__6 * 
		d__6)) * pow_di(&d__8, &i__3);
    }
    return ret_val;
} /* dg_ */

