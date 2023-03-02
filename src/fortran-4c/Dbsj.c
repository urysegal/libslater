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
/* CC   Function: Dbsj                                                 CCC */
/* CC                                                                  CCC */
/* CC   Evaluation of the kth derivative of the spherical Bessel       CCC */
/* CC   function at the point x.                                       CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - n : order of the spherical Bessel function                CCC */
/* CC      - v : norm of the vector v                                  CCC */
/* CC      - x : point at which we evaluate the derivative             CCC */
/* CC      - k : order of the derivative                               CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal dbsj_(integer *n, doublereal *x, doublereal *v, integer *k)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, 
	    d__10;

    /* Local variables */
    static doublereal dj[1006];
    static integer nm;
    static doublereal sj[1006];
    extern /* Subroutine */ int sphj_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);

    i__1 = *n + 5;
    d__1 = *v * *x;
    sphj_(&i__1, &d__1, &nm, sj, dj);
    if (*k == 0) {
	ret_val = sj[*n];
    } else if (*k == 1) {
	ret_val = (doublereal) (*n) / *x * sj[*n] - *v * sj[*n + 1];
    } else if (*k == 2) {
/* Computing 2nd power */
	d__1 = *x;
/* Computing 2nd power */
	d__2 = *v;
	ret_val = (doublereal) (*n * (*n - 1)) / (d__1 * d__1) * sj[*n] - (
		doublereal) ((*n << 1) + 1) * *v / *x * sj[*n + 1] + d__2 * 
		d__2 * sj[*n + 2];
    } else if (*k == 3) {
/* Computing 3rd power */
	d__1 = *x;
/* Computing 2nd power */
	d__2 = *x;
/* Computing 2nd power */
	d__3 = *v;
/* Computing 3rd power */
	d__4 = *v;
	ret_val = (doublereal) (*n * (*n - 1) * (*n - 2)) / (d__1 * (d__1 * 
		d__1)) * sj[*n] - (doublereal) (*n * 3 * *n) * *v / (d__2 * 
		d__2) * sj[*n + 1] + (doublereal) (*n * 3 + 3) * (d__3 * d__3)
		 / *x * sj[*n + 2] - d__4 * (d__4 * d__4) * sj[*n + 3];
    } else if (*k == 4) {
/* Computing 4th power */
	d__1 = *x, d__1 *= d__1;
/* Computing 3rd power */
	d__2 = *x;
/* Computing 2nd power */
	d__3 = *v;
/* Computing 2nd power */
	d__4 = *x;
/* Computing 3rd power */
	d__5 = *v;
/* Computing 4th power */
	d__6 = *v, d__6 *= d__6;
	ret_val = (doublereal) (*n * (*n - 1) * (*n - 2) * (*n - 3)) / (d__1 *
		 d__1) * sj[*n] - (doublereal) ((*n << 1) * (*n - 1) * ((*n <<
		 1) - 1)) * *v / (d__2 * (d__2 * d__2)) * sj[*n + 1] + (
		doublereal) (((*n << 1) * (*n + 1) + 1) * 3) * (d__3 * d__3) /
		 (d__4 * d__4) * sj[*n + 2] - (doublereal) (((*n << 1) + 3 )<< 1) * (d__5 * (d__5 * d__5)) / *x * sj[*n + 3] + d__6 * d__6 *
		sj[*n + 4];
    } else if (*k == 5) {
/* Computing 5th power */
	d__1 = *x, d__2 = d__1, d__1 *= d__1;
/* Computing 4th power */
	d__3 = *x, d__3 *= d__3;
/* Computing 2nd power */
	d__4 = *v;
/* Computing 3rd power */
	d__5 = *x;
/* Computing 3rd power */
	d__6 = *v;
/* Computing 2nd power */
	d__7 = *x;
/* Computing 4th power */
	d__8 = *v, d__8 *= d__8;
/* Computing 5th power */
	d__9 = *v, d__10 = d__9, d__9 *= d__9;
	ret_val = (doublereal) (*n * (*n - 1) * (*n - 2) * (*n - 3) * (*n - 4)
		) / (d__2 * (d__1 * d__1)) * sj[*n] - (doublereal) (*n * 5 * (
		*n - 1) * (*n - 1) * (*n - 2)) * *v / (d__3 * d__3) * sj[*n + 
		1] + (doublereal) (*n * 5 * ((*n << 1) * *n + 1)) * (d__4 * 
		d__4) / (d__5 * (d__5 * d__5)) * sj[*n + 2] - (doublereal) (((
		*n << 1) * *n + (*n << 2) + 3) * 5) * (d__6 * (d__6 * d__6)) /
		 (d__7 * d__7) * sj[*n + 3] + (doublereal) ((*n + 2) * 5) * (
		d__8 * d__8) / *x * sj[*n + 4] - d__10 * (d__9 * d__9) * sj[*
		n + 5];
    }
    return ret_val;
} /* dbsj_ */

