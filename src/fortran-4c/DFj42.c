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

/*     Last change:  S     4 May 2007   11:21 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC   Programmer: Hassan Safouhi                                     CCC */
/* CC                                                                  CCC */
/* CC   Function: DFj                                                  CCC */
/* CC                                                                  CCC */
/* CC   Evaluation of the kth derivative of Fj   at the  point x.      CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal dfj_(integer *nx, integer *n1, integer *n2, doublereal *r1, 
	doublereal *r2, doublereal *v, doublereal *a1, doublereal *b1, 
	doublereal *a2, doublereal *b2, integer *ng1, integer *ng2, integer *
	nj, doublereal *x, integer *k)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j, l, m;
    extern doublereal dg_(doublereal *, doublereal *, doublereal *, integer *,
	     integer *);
    static integer ie;
    extern doublereal dx_(integer *, doublereal *, integer *);
    static integer nn1, nn2;
    extern doublereal dtk_(integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *);
    static doublereal sum1, sum2, sum3, sum4;
    extern doublereal dbsj_(integer *, doublereal *, doublereal *, integer *),
	     binom_(integer *, integer *);

    nn1 = (*n1 << 1) + 1 - *ng1;
    nn2 = (*n2 << 1) + 1 - *ng2;
    ret_val = 0.;
    i__1 = *k;
    for (j = 0; j <= i__1; ++j) {
	sum1 = 0.;
	i__2 = *k - j;
	for (i__ = 0; i__ <= i__2; ++i__) {
	    sum2 = 0.;
	    i__3 = *k - j - i__;
	    for (l = 0; l <= i__3; ++l) {
		sum3 = 0.;
		i__4 = *k - j - i__ - l;
		for (m = 0; m <= i__4; ++m) {
		    sum4 = 0.;
		    i__5 = *k - j - i__ - l - m;
		    for (ie = 0; ie <= i__5; ++ie) {
			i__6 = *k - j - i__ - l - m;
			i__7 = *k - j - i__ - l - m - ie;
			sum4 += binom_(&i__6, &ie) * dtk_(n2, r2, a2, b2, x, &
				ie) * dbsj_(nj, x, v, &i__7);
/* L25: */
		    }
		    i__5 = *k - j - i__ - l;
		    sum3 += binom_(&i__5, &m) * dtk_(n1, r1, a1, b1, x, &m) * 
			    sum4;
/* L20: */
		}
		i__4 = *k - j - i__;
		sum2 += binom_(&i__4, &l) * dg_(a2, b2, x, &nn2, &l) * sum3;
/* L15: */
	    }
	    i__3 = *k - j;
	    sum1 += binom_(&i__3, &i__) * dg_(a1, b1, x, &nn1, &i__) * sum2;
/* L10: */
	}
	ret_val += binom_(k, &j) * dx_(nx, x, &j) * sum1;
/* L5: */
    }
    return ret_val;
} /* dfj_ */

