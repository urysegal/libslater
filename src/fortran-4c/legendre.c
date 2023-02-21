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
/* CC   Subroutine: legendre                                           CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the associated Legendre polynomial   CCC */
/* CC   at the point x.                                                CCC */
/* CC                                                                  CCC */
/* CC   This subroutine is based on the following 3 relations:         CCC */
/* CC            P(l=m,m) =   (2m-1)!! * (1-x^2)**m/2                  CCC */
/* CC          P(l=m+1,m) = x * (2m+1) * P(l=m,m)                      CCC */
/* CC       (l -m) P(l,m) = x * (2l-1) * P(l-1,m) - (l+m-1) * P(l-2,m) CCC */
/* CC                                                                  CCC */
/* CC   The result is stored in the array Plm such that:               CCC */
/* CC              P(l,m) = Plm (l*(l+1)/2 + m)                        CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int legendre_(doublereal *x, integer *lmax, doublereal *plm)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *), sqrt(doublereal);

    /* Local variables */
    static integer l, lj, mj;
    static doublereal dfact;

    dfact = 1.;
    i__1 = *lmax;
    for (l = 0; l <= i__1; ++l) {
/*     Evaluation of P(l=m,m) */
	d__1 = 1. - *x * *x;
	plm[l * (l + 1) / 2 + l] = dfact * sqrt(pow_di(&d__1, &l));
/*     Evaluation of P(l=m+1,m) */
	if (l < *lmax) {
	    plm[(l + 1) * (l + 2) / 2 + l] = (doublereal) ((l << 1) + 1) * *x 
		    * plm[l * (l + 1) / 2 + l];
	}
/*     Use of the third relation: Evaluation of P(l+2, m) */
	i__2 = *lmax;
	for (lj = l + 2; lj <= i__2; ++lj) {
	    mj = l;
	    plm[lj * (lj + 1) / 2 + mj] = ((doublereal) ((lj << 1) - 1) * *x *
		     plm[lj * (lj - 1) / 2 + mj] - (doublereal) (lj + mj - 1) 
		    * plm[(lj - 1) * (lj - 2) / 2 + mj]) / (doublereal) (lj - 
		    mj);
/* L10: */
	}
	dfact *= (doublereal) ((l << 1) + 1);
/* L5: */
    }
    return 0;
} /* legendre_ */

