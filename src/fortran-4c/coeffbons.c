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
static doublereal c_b3 = 2.;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: coeffbons                                          CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the coefficients which appear        CCC */
/* CC   when we express Slater-type orbitals |nlm> as finite linear    CCC */
/* CC   combinations of B functions |plm>.  The result is stored in    CCC */
/* CC   the array Ap[].                                                CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int coeffbons_(integer *n, integer *l, integer *pmin, 
	integer *pmax, doublereal *fact, doublereal *dfact, doublereal *ap)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static integer p;

    *pmin = (*n - *l + 1) / 2;
    *pmax = *n - *l;
    i__1 = *pmax;
    for (p = *pmin; p <= i__1; ++p) {
	i__2 = *pmax - p;
	i__3 = *l + p;
	ap[p] = pow_di(&c_b2, &i__2) * pow_di(&c_b3, &i__3) * fact[*pmax] * 
		fact[*l + p] / (fact[(p << 1) - *n + *l] * dfact[(*pmax - p) *
		 2]);
    }
    return 0;
} /* coeffbons_ */

