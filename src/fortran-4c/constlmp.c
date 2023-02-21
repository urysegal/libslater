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
/* CC   This subroutine evaluates the coefficients appearing in the    CCC */
/* CC   sum over l' and m' of equation (43) of the three-center        CCC */
/* CC   nuclear attraction integral.                                   CCC */
/* CC                                                                  CCC */
/* CC   This series of coefficients is common to all of the B          CCC */
/* CC   functions arising in the decomposition of the base of Slater   CCC */
/* CC   functions.  Consequently, this function can be called before   CCC */
/* CC   the function which calculates the integrals.                   CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC        - l, m : azimuthal and magnetic quantum numbers           CCC */
/* CC                                                                  CCC */
/* CC   output :                                                       CCC */
/* CC       - Cstlmp[] : array of coefficients                         CCC */
/* CC                 Cstlmp[lp*(lp+1)+mp] contains the term (lp, mp). CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int constlmp_(integer *l, integer *m, doublereal *fact, 
	doublereal *dfact, doublereal *cstlmp)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer lp, mp;
    extern doublereal gaunt1_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *);

    i__1 = *l;
    for (lp = 0; lp <= i__1; ++lp) {
	i__2 = lp;
	for (mp = -lp; mp <= i__2; ++mp) {
	    i__3 = *l - lp;
	    i__4 = *m - mp;
	    cstlmp[lp * (lp + 1) + mp] = gaunt1_(l, m, &lp, &mp, &i__3, &i__4,
		     fact) / (dfact[(lp << 1) + 1] * dfact[((*l - lp) << 1) + 1]
		    );
	}
    }
    return 0;
} /* constlmp_ */

