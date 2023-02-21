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

static integer c__9 = 9;
static integer c__1 = 1;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Subroutine: dfactorial                                         CCC */
/* CC                                                                  CCC */
/* CC   This subroutine calculates n!! (double factorial) for n=0 up   CCC */
/* CC   to nmax. The result is stored in the array Dfact such that     CCC */
/* CC   Dfact[n] = n!!                                                 CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int dfactorial_(integer *imax, doublereal *dfact)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };


    if (*imax < 0) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "error in the function dfactorial: nmax<0", (
		ftnlen)40);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    dfact[0] = 1.;
    if (*imax > 0) {
	dfact[1] = 1.;
	i__1 = *imax;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    dfact[i__] = (doublereal) i__ * dfact[i__ - 2];
	}
    }
    return 0;
} /* dfactorial_ */

