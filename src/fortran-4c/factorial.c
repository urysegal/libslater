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
/* CC   Subroutine: factorial                                          CCC */
/* CC                                                                  CCC */
/* CC   This function calculates n! for n=0 up to n=nmax. The result   CCC */
/* CC   is stored in the array fact such that:  Fact[n] = n!           CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int factorial_(integer *nmax, doublereal *fact)
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


    if (*nmax < 0) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "error in the function factorial: nmax < 0", (
		ftnlen)41);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    fact[0] = 1.;
    i__1 = *nmax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fact[i__] = (doublereal) i__ * fact[i__ - 1];
    }
    return 0;
} /* factorial_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC   Function: fact                                                 CCC */
/* CC                                                                  CCC */
/* CC   This function evaluates the factorial n!                       CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal fact_(integer *n)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__;

    ret_val = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val *= (doublereal) i__;
/* L5: */
    }
    return ret_val;
} /* fact_ */

