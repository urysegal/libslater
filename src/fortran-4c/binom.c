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

//static integer c__9 = 9;
//static integer c__1 = 1;
extern double boost_factorial(unsigned  n);
extern double boost_choose(unsigned n, unsigned k);

/*     Last change:  S    14 Apr 2007   10:54 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC   Programmer: Hassan Safouhi                                     CCC */
/* CC                                                                  CCC */
/* CC   Function: binom                                                CCC */
/* CC                                                                  CCC */
/* CC   This function evaluates the binomial coefficient               CCC */
/* CC                  nCm = n!/(m!(n-m)!) for n > m                   CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal binom_(integer *n, integer *m)
{
    return boost_choose(*n, *m);
#if 0
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;
    extern doublereal fact_(integer *);
    static doublereal prod;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };


    if (*n < *m) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "binomial is impossible", (ftnlen)22);
	e_wsle();
	s_stop("", (ftnlen)0);
    }
    if (*m == *n) {
	ret_val = 1.;
    } else {
	prod = 1.;
	i__1 = *n;
	for (i__ = *m + 1; i__ <= i__1; ++i__) {
	    prod *= (doublereal) i__;
/* L5: */
	}
	i__1 = *n - *m;
	ret_val = prod / boost_factorial(i__1);
    }
    return ret_val;
#endif
} /* binom_ */

