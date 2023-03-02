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

/*     Last change:  S    14 Apr 2007   10:53 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC   Programmer: Hassan Safouhi                                     CCC */
/* CC                                                                  CCC */
/* CC   Subroutine: bessel_red                                         CCC */
/* CC                                                                  CCC */
/* CC   This subroutine evaluates the reduced Bessel functions of      CCC */
/* CC   half-integral orders up to nmax+1/2, at the point z. They are  CCC */
/* CC   evaluated using the following recurrence relation:             CCC */
/* CC                                                                  CCC */
/* CC     k_{n+1/2}(z) = z^2 k_{n-2+1/2}(z) + (2n-1) k_{n-1+1/2}(z),   CCC */
/* CC                                                                  CCC */
/* CC   knowing that k_{1/2}(z) = e^(-z)                               CCC */
/* CC                                                                  CCC */
/* CC   The values are stored in the array Bred, such that:            CCC */
/* CC                                                                  CCC */
/* CC                    Bred[i] = k_{i+1/2}(z)                        CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int bessel_red__(integer *nmax, doublereal *z__, doublereal *
	bred)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer i__;

/* .... initialization of the series */
    bred[0] = exp(-(*z__));
    if (*nmax > 0) {
	bred[1] = (*z__ + 1.) * exp(-(*z__));
	i__1 = *nmax;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    bred[i__] = (doublereal) ((i__ << 1) - 1) * bred[i__ - 1] + *z__ *
		     *z__ * bred[i__ - 2];
	}
    }
    return 0;
} /* bessel_red__ */

