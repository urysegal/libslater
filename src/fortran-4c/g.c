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
/* CC   Function: G                                                    CCC */
/* CC                                                                  CCC */
/* CC   This function evaluates Gamma(s,x)                             CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - a = (1-s)*zeta1*zeta1 + s * zeta2*zeta2                   CCC */
/* CC      - b = s*(1-s)                                               CCC */
/* CC      - x = point at which we evaluate the integral               CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal g_(doublereal *a, doublereal *b, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    ret_val = sqrt(*a + *b * *x * *x);
    return ret_val;
} /* g_ */

