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
/* CC   This function evaluates the integrand of the semi-infinte      CCC */
/* CC   integrals appearing in the analytical expression of the four-  CCC */
/* CC   center Coulomb integrals.                                      CCC */
/* CC                                                                  CCC */
/* CC   The evaluation of these integrals is based on an infinite      CCC */
/* CC   series of definite integrals on [x_i, x_{i+1}] where  x_i are  CCC */
/* CC   the successive roots of the spherical Bessel function.         CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - nx, nu, ng, lambda : parameters of the integral           CCC */
/* CC      - r : distance between the centers.                         CCC */
/* CC      - v : norm of the vector v.                                 CCC */
/* CC      - a, b : parameters for gamma.                              CCC */
/* CC      - x : point at which we evaluate the integrand.             CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal fa4_(integer *nx, integer *lambda, doublereal *v, integer *nu12, 
	integer *ng12, doublereal *r12, doublereal *a12, doublereal *b12, 
	integer *nu34, integer *ng34, doublereal *r34, doublereal *a34, 
	doublereal *b34, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    extern doublereal g_(doublereal *, doublereal *, doublereal *);
    static doublereal dj[1001];
    static integer nm;
    static doublereal sj[1001], z12, z34;
    extern doublereal hatk_(integer *, doublereal *);
    extern /* Subroutine */ int sphj_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal xpnx, z12png, z34png;

    z12 = g_(a12, b12, x);
    z34 = g_(a34, b34, x);
    z12png = pow_di(&z12, ng12);
    z34png = pow_di(&z34, ng34);
    xpnx = pow_di(x, nx);
    d__1 = *v * *x;
    sphj_(lambda, &d__1, &nm, sj, dj);
    d__1 = *r12 * z12;
    d__2 = *r34 * z34;
    ret_val = xpnx * hatk_(nu12, &d__1) / z12png * hatk_(nu34, &d__2) / 
	    z34png * sj[*lambda];
    return ret_val;
} /* fa4_ */

