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

/*     Last change:  H    14 Mar 2009   11:23 am */
doublereal fscj_(integer *nx, integer *nu12, integer *nu34, doublereal *ab, 
	doublereal *cd, doublereal *a12, doublereal *b12, doublereal *a34, 
	doublereal *b34, integer *ng12, integer *ng34, doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    extern doublereal g_(doublereal *, doublereal *, doublereal *);
    static doublereal z12, z34;
    extern doublereal hatk_(integer *, doublereal *);
    static doublereal xpnx, zpng12, zpng34;

    z12 = g_(a12, b12, x);
    z34 = g_(a34, b34, x);
    zpng12 = pow_di(&z12, ng12);
    zpng34 = pow_di(&z34, ng34);
    xpnx = pow_di(x, nx);
    d__1 = *ab * z12;
    d__2 = *cd * z34;
    ret_val = xpnx * hatk_(nu12, &d__1) / zpng12 * hatk_(nu34, &d__2) / 
	    zpng34;
    return ret_val;
} /* fscj_ */

