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
/* CC   Pour l'application de la transformation D avec la reduction de CCC */
/* CC   l'ordre de l'equation differentielle, on a besoin de la        CCC */
/* CC   fonction HFj(x)                                                CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal hfa4_(integer *nx, integer *lambda, doublereal *v, integer *nu12, 
	integer *ng12, doublereal *r12, doublereal *a12, doublereal *b12, 
	integer *nu34, integer *ng34, doublereal *r34, doublereal *a34, 
	doublereal *b34, doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    extern doublereal g_(doublereal *, doublereal *, doublereal *);
    static doublereal dj[1002];
    static integer nm;
    static doublereal sj[1002], z12, z34;
    extern doublereal hatk_(integer *, doublereal *);
    extern /* Subroutine */ int sphj_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal xpnx, z12png, z34png;

    z12 = g_(a12, b12, x);
    z34 = g_(a34, b34, x);
    z12png = pow_di(&z12, ng12);
    z34png = pow_di(&z34, ng34);
    xpnx = pow_di(x, nx);
    i__1 = *lambda + 1;
    d__1 = *v * *x;
    sphj_(&i__1, &d__1, &nm, sj, dj);
    d__1 = *r12 * z12;
    d__2 = *r34 * z34;
    ret_val = xpnx * hatk_(nu12, &d__1) / z12png * hatk_(nu34, &d__2) / 
	    z34png * (-(*v)) * sj[*lambda + 1];
    return ret_val;
} /* hfa4_ */

doublereal hfa4w_(integer *nx, integer *lambda, doublereal *v, integer *nu12, 
	integer *ng12, doublereal *r12, doublereal *a12, doublereal *b12, 
	integer *nu34, integer *ng34, doublereal *r34, doublereal *a34, 
	doublereal *b34, doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    extern doublereal g_(doublereal *, doublereal *, doublereal *);
    static doublereal dj[1002];
    static integer nm;
    static doublereal sj[1002], z12, z34;
    extern doublereal hatk_(integer *, doublereal *);
    extern /* Subroutine */ int sphj_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal xpnx, z12png, z34png;

    z12 = g_(a12, b12, x);
    z34 = g_(a34, b34, x);
    z12png = pow_di(&z12, ng12);
    z34png = pow_di(&z34, ng34);
    xpnx = pow_di(x, nx);
    i__1 = *lambda + 1;
    d__1 = *v * *x;
    sphj_(&i__1, &d__1, &nm, sj, dj);
    d__1 = *r12 * z12;
    d__2 = *r34 * z34;
    ret_val = xpnx * hatk_(nu12, &d__1) / z12png * hatk_(nu34, &d__2) / 
	    z34png * (-(*v)) * sj[*lambda + 1];
    return ret_val;
} /* hfa4w_ */

doublereal hfa4wx_(integer *nx, integer *lambda, doublereal *v, integer *nu12,
	 integer *ng12, doublereal *r12, doublereal *a12, doublereal *b12, 
	integer *nu34, integer *ng34, doublereal *r34, doublereal *a34, 
	doublereal *b34, doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double pow_di(doublereal *, integer *);

    /* Local variables */
    static doublereal dj[1002];
    static integer nm;
    static doublereal sj[1002];
    extern /* Subroutine */ int sphj_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal xpnx;

    xpnx = pow_di(x, nx);
    i__1 = *lambda + 1;
    d__1 = *v * *x;
    sphj_(&i__1, &d__1, &nm, sj, dj);
    ret_val = xpnx * (-(*v)) * sj[*lambda + 1];
    return ret_val;
} /* hfa4wx_ */

