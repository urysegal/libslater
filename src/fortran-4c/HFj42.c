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
doublereal hfj_(integer *nx, integer *nu12, integer *nu34, doublereal *ab, 
	doublereal *cd, doublereal *v, doublereal *a12, doublereal *b12, 
	doublereal *a34, doublereal *b34, integer *ng12, integer *ng34, 
	integer *nj, doublereal *x)
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
    static doublereal xpnx, zpng12, zpng34;

    z12 = g_(a12, b12, x);
    z34 = g_(a34, b34, x);
    zpng12 = pow_di(&z12, ng12);
    zpng34 = pow_di(&z34, ng34);
    xpnx = pow_di(x, nx);
    i__1 = *nj + 1;
    d__1 = *v * *x;
    sphj_(&i__1, &d__1, &nm, sj, dj);
    d__1 = *ab * z12;
    d__2 = *cd * z34;
    ret_val = xpnx * hatk_(nu12, &d__1) / zpng12 * hatk_(nu34, &d__2) / 
	    zpng34 * (-(*v)) * sj[*nj + 1];
/*      nn1=2*n1+1-ng1 */
/*      nn2=2*n2+1-ng2 */
/*      z1 = g(a1,b1,x) */
/*      z2 = g(a2,b2,x) */
/*      HFj = x**nx * z1**nn1 * z2**nn2  * tk(n1,r1,z1) * */
/*     $      tk(n2,r2,z2) * (- v) * SJ(nj+1) */
/*      Dbsj(nj,x,v,1) */
    return ret_val;
} /* hfj_ */

