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

/*     Last change:  S    14 Apr 2007   11:19 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC   Programmer: Lilian Berlu                                       CCC */
/* CC   Function: gaunt1                                               CCC */
/* CC                                                                  CCC */
/* CC   This function evaluates the coefficient <L3M3 | L2 M2 | L1 M1> CCC */
/* CC                                                                  CCC */
/* CC   It is based on relation 5.83 p 55 of "Angular Momentum         CCC */
/* CC   Techniques in Quantum Mechanics" by V. Devanathan.             CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal gaunt1_(integer *l3, integer *m3, integer *l2, integer *m2, 
	integer *l1, integer *m1, doublereal *fact)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b, c__;
    static integer nu;
    static doublereal cg1, cg2, cnu, dnu, cste;
    static integer nu_max__;
    static doublereal fact_nu__;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*  Programmer : Hassan Safouhi                  sept  25 2002 */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* .... pi */
/* .... Tiny   => the smallest number */
/* .... Huge   => the largest number */
/* .... Eps    => Epsilon used in general */
/* .... EpsV   => limit to consider v as zero */
/* .... EpsLu  => Epsilon used in LU subroutine */
/* .... EpsSS  =>  limit used for the pre-determined accuracy for the sum using infinite series */
/* .... EpsDB  =>  limit used for the pre-determined accuracy for \bar{D} */
/* .... EpsEP  =>  limit used for the pre-determined accuracy for epsilon algorithm */
/* .... EpsUL  =>  limit used for the pre-determined accuracy for Levin's u */
/* .... EpsHD  =>  limit used for the pre-determined accuracy for H\bar{D} */
/* .... EpsER  =>  limit used for the pre-determined accuracy for epsilon RECURSIVE */
/* .... EpsUR  =>  limit used for the pre-determined accuracy for Levin's u RECURSIVE */
/* .... EpsSD  =>  limit used for the pre-determined accuracy for S\bar{D} */
/* .... EpsWD  =>  limit used for the pre-determined accuracy for W */
/* .... EpsSC  =>  limit used for the pre-determined accuracy for SD Sin Cos */
/* .....FOR THE CALCULATION IN bESSEL SIN-COS, USE THE FOLLOWING */
/* ......... data EpsSC, EpsSD, EpsWD / 1.0d-13, 1.0d-14, 1.0d-14/ */
/* .... Cutoff */
/*      data Cutoff/1.0d-12/ */
/* .... RADEG => Convert degrees to radians */
/*      data RADEG /1.74532925199432958d-02/ */
/* .... nombre maximal d'iteration dans le cycle SCF */
/*      data MAX_ITER / 2000/ */
    if ((*l1 + *l2 + *l3 + 1) / 2 == (*l1 + *l2 + *l3) / 2 && *m1 + *m2 == *
	    m3 && abs(*m3) <= *l3 && abs(*m2) <= *l2 && abs(*m1) <= *l1 && *
	    l1 + *l2 - *l3 >= 0 && *l1 + *l3 - *l2 >= 0 && *l3 + *l2 - *l1 >= 
	    0) {
/* .... normalization constant of the integral */
	cste = sqrt((*l1 * 2. + 1.) * (*l2 * 2. + 1.) / (pi * 4. * (*l3 * 2. 
		+ 1.)));
/* .... first Clebsch-Gordan coefficient [l1 m1, l2 m2, l3 m3] */
	a = fact[*l1 + *l2 - *l3] * fact[*l1 - *l2 + *l3] * fact[-(*l1) + *l2 
		+ *l3] / fact[*l1 + *l2 + *l3 + 1];
	b = fact[*l1 + *m1] * fact[*l2 + *m2] * fact[*l3 + *m3] * fact[*l1 - *
		m1] * fact[*l2 - *m2] * fact[*l3 - *m3];
/* .... the highest order for nu is determined by the smallest factorial */
/* Computing MIN */
/* Computing MIN */
	i__3 = *l1 - *m1, i__4 = *l2 + *m2;
	i__1 = *l1 + *l2 - *l3, i__2 = min(i__3,i__4);
	nu_max__ = min(i__1,i__2);
	fact_nu__ = 1.;
	c__ = 0.;
	i__1 = nu_max__;
	for (nu = 0; nu <= i__1; ++nu) {
	    if (*l1 + *l2 - *l3 - nu >= 0 && *l1 - *m1 - nu >= 0 && *l2 + *m2 
		    - nu >= 0 && *l3 - *l2 + *m1 + nu >= 0 && *l3 - *l1 - *m2 
		    + nu >= 0) {
		cnu = fact[*l1 + *l2 - *l3 - nu] * fact[*l1 - *m1 - nu] * 
			fact[*l2 + *m2 - nu] * fact[*l3 - *l2 + *m1 + nu] * 
			fact[*l3 - *l1 - *m2 + nu];
		c__ += fact_nu__ / cnu;
	    }
	    dnu = (doublereal) (nu + 1);
	    fact_nu__ /= dnu;
	}
	cg1 = sqrt((*l3 * 2. + 1.) * a * b) * c__;
/* .... second Clebsch-Gordan coefficient [l1 0, l2 0, l3 0] */
	b = fact[*l1] * fact[*l1] * fact[*l2] * fact[*l2] * fact[*l3] * fact[*
		l3];
/* .... the highest order for nu is determined by the smallest factorial */
/* Computing MIN */
	i__1 = *l1 + *l2 - *l3, i__2 = min(*l1,*l2);
	nu_max__ = min(i__1,i__2);
	c__ = 0.;
	fact_nu__ = 1.;
	i__1 = nu_max__;
	for (nu = 0; nu <= i__1; ++nu) {
	    if (*l1 + *l2 - *l3 - nu >= 0 && *l1 - nu >= 0 && *l2 - nu >= 0 &&
		     *l3 - *l2 + nu >= 0 && *l3 - *l1 + nu >= 0) {
		cnu = fact[*l1 + *l2 - *l3 - nu] * fact[*l1 - nu] * fact[*l2 
			- nu] * fact[*l3 - *l2 + nu] * fact[*l3 - *l1 + nu];
		c__ += fact_nu__ / cnu;
	    }
	    dnu = (doublereal) (nu + 1);
	    fact_nu__ /= dnu;
	}
	cg2 = sqrt((*l3 * 2. + 1.) * a * b) * c__;
/* .... Gaunt coefficient */
	ret_val = cste * cg1 * cg2;
    } else {
	ret_val = 0.;
    }
    return ret_val;
} /* gaunt1_ */

