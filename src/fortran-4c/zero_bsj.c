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
static integer c__3 = 3;
static integer c__5 = 5;

/*     Last change:  H     4 Jul 2010    2:05 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC   Programmer: Hassan Safouhi                                     CCC */
/* CC   Subroutine : zero_bsj                                          CCC */
/* CC                                                                  CCC */
/* CC   This subroutine finds the first nmax roots of the spherical    CCC */
/* CC   Bessel functions j_{l+1/2} where l=0,1,2,...,lmax.             CCC */
/* CC                                                                  CCC */
/* CC   Input :                                                        CCC */
/* CC      - nmax  : number of roots required                          CCC */
/* CC      - lmax  : highest order of Bessel functions for which we    CCC */
/* CC                will find the roots.                              CCC */
/* CC                                                                  CCC */
/* CC   Output :                                                       CCC */
/* CC      - rac[] : array containing the roots                        CCC */
/* CC        rac[nmax*(l-1)+i] == the ith zero of j_{l+1/2}            CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int zero_bsj__(integer *lmax, integer *nmax, doublereal *rac)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double acos(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static doublereal a, b;
    static integer i__, l;
    static doublereal pi, xn;
    static integer ier;
    static doublereal diff, xbar[60]	/* was [30][2] */;
    static integer kmax;
    extern /* Subroutine */ int racdi_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    static doublereal terme;

    /* Fortran I/O blocks */
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };


    /* Parameter adjustments */
    --rac;

    /* Function Body */
    kmax = 200;
    pi = acos(-1.);
/* .... roots of j_{1+1/2} */
    xbar[0] = 4.4934094579090642085;
    xbar[30] = 7.7252518369377067842;
/* .... roots of j_{2+1/2} */
    xbar[1] = 5.7634591968945496632;
    xbar[31] = 9.095011330476355127;
/* .... roots of j_{3+1/2} */
    xbar[2] = 6.9879320005005194361;
    xbar[32] = 10.417118547379365268;
/* .... roots of j_{4+1/2} */
    xbar[3] = 8.1825614525712424552;
    xbar[33] = 11.704907154570390659;
/* .... roots of j_{5+1/2} */
    xbar[4] = 9.3558121110427467926;
    xbar[34] = 12.96653017277434472;
/* .... roots of j_{6+1/2} */
    xbar[5] = 10.512835408093998524;
    xbar[35] = 14.207392458842459604;
/* .... roots of j_{7+1/2} */
    xbar[6] = 11.657032192516371794;
    xbar[36] = 15.431289210268378298;
/* .... roots of j_{8+1/2} */
    xbar[7] = 12.790781711972119439;
    xbar[37] = 16.641002881512189759;
/* .... roots of j_{9+1/2} */
    xbar[8] = 13.915822610504896772;
    xbar[38] = 17.838643199205325374;
/* .... roots of j_{10+1/2} */
    xbar[9] = 15.033469303743437706;
    xbar[39] = 19.025853536127758758;
/* .... roots of j_{11+1/2} */
    xbar[10] = 16.14474294230134177;
    xbar[40] = 20.203942632811727975;
/* .... roots of j_{12+1/2} */
    xbar[11] = 17.25045478412596367;
    xbar[41] = 21.373972181162748996;
/* .... roots of j_{13+1/2} */
    xbar[12] = 18.351261495947273517;
    xbar[42] = 22.536817071119838118;
/* .... roots of j_{14+1/2} */
    xbar[13] = 19.447703108094735569;
    xbar[43] = 23.693208037471379157;
/* .... roots of j_{15+1/2} */
    xbar[14] = 20.540229825048211154;
    xbar[44] = 24.843762597586348306;
/* .... roots of j_{16+1/2} */
    xbar[15] = 21.62922143659035612;
    xbar[45] = 25.989007976413589063;
/* .... roots of j_{17+1/2} */
    xbar[16] = 22.71500167497224254;
    xbar[46] = 27.129398412300741228;
/* .... roots of j_{18+1/2} */
    xbar[17] = 23.797849034128304879;
    xbar[47] = 28.265328436642370491;
/* .... roots of j_{19+1/2} */
    xbar[18] = 24.878005058207190103;
    xbar[48] = 29.397143213440884324;
/* .... roots of j_{20+1/2} */
    xbar[19] = 25.955680785040136982;
    xbar[49] = 30.525146695246320405;
/* .... roots of j_{21+1/2} */
    xbar[20] = 27.031061821350016904;
    xbar[50] = 31.649608132509527699;
/* .... roots of j_{22+1/2} */
    xbar[21] = 28.10431238769707582;
    xbar[51] = 32.770767324195762804;
/* .... roots of j_{23+1/2} */
    xbar[22] = 29.175578576919045296;
    xbar[52] = 33.888838894135169113;
/* .... roots of j_{24+1/2} */
    xbar[23] = 30.24499100461545212;
    xbar[53] = 35.004015804722693939;
/* .... roots of j_{25+1/2} */
    xbar[24] = 31.312666984322407713;
    xbar[54] = 36.116472267412190614;
/* .... roots of j_{26+1/2} */
    xbar[25] = 32.378712327200142341;
    xbar[55] = 37.226366171563002225;
/* .... roots of j_{27+1/2} */
    xbar[26] = 33.443222842246179495;
    xbar[56] = 38.333841125320523187;
/* .... roots of j_{28+1/2} */
    xbar[27] = 34.506285595548477754;
    xbar[57] = 39.439028181452133264;
/* .... roots of j_{29+1/2} */
    xbar[28] = 35.567979974076099089;
    xbar[58] = 40.542047305426805018;
/* .... roots of j_{30+1/2} */
    xbar[29] = 36.628378589713436836;
    xbar[59] = 41.643008631132493502;
    i__1 = *lmax;
    for (l = 1; l <= i__1; ++l) {
	diff = xbar[l + 29] - xbar[l - 1];
	terme = xbar[l + 29];
	rac[(l - 1) * *nmax + 1] = xbar[l - 1];
	rac[(l - 1) * *nmax + 2] = xbar[l + 29];
	i__2 = *nmax;
	for (i__ = 3; i__ <= i__2; ++i__) {
	    a = terme + pi * .5;
	    terme += diff;
	    b = terme;
	    racdi_(&l, &a, &b, &xn, &kmax, &ier);
	    rac[(l - 1) * *nmax + i__] = xn;
	    diff = rac[(l - 1) * *nmax + i__] - rac[(l - 1) * *nmax + i__ - 1]
		    ;
/* .... Test to ensure that xn != NAN */
	    if (xn != xn) {
		s_wsle(&io___12);
		do_lio(&c__9, &c__1, "in zer_bsj", (ftnlen)10);
		e_wsle();
		s_wsle(&io___13);
		do_lio(&c__9, &c__1, "rac(", (ftnlen)4);
		i__3 = (l - 1) * *nmax + i__;
		do_lio(&c__3, &c__1, (char *)&i__3, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, ") =", (ftnlen)3);
		do_lio(&c__5, &c__1, (char *)&xn, (ftnlen)sizeof(doublereal));
		e_wsle();
		s_wsle(&io___14);
		do_lio(&c__9, &c__1, "lamda = ", (ftnlen)8);
		do_lio(&c__3, &c__1, (char *)&l, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " ", (ftnlen)1);
		do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " ith root", (ftnlen)9);
		e_wsle();
		s_stop("", (ftnlen)0);
	    }
	}
    }
    return 0;
} /* zero_bsj__ */

