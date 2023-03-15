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

#include <f2c.h>


/* Table of constant values */

//static integer c__200 = 200;
//static integer c__100 = 100;
static integer c__48 = 48;
static integer c__20 = 20;
static integer c__96 = 96;
static integer c__30 = 30;
static integer c__10000 = 10000;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;
static doublereal c_b91 = 2.;
static doublereal c_b101 = -1.;
static doublereal c_b104 = -2.;

//static doublereal cnp[20302];
static doublereal /*fact[101], */ rlag[256], xrac[3000000];
static doublereal /* dfact[101], */ xbar[64];
static doublereal xh2[64], xh3[64], wlag[256];
static doublereal  xh[64], xbar2[64], xbar3[64];




extern double boost_choose(unsigned n, unsigned k);

extern double boost_factorial(unsigned  n);
extern double boost_double_factorial(unsigned  n);


void
initarrays()
{
    static int done = 0;
  //  extern /* Subroutine */ int combinatorial_(integer *, doublereal *);
   // extern /* Subroutine */ int factorial_(integer *, doublereal *);
  //  extern /* Subroutine */ int dfactorial_(integer *, doublereal *);
    extern /* Subroutine */ int gaussleg_(integer *, doublereal *, doublereal
    *), zero_bsj__(integer *, integer *, doublereal *);
    extern /* Subroutine */ int  gausslag_(integer *, doublereal *, doublereal *);


    if ( !done ) {
    //    combinatorial_(&c__200, cnp);
      //  factorial_(&c__100, fact);
    //    dfactorial_(&c__100, dfact);
        gaussleg_(&c__48, xbar, xh);
        gaussleg_(&c__20, xbar2, xh2);
        gaussleg_(&c__48, xbar3, xh3);
        gausslag_(&c__96, rlag, wlag);
        zero_bsj__(&c__30, &c__10000, xrac);
        done = 1;
    }
}


/*     Last change:  H     2 Mar 2009    2:35 pm */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC    Programmer: Hassan Safouhi                                    CCC */
/* CC                                                                  CCC */
/* CC   This program evaluates the complete expression of four-center  CCC */
/* CC   two-electron Coulomb integral over STF functions:              CCC */
/* CC                                                                  CCC */
/* CC     <X1(r1-OA) X3(r2-OC) | 1/|r1 - r2| |X2(r1-OB) X4(r2-OD)>     CCC */
/* CC                                                                  CCC */
/* CC   The approaches used in this program are :                      CCC */
/* CC   - The infinite series with spherical Bessel function           CCC */
/* CC   - The infinite series with the sine function                   CCC */
/* CC      (S\overline{D} method)                                      CCC */
/* CC   - \overline{D} transformation with 4th order differential      CCC */
/* CC      equation                                                    CCC */
/* CC   - \overline{D} transformation with 2nd order differential      CCC */
/* CC      equation                                                    CCC */
/* CC   - \overline{D} transformation with 2nd order differential      CCC */
/* CC      equation using the W-algorithm of Sidi                      CCC */
/* CC   - S\overline{D} approach                                       CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
int fortran_fourc(    integer np1, integer l1, integer m1, doublereal zeta1,
                                  integer np2, integer l2, integer m2, doublereal zeta2,
                                  integer np3, integer l3, integer m3, doublereal zeta3,
                                  integer np4, integer l4, integer m4, doublereal zeta4,
                                  doublereal xa, doublereal ya, doublereal za,
                                  doublereal xb, doublereal yb, doublereal zb,
                                  doublereal xc, doublereal yc, doublereal zc,
                                  doublereal xd, doublereal yd, doublereal zd,
                                  doublereal *sdbar_res, doublereal *wgrep_res /* out vars */
)
{
    /* Initialized data */

    static doublereal pi = 3.14159265358979324;
    static doublereal tiny = 1e-150;
    static doublereal eps = 1e-15;

    /* Format strings */
    static char fmt_42[] = "(\002 SDbar32_lrzv: \002,d20.10,3(1x,d10.4),1x,d"
	    "10.4)";
    static char fmt_43[] = "(\002 GREP-W32    : \002,d20.10,3(1x,d10.4),1x,d"
	    "10.4)";
    static char fmt_51[] = "(4(3(i2,\002&\002),f5.2,\002&\002),3(f10.5,\002"
	    "&\002,f8.3,\002&\002,f8.3,\002&\002),d21.10,\002&\002,d10.4,\002&"
	    "\002,f8.3,\002\\ \002)";
    static char fmt_24[] = "(75(\002-\002))";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, 
	    i__11, i__12, i__13, i__14, i__15, i__16, i__17, i__18, i__19, 
	    i__20, i__21, i__22, i__23, i__24;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer /*s_rsle(cilist *),*/ do_lio(integer *, integer *, char *, ftnlen),
	    e_rsle(void);
    double sqrt(doublereal), pow_di(doublereal *, integer *), cos(doublereal),
	     sin(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsle(cilist *), e_wsle(void);

    /* Local variables */
    static integer delta_l12__, delta_l34__;
    static doublereal cste_l1234__, timecode;
    extern /* Subroutine */ int cpu_time__(doublereal *), gausslag_(integer *,
	     doublereal *, doublereal *);
    static doublereal cste_llp__;
    extern /* Subroutine */ int gaussleg_(integer *, doublereal *, doublereal 
	    *), zero_bsj__(integer *, integer *, doublereal *), constlmp_(
	    integer *, integer *,  doublereal *);
    static integer k, l, m;
    doublereal s, t, v[4096];
    extern /* Subroutine */ int
	    coeffbons_(integer *, integer *, integer *, integer *, doublereal *);
    static integer n1, n2, n3, n4;
    static doublereal cstmprod12, cstmprod34, ab, a12, ad, b12, cd, a34, b34;
    static integer j12, jc, l12, m12;
    static doublereal dk;
    static integer j34, l34, m34, ni;
    static integer mp, lambda_min__, lambda_max__, mu, lp, nx, ns, nt;
    //extern /* Subroutine */ int dfactorial_(integer *, doublereal *);
    doublereal ap1[22], ap2[22], ap3[22], ap4[22];
    static integer lp1, mp1;

    extern /* Subroutine */ int harmonique_(doublereal *, integer *,
	    doublereal *);
    static integer lp2, mp2, lp3;
    static doublereal tstartcode;
    static integer mp3, lp4, mp4;
    static doublereal rac[100001];
    static integer ng12, iii, ng34;
    static integer nu12, nsd, nwd, nu34;
    static doublereal tsd, twd;
    static integer nst;
    static doublereal cnp1;
    static integer ngt1, ngt2, ngt3, ngt4, ngt5;
    extern /* Subroutine */ int sdbar42_lrzv__(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, /* doublereal * ,*/  doublereal *,
	     integer *, doublereal *);
    static doublereal cste;
    static integer lmin;
    static doublereal cosi;
    static integer lmax;
    doublereal sine, phiv[4096], ttsd;
    extern /* Subroutine */ int vharmonique2_(integer *, integer *, 
	    doublereal *, doublereal *);
    static doublereal pmun, ttwd, ylmv[26873857],
	    tfinishtcode,  phiab,
	    phiad, phicd, facts, ylmab[6561];
    static integer l12min, l12max;
    doublereal ylmcd[6561];
    static integer l34min, l34max;
    static doublereal valsd;
    extern /* Subroutine */ int gaunt_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *);
    static integer lpmin, lpmax;
    static doublereal valwd, const__;
    static integer nylmv;
    doublereal coulbstf42_isd__, coulbstf42_iwd__, argnt1[41], argnt2[
	    41], argnt3[41], argnt4[41], argnt5[81], cnorm1;
    integer np1min, np2min, np1max, np2max, np3min, np3max, np4min,
	    np4max;
    doublereal cnorm2, cnorm3, cnorm4, cstmp1, cstmp2, cstmp3, cstmp4;
    integer lambda;
    doublereal coulbstf42_rsd__=0, coulbstf42_rwd__=0, coulbstf42_rss__=0,
	    coulbstf42_iss__=0, thetab, thetad;
    extern /* Subroutine */ int cartco_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal thetcd;
    static doublereal cste_j__;
    static integer nylmab, nylmcd;
    doublereal thetav[4096], cnormp;
    extern /* Subroutine */ int vcoord3_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    doublereal cstlmp1[442], cstlmp2[442], cstlmp3[442], cstlmp4[442],
	    cste_lb__, cste_ap__, err_isd__, err_iwd__;
    extern /* Subroutine */ int grepwd42_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    ;
    static doublereal err_rsd__, err_rwd__;

    /* Fortran I/O blocks */
    static cilist io___212 = { 0, 6, 0, fmt_42, 0 };
    static cilist io___213 = { 0, 6, 0, fmt_43, 0 };
    static cilist io___214 = { 0, 12, 0, fmt_51, 0 };
    static cilist io___215 = { 0, 13, 0, fmt_51, 0 };
    static cilist io___216 = { 0, 6, 0, 0, 0 };
    static cilist io___217 = { 0, 6, 0, fmt_24, 0 };
    static cilist io___220 = { 0, 6, 0, 0, 0 };


/* .... MAX_REP_SS : Maximal number of iterations for the subroutine using the infinite series */
/* .... MAX_REP_UL :   Maximal number of iterations for the subroutine using the Levin's u transform */
/* .... MAX_REP_EP : Maximal number of iterations for the subroutine using the epsilon algorithm of Wynn */
/* .... MAX_REP_UR :   Maximal number of iterations for the subroutine using the Levin's u recurrence */
/* .... MAX_REP_ER : Maximal number of iterations for the subroutine using the epsilon recurrence */
/* .... MAX_REP_HD :  Maximal number of iterations for the subroutine using H\bar{D} */
/* .... MAX_REP_DB :   Maximal number of iterations for the subroutine using \bar{D} */
/* .... MAX_REP_SD :  Maximal number of iterations for the subroutine using S\bar{D} */
/* .... MAX_REP_WD : Maximal number of iterations for the subroutine using H\bar{D} and the W-algorithm of Sid
i */
/* .... MAX_REP_SC : Maximal number of iterations for the subroutine using Sin-Cos-S\bar{D} */
/* .... MAX_REP_TC : Maximal number of iterations for the subroutine using Sin-Cos-S\bar{D} */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* CC                                                                  CCC */
/* CC  Here we define the maxima of all parameters used in the code.   CCC */
/* CC                                                                  CCC */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* .... MAX_N : maximal value of the principal quantum number n */
/* .... MAX_L : maximal value of the quantum number l */
/* .... For the Four-center Two-electron Coulomb integral */
/* .... MAX_bsj : Maximum order of spherical Bessel functions in zero_bsj (lambda = 30) */
/* .... MAX_Fact : Size of arrays Fact, Dfact */
/* .... MAX_Cnp : Size of arrays Cnp */
/* .... MAX_LEG : Maximum order of the Gauss-Legendre quadrature */
/* .... MAX_RAC : Maximum order of the Gauss-Legendre quadrature */
/* .... MAX_LAG : Maximum order of the Gauss-Laguerre quadrature */
/* .... For the Four-center Two-electron Coulomb integral */
/* .... MAX_DEV : Maximum number of terms in infinite series */
/* .... MAX_DIV : Maximal number of subdivisions of the finite intervals when using Gauss-Legendre */
/* .... MAX_DB : Maximum order for Dbar */
/* .... MAX_HD : Maximum order for HDbar equal to MAX_DB for LU to operate correctly */
/* .... MAX_UL : maximal order of Levin u */
/* .... MAX_EP : maximal order of epsilon algorithm */
/* .... MAX_UR : maximal order of Levin u using recurrence formulae */
/* .... MAX_ER : maximal order of Epsilon using recurrence formulae */
/* .... MAX_SD : Maximum order for SDbar */
/* .... MAX_WD  : Maximum order for GREP-W */
/* .... Here starts Coulomb42 */
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
    cpu_time__(&tstartcode);
/* .....initialization of arrays Cnp, Fact, Dfact, zero_bsj */
    initarrays();

    i__1 = 1;

    for (iii = 1; iii <= i__1; ++iii) {
/* .... Read the parameters for the integral */


	    cartco_(&xa, &ya, &za, &xb, &yb, &zb, &ab, &thetab, &phiab);
	    cartco_(&xa, &ya, &za, &xd, &yd, &zd, &ad, &thetad, &phiad);
	    cartco_(&xc, &yc, &zc, &xd, &yd, &zd, &cd, &thetcd, &phicd);
/* .... initialization of the series */
	coulbstf42_isd__ = 0.;
	coulbstf42_iwd__ = 0.;
/* .... Initialization of counters */
	ttsd = 0.;
	ttwd = 0.;
/* .... Decomposition of STFs in terms of B functions. */
	coeffbons_(&np1, &l1, &np1min, &np1max,  ap1);
	coeffbons_(&np2, &l2, &np2min, &np2max,  ap2);
	coeffbons_(&np3, &l3, &np3min, &np3max,  ap3);
	coeffbons_(&np4, &l4, &np4min, &np4max,  ap4);
/* .... evaluation of factors from the multiplication theorem of the */
/* .... solid spherical harmonics. These coefficients are common */
/* .... to all the B functions present in the integrals over Slater functions. */
	constlmp_(&l1, &m1,  cstlmp1);
	constlmp_(&l2, &m2,  cstlmp2);
	constlmp_(&l3, &m3,  cstlmp3);
	constlmp_(&l4, &m4,  cstlmp4);
/* .... Spherical Haromincs */
	i__2 = l1 + l2 + l3 + l4;
	harmonique_(&thetab, &i__2, ylmab);
	i__2 = l1 + l2 + l3 + l4;
	harmonique_(&thetcd, &i__2, ylmcd);
/* .... array containing the coordinates of \vect{v} */
	vcoord3_(&ab, &thetab, &phiab, &cd, &thetcd, &phicd, &ad, &thetad, &
		phiad, &c__48, xbar3, v, thetav, phiv);
/* .... array of spherical harmonics Ylmv. */
/* .... the maximum value for lambda is l1+l2+l3+l4 */
	i__2 = l1 + l2 + l3 + l4;
	vharmonique2_(&i__2, &c__48, thetav, ylmv);
/* .... Normalization constants */
	cnorm1 = pow_di(&c_b91, &np1) * sqrt(zeta1 * 2. / boost_factorial(np1 * 2)) *
		zeta1;
	cnorm2 = pow_di(&c_b91, &np2) * sqrt(zeta2 * 2. / boost_factorial(np2 * 2)) *
		zeta2;
	cnorm3 = pow_di(&c_b91, &np3) * sqrt(zeta3 * 2. / boost_factorial(np3 * 2)) *
		zeta3;
	cnorm4 = pow_di(&c_b91, &np4) * sqrt(zeta4 * 2. / boost_factorial(np4 * 2)) *
		zeta4;
	cnormp = cnorm1 * cnorm2 * cnorm3 * cnorm4;
/* ...  Constante d'integration */
/* Computing 5th power */
	d__1 = pi, d__2 = d__1, d__1 *= d__1;
	i__2 = l1 + l2 + l3 + l4 + 2;
	i__3 = l1 - 1;
	i__4 = l2 - 1;
	i__5 = l3 - 1;
	i__6 = l4 - 1;
	const__ = d__2 * (d__1 * d__1) * 8192. / pow_di(&c_b91, &i__2) *
            boost_double_factorial((l1 << 1) + 1) * boost_double_factorial((l2 << 1) + 1) * boost_double_factorial((l3 << 1)
		+ 1) * boost_double_factorial((l4 << 1) + 1) * pow_di(&zeta1, &i__3) * pow_di(&
		zeta2, &i__4) * pow_di(&zeta3, &i__5) * pow_di(&zeta4, &i__6);
	i__2 = l1;
	for (lp1 = 0; lp1 <= i__2; ++lp1) {
	    i__3 = lp1;
	    for (mp1 = -lp1; mp1 <= i__3; ++mp1) {
		i__4 = l2;
		for (lp2 = 0; lp2 <= i__4; ++lp2) {
		    i__5 = lp2;
		    for (mp2 = -lp2; mp2 <= i__5; ++mp2) {
/* .... evaluation of: Constlmp */
			cstmp1 = cstlmp1[lp1 * (lp1 + 1) + mp1];
			cstmp2 = cstlmp2[lp2 * (lp2 + 1) + mp2];
			cstmprod12 = cstmp1 * cstmp2;
			if (abs(cstmprod12) > tiny) {
/* .... evaluation of: <l2-lp2 m2-mp2| l1-lp1 m1-mp1 | l34 m34> */
			    i__6 = l1 - lp1;
			    i__7 = m1 - mp1;
			    i__8 = l2 - lp2;
			    i__9 = m2 - mp2;
			    gaunt_(&i__6, &i__7, &i__8, &i__9, &l12min, &
				    l12max, &m12, &ngt1, argnt1);
/* .... evaluation of: lp2 mp2|lp1 mp1|l m> */
			    gaunt_(&lp1, &mp1, &lp2, &mp2, &lmin, &lmax, &m, &
				    ngt2, argnt2);
			    i__6 = l3;
			    for (lp3 = 0; lp3 <= i__6; ++lp3) {
				i__7 = lp3;
				for (mp3 = -lp3; mp3 <= i__7; ++mp3) {
				    i__8 = l4;
				    for (lp4 = 0; lp4 <= i__8; ++lp4) {
					i__9 = lp4;
					for (mp4 = -lp4; mp4 <= i__9; ++mp4) {
/*     factors evaluated using Constlmp */
					    cstmp3 = cstlmp3[lp3 * (lp3 + 1) 
						    + mp3];
					    cstmp4 = cstlmp4[lp4 * (lp4 + 1) 
						    + mp4];
					    cstmprod34 = cstmp3 * cstmp4;
					    if (abs(cstmprod34) > tiny) {
/* .... evaluation of:  <l3-lp3 m3-mp3| l4-lp4 m4-mp4 | l34 m34> */
			  i__10 = l4 - lp4;
			  i__11 = m4 - mp4;
			  i__12 = l3 - lp3;
			  i__13 = m3 - mp3;
			  gaunt_(&i__10, &i__11, &i__12, &i__13, &l34min, &
				  l34max, &m34, &ngt3, argnt3);
/* .... evaluation of:  <lp4 mp4|lp3 mp3|lp mp> */
			  gaunt_(&lp3, &mp3, &lp4, &mp4, &lpmin, &lpmax, &mp, 
				  &ngt4, argnt4);
			  i__10 = l12max;
			  for (l12 = l12min; l12 <= i__10; l12 += 2) {
			      i__11 = l34max;
			      for (l34 = l34min; l34 <= i__11; l34 += 2) {
				  ngt1 = (l12 - l12min) / 2 + 1;
				  ngt3 = (l34 - l34min) / 2 + 1;
				  cste_l1234__ = argnt1[ngt1 - 1] * argnt3[
					  ngt3 - 1];
/* .... evaluation of:  <l12 m12 | l34 m34 | lambda mu> */
				  gaunt_(&l34, &m34, &l12, &m12, &
					  lambda_min__, &lambda_max__, &mu, &
					  ngt5, argnt5);
				  i__12 = lmax;
				  for (l = lmin; l <= i__12; l += 2) {
				      i__13 = lpmax;
				      for (lp = lpmin; lp <= i__13; lp += 2) {
					  ngt2 = (l - lmin) / 2 + 1;
					  ngt4 = (lp - lpmin) / 2 + 1;
					  nylmab = l * (l + 1) + m;
					  nylmcd = lp * (lp + 1) + mp;
					  d__1 = -cd;
					  cste_llp__ = pow_di(&ab, &l) * 
						  argnt2[ngt2 - 1] * ylmab[
						  nylmab] * pow_di(&d__1, &lp)
						   * argnt4[ngt4 - 1] * ylmcd[
						  nylmcd];
					  if (abs(cste_llp__) > tiny) {
			i__14 = lambda_max__;
			for (lambda = lambda_min__; lambda <= i__14; lambda +=
				 2) {
			    ngt5 = (lambda - lambda_min__) / 2 + 1;
			    ni = l1 + lp1 + l2 + lp2 + l3 + lp3 + l4 + lp4 + 
				    lambda;
			    i__15 = ni / 2 + l1 + l2 + lp2 + lp4 + lambda;
			    pmun = pow_di(&c_b101, &i__15);
			    cste_lb__ = pmun * argnt5[ngt5 - 1];
/* .... Zeros of the spherical Bessel function */
			    rac[0] = 0.;
			    if (lambda == 0) {
				for (k = 1; k <= 10000; ++k) {
				    dk = (doublereal) k;
				    rac[k] = dk * pi;
				}
			    } else {
				for (k = 1; k <= 10000; ++k) {
				    rac[k] = xrac[(lambda - 1) * 10000 + k - 
					    1];
				}
			    }
			    i__15 = np1max;
			    for (n1 = np1min; n1 <= i__15; ++n1) {
				i__16 = np2max;
				for (n2 = np2min; n2 <= i__16; ++n2) {
				    i__17 = np3max;
				    for (n3 = np3min; n3 <= i__17; ++n3) {
					i__18 = np4max;
					for (n4 = np4min; n4 <= i__18; ++n4) {
		      cste_ap__ = ap1[n1] * ap2[n2] * ap3[n3] * ap4[n4];
		      i__19 = n1 << 1;
		      i__20 = n2 << 1;
		      i__21 = n3 << 1;
		      i__22 = n4 << 1;
		      i__23 = n1 + n2 + n3 + n4;
		      cste = boost_factorial(n1 + l1 + n2 + l2 + 1) * boost_factorial(n3 + l3 + n4
			      + l4 + 1) / (boost_factorial(n1 + l1) * boost_factorial(n2 + l2) *
                                boost_factorial(n3 + l3)* boost_factorial(n4 + l4)) * pow_di(&zeta1,
			      &i__19) * pow_di(&zeta2, &i__20) * pow_di(&
			      zeta3, &i__21) * pow_di(&zeta4, &i__22) / 
			      pow_di(&c_b91, &i__23);
		      delta_l12__ = (lp1 + lp2 - l) / 2;
		      delta_l34__ = (lp3 + lp4 - lp) / 2;
		      i__19 = delta_l12__;
		      for (j12 = 0; j12 <= i__19; ++j12) {
			  i__20 = delta_l34__;
			  for (j34 = 0; j34 <= i__20; ++j34) {
			   //   cnp1 = cnp[delta_l12__ * (delta_l12__ + 1) / 2
				 //     + j12] * cnp[delta_l34__ * (delta_l34__
				   //   + 1) / 2 + j34];
                   cnp1 = boost_choose( delta_l12__ , j12) * boost_choose(delta_l34__, j34) ;
			      i__21 = j12 + j34;
			      cste_j__ = pow_di(&c_b104, &i__21) * cnp1 / (
                          boost_factorial(n1 + n2 + l1 + l2 - j12 + 1) *
                          boost_factorial(n3 + n4 + l3 + l4 - j34 + 1));
/* .... Determinationof the integration parameters: */
			      nx = l1 - lp1 + l2 - lp2 + l3 - lp3 + l4 - lp4;
			      nu12 = n1 + n2 + l1 + l2 - l - j12;
			      nu34 = n3 + n4 + l3 + l4 - lp - j34;
			      ng12 = ((n1 + l1 + n2 + l2) << 1) - (lp1 + lp2 +
				      l) + 1;
			      ng34 = ((n3 + l3 + n4 + l4) << 1) - (lp3 + lp4 +
				      lp) + 1;
			      for (ns = 1; ns <= 48; ++ns) {
				  for (nt = 1; nt <= 48; ++nt) {
				      s = (xbar3[ns - 1] + 1.) * .5;
/* Computing 2nd power */
				      d__1 = zeta1;
/* Computing 2nd power */
				      d__2 = zeta2;
				      a12 = (1. - s) * (d__1 * d__1) + s * (
					      d__2 * d__2);
				      b12 = s * (1. - s);
				      t = (xbar3[nt - 1] + 1.) * .5;
/* Computing 2nd power */
				      d__1 = zeta4;
/* Computing 2nd power */
				      d__2 = zeta3;
				      a34 = (1. - t) * (d__1 * d__1) + t * (
					      d__2 * d__2);
				      b34 = t * (1. - t);
/*     R&I = e^{i (m*phiab + mp*phicd + mu*phiv) } */
				      nst = (ns - 1) * 48 + nt;
				      cosi = cos(m * phiab + mp * phicd + mu *
					       phiv[nst - 1]);
				      sine = sin(m * phiab + mp * phicd + mu *
					       phiv[nst - 1]);
/*     Control of accuracy */
				      if (abs(cosi) < eps) {
					  cosi = 0.;
				      }
				      if (abs(sine) < eps) {
					  sine = 0.;
				      }
				      nylmv = nst - 1 + (lambda * (lambda + 1)
					       + mu) * 2304;
				      if (v[nst - 1] > tiny && ylmv[nylmv] != 
					      0.) {
					  i__21 = n2 + l2 + l1 - lp1;
					  d__1 = 1. - s;
					  i__22 = n1 + l1 + l2 - lp2;
					  i__23 = n3 + l3 + l4 - lp4;
					  d__2 = 1. - t;
					  i__24 = n4 + l4 + l3 - lp3;
					  facts = xh3[ns - 1] * .25 * xh3[nt 
						  - 1] * pow_di(&s, &i__21) * 
						  pow_di(&d__1, &i__22) * 
						  pow_di(&t, &i__23) * pow_di(
						  &d__2, &i__24) * ylmv[nylmv]
						   * cste_l1234__ * 
						  cstmprod12 * cstmprod34 * 
						  cste_llp__ * cste_lb__ * 
						  cste_j__ * cste * cste_ap__ 
						  * const__ * cnormp;
/* .... S\bar{D}: rec lim */
					  jc = 0;
					  sdbar42_lrzv__(&nx, &lambda, &v[nst 
						  - 1], &nu12, &ng12, &ab, &
						  a12, &b12, &nu34, &ng34, &
						  cd, &a34, &b34, &jc, &c__20,
						   xbar2, xh2, &c__96, rlag, 
						  wlag,
                          //cnp,
                          &valsd, &nsd, &tsd);
					  coulbstf42_rsd__ += facts * valsd * 
						  cosi;
					  coulbstf42_isd__ += facts * valsd * 
						  sine;
					  ttsd += tsd;
/* .... GREP-WD */
					  grepwd42_(&nx, &lambda, &v[nst - 1],
						   &nu12, &ng12, &ab, &a12, &
						  b12, &nu34, &ng34, &cd, &
						  a34, &b34, &c__20, xbar2, 
						  xh2, &c__96, rlag, wlag, 
						  rac, &valwd, &nwd, &twd);
					  coulbstf42_rwd__ += facts * valwd * 
						  cosi;
					  coulbstf42_iwd__ += facts * valwd * 
						  sine;
					  ttwd += twd;
				      }
/* L40: */
				  }
			      }
/* L35: */
			  }
		      }
/* L30: */
					}
				    }
				}
			    }
			}
					  }
/* L25: */
				      }
				  }
/* L20: */
			      }
			  }
					    }
/* L15: */
					}
				    }
				}
			    }
			}
/* L10: */
		    }
		}
	    }
	}

    *sdbar_res =   coulbstf42_rsd__;
    *wgrep_res = coulbstf42_rwd__;

/* .... Erreurs */
	err_rsd__ = (d__1 = coulbstf42_rsd__ - coulbstf42_rss__, abs(d__1));
	err_rwd__ = (d__1 = coulbstf42_rwd__ - coulbstf42_rss__, abs(d__1));
	err_isd__ = (d__1 = coulbstf42_isd__ - coulbstf42_iss__, abs(d__1));
	err_iwd__ = (d__1 = coulbstf42_iwd__ - coulbstf42_iss__, abs(d__1));
/* .... Results */
	s_wsfe(&io___212);
	do_fio(&c__1, (char *)&coulbstf42_rsd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&err_rsd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&coulbstf42_isd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&err_isd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ttsd, (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___213);
	do_fio(&c__1, (char *)&coulbstf42_rwd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&err_rwd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&coulbstf42_iwd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&err_iwd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ttwd, (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___214);
	do_fio(&c__1, (char *)&n1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&l1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&m1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&zeta1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&n2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&l2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&m2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&zeta2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&n3, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&l3, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&m3, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&zeta3, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&n4, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&l4, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&m4, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&zeta4, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ab, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&thetab, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&phiab, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ad, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&thetad, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&phiad, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&thetcd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&phicd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&coulbstf42_rsd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&err_rsd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ttsd, (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___215);
	do_fio(&c__1, (char *)&n1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&l1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&m1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&zeta1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&n2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&l2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&m2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&zeta2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&n3, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&l3, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&m3, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&zeta3, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&n4, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&l4, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&m4, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&zeta4, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ab, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&thetab, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&phiab, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ad, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&thetad, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&phiad, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&cd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&thetcd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&phicd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&coulbstf42_rwd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&err_rwd__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&ttwd, (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsle(&io___216);
	do_lio(&c__3, &c__1, (char *)&iii, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsfe(&io___217);
	e_wsfe();
/* L666: */
    }
/* .... FORMATS      FORMATS       FORMATS       FORMATS      FORMATS */
    cpu_time__(&tfinishtcode);
    timecode = tfinishtcode - tstartcode;
    s_wsle(&io___220);
    do_lio(&c__5, &c__1, (char *)&timecode, (ftnlen)sizeof(doublereal));
    e_wsle();
    return 0;
} /* MAIN__ */
