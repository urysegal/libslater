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

static integer c_n1 = -1;
static integer c__1 = 1;

/*     Last change:  S    14 Apr 2007   11:19 pm */
/* Subroutine */ int gaunt_(integer *l2, integer *m2, integer *l3, integer *
	m3, integer *l1min, integer *l1max, integer *m1, integer *ngaunt, 
	doublereal *argnt)
{
    /* Initialized data */

    static doublereal tiny = 1e-60;
    static doublereal srtiny = 1e-30;
    static doublereal srhuge = 1e30;
    static doublereal pi4 = 12.56637061435917295385057;

    /* Format strings */
    static char fmt_1000[] = "(\0020\002,\002 *** ERROR IN SUBROUTINE GAUNT "
	    "***\002/\002  \002,33(\002-\002))";
    static char fmt_1010[] = "(\0020\002,\002ILLEGAL ANGULAR MOMENTUM QUANTU"
	    "M NUMBERS WERE USED  S INPUT PARAMETERS\002/\0020\002,\002L2 ="
	    " \002,i4,\002  M2 = \002,i5,\002  L3 = \002,i4,\002  M3 = \002,i"
	    "5)";
    static char fmt_1020[] = "(\0020\002,\002  DIMENSION OF ARGNT NOT LARGE "
	    "ENOUGH TO STORE ALL  HE GAUNT COEFFICIENTS REQUIRED\002/\0020"
	    "\002,\002LARRAY = \002,i4,\002   NGAUNT = \002,i4)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer pow_ii(integer *, integer *), s_wsfe(cilist *), e_wsfe(void), 
	    do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;
    static doublereal c1, c2;
    static integer l1;
    static doublereal x1, x2, x3, y2, y3, y1;
    static integer la, lb, l12;
    static doublereal dl;
    static integer mm, lp1, l1p2, m1p2;
    static doublereal sum;
    static integer labm, lold, npar, lnew;
    static doublereal sump, c1old, c1new, ratio, cnorm;
    static integer min3jm;
    static doublereal facold, facnew;
    static integer lmatch;
    static doublereal sumbac;
    static integer minpar, igaunt, kgaunt;
    static doublereal sumfor;

    /* Fortran I/O blocks */
    static cilist io___43 = { 0, 4, 0, fmt_1000, 0 };
    static cilist io___44 = { 0, 4, 0, fmt_1010, 0 };
    static cilist io___45 = { 0, 4, 0, fmt_1000, 0 };
    static cilist io___46 = { 0, 4, 0, fmt_1020, 0 };


/* ********************************************************************** */

/*  PROGRAM:  SUBROUTINE GAUNT. */

/*  PROGRAMMER:  ERNST JOACHIM WENIGER. */

/*  DATE:  REGENSBURG, 28 04 1981. */


/*  DESCRIPTION: */
/*  ----------- */

/*  THIS SUBROUTINE CALCULATES A STRING OF GAUNT-COEFFICIENTS */

/*              ( L3 M3 / L2 M2 / L1 M1 ) */

/*  FOR ALL ALLOWED L1-VALUES, WHERE L2, M2, L3 AND M3 ARE FIXED INPUT */
/*  QUANTITIES. */

/*  THIS YIELDS THE FOLLOWING CONDITIONS FOR L1 AND M1: */

/*  M1 = M3 - M2. */

/*  L1 = L1MIN, L1MIN + 2, ... , L1MAX - 2 , L1MAX. */

/*  L1MAX = L2 + L3. */

/*  L1MIN = MAX0(IABS(L2-L3),IABS(M1)) */
/*          IF MAX0(IABS(L2-L3),IABS(M1)) + L1MAX IS EVEN AND */
/*  L1MIN = MAX0(IABS(L2-L3),IABS(M1)) + 1 */
/*          IF MAX0(IABS(L2-L3),IABS(M1)) + L1MAX IS ODD. */

/*  THE RESULTING GAUNT-COEFFICIENTS ( L3 M3 / L2 M2 / L1 M3-M2 ) ARE */
/*  STORED IN THE 1-DIMENSIONAL ARRAY ARGNT AS ARGNT((L1-L1MIN)/2 + 1), */
/*  I.E. THE ARRAY ARGNT HAS TO BE AT LEAST OF LENGTH */
/*  NGAUNT = (L1MAX - L1MIN)/2 + 1. */


/*  IN THIS PROGRAM THE GAUNT-COEFFICIENTS ARE CALCULATED USING THEIR */
/*  REPRESENTATION IN TERMS 3JM-SYMBOLS ( M.ROTENBERG ET AL., THE 3J- */
/*  AND 6J-SYMBOLS, P.5, EQ.(1.19) ), I.E. THE GAUNT COEFFICIENTS ARE */
/*  REPRESENTED AS THE PRODUCT OF TWO 3JM-SYMBOLS. ONE - USUALLY */
/*  CALLED THE PARITY COEFFICIENT - POSSESSES ONLY MAGNETIC QUANTUM */
/*  NUMBERS M1, M2 AND M3 WITH M1 = M2 = M3 = 0. THE OTHER CONTAINS */
/*  MAGNETIC QUANTUM NUMBERS WHERE IN THE GENERAL CASE AT LEAST TWO */
/*  M-VALUES ARE DIFFERENT FROM ZERO. */
/*  ACCORDINGLY TWO STRINGS OF 3JM-SYMBOLS HAVE TO BE CALCULATED. */
/*  IN EITHER CASE THE 3JM-SYMBOLS WILL BE CALCULATED RECURSIVELY. */

/*  THE 3JM-SYMBOLS WITH NON-VANISHING MAGNETIC QUANTUM NUMBERS ARE */
/*  CALCULATED WITH THE HELP OF A 3-TERM RECURSION EQUATION IN L1. FOR */
/*  A DISCUSSION OF THIS RECURSION EQUATION AND THE STABILITY PROBLEMS */
/*  ASSOCIATED WITH IT SEE K.SCHULTEN AND R.G.GORDON, J.MATH.PHYS.16, */
/*  1961-1970 (1975) AND IBID. 16,01971-1988 (1975). */

/*  FOR THE SAKE OF THE NUMERICAL STABILITY THE RECURSION WILL PROCEED */
/*  SIMULTANEOUSLY FORWARDS AND BACKWARDS, STARTING FROM L1MIN = */
/*  MAX0(IABS(L2-L3),IABS(M1)) AND L1MAX = L2 + L3, RESPECTIVELY. */
/*  THE TWO RECURSIONS PROCEED IN EITHER DIRECTION AS LONG AS NUMERICAL */
/*  STABILITY IS GUARANTEED. */

/*  AT THE UPPER AND THE LOWER BOUNDARY OF THE DOMAIN OF ALLOWED L1- */
/*  VALUES THE 3-TERM RECURSION EQUATION REDUCES TO A 2-TERM RECURSION */
/*  EQUATION. THEREFORE, FOR THE TWO 3JM-SYMBOLS AT THE BOUNDARY OF */
/*  THE L1-DOMAIN WE CAN CHOSE  ARBITRARY STARTING VALUES. THE TWO */
/*  RECURSION SERIES ARE MATCHED AT THREE INTERMEDIATE POINTS BY LEAST */
/*  SQUARE FIT. THE REMAINING INDETERMINACY CAN THEN BE REMOVED WITH */
/*  THE HELP OF THE NORMALIZING RELATIONSHIP OF THE 3JM-SYMBOLS. */

/*  IN ORDER TO MINIMIZE THE NUMBER OF SQUARE ROOT EVALUATIONS WE do NOT */
/*  COMPUTE THE PARITY COEFFICIENTS DIRECTLY BUT INSTEAD THE PARITY */
/*  COEFFICIENTS MULTIPLIED BY DSQRT(2*L1+1). THESE QUANTITIES ARE */
/*  CALCULATED USING A HOMOGENEOUS 2-TERM RECURSION EQUATION WHICH CAN */
/*  BE DERIVED FROM THE HOMOGENEOUS 3-TERM RECURSION EQUATION USED FOR */
/*  THE COMPUTATION OF THE GENERAL 3JM-SYMBOLS. */
/*  IT SHOULD BE NOTED THAT UNLIKE AS IN THE GENERAL CASE THE 2-TERM */
/*  RECURSION IS NUMMERICALLY STABLE IN EITHER DIRECTION. */
/*  THE RECURSION STARTS AT L1 = IABS(L2-L3) WITH AN ARBITRARY STARTING */
/*  VALUE. THE REMAINING INDETERMINACY IS THEN REMOVED WITH THE HELP OF */
/*  THE NORMALIZATION RELATIONSHIP OF THE PARITY COEFFICIENTS. */

/*  THE GAUNT COEFFICIENTS ARE THEN CALCULATED ACCORDING TO EQ.(1.19) */
/*  IN M.ROTENGERG ET AL. , THE 3J- AND 6J-SYMBOLS. */

/*  IF GAUNT-COEFFICIENTS OF THE TYPE ( L3 0 / L2 0 / L1 0 ) HAVE TO BE */
/*  CALCULATED WE USE THE SQUARED FORM OF THE 2-TERM RECURSION EQUATION */
/*  TO CALCULATE THE SQUARES OF THE 3JM-SYMBOLS WITH M1 = M2 = M3 = 0 */
/*  DIRECTLY. */


/*  LIST OF PARAMETERS: */
/*  ------------------ */

/*  INPUT PARAMETERS: */

/*  L2    : ANGULAR MOMENTUM QUANTUM NUMBER. */
/*  M2    : MAGNETIC QUANTUM NUMBER. */
/*  L3    : ANGULAR MOMENTUM QUANTUM NUMBER. */
/*  M3    : MAGNETIC QUANTUM NUMBER. */
/*  LARRAY: LENGTH OF THE ARRAY ARGNT. */

/*  OUTPUT PARAMETERS: */

/*  L1MIN : MINIMUM VALUE OF THE ANGULAR MOMENTUM QUANTUM NUMBER L1. */
/*  L1MAX : MAXIMUM VALUE OF THE ANGULAR MOMENTUM QUANTUM NUMBER L1. */
/*  M1    : MAGNETIC QUANTUM NUMBER. */
/*  NGAUNT: NUMBER OF GAUNT COEFFICIENTS CALCULATED AND STORED. */
/*  ARGNT : 1-DIMENSIONAL ARRAY CONTAINING THE GAUNT-COEFFICIENTS. */

/*  MACHINE-DEPendENT QUANTITIES: */

/*  HUGE,TINY,SRHUGE,SRTINY,PI4. */

/*  TINY SHOULD BE SET CLOSE TO BUT NOT IDENTICAL WITH THE SMALLEST */
/*  FLOATING POINT NUMBER THAT IS REPRESENTABLE ON THE COMPUTER */
/*  SRTINY IS THE SQUARE ROOT OF TINY. */

/*  HUGE SHOULD BE SET CLOSE TO BUT NOT IDENTICAL WITH THE LARGEST */
/*  FLOATING POINT NUMBER THAT IS REPRESENTABLE ON THE COMPUTER */
/*  SRHUGE IS THE SQUARE ROOT OF HUGE. */

/*  IN ADDITION IT IS ASSUMED THAT HUGE = TINY**(-1) HOLDS. */

/*  PI4 ( CORRESPONDING TO 4.D0*PI ). */

/* ********************************************************************** */
/* .... DIMENSION ARGNT(LARRAY) */

    /* Parameter adjustments */
    --argnt;

    /* Function Body */


/*  CHECK RELATIVE MAGNITUDES OF L- AND M-VALUES. */

    if (*l2 < abs(*m2)) {
	goto L360;
    }
    if (*l3 < abs(*m3)) {
	goto L360;
    }

/*  CALCULATE M1 AND DETERMINE THE LIMITS FOR L1. */

    *m1 = *m3 - *m2;
    *l1max = *l2 + *l3;
/* Computing MAX */
    i__2 = (i__1 = *l2 - *l3, abs(i__1)), i__3 = abs(*m1);
    *l1min = ((*l1max + max(i__2,i__3) + 1) / 2 << 1) - *l1max;

/*  DETERMINE THE NUMBER OF NON-VANISHING GAUNT-COEFFICIENTS. */

    *ngaunt = (*l1max - *l1min) / 2 + 1;

/*  CHECK THE WHETHER DIMENSION OF ARGNT IS LARGE ENOUGH TO STORE ALL TH */
/*  GAUNT-COEFFICIENTS REQUIRED. */

/* .... if(LARRAY.LT.NGAUNT)  GO TO 370 */

/*  CHECK WHETHER ALL ANGULAR MOMENTUM QUANTUM NUMBERS VANISH, */
/*  I.E. L1 = L2 = L3 = 0. */

    if (*l1max == 0) {
	goto L340;
    }

/*  CALCULATION OF SOME VARIABLES WHICH ASSUME CONSTANT VALUES DURING */
/*  THE RECURSIVE CALCULATION OF THE UNNORMALIZED 3JM-SYMBOLS. */

    la = (*l2 - *l3) * (*l2 - *l3);
    lb = (*l2 + *l3 + 1) * (*l2 + *l3 + 1);

/*  CHECK WHETER M1 = M2 = M3 = 0 HOLDS. */

    if (*m2 != 0) {
	goto L30;
    }
    if (*m3 != 0) {
	goto L30;
    }

/* ---------------------------------------------------------------------- */

/*  THIS IS REACHED IF M1 = M2 = M3 = 0 HOLDS. */

/*  THE RECURSION EQUATION USED HERE IS A 2-TERM RECURSION EQUATION WHIC */
/*  CAN BE OBTAINED BY SQUARING THE 2-TERM RECURSION SATISFIED BY THE */
/*  PARITY COEFFICIENTS. */

/*  INITIALIZATION OF THE ANGULAR MOMENTUM QUANTUM NUMBER L1. */

    l1 = *l1min;
    x1 = srtiny;
    dl = (doublereal) (l1 + l1 + 1);
    sum = x1 * dl;

/*  FOR THE SAKE OF EFFICIENCY THE SQUARES OF THE UNNORMALIZED PARITY */
/*  COEFFICIENTS ARE NOW MULTIPLIED BY THE SQUARE ROOT OF 2*L1 + 1. */

    argnt[1] = x1 * sqrt(dl);

/*  RECURSIVE CALCULATION OF THE SQUARES OF THE OTHER */
/*  UNNORMALIZED PARITY COEFFICIENTS. */

    i__1 = *ngaunt;
    for (igaunt = 2; igaunt <= i__1; ++igaunt) {

/*  CALCULATION OF SOME VARIABLES THAT ARE REQUIRED FOR THE RECURSIVE */
/*  CALCULATION OF THE NEXT SQUARED PARITY COEFFICIENT. */

	l1 += 2;
	l1p2 = l1 * l1;
	lp1 = (l1 - 1) * (l1 - 1);

/*  RECURSIVE CALCULATION OF THE NEXT SQUARED PARITY COEFFICIENT. */

	x1 = (doublereal) ((lp1 - la) * (lb - lp1)) / (doublereal) ((l1p2 - 
		la) * (lb - l1p2)) * x1;
	dl = (doublereal) (l1 + l1 + 1);
	sum += dl * x1;
	argnt[igaunt] = x1 * sqrt(dl);
/* L10: */
    }

/*  CALCULATION OF THE NORMALIZATION CONSTANT. */

    sum = sqrt((doublereal) ((*l2 + *l2 + 1) * (*l3 + *l3 + 1)) / pi4) / sum;

/*  NORMALIZE THE STRING OF GAUNT COEFFICIENTS ( L3 0 / L2 0 / L1 0 ) */

    i__1 = *ngaunt;
    for (igaunt = 1; igaunt <= i__1; ++igaunt) {
	argnt[igaunt] *= sum;
/* L20: */
    }
    return 0;

/*  end OF THE RECURSIVE CALCULATION OF THOSE GAUNT COEFFICIENTS WHERE */
/*  M1 = M2 = M3 = 0 HOLDS. */

/* ---------------------------------------------------------------------- */

/*  THIS IS REACHED IF NOT ALL OF THE MAGNETIC QUANTUM NUMBERS M1, M2 AN */
/*  M3 VANISH. */

L30:

/*  DETERMINE THE MINIMAL VALUE WHICH L1 CAN ASSUME DURING THE RECURSIVE */
/*  CALCULATION OF THE UNNORMALIZED 3JM-SYMBOLS WITH NON-VANISHING */
/*  MAGNETIC QUANTUM NUMBERS. */

/* Computing MAX */
    i__2 = (i__1 = *l2 - *l3, abs(i__1)), i__3 = abs(*m1);
    min3jm = max(i__2,i__3);

/*  CHECK WHETHER MIN3JM = L1MAX HOLDS, I.E. ONLY ONE 3JM-SYMBOL WITH */
/*  NON-VANISHING MAGNETIC QUANTUM NUMBERS HAS TO BE CALCULATED. */

    if (min3jm == *l1max) {
	goto L270;
    }

/*  CALCULATION OF SOME VARIABLES WHICH ASSUME CONSTANT VALUES DURING TH */
/*  RECURSIVE CALCULATION OF THE UNNORMALIZED 3JM-SYMBBOLS WITH NON- */
/*  VANISHING MAGNETIC QUANTUM NUMBERS. */

    labm = *m1 * (*l2 * (*l2 + 1) - *l3 * (*l3 + 1));
    mm = *m2 + *m3;
    m1p2 = *m1 * *m1;

/*  INITIALIZATION OF THE DUMMY VARIABLE L1. */

    l1 = min3jm;

/*  CHECK WHETHER MIN3JM = L1MIN HOLDS. */

    if (min3jm == *l1min) {
	goto L50;
    }

/*  THIS IS REACHED IF L1MIN = MIN3JM + 1 HOLDS. IN THAT CASE THE FIRST */
/*  UNNORMALIZED 3JM-SYMBOL CALCULATED WILL NOT BE STORED IN ARGNT. */

/*  CALCULATION OF THE FIRST UNNORMALIZED 3JM-SYMBOL */

    x3 = srtiny;
    l12 = l1 + l1 + 1;
    sumfor = (doublereal) l12 * tiny;

/*  CALCULATION OF SOME VARIABLES REQUIRED FOR THE RECURSIVE */
/*  CALCULATION OF THE SECOND UNNORMALIZED 3JM-SYMBOL. */

    ++l1;
    l1p2 = l1 * l1;
    facold = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));
    c2 = 1. / ((doublereal) (l1 - 1) * facold);
    c1 = c2 * (doublereal) (l12 * (labm + mm * l1 * (l1 - 1)));

/*  RECURSIVE CALCULATION OF THE SECOND UNNORMALIZED 3JM-SYMBOL. */

    x2 = c1 * x3;
    c1old = abs(c1);
    argnt[1] = x2;
    l12 += 2;
    sumfor += (doublereal) l12 * x2 * x2;

/*  CHECK WHETHER L1 = L1MAX HOLDS. */

    if (l1 == *l1max) {
	goto L260;
    }

/*  RECURSIVE CALCULATION OF THE NEXT 3JM-SYMBOL. */


/*  CALCULATION OF SOME COEFFICIENTS THAT ARE REQUIRED FOR THE RECURSIVE */
/*  CALCULATION OF THE NEXT 3JM-SYMBOL. */

    ++l1;
    l1p2 = l1 * l1;
    facnew = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));
    c2 = 1. / ((doublereal) (l1 - 1) * facnew);
    c1 = c2 * (doublereal) (l12 * (labm + mm * l1 * (l1 - 1)));

/*  RECURSIVE CALCULATION OF THE THIRD UNNORMALIZED 3JM-SYMBOL. */

    x1 = c1 * x2 - (doublereal) l1 * facold * c2 * x3;

/*  CHECK WHETHER THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */

    if (abs(x1) < srhuge) {
	goto L40;
    }

/*  THIS IS REACHED IF THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */
/*  IN THAT CASE THE 3JM-SYMBOLS ALREADY CALCULATED AND SUMFOR HAVE TO B */
/*  RESCALED IN ORDER TO PREVENT OVERFLOW. */

    x1 *= srtiny;
    x2 *= srtiny;
    argnt[1] *= srtiny;
    sumfor *= tiny;
L40:

/*  CHECK WHETHER DABS(C1) IS DECREASING. */
/*  AS LONG AS DABS(C1) IS DECREASING THE RECURSION PROCEEDS TOWARDS */
/*  INCREASING 3JM-VALUES AND IS NUMERICALLY STABLE. ONCE AN INCREASE */
/*  OF DABS(C1) IS DETECTED THE RECURSION DIRECTION IS REVERSED. */

    c1new = abs(c1);
    if (c1new >= c1old) {
	goto L130;
    }

/*  CALCULATION OF SUMFOR. */

    l12 += 2;
    sumfor += (doublereal) l12 * x1 * x1;
    goto L70;

/*  THIS REACHED IF L1MIN = MIN3JM HOLDS. */
/*  IN THAT CASE THE FIRST 3JM-SYMBOL CALCULATED IS ALSO NEEDED FOR THE */
/*  COMPUTATION OF THE FIRST GAUNT COEFFICIENT AND HAS TO BE STORED IN */
/*  THE ARRAY ARGNT. */

L50:

/*  CALCULATION OF THE FIRST UNNORMALIZED 3JM-SYMBOL. */

    x2 = srtiny;
    l12 = l1 + l1 + 1;
    sumfor = (doublereal) l12 * tiny;
    argnt[1] = x2;

/*  CHECK WHETHER THE RECURSION STARTS AT L1 = 0. */

    if (*l1min == 0) {
	goto L60;
    }

/*  RECURSIVE CALCULATION OF THE SECOND UNNORMALIZED 3JM-SYMBOL. */

    ++l1;
    l1p2 = l1 * l1;
    facnew = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));
    c1 = (doublereal) (l12 * (labm + mm * l1 * (l1 - 1))) / ((doublereal) (l1 
	    - 1) * facnew);
    x1 = c1 * x2;
    l12 += 2;
    sumfor += (doublereal) l12 * x1 * x1;
    c1new = abs(c1);
    goto L70;

/*  THIS IS REACHED IF L1MIN = 0 HOLDS. IN THAT CASE THE 2-TERM-RECURSIO */
/*  EQUATION WHICH HOLDS AT THE LOWER BOUNDARY OF THE ALLOWED L1-DOMAIN */
/*  HAS TO BE MODIFIED. */

L60:

/*  CALCULATION OF THE SECOND 3JM-SYMBOL. */

    l1 = 1;
    l12 = 3;
    c1 = (doublereal) (*m2) / sqrt((doublereal) (*l2 * (*l2 + 1)));
    x1 = c1 * x2;
    sumfor += x1 * 3. * x1;
    facnew = sqrt((doublereal) (lb - 1));
    c1new = abs(c1);
L70:

/*  INITIALIZATION OF THE DUMMY VARIABLE IGAUNT. */

    igaunt = 1;

/*  RECURSIVE CALCULATION OF THE REMAINING UNNORMALIZED 3JM-SYMBOLS. */

L80:

/*  CALCULATION OF SOME VARIABLES AND COEFFICIENTS THAT ARE REQUIRED FOR */
/*  THE RECURSIVE CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    facold = facnew;
    c1old = c1new;
    x3 = x2;
    x2 = x1;
    ++l1;
    l1p2 = l1 * l1;
    facnew = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));
    c2 = 1. / ((doublereal) (l1 - 1) * facnew);
    c1 = c2 * (doublereal) (l12 * (labm + mm * l1 * (l1 - 1)));
    l12 += 2;

/*  RECURSIVE CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    x1 = c1 * x2 - (doublereal) l1 * facold * c2 * x3;

/*  CHECK WHETHER THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */

    if (abs(x1) < srhuge) {
	goto L100;
    }

/*  THIS IS REACHED IF THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */
/*  IN THAT CASE THE 3JM-SYMBOLS ALREADY CALCULATED AND SUMFOR HAVE TO B */
/*  RESCALED IN ORDER TO PREVENT OVERFLOW. */

    x1 *= srtiny;
    x2 *= srtiny;
    x3 *= srtiny;
    i__1 = igaunt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	argnt[i__] *= srtiny;
/* L90: */
    }
    sumfor *= tiny;
L100:

/*  CHECK WHETHER L1 = L1MAX HOLDS. */

    if (l1 == *l1max) {
	goto L250;
    }

/*  CHECK WHETHER DABS(C1) IS DECREASING. */
/*  AS LONG AS DABS(C1) IS DECREASING THE RECURSION PROCEEDS TOWARDS */
/*  INCREASING 3JM-VALUES AND IS NUMERICALLY STABLE. ONCE AN INCREASE */
/*  OF DABS(C1) IS DETECTED THE RECURSION DIRECTION IS REVERSED. */

    c1new = abs(c1);
    if (c1new >= c1old) {
	goto L130;
    }
    ++igaunt;
    sumfor += (doublereal) l12 * x1 * x1;

/*  STORAGE OF THOSE UNNORMALIZED 3JM-SYMBOL WITH NON-VANISHING MAGNETIC */
/*  QUANTUM NUMBERS M1, M2 AND M3 THAT ARE REQUIRED FOR THE CALCULATION */
/*  OF THE GAUNT COEFFICIENTS IN THE ARRAY ARGNT. */

    argnt[igaunt] = x1;

/*  CALCULATION OF SOME VARIABLES AND COEFFICIENTS THAT ARE REQUIRED FOR */
/*  THE RECURSIVE CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    facold = facnew;
    c1old = c1new;
    x3 = x2;
    x2 = x1;
    ++l1;
    l1p2 = l1 * l1;
    facnew = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));
    c2 = 1. / ((doublereal) (l1 - 1) * facnew);
    c1 = c2 * (doublereal) (l12 * (labm + mm * l1 * (l1 - 1)));

/*  RECURSIVE CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    x1 = c1 * x2 - (doublereal) l1 * facold * c2 * x3;

/*  CHECK WHETHER THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */

    if (abs(x1) < srhuge) {
	goto L120;
    }

/*  THIS IS REACHED IF THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */
/*  IN THAT CASE THE 3JM-SYMBOLS ALREADY CALCULATED AND SUMFOR HAVE TO B */
/*  RESCALED IN ORDER TO PREVENT OVERFLOW. */

    x1 *= srtiny;
    x2 *= srtiny;
    x3 *= srtiny;
    i__1 = igaunt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	argnt[i__] *= srtiny;
/* L110: */
    }
    sumfor *= tiny;
L120:

/*  CHECK WHETHER DABS(C1) IS DECREASING. */
/*  AS LONG AS DABS(C1) IS DECREASING THE RECURSION PROCEEDS TOWARDS */
/*  INCREASING 3JM-VALUES AND IS NUMERICALLY STABLE. ONCE AN INCREASE */
/*  OF DABS(C1) IS DETECTED THE RECURSION DIRECTION IS REVERSED. */

    c1new = abs(c1);
    if (c1new >= c1old) {
	goto L130;
    }
    l12 += 2;
    sumfor += (doublereal) l12 * x1 * x1;
    goto L80;

/*  THIS IS REACHED IF THE STABILITY OF THE RECURSION IS REVERSED. */

L130:
    lmatch = l1 - 1;

/*  end OF FORWARD RECURSION. */

/* --------------------------------------------------------------------- */

/*  START OF BACKWARD RECURSION. */

/*  RECURSIVE CALCULATION OF THE UNNORMALIZED 3JM-SYMBOLS STARTING AT */
/*  L1 = L1MAX IN DIRECTION OF DECREASING L1-VALUES UNTIL FORWARD AND */
/*  BACKWARD RECURSION OVERLAP AT THREE SUCCESSIVE POINTS: */
/*  L1 = LMATCH + 1, L1 = LMATCH AND L1 = LMATCH - 1. */

/*  INITIALIZATION OF THE DUMMY VARIABLE L1. */

    l1 = *l1max;

/*  CALCULATION OF THE FIRST 3JM-SYMBOL. */

    y2 = srtiny;
    l12 = l1 + l1 + 1;
    sumbac = (doublereal) l12 * tiny;

/*  STORAGE OF THE FIRST UNNORMALIZED 3JM-SYMBOL IN THE ARRAY ARGNT. */

    argnt[*ngaunt] = y2;

/*  CALCULATION OF THE SECOND UNNORMALIZED 3JM-SYMBOL WITH THE HELP OF */
/*  THE SIMPLIFIED 2-TERM RECURSION EQUATION WHICH HOLDS AT THE UPPER */
/*  BOUNDARY OF THE DOMAIN OF ALLOWED L1-VALUES. */

    l1p2 = l1 * l1;
    facnew = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));
    y3 = (doublereal) (l12 * (labm + mm * l1 * (l1 + 1))) / ((doublereal) (l1 
	    + 1) * facnew) * y2;
    --l1;

/*  INITIALIZATION OF THE DUMMY VARIABLE KGAUNT. */

    kgaunt = *ngaunt;

/*  CALCULATION OF THE OTHER 3JM-SYMBOLS REQUIRED USING BACKWARD */
/*  RECURSION. */

L140:
    l12 += -2;
    sumbac += (doublereal) l12 * y3 * y3;

/*  CALCULATION OF SOME VARIABLES AND COEFFICIENTS THAT ARE REQUIRED FOR */
/*  THE RECURSIVE CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    l1p2 = l1 * l1;
    y1 = y2;
    y2 = y3;
    facold = facnew;
    facnew = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));

/*  CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    y3 = ((doublereal) (l12 * (labm + mm * l1 * (l1 + 1))) * y2 - facold * (
	    doublereal) l1 * y1) / ((doublereal) (l1 + 1) * facnew);
    --l1;

/*  CHECK WHETHER THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */

    if (abs(y3) < srhuge) {
	goto L160;
    }

/*  THIS IS REACHED IF THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */
/*  IN THAT CASE THE 3JM-SYMBOLS ALREADY CALCULATED AND SUMBAC HAVE TO B */
/*  RESCALED IN ORDER TO PREVENT OVERFLOW. */

    y3 *= srtiny;
    y2 *= srtiny;
    i__1 = *ngaunt;
    for (i__ = kgaunt; i__ <= i__1; ++i__) {
	argnt[i__] *= srtiny;
/* L150: */
    }
    sumbac *= tiny;
L160:

/*  CHECK WHETHER L1 = LMATCH HOLDS. */

    if (l1 == lmatch) {
	goto L190;
    }

/*  STORAGE OF THE UNNORMALIZED 3JM-SYMBOL WITH NON-VANISHING MAGNETIC */
/*  QUANTUM NUMBERS M1, M2 AND M3 REQUIRED FOR THE CALCULATION OF THE */
/*  GAUNT COEFFICIENT IN THE ARRAY ARGNT. */

    --kgaunt;
    argnt[kgaunt] = y3;

/*  CALCULATION OF SUMBAC. */

    l12 += -2;
    sumbac += (doublereal) l12 * y3 * y3;

/*  CALCULATION OF SOME VARIABLES AND COEFFICIENTS THAT ARE REQUIRED FOR */
/*  THE RECURSIVE CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    l1p2 = l1 * l1;
    y1 = y2;
    y2 = y3;
    facold = facnew;
    facnew = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));

/*  CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    y3 = ((doublereal) (l12 * (labm + mm * l1 * (l1 + 1))) * y2 - facold * (
	    doublereal) l1 * y1) / ((doublereal) (l1 + 1) * facnew);
    --l1;

/*  CHECK WHETHER THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */

    if (abs(y3) < srhuge) {
	goto L180;
    }

/*  THIS IS REACHED IF THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */
/*  IN THAT CASE THE 3JM-SYMBOLS ALREADY CALCULATED AND SUMBAC HAVE TO B */
/*  RESCALED IN ORDER TO PREVENT OVERFLOW. */

    y3 *= srtiny;
    y2 *= srtiny;
    i__1 = *ngaunt;
    for (i__ = kgaunt; i__ <= i__1; ++i__) {
	argnt[i__] *= srtiny;
/* L170: */
    }
    sumbac *= tiny;
L180:

/*  CHECK WHETHER L1 = LMATCH HOLDS. */

    if (l1 == lmatch) {
	goto L190;
    }

/*  THIS IS REACHED IF NOT ALL 3JM-SYMBOLS REQUIRED HAVE BEEN CALCULATED */
/*  USING BACKWARD RECURSION. */

    goto L140;

/*  THIS IS REACHED IF L1 = LMATCH HOLDS AND IF ALL BUT ONE OF THE */
/*  UNNORMALIZED 3JM-SYMBOLS REQUIRED ARE CALCULATED WITH THE */
/*  HELP OF BACKWARD RECURSION. */

L190:

/*  CALCULATION OF SOME VARIABLES AND COEFFICIENTS THAT ARE REQUIRED FOR */
/*  THE RECURSIVE CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    l12 += -2;
    l1p2 = l1 * l1;
    y1 = y2;
    y2 = y3;
    facold = facnew;
    facnew = sqrt((doublereal) ((l1p2 - la) * (lb - l1p2) * (l1p2 - m1p2)));

/*  CALCULATION OF THE NEXT UNNORMALIZED 3JM-SYMBOL. */

    y3 = ((doublereal) (l12 * (labm + mm * l1 * (l1 + 1))) * y2 - facold * (
	    doublereal) l1 * y1) / ((doublereal) (l1 + 1) * facnew);

/*  CHECK WHETHER THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */

    if (abs(y3) < srhuge) {
	goto L210;
    }

/*  THIS IS REACHED IF THE LAST UNNORMALIZED 3JM-SYMBOL EXCEEDS SRHUGE. */
/*  IN THAT CASE THE 3JM-SYMBOLS ALREADY CALCULATED AND SUMBAC HAVE TO B */
/*  RESCALED IN ORDER TO PREVENT OVERFLOW. */

    y3 *= srtiny;
    y2 *= srtiny;
    i__1 = *ngaunt;
    for (i__ = kgaunt; i__ <= i__1; ++i__) {
	argnt[i__] *= srtiny;
/* L200: */
    }
    sumbac *= tiny;
L210:

/*  end OF BACKWARD RECURSION. */

/* ---------------------------------------------------------------------- */

/*  DETRMINATION OF THE NORMALIZATION CONSTANT. */

/*  DETERMINE NOW RATIO SUCH THAT YI = XI*RATIO (I = 1,2,3) HOLDS WITH */
/*  MINIMAL ERROR. */

    ratio = (x1 * y1 + x2 * y2 + x3 * y3) / (x1 * x1 + x2 * x2 + x3 * x3);

/*  CHECK WHETHER DABS(RATIO) IS LESS THAN ONE. */

    if (abs(ratio) < 1.) {
	goto L230;
    }

/*  THIS IS REACHED IF DABS(RATIO) IS GREATER THAN ONE. */

    igaunt = kgaunt - 1;
    i__1 = igaunt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	argnt[i__] *= ratio;
/* L220: */
    }
    sum = ratio * ratio * sumfor + sumbac;
    goto L280;

/*  THIS IS REACHED IF DABS(RATIO) IS LESS THAN ONE. */

L230:
    ratio = 1. / ratio;
    i__1 = *ngaunt;
    for (i__ = kgaunt; i__ <= i__1; ++i__) {
	argnt[i__] *= ratio;
/* L240: */
    }
    sum = sumfor + ratio * ratio * sumbac;
    goto L280;

/*  THIS IS REACHED IF ALL UNNORMALIZED 3JM-SYMBOLS WITH NON-VANISHING */
/*  MAGNETIC QUANTUM NUMBERS WERE CALCULATED WITH UPWARD REDURSION AND I */
/*  THE LAST 3JM-SYMBOL HAS TO BE STORED IN THE ARRAY ARGNT. */

L250:
    sum = sumfor + (doublereal) l12 * x1 * x1;
    argnt[*ngaunt] = x1;
    goto L280;

/*  THIS IS REACHED IF ALL UNNORMALIZED 3JM-SYMBOLS WITH NON-VANISHING */
/*  MAGNETIC QUANTUM NUMBERS M1, M2 AND M3 WERE CALCULATED WITH UPWARD */
/*  RECURSION. */
L260:
    sum = sumfor;
    goto L280;

/*  THIS IS REACHED IF ONLY ONE 3JM-SYMBOL WITH NON-VANISHING MAGNETIC */
/*  QUANTUM NUMBERS M1, M2 AND M3 HAS TO BE CALCULATED, BUT MORE THAN */
/*  ONE PARITY COEFFICIENT HAS TO BE CALCULATED. */

L270:

/*  CHECK WHETHER L2 = 0 HOLDS. */

    if (*l2 == 0) {
	goto L340;
    }

/*  CHECK WHETHER L3 = 0 HOLDS. */

    if (*l3 == 0) {
	goto L350;
    }
    argnt[1] = 1. / sqrt((doublereal) (*l1max + *l1max + 1));
    sum = 1.;
/* ---------------------------------------------------------------------- */

/*  START OF THE RECURSIVE CALCULATION OF THE SO-CALLED PARITY */
/*  COEFFICIENTS, I.E. 3JM-SYMBOLS WITH M1 = M2 = M3 = 0. */

L280:

/*  DETERMINE LOWER LIMIT OF THE L1-DOMAIN OF THE PARITY COEFFICIENTS. */

    minpar = (i__1 = *l2 - *l3, abs(i__1));

/*  INITIALIZATION OF THE DUMMY VARIABLES L1 AND KGAUNT. */

    l1 = minpar;
    kgaunt = 1;

/*  CALCULATION OF THE FIRST UNNORMALIZED PARITY COEFFICIENT. */

    x1 = srtiny;
    sump = tiny;
    lnew = l1 + l1 + 1;

/*  CHECK WHETHER MINPAR = L1MIN HOLDS. */

    if (minpar < *l1min) {
	goto L290;
    }

/*  THIS IS REACHED IF MINPAR = L1MIN HOLDS. IN THAT CASE THE FIRST */
/*  UNNORMALIZED PARITY COEFFICIENT IS REQUIRED FOR THE CALCULATION OF */
/*  THE FIRST GAUNT COEFFICIENT AND HAS TO BE STORED IN THE ARRAY ARGNT. */

    argnt[1] *= x1;
    ++kgaunt;
L290:

/*  CHECK WHETHER THERE ARE SOME MORE UNNORMALIZED PARITY COEFFICIENTS T */
/*  BE CALCULATED THAT ARE ONLY REQUIRED FOR THE CALCULATION OF THE */
/*  NORMALIZATION CONSTANT OR WHETHER THE NEXT UNNORMALIZED PARITY */
/*  COEFFICIENT WILL ALSO BE REQUIRED FOR THE CALCULATION OF THE FIRST */
/*  GAUNT COEFFICIENT OR WHETHER ALL UNNORMALIZED PARITY COEFFICIENT WIL */
/*  ALSO BE REQUIRED FOR THE CALCULATION OF THE GAUNT COEFFICIENTS. */

    if (minpar + 2 >= *l1min) {
	goto L310;
    }

/*  THIS IS REACHED IF THERE ARE STILL SOME MORE UNNORMALIZED PARITY */
/*  COEFFICIENTS TO BE CALCULATED THAT WILL ONLY BE REQUIRED FOR THE */
/*  CALCULATION OF THE NORMALIZATION CONSTANT. */

/*  DETERMINE THE NUMBER OF UNNORMALIZED PARITY COEFFICIENTS THAT ARE TO */
/*  BE CALCULATED AND THAT ARE ONLY REQUIRED FOR THE CALCULATION OF THE */
/*  NORMALIZATION CONSTANT OF THE RECURSION STRING. */

    npar = (*l1min - minpar) / 2;

/*  CALCULATION OF THOSE UNNORMALIZED PARITY COEFFICIENTS THAT ARE ONLY */
/*  REQUIRED FOR THE CALCULATION OF THE NORMALIZATION CONSTANT. */

    i__1 = npar;
    for (i__ = 2; i__ <= i__1; ++i__) {
	l1 += 2;
	lp1 = (l1 - 1) * (l1 - 1);
	l1p2 = l1 * l1;
	lold = lnew;
	lnew += 4;
	x1 = -sqrt((doublereal) (lnew * (lp1 - la) * (lb - lp1)) / (
		doublereal) (lold * (l1p2 - la) * (lb - l1p2))) * x1;
	sump += x1 * x1;
/* L300: */
    }

/*  THIS IS REACHED IF THE NEXT UNNORMALIZED PARITY COEFFICIENT WILL */
/*  BE REQUIRED FOR THE CALCULATION OF THE GAUNT COEFFICIENTS. */

L310:
    i__1 = *ngaunt;
    for (igaunt = kgaunt; igaunt <= i__1; ++igaunt) {
	l1 += 2;
	lp1 = (l1 - 1) * (l1 - 1);
	l1p2 = l1 * l1;
	lold = lnew;
	lnew += 4;
	x1 = -sqrt((doublereal) (lnew * (lp1 - la) * (lb - lp1)) / (
		doublereal) (lold * (l1p2 - la) * (lb - l1p2))) * x1;
	sump += x1 * x1;
	argnt[igaunt] *= x1;
/* L320: */
    }

/*  end OF THE RECURSIVE CALCULATION OF THE SO-CALLED PARITY COEFFICIENT */

/* ---------------------------------------------------------------------- */

/*  CALCULATION OF THE NORMALIZATION CONSTANT OF THE STRING OF GAUNT */
/*  COEFFICIENTS. */

    cnorm = sqrt((doublereal) ((*l2 + *l2 + 1) * (*l3 + *l3 + 1)) / pi4) / (
	    sqrt(sum) * sqrt(sump));

/*  CHECK WHETHER THE SIGN OF THE HIGHEST UNNORMALIZED GAUNT COEFFICIENT */
/*  IS CORRECT. */

    if (argnt[*ngaunt] != pow_ii(&c_n1, m2) * (d__1 = argnt[*ngaunt], abs(
	    d__1))) {
	cnorm = -cnorm;
    }

/*  NORMALIZE THE STRING OF GAUNT COEFFICIENTS ( L3 M3 / L2 M2 / L1 M1 ) */

    i__1 = *ngaunt;
    for (igaunt = 1; igaunt <= i__1; ++igaunt) {
	argnt[igaunt] *= cnorm;
/* L330: */
    }
    return 0;

/*  THIS IS REACHED IF EITHER L1 = L2 = L3 = 0 OR L2 = 0 BUT NOT L3 = 0 */
/*  HOLDS. IN THE FIRST CASE WE HAVE */
/*  ( 0 0 / 0 0 / 0 0 ) = 1./DSQRT(4*PI) */
/*  AND IN THE SECOND CASE WE HAVE */
/*  ( L3 M3 / 0 0 / L3 M3 ) = 1./DSQRT(4*PI). */

L340:
    argnt[1] = 1. / sqrt(pi4);
    return 0;

/*  THIS IS REACHED IF L3 = 0 HOLDS. */
/*  IN THAT CASE WE HAVE */
/*  ( L2 -M2 / 0 0 / L2 -M2 ) = (-1)**M2/DSQRT(PI4) */

L350:
    argnt[1] = pow_ii(&c_n1, m2) / sqrt(pi4);
    return 0;

/* ---------------------------------------------------------------------- */

/*  ERROR SECTION */
/*  ------------- */

/*  THIS REACHED IF ILLEGAL ANGULAR MOMENTUM QUANTUM NUMBERS L2, M2, L3 */
/*  AND M3 WERE USED AS INPUT PARAMETERS. */

L360:
    s_wsfe(&io___43);
    e_wsfe();
    s_wsfe(&io___44);
    do_fio(&c__1, (char *)&(*l2), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*m2), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*l3), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*m3), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);

/*  THIS IS REACHED IF DIMENSION OF ARGNT IS NOT LARGE ENOUGH TO STORE */
/*  ALL THE GAUNT COEFFICIENTS REQUIRED, I.E. NGAUNT IS GREATER */
/*  THAN LARRAY. */

/* L370: */
    s_wsfe(&io___45);
    e_wsfe();
    s_wsfe(&io___46);
    do_fio(&c__1, (char *)&(*ngaunt), (ftnlen)sizeof(integer));
    e_wsfe();
/*     WRITE(4,1020) LARRAY,NGAUNT */
    s_stop("", (ftnlen)0);
    return 0;
} /* gaunt_ */

