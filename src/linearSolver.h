/*
 * linearSolver.h
 *
 * Created: Sat 28 Mar 2020 03:09:46 PM CET
 * Author : hhh
 */

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <stdbool.h>

extern int nKdim;
extern int nNewtonIter;
extern int nNewtonIterGlobal;

extern int nGMRESiterGlobal;

extern int nInnerNewton;

extern int nInnerGMRES;

extern int iterGlobal;
extern bool usePrecond;

extern double rEps0, srEps0;

extern double eps2newton;

extern double epsGMRES;
extern double gammaEW;

extern double **XK, **R_XK;
extern double ***Dinv;
extern double ****lowerUpper;
extern long ***elemToElem;

//extern bool implicitInitIsDone;

void initLinearSolver(void);

#endif
