/*
 * linearSolver.h
 *
 * Created: Sat 28 Mar 2020 03:09:46 PM CET
 * Author : hhh
 */

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <stdbool.h>

#include "main.h"
#include "mesh.h"

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

void initLinearSolver(void);
double vectorDotProduct(double A[NVAR][nElems], double B[NVAR][nElems]);
void buildMatrix(double time, double dt);
void GMRES_M(double time, double dt, double alpha, double beta, double B[NVAR][nElems],
		double normB, double *abortCrit, double deltaX[NVAR][nElems]);

#endif
