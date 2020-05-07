/** \file
 *
 * \author hhh
 * \date Sat 28 Mar 2020 03:09:46 PM CET
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

extern bool usePrecond;

extern double rEps0;
extern double srEps0;

extern double eps2newton;
extern double eps2newton_sq;

extern double epsGMRES;
extern double gamEW;

extern double **XK;
extern double **R_XK;

void initLinearSolver(void);
double vectorDotProduct(double **A, double **B);
void GMRES_M(double time, double dt, double alpha, double **B,
		double normB, double *abortCrit, double **deltaX);
void freeLinearSolver(void);

#endif
