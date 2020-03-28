/*
 * linearSolver.c
 *
 * Created: Sat 28 Mar 2020 03:41:56 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "main.h"
#include "linearSolver.h"
#include "readInTools.h"
#include "mesh.h"
#include "timeDiscretization.h"
#include "memTools.h"

/* extern variables */
int nKdim;
int nNewtonIter;
int nNewtonIterGlobal;

int nGMRESiterGlobal;

int nInnerNewton;

int nInnerGMRES;

int iterGlobal;
bool usePrecond;

double rEps0, srEps0;

double eps2newton;

double epsGMRES;
double gammaEW;

double **XK, **R_XK;
double ***Dinv;
double ****lowerUpper;
long ***elemToElem;

//bool implicitInitIsDone = false;

/*
 * initialize linear solver
 */
void initLinearSolver(void)
{
	if (isImplicit) {
		nKdim = getInt("nKdim", "5");

		eps2newton = getDbl("epsNewton", "0.001");
		eps2newton *= eps2newton;

		epsGMRES = getDbl("epsGMRES", "0.001");

		rEps0 = sqrt(DBL_EPSILON);
		srEps0 = 1.0 / rEps0;

		nNewtonIter = getInt("nNewtonIter", "20");

		nNewtonIterGlobal = 0;
		nGMRESiterGlobal = 0;
		nInnerNewton = 0;
		nInnerGMRES = 0;

		gammaEW = getDbl("gammaEW", "0.9");

		XK = dyn2DdblArray(NVAR, nElems);
		R_XK = dyn2DdblArray(NVAR, nElems);

		usePrecond = getBool("precond", "F");
		if (usePrecond) {
			Dinv = dyn3DdblArray(NVAR, NVAR, nElems);
			lowerUpper = dyn4DdblArray(NVAR, NVAR, 4, nElems);
			elemToElem = dyn3DintArray(2, 4, nElems);

			for (long iElem = 0; iElem < nElems; ++iElem) {
				elem_t *aElem = elem[iElem];
				side_t *aSide = aElem->firstSide;
				long iSide1 = 0;
				while (aSide) {
					long NBelemId = aSide->connection->elem->id;
					if ((NBelemId > 0) && (NBelemId < nElems)) {
						side_t *bSide = aSide->connection->elem->firstSide;
						long iSide2 = 0;
						while (bSide) {
							long elemId = aSide->connection->elem->id;
							if ((elemId > 0) && (elemId < nElems) && (bSide->connection->elem->id == iElem)) {
								elemToElem[0][iSide2][NBelemId] = iElem;
								elemToElem[1][iSide2][NBelemId] = iSide1;
							}
							iSide2++;
							bSide = bSide->nextElemSide;
						}
					}
					iSide1++;
					aSide = aSide->nextElemSide;
				}
			}
		}
	}
}
