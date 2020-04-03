/*
 * linearSolver.c
 *
 * Created: Sat 28 Mar 2020 03:41:56 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "main.h"
#include "linearSolver.h"
#include "readInTools.h"
#include "mesh.h"
#include "timeDiscretization.h"
#include "memTools.h"
#include "fluxCalculation.h"

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
			Dinv = dyn3DdblArray(nElems, NVAR, NVAR);
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

/*
 * compute dot product for vectors a and b: result = a * b
 */
double vectorDotProduct(double A[NVAR][nElems], double B[NVAR][nElems])
{
	double res = 0.0;

	#pragma omp parallel for reduction(+:res)
	for (long iElem = 0; iElem < nElems; ++iElem) {
		res += A[RHO][iElem] * B[RHO][iElem];
		res += A[MX][iElem]  * B[MX][iElem];
		res += A[MY][iElem]  * B[MY][iElem];
		res += A[E][iElem]   * B[E][iElem];
	}

	return res;
}

/*
 * compute inverse of a 4x4 matrix
 */
bool calcDinv(double **A, double **Ainv)
{
	double det = A[0][0]*(A[1][1]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])
		+A[1][2]*(A[2][3]*A[3][1]-A[2][1]*A[3][3])+A[1][3]*(A[2][1]
		*A[3][2]-A[2][2]*A[3][1]))-A[0][1]*(A[1][0]*(A[2][2]*A[3][3]
		-A[2][3]*A[3][2])+A[1][2]*(A[2][3]*A[3][0]-A[2][0]*A[3][3])
		+A[1][3]*(A[2][0]*A[3][2]-A[2][2]*A[3][0]))+A[0][2]*(A[1][0]
		*(A[2][1]*A[3][3]-A[2][3]*A[3][1])+A[1][1]*(A[2][3]*A[3][0]
		-A[2][0]*A[3][3])+A[1][3]*(A[2][0]*A[3][1]-A[2][1]*A[3][0]))
		-A[0][3]*(A[1][0]*(A[2][1]*A[3][2]-A[2][2]*A[3][1])+A[1][1]
		*(A[2][2]*A[3][0]-A[2][0]*A[3][2])+A[1][2]*(A[2][0]*A[3][1]
		-A[2][1]*A[3][0]));

	if (fabs(det) < DBL_EPSILON) {
		return false;
	}

	double coFac[4][4];

	coFac[0][0] = A[1][1]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])+A[1][2]*(A[2][3]
		*A[3][1]-A[2][1]*A[3][3])+A[1][3]*(A[2][1]*A[3][2]-A[2][2]*A[3][1]);
	coFac[0][1] = A[1][0]*(A[2][3]*A[3][2]-A[2][2]*A[3][3])+A[1][2]*(A[2][0]
		*A[3][3]-A[2][3]*A[3][0])+A[1][3]*(A[2][2]*A[3][0]-A[2][0]*A[3][2]);
	coFac[0][2] = A[1][0]*(A[2][1]*A[3][3]-A[2][3]*A[3][1])+A[1][1]*(A[2][3]
		*A[3][0]-A[2][0]*A[3][3])+A[1][3]*(A[2][0]*A[3][1]-A[2][1]*A[3][0]);
	coFac[0][3] = A[1][0]*(A[2][2]*A[3][1]-A[2][1]*A[3][2])+A[1][1]*(A[2][0]
		*A[3][2]-A[2][2]*A[3][0])+A[1][2]*(A[2][1]*A[3][0]-A[2][0]*A[3][1]);
	coFac[1][0] = A[0][1]*(A[2][3]*A[3][2]-A[2][2]*A[3][3])+A[0][2]*(A[2][1]
		*A[3][3]-A[2][3]*A[3][1])+A[0][3]*(A[2][2]*A[3][1]-A[2][1]*A[3][2]);
	coFac[1][1] = A[0][0]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])+A[0][2]*(A[2][3]
		*A[3][0]-A[2][0]*A[3][3])+A[0][3]*(A[2][0]*A[3][2]-A[2][2]*A[3][0]);
	coFac[1][2] = A[0][0]*(A[2][3]*A[3][1]-A[2][1]*A[3][3])+A[0][1]*(A[2][0]
		*A[3][3]-A[2][3]*A[3][0])+A[0][3]*(A[2][1]*A[3][0]-A[2][0]*A[3][1]);
	coFac[1][3] = A[0][0]*(A[2][1]*A[3][2]-A[2][2]*A[3][1])+A[0][1]*(A[2][2]
		*A[3][0]-A[2][0]*A[3][2])+A[0][2]*(A[2][0]*A[3][1]-A[2][1]*A[3][0]);
	coFac[2][0] = A[0][1]*(A[1][2]*A[3][3]-A[1][3]*A[3][2])+A[0][2]*(A[1][3]
		*A[3][1]-A[1][1]*A[3][3])+A[0][3]*(A[1][1]*A[3][2]-A[1][2]*A[3][1]);
	coFac[2][1] = A[0][0]*(A[1][3]*A[3][2]-A[1][2]*A[3][3])+A[0][2]*(A[1][0]
		*A[3][3]-A[1][3]*A[3][0])+A[0][3]*(A[1][2]*A[3][0]-A[1][0]*A[3][2]);
	coFac[2][2] = A[0][0]*(A[1][1]*A[3][3]-A[1][3]*A[3][1])+A[0][1]*(A[1][3]
		*A[3][0]-A[1][0]*A[3][3])+A[0][3]*(A[1][0]*A[3][1]-A[1][1]*A[3][0]);
	coFac[2][3] = A[0][0]*(A[1][2]*A[3][1]-A[1][1]*A[3][2])+A[0][1]*(A[1][0]
		*A[3][2]-A[1][2]*A[3][0])+A[0][2]*(A[1][1]*A[3][0]-A[1][0]*A[3][1]);
	coFac[3][0] = A[0][1]*(A[1][3]*A[2][2]-A[1][2]*A[2][3])+A[0][2]*(A[1][1]
		*A[2][3]-A[1][3]*A[2][1])+A[0][3]*(A[1][2]*A[2][1]-A[1][1]*A[2][2]);
	coFac[3][1] = A[0][0]*(A[1][2]*A[2][3]-A[1][3]*A[2][2])+A[0][2]*(A[1][3]
		*A[2][0]-A[1][0]*A[2][3])+A[0][3]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]);
	coFac[3][2] = A[0][0]*(A[1][3]*A[2][1]-A[1][1]*A[2][3])+A[0][1]*(A[1][0]
		*A[2][3]-A[1][3]*A[2][0])+A[0][3]*(A[1][1]*A[2][0]-A[1][0]*A[2][1]);
	coFac[3][3] = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])+A[0][1]*(A[1][2]
		*A[2][0]-A[1][0]*A[2][2])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			Ainv[i][j] = coFac[j][i] / det;
		}
	}

	return true;
}

/*
 * compute the global jacobian matrix by use of finite differences
 */
void buildMatrix(double time, double dt)
{
	double ***D = dyn3DdblArray(nElems, NVAR, NVAR);

	/* fd for diagonal element */
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		for (int iVar = 0; iVar < NVAR; ++iVar) {
			fluxJacobianFD(time, aElem, iVar);

			D[iElem][RHO][iVar] += aElem->u_t[RHO];
			D[iElem][MX][iVar]  += aElem->u_t[MX];
			D[iElem][MY][iVar]  += aElem->u_t[MY];
			D[iElem][E][iVar]   += aElem->u_t[E];

			side_t *aSide = aElem->firstSide;
			long iSide = 0;
			while (aSide) {
				long NBelemId = aSide->connection->elem->id;

				if ((NBelemId >= 0) && (NBelemId < nElems)) {
					long NBsideId = elemToElem[1][iSide][iElem];

					lowerUpper[RHO][iVar][NBsideId][NBelemId] += aSide->connection->elem->u_t[RHO];
					lowerUpper[MX][iVar][NBsideId][NBelemId]  += aSide->connection->elem->u_t[MX];
					lowerUpper[MY][iVar][NBsideId][NBelemId]  += aSide->connection->elem->u_t[MY];
					lowerUpper[E][iVar][NBsideId][NBelemId]   += aSide->connection->elem->u_t[E];
				}

				aSide = aSide->nextElemSide;
				iSide++;
			}
		}

		for (int iVar1 = 0; iVar1 < NVAR; ++iVar1) {
			for (int iVar2 = 0; iVar2 < NVAR; ++iVar2) {
				D[iElem][iVar1][iVar2] *= - dt;
			}

			D[iElem][iVar1][iVar1] += 1.0;
		}

		/* backup */
		bool isOK = calcDinv(D[iElem], Dinv[iElem]);
		if (!isOK) {
			printf("| LUSGS D-Matrix is singular at Element %ld\n", iElem);
			exit(1);
		}
	}

	// TODO: maybe dt can be multiplied above
	for (int i = 0; i < NVAR; ++i) {
		for (int j = 0; j < NVAR; ++j) {
			for (int k = 0; k < 4; ++k) {
				for (int l = 0; l < nElems; ++l) {
					lowerUpper[i][j][k][l] *= - dt;
				}
			}
		}
	}

	free(D);
}

/*
 * uses matrix free to solve the linear system, deltaX=0 is the initial guess
 * X0 is already stored in U
 */
void GMRES_M(double time, double dt, double alpha, double beta, double B[NVAR][nElems],
		double normB, double *abortCrit, double deltaX[NVAR][nElems])
{
	*abortCrit = epsGMRES * normB;
}
