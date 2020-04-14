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
#include "equationOfState.h"
#include "finiteVolume.h"

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
double eps2newton_sq;

double epsGMRES;
double gammaEW;

double **XK, **R_XK;
double ***Dinv;
double ****lowerUpper;
long ***elemToElem;

double ***V;
double ***Z;

/*
 * initialize linear solver
 */
void initLinearSolver(void)
{
	if (isImplicit) {
		nKdim = getInt("nKdim", "5");

		eps2newton = getDbl("epsNewton", "0.001");
		eps2newton *= eps2newton;
		eps2newton_sq = sqrt(eps2newton);

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
					if ((NBelemId >= 0) && (NBelemId < nElems)) {
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

		V = dyn3DdblArray(nKdim, NVAR, nElems);
		Z = dyn3DdblArray(nKdim, NVAR, nElems);

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

	coFac[0][0] = A[1][1]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])+A[1][2]*(A[2][3]*A[3][1]-A[2][1]*A[3][3])+A[1][3]*(A[2][1]*A[3][2]-A[2][2]*A[3][1]);
	coFac[0][1] = A[1][0]*(A[2][3]*A[3][2]-A[2][2]*A[3][3])+A[1][2]*(A[2][0]*A[3][3]-A[2][3]*A[3][0])+A[1][3]*(A[2][2]*A[3][0]-A[2][0]*A[3][2]);
	coFac[0][2] = A[1][0]*(A[2][1]*A[3][3]-A[2][3]*A[3][1])+A[1][1]*(A[2][3]*A[3][0]-A[2][0]*A[3][3])+A[1][3]*(A[2][0]*A[3][1]-A[2][1]*A[3][0]);
	coFac[0][3] = A[1][0]*(A[2][2]*A[3][1]-A[2][1]*A[3][2])+A[1][1]*(A[2][0]*A[3][2]-A[2][2]*A[3][0])+A[1][2]*(A[2][1]*A[3][0]-A[2][0]*A[3][1]);
	coFac[1][0] = A[0][1]*(A[2][3]*A[3][2]-A[2][2]*A[3][3])+A[0][2]*(A[2][1]*A[3][3]-A[2][3]*A[3][1])+A[0][3]*(A[2][2]*A[3][1]-A[2][1]*A[3][2]);
	coFac[1][1] = A[0][0]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])+A[0][2]*(A[2][3]*A[3][0]-A[2][0]*A[3][3])+A[0][3]*(A[2][0]*A[3][2]-A[2][2]*A[3][0]);
	coFac[1][2] = A[0][0]*(A[2][3]*A[3][1]-A[2][1]*A[3][3])+A[0][1]*(A[2][0]*A[3][3]-A[2][3]*A[3][0])+A[0][3]*(A[2][1]*A[3][0]-A[2][0]*A[3][1]);
	coFac[1][3] = A[0][0]*(A[2][1]*A[3][2]-A[2][2]*A[3][1])+A[0][1]*(A[2][2]*A[3][0]-A[2][0]*A[3][2])+A[0][2]*(A[2][0]*A[3][1]-A[2][1]*A[3][0]);
	coFac[2][0] = A[0][1]*(A[1][2]*A[3][3]-A[1][3]*A[3][2])+A[0][2]*(A[1][3]*A[3][1]-A[1][1]*A[3][3])+A[0][3]*(A[1][1]*A[3][2]-A[1][2]*A[3][1]);
	coFac[2][1] = A[0][0]*(A[1][3]*A[3][2]-A[1][2]*A[3][3])+A[0][2]*(A[1][0]*A[3][3]-A[1][3]*A[3][0])+A[0][3]*(A[1][2]*A[3][0]-A[1][0]*A[3][2]);
	coFac[2][2] = A[0][0]*(A[1][1]*A[3][3]-A[1][3]*A[3][1])+A[0][1]*(A[1][3]*A[3][0]-A[1][0]*A[3][3])+A[0][3]*(A[1][0]*A[3][1]-A[1][1]*A[3][0]);
	coFac[2][3] = A[0][0]*(A[1][2]*A[3][1]-A[1][1]*A[3][2])+A[0][1]*(A[1][0]*A[3][2]-A[1][2]*A[3][0])+A[0][2]*(A[1][1]*A[3][0]-A[1][0]*A[3][1]);
	coFac[3][0] = A[0][1]*(A[1][3]*A[2][2]-A[1][2]*A[2][3])+A[0][2]*(A[1][1]*A[2][3]-A[1][3]*A[2][1])+A[0][3]*(A[1][2]*A[2][1]-A[1][1]*A[2][2]);
	coFac[3][1] = A[0][0]*(A[1][2]*A[2][3]-A[1][3]*A[2][2])+A[0][2]*(A[1][3]*A[2][0]-A[1][0]*A[2][3])+A[0][3]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]);
	coFac[3][2] = A[0][0]*(A[1][3]*A[2][1]-A[1][1]*A[2][3])+A[0][1]*(A[1][0]*A[2][3]-A[1][3]*A[2][0])+A[0][3]*(A[1][1]*A[2][0]-A[1][0]*A[2][1]);
	coFac[3][3] = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])+A[0][1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

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

					lowerUpper[RHO][iVar][NBsideId][NBelemId]
						-= aSide->connection->elem->u_t[RHO] * dt;

					lowerUpper[MX][iVar][NBsideId][NBelemId]
						-= aSide->connection->elem->u_t[MX] * dt;

					lowerUpper[MY][iVar][NBsideId][NBelemId]
						-= aSide->connection->elem->u_t[MY] * dt;

					lowerUpper[E][iVar][NBsideId][NBelemId]
						-= aSide->connection->elem->u_t[E] * dt;
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

	free(D);

	/* TODO: get rid of */
	//for (int i = 0; i < NVAR; ++i) {
	//	for (int j = 0; j < NVAR; ++j) {
	//		for (int k = 0; k < 4; ++k) {
	//			for (int l = 0; l < nElems; ++l) {
	//				lowerUpper[i][j][k][l] *= - dt;
	//			}
	//		}
	//	}
	//}

}

/*
 * LUSGS preconditioner, uses precomputed Block-LUSGS, B is the old vector and
 * deltaX is preconditioned vector
 */
void LUSGS_FD(double time, double dt, double **B, double **deltaX)
{
	/* compute LU-SGS with FDs */
	double deltaXstar[NVAR][nElems];
	memset(deltaXstar, 0, NVAR * nElems * sizeof(double));

	/* forward sweep */
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		side_t *aSide = aElem->firstSide;
		long iSide = 0;

		while (aSide) {
			long NBelemId = aSide->connection->elem->id;
			if ((NBelemId >= 0) && (NBelemId < iElem)) {
				/* rotate state into normal direction */
				double tmp[4] = {0.0};
				for (int i = 0; i < NVAR; ++i) {
					for (int j = 0; j < NVAR; ++j) {
						tmp[i] += lowerUpper[i][j][iSide][iElem] * deltaXstar[j][NBelemId];
					}
				}

				deltaXstar[0][iElem] += tmp[0];
				deltaXstar[1][iElem] += tmp[1];
				deltaXstar[2][iElem] += tmp[2];
				deltaXstar[3][iElem] += tmp[3];
			}

			aSide = aSide->nextElemSide;
			iSide++;
		}

		/* calculate deltaXstar */
		double tmp[4] = {0.0};
		for (int i = 0; i < NVAR; ++i) {
			for (int j = 0; j < NVAR; ++j) {
				tmp[i] += Dinv[iElem][i][j] * (B[j][iElem] - deltaXstar[j][iElem]);
			}
		}
		deltaXstar[0][iElem] = tmp[0];
		deltaXstar[1][iElem] = tmp[1];
		deltaXstar[2][iElem] = tmp[2];
		deltaXstar[3][iElem] = tmp[3];
	}

	/* backward sweep */
	#pragma omp parallel for
	for (long iElem = nElems - 1; iElem >= 0; --iElem) {
		elem_t *aElem = elem[iElem];
		side_t *aSide = aElem->firstSide;
		long iSide = 0;
		while (aSide) {
			long NBelemId = aSide->connection->elem->id;
			if ((NBelemId > iElem) && (NBelemId < nElems)) {
				/* rotate state into normal direction */
				double tmp[4] = {0.0};
				for (int i = 0; i < NVAR; ++i) {
					for (int j = 0; j < NVAR; ++j) {
						tmp[i] += lowerUpper[i][j][iSide][iElem] * deltaX[j][NBelemId];
					}
				}

				deltaX[0][iElem] += tmp[0];
				deltaX[1][iElem] += tmp[1];
				deltaX[2][iElem] += tmp[2];
				deltaX[3][iElem] += tmp[3];
			}

			aSide = aSide->nextElemSide;
			iSide++;
		}

		/* calculate deltaXstar */
		double tmp[4] = {0.0};
		for (int i = 0; i < NVAR; ++i) {
			for (int j = 0; j < NVAR; ++j) {
				tmp[i] += Dinv[iElem][i][j] * deltaX[j][iElem];
			}
		}

		deltaX[0][iElem] = deltaXstar[0][iElem] - tmp[0];
		deltaX[1][iElem] = deltaXstar[1][iElem] - tmp[1];
		deltaX[2][iElem] = deltaXstar[2][iElem] - tmp[2];
		deltaX[3][iElem] = deltaXstar[3][iElem] - tmp[3];
	}
}

/*
 * computes matrix vector product using spatial operator and finite difference
 * approach, A is operator at linearization state xk (Newton iteration)
 */
void matrixVector(double time, double dt, double alpha, double beta, double **V,
		double res[NVAR][nElems])
{
	/* prerequisites for FD matrix vector approximation */
	double epsFD = 0.0;
	#pragma omp parallel for reduction(+:epsFD)
	for (long iElem = 0; iElem < nElems; ++iElem) {
		epsFD += V[RHO][iElem] * V[RHO][iElem];
		epsFD += V[MX][iElem]  * V[MX][iElem];
		epsFD += V[MY][iElem]  * V[MY][iElem];
		epsFD += V[E][iElem]   * V[E][iElem];
	}
	epsFD = rEps0 / sqrt(epsFD);

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		aElem->cVar[RHO] = XK[RHO][iElem] + epsFD * V[RHO][iElem];
		aElem->cVar[MX]  = XK[MX][iElem]  + epsFD * V[MX][iElem];
		aElem->cVar[MY]  = XK[MY][iElem]  + epsFD * V[MY][iElem];
		aElem->cVar[E]   = XK[E][iElem]   + epsFD * V[E][iElem];

		consPrim(aElem->cVar, aElem->pVar);
	}

	fvTimeDerivative(time, iterGlobal);

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		res[RHO][iElem] = V[RHO][iElem] - alpha * dt * (aElem->u_t[RHO] - R_XK[RHO][iElem]) / epsFD;
		res[MX][iElem]  = V[MX][iElem]  - alpha * dt * (aElem->u_t[MX]  - R_XK[MX][iElem])  / epsFD;
		res[MY][iElem]  = V[MY][iElem]  - alpha * dt * (aElem->u_t[MY]  - R_XK[MY][iElem])  / epsFD;
		res[E][iElem]   = V[E][iElem]   - alpha * dt * (aElem->u_t[E]   - R_XK[E][iElem])   / epsFD;
	}
}

/*
 * uses matrix free to solve the linear system, deltaX=0 is the initial guess
 * X0 is already stored in U
 */
void GMRES_M(double time, double dt, double alpha, double beta, double B[NVAR][nElems],
		double normB, double *abortCrit, double deltaX[NVAR][nElems])
{
	*abortCrit = epsGMRES * normB;

	double R0[NVAR][nElems];
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		R0[RHO][iElem] = - B[RHO][iElem];
		R0[MX][iElem]  = - B[MX][iElem];
		R0[MY][iElem]  = - B[MY][iElem];
		R0[E][iElem]   = - B[E][iElem];
	}

	memset(deltaX, 0, NVAR * nElems * sizeof(double));

	double normR0 = normB;

	nInnerGMRES = 0;

	/* GMRES(m) */
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		V[0][RHO][iElem] = R0[RHO][iElem] / normR0;
		V[0][MX][iElem]  = R0[MX][iElem]  / normR0;
		V[0][MY][iElem]  = R0[MY][iElem]  / normR0;
		V[0][E][iElem]   = R0[E][iElem]   / normR0;
	}


	double gam[nKdim + 1];
	gam[0] = normR0;

	int m;
	double H[nKdim + 1][nKdim + 1], C[nKdim], S[nKdim];
	memset(H, 0, (nKdim + 1) * (nKdim + 1) * sizeof(double));
	for (m = 0; m < nKdim; ++m) {
		nInnerGMRES++;

		if (usePrecond) {
			LUSGS_FD(time, dt, V[m], Z[m]);
		} else {
			memcpy(Z[m], V[m], NVAR * nElems * sizeof(double));
		}

		double W[NVAR][nElems];
		matrixVector(time, dt, alpha, beta, Z[m], W);

		/* Gram-Schmidt */
		for (int nn = 0; nn <= m; ++nn) {
			double res = 0.0;
			#pragma omp parallel for reduction(+:res)
			for (long iElem = 0; iElem < nElems; ++iElem) {
				res += V[nn][RHO][iElem] * W[RHO][iElem];
				res += V[nn][MX][iElem]  * W[MX][iElem];
				res += V[nn][MY][iElem]  * W[MY][iElem];
				res += V[nn][E][iElem]   * W[E][iElem];
			}
			H[nn][m] = res;

			#pragma omp parallel for
			for (int iElem = 0; iElem < nElems; ++iElem) {
				W[RHO][iElem] -= H[nn][m] * V[nn][RHO][iElem];
				W[MX][iElem]  -= H[nn][m] * V[nn][MX][iElem];
				W[MY][iElem]  -= H[nn][m] * V[nn][MY][iElem];
				W[E][iElem]   -= H[nn][m] * V[nn][E][iElem];
			}
		}

		double res = vectorDotProduct(W, W);

		H[m + 1][m] = sqrt(res);

		/* Givens rotation */
		for (int nn = 0; nn <= m - 1; ++nn) {
			double tmp   = C[nn] * H[nn][m] + S[nn] * H[nn + 1][m];
			H[nn + 1][m] = - S[nn] * H[nn][m] + C[nn] * H[nn + 1][m];
			H[nn][m]     = tmp;
		}

		double bet = sqrt(H[m][m] * H[m][m] + H[m + 1][m] * H[m + 1][m]);
		S[m] = H[m + 1][m] / bet;
		C[m] = H[m][m] / bet;
		H[m][m] = bet;
		gam[m + 1] = - S[m] * gam[m];
		gam[m] = C[m] * gam[m];

		if ((fabs(gam[m + 1]) <= *abortCrit) || (m == nKdim - 1)) {
			double alp[nKdim];
			for (int nn = m; nn >= 0; --nn) {
				alp[nn] = gam[nn];

				for (int o = nn + 1; o <= m; ++o) {
					alp[nn] -= H[nn][o] * alp[o];
				}

				alp[nn] /= H[nn][nn];
			}

			for (int nn = 0; nn <= m; ++nn) {
				#pragma omp parallel for
				for (long iElem = 0; iElem < nElems; ++iElem) {
					deltaX[RHO][iElem] += alp[nn] * Z[nn][RHO][iElem];
					deltaX[MX][iElem]  += alp[nn] * Z[nn][MX][iElem];
					deltaX[MY][iElem]  += alp[nn] * Z[nn][MY][iElem];
					deltaX[E][iElem]   += alp[nn] * Z[nn][E][iElem];
				}
			}
			nGMRESiterGlobal += nInnerGMRES;

			return;
		} else {
			/* no convergence, next iteration */
			#pragma omp parallel for
			for (long iElem = 0; iElem < nElems; ++iElem) {
				V[m + 1][RHO][iElem] = W[RHO][iElem] / H[m + 1][m];
				V[m + 1][MX][iElem]  = W[MX][iElem]  / H[m + 1][m];
				V[m + 1][MY][iElem]  = W[MY][iElem]  / H[m + 1][m];
				V[m + 1][E][iElem]   = W[E][iElem]   / H[m + 1][m];
			}
		}
	}

	free(V);
	free(Z);

	printf("| GMRES not converged with %d iterations:\n", nInnerGMRES);
	printf("| Norm_R0: %g\n", fabs(gam[0]));
	printf("| Norm_R : %g\n", fabs(gam[m]));
	exit(1);
}
