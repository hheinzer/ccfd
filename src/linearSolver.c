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
double ***D;
double ***Dinv;
double **dRdU;

double ***V;
double ***Z;
double **R0;
double **W;
double **deltaXstar;

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
			D = dyn3DdblArray(nElems, NVAR, NVAR);
			Dinv = dyn3DdblArray(nElems, NVAR, NVAR);
			deltaXstar = dyn2DdblArray(NVAR, nElems);
			dRdU = dyn2DdblArray(NVAR * nElems, NVAR * nElems);
		}

		V = dyn3DdblArray(nKdim, NVAR, nElems);
		Z = dyn3DdblArray(nKdim, NVAR, nElems);
		R0 = dyn2DdblArray(NVAR, nElems);
		W = dyn2DdblArray(NVAR, nElems);
	}
}

/*
 * compute dot product for vectors a and b: result = a * b
 */
double vectorDotProduct(double **A, double **B)
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

	if (fabs(det) <= 1e-10) {
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
	#pragma omp parallel for
	for (long iElem = 0; iElem < NVAR * nElems; ++iElem) {
		for (long jElem = 0; jElem < NVAR * nElems; ++jElem) {
			dRdU[iElem][jElem] = 0.0;
		}
	}

	long s = 0;
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		for (int iVar = 0; iVar < NVAR; ++iVar) {
			aElem->pVar[iVar] += rEps0;
			fvTimeDerivative(t, 0);
			aElem->pVar[iVar] -= rEps0;

			long r = iElem * NVAR;
			for (int jVar = 0; jVar < NVAR; ++jVar) {
				dRdU[r++][s] += (aElem->u_t[jVar]
						- R_XK[jVar][iElem]) * srEps0;

			}

			side_t *aSide = aElem->firstSide;
			while (aSide) {
				if ((aSide->connection->elem->id >=0)
						&& (aSide->connection->elem->id) < nElems) {
					long jElem = aSide->connection->elem->id;
					elem_t *bElem = elem[jElem];

					r = jElem * NVAR;
					for (int jVar = 0; jVar < NVAR; ++jVar) {
						dRdU[r++][s] += (bElem->u_t[jVar]
								- R_XK[jVar][jElem]) * srEps0;

					}
				}

				aSide = aSide->nextElemSide;
			}

			++s;
		}
	}

	#pragma omp parallel for
	for (long iElem = 0; iElem < NVAR * nElems; ++iElem) {
		for (long jElem = 0; jElem < NVAR * nElems; ++jElem) {
			dRdU[iElem][jElem] *= - dt;
		}

		dRdU[iElem][iElem] += 1.0;
	}

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		long r = iElem * NVAR;

		for (int iVar = 0; iVar < NVAR; ++iVar) {
			for (int jVar = 0; jVar < NVAR; ++jVar) {
				D[iElem][iVar][jVar] = dRdU[r + iVar][r + jVar];
			}
		}

		bool isOK = calcDinv(D[iElem], Dinv[iElem]);
		if (!isOK) {
			printf("| LUSGS D-Matrix is singular at Element %ld\n", iElem);
			exit(1);
		}

	}
}

/*
 * LUSGS preconditioner
 */
void LUSGS(double time, double dt, double **B, double **delX)
{
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		for (int iVar = 0; iVar < NVAR; ++iVar) {
			delX[iVar][iElem] = 0.0;
			deltaXstar[iVar][iElem] = 0.0;
		}
	}

	/* forward sweep */
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		side_t *aSide = aElem->firstSide;

		double tmp1[NVAR] = {0.0};
		while (aSide) {
			long NBelemID = aSide->connection->elem->id;
			if ((NBelemID >= 0) && (NBelemID < iElem)) {
				long r = iElem * NVAR;
				long s = NBelemID * NVAR;
				for (int iVar = 0; iVar < NVAR; ++iVar) {
					for (int jVar = 0; jVar < NVAR; ++jVar) {
						tmp1[iVar] += dRdU[r + iVar][s + jVar]
							* deltaXstar[jVar][NBelemID];
					}
				}
			}

			aSide = aSide->nextElemSide;
		}

		double tmp2[NVAR] = {0.0};
		for (int iVar = 0; iVar < NVAR; ++iVar) {
			for (int jVar = 0; jVar < NVAR; ++jVar) {
				tmp2[iVar] += Dinv[iElem][iVar][jVar]
					* (B[jVar][iElem] - tmp1[jVar]);
			}
			deltaXstar[iVar][iElem] = tmp2[iVar];
		}
	}

	/* backwards sweep */
	for (long iElem = nElems - 1; iElem >= 0; --iElem) {
		elem_t *aElem = elem[iElem];
		side_t *aSide = aElem->firstSide;

		double tmp1[NVAR] = {0.0};
		while (aSide) {
			long NBelemID = aSide->connection->elem->id;
			if ((NBelemID > iElem) && (NBelemID < nElems)) {
				long r = iElem * NVAR;
				long s = NBelemID * NVAR;
				for (int iVar = 0; iVar < NVAR; ++iVar) {
					for (int jVar = 0; jVar < NVAR; ++jVar) {
						tmp1[iVar] += dRdU[r + iVar][s + jVar]
							* delX[jVar][NBelemID];
					}
				}
			}

			aSide = aSide->nextElemSide;
		}

		double tmp2[NVAR] = {0.0};
		for (int iVar = 0; iVar < NVAR; ++iVar) {
			for (int jVar = 0; jVar < NVAR; ++jVar) {
				tmp2[iVar] += Dinv[iElem][iVar][jVar]* tmp1[jVar];
			}
			delX[iVar][iElem] = deltaXstar[iVar][iElem] - tmp2[iVar];
		}
	}
}

/*
 * computes matrix vector product using spatial operator and finite difference
 * approach, A is operator at linearization state xk (Newton iteration)
 */
void matrixVector(double time, double dt, double alpha, double beta, double **v,
		double **res)
{
	/* prerequisites for FD matrix vector approximation */
	double epsFD = vectorDotProduct(v, v);
	epsFD = rEps0 / sqrt(epsFD);

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		aElem->cVar[RHO] = XK[RHO][iElem] + epsFD * v[RHO][iElem];
		aElem->cVar[MX]  = XK[MX][iElem]  + epsFD * v[MX][iElem];
		aElem->cVar[MY]  = XK[MY][iElem]  + epsFD * v[MY][iElem];
		aElem->cVar[E]   = XK[E][iElem]   + epsFD * v[E][iElem];

		consPrim(aElem->cVar, aElem->pVar);
	}

	fvTimeDerivative(time, iterGlobal);

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		res[RHO][iElem] = v[RHO][iElem] - alpha * dt * (aElem->u_t[RHO] - R_XK[RHO][iElem]) / epsFD;
		res[MX][iElem]  = v[MX][iElem]  - alpha * dt * (aElem->u_t[MX]  - R_XK[MX][iElem])  / epsFD;
		res[MY][iElem]  = v[MY][iElem]  - alpha * dt * (aElem->u_t[MY]  - R_XK[MY][iElem])  / epsFD;
		res[E][iElem]   = v[E][iElem]   - alpha * dt * (aElem->u_t[E]   - R_XK[E][iElem])   / epsFD;
	}
}

/*
 * uses matrix free to solve the linear system, deltaX=0 is the initial guess
 * X0 is already stored in U
 */
void GMRES_M(double time, double dt, double alpha, double beta, double **B,
		double normB, double *abortCrit, double **delX)
{
	*abortCrit = epsGMRES * normB;

	double normR0 = normB;

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		R0[RHO][iElem] = - B[RHO][iElem];
		R0[MX][iElem]  = - B[MX][iElem];
		R0[MY][iElem]  = - B[MY][iElem];
		R0[E][iElem]   = - B[E][iElem];

		delX[RHO][iElem] = 0.0;
		delX[MX][iElem]  = 0.0;
		delX[MY][iElem]  = 0.0;
		delX[E][iElem]   = 0.0;

		V[0][RHO][iElem] = R0[RHO][iElem] / normR0;
		V[0][MX][iElem]  = R0[MX][iElem]  / normR0;
		V[0][MY][iElem]  = R0[MY][iElem]  / normR0;
		V[0][E][iElem]   = R0[E][iElem]   / normR0;
	}

	nInnerGMRES = 0;

	double gam[nKdim + 1];
	gam[0] = normR0;

	int m;
	double H[nKdim + 1][nKdim + 1], C[nKdim], S[nKdim];

	if (usePrecond) {
		buildMatrix(t, dt);
	}

	for (m = 0; m < nKdim; ++m) {
		nInnerGMRES++;

		if (usePrecond) {
			LUSGS(time, dt, V[m], Z[m]);
		} else {
			#pragma omp parallel for
			for (long iElem = 0; iElem < nElems; ++iElem) {
				Z[m][RHO][iElem] = V[m][RHO][iElem];
				Z[m][MX][iElem]  = V[m][MX][iElem];
				Z[m][MY][iElem]  = V[m][MY][iElem];
				Z[m][E][iElem]   = V[m][E][iElem];
			}
		}

		matrixVector(time, dt, alpha, beta, Z[m], W);

		/* Gram-Schmidt */
		for (int nn = 0; nn <= m; ++nn) {
			H[nn][m] = vectorDotProduct(V[nn], W);

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
					delX[RHO][iElem] += alp[nn] * Z[nn][RHO][iElem];
					delX[MX][iElem]  += alp[nn] * Z[nn][MX][iElem];
					delX[MY][iElem]  += alp[nn] * Z[nn][MY][iElem];
					delX[E][iElem]   += alp[nn] * Z[nn][E][iElem];
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

	printf("| GMRES not converged with %d iterations:\n", nInnerGMRES);
	printf("| Norm_R0: %g\n", fabs(gam[0]));
	printf("| Norm_R : %g\n", fabs(gam[m]));
	exit(1);
}

/*
 * free all memory that was allocated for
 */
void freeLinearSolver(void)
{
	if (isImplicit) {
		free(XK);
		free(R_XK);
		free(V);
		free(Z);
		free(R0);
		free(W);

		if (usePrecond) {
			free(deltaXstar);
			free(D);
			free(Dinv);
			free(dRdU);
		}
	}
}
