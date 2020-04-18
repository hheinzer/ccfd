/** \file
 *
 * \brief Contains the reconstruction and limiter functions
 *
 * \author hhh
 * \date Sat 28 Mar 2020 10:17:02 AM CET
 */

#include <math.h>

#include "main.h"
#include "reconstruction.h"
#include "finiteVolume.h"
#include "mesh.h"

/* extern variables */
int limiter;				/**< limiter selection */
double venk_k;				/**< constant for Venkatakrishnan limiter */

/**
 * \brief Limiter after Barth & Jespersen
 * \note 2D, unstructured limiter
 * \param[in] *aElem Pointer to an element
 */
void limiterBarthJespersen(elem_t *aElem)
{
	/* determine uMin and uMax */
	double uMin[NVAR], uMax[NVAR];
	uMin[RHO] = uMax[RHO] = aElem->pVar[RHO];
	uMin[VX]  = uMax[VX]  = aElem->pVar[VX];
	uMin[VY]  = uMax[VY]  = aElem->pVar[VY];
	uMin[P]   = uMax[P]   = aElem->pVar[P];

	side_t *aSide = aElem->firstSide;
	while (aSide) {
		uMin[RHO] = fmin(uMin[RHO], aSide->connection->elem->pVar[RHO]);
		uMin[VX]  = fmin(uMin[VX],  aSide->connection->elem->pVar[VX]);
		uMin[VY]  = fmin(uMin[VY],  aSide->connection->elem->pVar[VY]);
		uMin[P]   = fmin(uMin[P],   aSide->connection->elem->pVar[P]);

		uMax[RHO] = fmax(uMax[RHO], aSide->connection->elem->pVar[RHO]);
		uMax[VX]  = fmax(uMax[VX],  aSide->connection->elem->pVar[VX]);
		uMax[VY]  = fmax(uMax[VY],  aSide->connection->elem->pVar[VY]);
		uMax[P]   = fmax(uMax[P],   aSide->connection->elem->pVar[P]);

		aSide = aSide->nextElemSide;
	}

	double minDiff[NVAR], maxDiff[NVAR];
	minDiff[RHO] = uMin[RHO] - aElem->pVar[RHO];
	minDiff[VX]  = uMin[VX]  - aElem->pVar[VX];
	minDiff[VY]  = uMin[VY]  - aElem->pVar[VY];
	minDiff[P]   = uMin[P]   - aElem->pVar[P];

	maxDiff[RHO] = uMax[RHO] - aElem->pVar[RHO];
	maxDiff[VX]  = uMax[VX]  - aElem->pVar[VX];
	maxDiff[VY]  = uMax[VY]  - aElem->pVar[VY];
	maxDiff[P]   = uMax[P]   - aElem->pVar[P];

	/* loop over all edges: determine phi */
	double phi[NVAR] = {1.0, 1.0, 1.0, 1.0}, phiLoc[NVAR], uDiff[NVAR];
	aSide = aElem->firstSide;
	while (aSide) {
		phiLoc[RHO] = 1.0;
		phiLoc[VX]  = 1.0;
		phiLoc[VY]  = 1.0;
		phiLoc[P]   = 1.0;
		for (int iVar = 0; iVar < NVAR; ++iVar) {
			uDiff[iVar] = aElem->u_x[iVar] * aSide->GP[X]
				    + aElem->u_y[iVar] * aSide->GP[Y];
			if (uDiff[iVar] > 0.0) {
				phiLoc[iVar] = fmin(1.0, maxDiff[iVar] / uDiff[iVar]);
			} else if (uDiff[iVar] < 0.0) {
				phiLoc[iVar] = fmin(1.0, minDiff[iVar] / uDiff[iVar]);
			}
		}

		phi[RHO] = fmin(phi[RHO], phiLoc[RHO]);
		phi[VX]  = fmin(phi[VX],  phiLoc[VX]);
		phi[VY]  = fmin(phi[VY],  phiLoc[VY]);
		phi[P]   = fmin(phi[P],   phiLoc[P]);

		aSide = aSide->nextElemSide;
	}

	/* compute limited gradients */
	aElem->u_x[RHO] *= phi[RHO];
	aElem->u_x[VX]  *= phi[VX];
	aElem->u_x[VY]  *= phi[VY];
	aElem->u_x[P]   *= phi[P];

	aElem->u_y[RHO] *= phi[RHO];
	aElem->u_y[VX]  *= phi[VX];
	aElem->u_y[VY]  *= phi[VY];
	aElem->u_y[P]   *= phi[P];
}

/**
 * \brief Limiter after Venkatakrishnan, with additional limiting parameter k
 * \note 2D, unstructured limiter
 * \param[in] *aElem Pointer to an element
 */
void limiterVenkatakrishnan(elem_t *aElem)
{
	/* determine uMin and uMax */
	double uMin[NVAR], uMax[NVAR];
	uMin[RHO] = uMax[RHO] = aElem->pVar[RHO];
	uMin[VX]  = uMax[VX]  = aElem->pVar[VX];
	uMin[VY]  = uMax[VY]  = aElem->pVar[VY];
	uMin[P]   = uMax[P]   = aElem->pVar[P];

	side_t *aSide = aElem->firstSide;
	while (aSide) {
		uMin[RHO] = fmin(uMin[RHO], aSide->connection->elem->pVar[RHO]);
		uMin[VX]  = fmin(uMin[VX],  aSide->connection->elem->pVar[VX]);
		uMin[VY]  = fmin(uMin[VY],  aSide->connection->elem->pVar[VY]);
		uMin[P]   = fmin(uMin[P],   aSide->connection->elem->pVar[P]);

		uMax[RHO] = fmax(uMax[RHO], aSide->connection->elem->pVar[RHO]);
		uMax[VX]  = fmax(uMax[VX],  aSide->connection->elem->pVar[VX]);
		uMax[VY]  = fmax(uMax[VY],  aSide->connection->elem->pVar[VY]);
		uMax[P]   = fmax(uMax[P],   aSide->connection->elem->pVar[P]);

		aSide = aSide->nextElemSide;
	}

	double minDiff[NVAR], maxDiff[NVAR], minDiffsq[NVAR], maxDiffsq[NVAR];
	minDiff[RHO] = uMin[RHO] - aElem->pVar[RHO];
	minDiff[VX]  = uMin[VX]  - aElem->pVar[VX];
	minDiff[VY]  = uMin[VY]  - aElem->pVar[VY];
	minDiff[P]   = uMin[P]   - aElem->pVar[P];

	maxDiff[RHO] = uMax[RHO] - aElem->pVar[RHO];
	maxDiff[VX]  = uMax[VX]  - aElem->pVar[VX];
	maxDiff[VY]  = uMax[VY]  - aElem->pVar[VY];
	maxDiff[P]   = uMax[P]   - aElem->pVar[P];

	minDiffsq[RHO] = minDiff[RHO] * minDiff[RHO];
	minDiffsq[VX]  = minDiff[VX]  * minDiff[VX];
	minDiffsq[VY]  = minDiff[VY]  * minDiff[VY];
	minDiffsq[P]   = minDiff[P]   * minDiff[P];

	maxDiffsq[RHO] = maxDiff[RHO] * maxDiff[RHO];
	maxDiffsq[VX]  = maxDiff[VX]  * maxDiff[VX];
	maxDiffsq[VY]  = maxDiff[VY]  * maxDiff[VY];
	maxDiffsq[P]   = maxDiff[P]   * maxDiff[P];

	/* loop over all edges: determine phi */
	double phi[NVAR] = {1.0, 1.0, 1.0, 1.0}, phiLoc[NVAR], uDiff[NVAR], uDiffsq[NVAR];
	aSide = aElem->firstSide;
	while (aSide) {
		phiLoc[RHO] = 1.0;
		phiLoc[VX]  = 1.0;
		phiLoc[VY]  = 1.0;
		phiLoc[P]   = 1.0;
		for (int iVar = 0; iVar < NVAR; ++iVar) {
			uDiff[iVar] = aElem->u_x[iVar] * aSide->GP[X]
				    + aElem->u_y[iVar] * aSide->GP[Y];
			uDiffsq[iVar] = uDiff[iVar] * uDiff[iVar];

			if (uDiff[iVar] > 0.0) {
				phiLoc[iVar] = 1.0 / uDiff[iVar] * (((maxDiffsq[iVar] + aElem->venkEps_sq) * uDiff[iVar]
							+ 2.0 * uDiffsq[iVar] * maxDiff[iVar])
						/ (maxDiffsq[iVar] + 2.0 * uDiffsq[iVar] + uDiff[iVar]
							* maxDiff[iVar] + aElem->venkEps_sq));
			} else if (uDiff[iVar] < 0.0) {
				phiLoc[iVar] = 1.0 / uDiff[iVar] * (((minDiffsq[iVar] + aElem->venkEps_sq) * uDiff[iVar]
							+ 2.0 * uDiffsq[iVar] * minDiff[iVar])
						/ (minDiffsq[iVar] + 2.0 * uDiffsq[iVar] + uDiff[iVar]
							* minDiff[iVar] + aElem->venkEps_sq));
			}
		}

		phi[RHO] = fmin(phi[RHO], phiLoc[RHO]);
		phi[VX]  = fmin(phi[VX],  phiLoc[VX]);
		phi[VY]  = fmin(phi[VY],  phiLoc[VY]);
		phi[P]   = fmin(phi[P],   phiLoc[P]);

		aSide = aSide->nextElemSide;
	}

	/* compute limited gradients */
	aElem->u_x[RHO] *= phi[RHO];
	aElem->u_x[VX]  *= phi[VX];
	aElem->u_x[VY]  *= phi[VY];
	aElem->u_x[P]   *= phi[P];

	aElem->u_y[RHO] *= phi[RHO];
	aElem->u_y[VX]  *= phi[VX];
	aElem->u_y[VY]  *= phi[VY];
	aElem->u_y[P]   *= phi[P];
}

/**
 * \brief Compute the gradients of dU/dx
 * \param[in] time Calculation time at which to perform the spatial reconstruction
 */
void spatialReconstruction(double time)
{
	if (spatialOrder == 1) {
		/* set side states to be equal to mean value */
		#pragma omp parallel for
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];
			side_t *aSide = aElem->firstSide;
			aElem->u_x[RHO] = 0.0;
			aElem->u_x[VX]  = 0.0;
			aElem->u_x[VY]  = 0.0;
			aElem->u_x[P]   = 0.0;

			aElem->u_y[RHO] = 0.0;
			aElem->u_y[VX]  = 0.0;
			aElem->u_y[VY]  = 0.0;
			aElem->u_y[P]   = 0.0;

			aElem->u_t[RHO] = 0.0;
			aElem->u_t[VX]  = 0.0;
			aElem->u_t[VY]  = 0.0;
			aElem->u_t[P]   = 0.0;

			while (aSide) {
				aSide->pVar[RHO] = aElem->pVar[RHO];
				aSide->pVar[VX]  = aElem->pVar[VX];
				aSide->pVar[VY]  = aElem->pVar[VY];
				aSide->pVar[P]   = aElem->pVar[P];
				aSide = aSide->nextElemSide;
			}
		}
	} else {
		/* reconstruction of values at side GPs */
		#pragma omp parallel for
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];
			aElem->u_x[RHO] = 0.0;
			aElem->u_x[VX]  = 0.0;
			aElem->u_x[VY]  = 0.0;
			aElem->u_x[P]   = 0.0;

			aElem->u_y[RHO] = 0.0;
			aElem->u_y[VX]  = 0.0;
			aElem->u_y[VY]  = 0.0;
			aElem->u_y[P]   = 0.0;

			aElem->u_t[RHO] = 0.0;
			aElem->u_t[VX]  = 0.0;
			aElem->u_t[VY]  = 0.0;
			aElem->u_t[P]   = 0.0;
		}

		setBCatBarys(time);

		#pragma omp parallel for
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];
			side_t *aSide = aElem->firstSide;
			while (aSide) {
				double pDiff[NVAR];
				pDiff[RHO] = aSide->connection->elem->pVar[RHO] - aSide->elem->pVar[RHO];
				pDiff[VX]  = aSide->connection->elem->pVar[VX]  - aSide->elem->pVar[VX];
				pDiff[VY]  = aSide->connection->elem->pVar[VY]  - aSide->elem->pVar[VY];
				pDiff[P]   = aSide->connection->elem->pVar[P]   - aSide->elem->pVar[P];

				aElem->u_x[RHO] += aSide->w[X] * pDiff[RHO];
				aElem->u_x[VX]  += aSide->w[X] * pDiff[VX];
				aElem->u_x[VY]  += aSide->w[X] * pDiff[VY];
				aElem->u_x[P]   += aSide->w[X] * pDiff[P];

				aElem->u_y[RHO] += aSide->w[Y] * pDiff[RHO];
				aElem->u_y[VX]  += aSide->w[Y] * pDiff[VX];
				aElem->u_y[VY]  += aSide->w[Y] * pDiff[VY];
				aElem->u_y[P]   += aSide->w[Y] * pDiff[P];

				aSide = aSide->nextElemSide;
			}
		}

		/* limit gradients and reconstruct values at side GPs */
		#pragma omp parallel for
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];

			/* limit gradient */
			switch (limiter) {
			case BARTHJESPERSEN:
				limiterBarthJespersen(aElem);
				break;
			case VENKATAKRISHNAN:
				limiterVenkatakrishnan(aElem);
				break;
			}

			/* reconstruct values at side GPs */
			side_t *aSide = aElem->firstSide;
			while (aSide) {
				double dx = aSide->GP[X];
				double dy = aSide->GP[Y];

				aSide->pVar[RHO] = aElem->pVar[RHO]
					+ dx * aElem->u_x[RHO] + dy * aElem->u_y[RHO];

				aSide->pVar[VX]  = aElem->pVar[VX]
					+ dx * aElem->u_x[VX]  + dy * aElem->u_y[VX];

				aSide->pVar[VY]  = aElem->pVar[VY]
					+ dx * aElem->u_x[VY]  + dy * aElem->u_y[VY];

				aSide->pVar[P]   = aElem->pVar[P]
					+ dx * aElem->u_x[P]   + dy * aElem->u_y[P];

				aSide = aSide->nextElemSide;
			}
		}
	}
}
