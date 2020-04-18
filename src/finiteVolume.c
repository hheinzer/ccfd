/** \file
 *
 * \brief Finite volume time derivative functions
 *
 * \author hhh
 * \date Fri 27 Mar 2020 05:10:43 PM CET
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "finiteVolume.h"
#include "mesh.h"
#include "readInTools.h"
#include "reconstruction.h"
#include "timeDiscretization.h"
#include "fluxCalculation.h"
#include "equation.h"
#include "source.h"

/* extern variables */
int spatialOrder;			/**< the spacial order to be used */
int fluxFunction;			/**< the flux function to be used */

/**
 * \brief Initialize the finite volume method
 */
void initFV(void)
{
	printf("\nInitialize Spacial Discretization:\n");
	spatialOrder = getInt("spatialOrder", "1");
	if (spatialOrder > 2) {
		printf("| ERROR: Spatial Discretization Order must be 1 or 2\n");
		exit(1);
	} else if (spatialOrder == 2) {
		limiter = getInt("limiter", "1");
		switch (limiter) {
		case BARTHJESPERSEN:
			printf("| Limiter: Barth & Jesperson\n");
			break;
		case VENKATAKRISHNAN:
			printf("| Limiter: Venkatakrishnan\n");
			venk_k = getDbl("venk_k", "1");

			for (long iElem = 0; iElem < nElems; ++iElem) {
				elem[iElem]->venkEps_sq =
					(venk_k * sqrt(elem[iElem]->area)) *
					(venk_k * sqrt(elem[iElem]->area)) *
					(venk_k * sqrt(elem[iElem]->area));
			}
			break;
		default:
			printf("| ERROR: Limiter must be either 1 or 2\n");
			exit(1);
		}
	}

	for (long iElem = 0; iElem < nElems; ++iElem) {
		for (int iVar = 0; iVar < NVAR; ++iVar) {
			elem[iElem]->source[iVar] = 0.0;
		}
	}
}

/**
 * \brief Perform the spacial operator of the finite volume scheme
 *
 * First, the local time step is calculated, then spacial gradients inside of
 * the cells are reconstructed. Following that, the boundary conditions at the
 * sides are applied and the numerical flux is calculated, using the specified
 * flux function. Finally, the source term is evaluated and the time derivatives
 * of all the elements are calculated.
 *
 * \param[in] time Calculation time at which to perform the finite volume differentiation
 */
void fvTimeDerivative(double time)
{
	/* set dt for boundary condition calculation */
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		aElem->dtLoc = 0.5 * aElem->dt * (timeOrder - 1);
	}

	spatialReconstruction(time);
	setBCatSides(time);
	fluxCalculation();

	if (doCalcSource) {
		calcSource(time);
	}

	/* time update of the conservative variables */
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		side_t *aSide = aElem->firstSide;

		aElem->u_t[RHO] = 0.0;
		aElem->u_t[VX]  = 0.0;
		aElem->u_t[VY]  = 0.0;
		aElem->u_t[E]   = 0.0;

		while (aSide) {
			aElem->u_t[RHO] += aSide->flux[RHO];
			aElem->u_t[VX]  += aSide->flux[VX];
			aElem->u_t[VY]  += aSide->flux[VY];
			aElem->u_t[E]   += aSide->flux[E];

			aSide = aSide->nextElemSide;
		}

		/* source term contribution */
		aElem->u_t[RHO] = (aElem->source[RHO] - aElem->u_t[RHO]) * aElem->areaq;
		aElem->u_t[VX]  = (aElem->source[VX]  - aElem->u_t[VX])  * aElem->areaq;
		aElem->u_t[VY]  = (aElem->source[VY]  - aElem->u_t[VY])  * aElem->areaq;
		aElem->u_t[E]   = (aElem->source[E]   - aElem->u_t[E])   * aElem->areaq;
	}
}
