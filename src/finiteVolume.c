/*
 * finiteVolume.c
 *
 * Created: Fri 27 Mar 2020 05:10:43 PM CET
 * Author : hhh
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
#include "source.h"

/* extern variables */
int spatialOrder;
int fluxFunction;

/*
 * initialize the finite volume method
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
		for (int i = 0; i < NVAR; ++i) {
			elem[iElem]->source[i] = 0.0;
		}
	}
}

/*
 * calling all subroutines to perform the spacial operator of the FV scheme
 * Result: residual vecto
 */
void FVtimeDerivative(double time, long iter)
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
	calcSource(time);

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

		/* source term */
		aElem->u_t[RHO] = (aElem->source[RHO] - aElem->u_t[RHO]) * aElem->areaq;
		aElem->u_t[VX]  = (aElem->source[VX]  - aElem->u_t[VX])  * aElem->areaq;
		aElem->u_t[VY]  = (aElem->source[VY]  - aElem->u_t[VY])  * aElem->areaq;
		aElem->u_t[E]   = (aElem->source[E]   - aElem->u_t[E])   * aElem->areaq;

	}
}
