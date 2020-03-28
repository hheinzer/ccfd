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
