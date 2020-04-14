/*
 * boundary.c
 *
 * Created: Tue 24 Mar 2020 10:10:51 AM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "boundary.h"
#include "readInTools.h"
#include "equation.h"
#include "mesh.h"
#include "equationOfState.h"
#include "exactFunction.h"

/* extern variables */
boundary_t *firstBC;
int nBC;
bool isPeriodic;

/*
 * initialize boundary conditions
 */
void initBoundary(void)
{
	printf("\nInitializing Boundary Conditions:\n");

	nBC = getInt("nBC", NULL);

	firstBC = NULL;

	for (int iBC = 0; iBC < nBC; ++iBC) {
		boundary_t *aBC = malloc(sizeof(boundary_t));
		if (!aBC) {
			printf("| ERROR: could not allocate aBC\n");
			exit(1);
		}

		aBC->next = firstBC;
		firstBC = aBC;

		int intIn = getInt("BCtype", NULL);
		aBC->BCtype = intIn / 100;
		aBC->BCid = intIn % 100;

		double Ma, alpha, c, v;
		switch (aBC->BCtype) {
		case SLIPWALL:
			printf("| BC Type: Slip Wall\n");
			break;
		case WALL: {
			printf("| BC Type: No Slip Wall\n");

			bool isAdiabatic = getBool("adiabaticWall", NULL);
			if (isAdiabatic) {
				aBC->isAdiabatic = true;
				printf("| Adiabatic Wall\n");
			} else {
				aBC->isAdiabatic = false;

				/* temperature or heat flux prescibed ? */
				if (countKeys("wallTemperature", -1) > 1) {
					aBC->isTemperaturePrescribed = true;
					aBC->temperature = getDbl("wallTemperature", NULL);
				} else if (countKeys("wallHeatFlux", -1) > 1) {
					aBC->isTemperaturePrescribed = false;
					aBC->temperature = getDbl("wallHeatFlux", NULL);
				} else {
					printf("| ERROR: Not possible to prescirbe Wall Heat Flux and Wall Temperature\n");
					exit(1);
				}
			}
			break;
		}
		case INFLOW:
			printf("| BC Type: Inflow\n");

			aBC->pVar[RHO] = getDbl("rho", NULL);

			Ma = getDbl("mach", NULL);

			alpha = getDbl("alpha", NULL);

			aBC->pVar[P] = getDbl("pressure", NULL);

			c = sqrt(gamma * aBC->pVar[P] / aBC->pVar[RHO]);
			v = Ma * c;
			aBC->pVar[VX] = v * cos(alpha * pi / 180.0);
			aBC->pVar[VY] = v * sin(alpha * pi / 180.0);
			break;
		case CHARACTERISTIC:
			printf("| BC Type: Characteristic\n");

			aBC->pVar[RHO] = getDbl("rho", NULL);

			Ma = getDbl("mach", NULL);

			alpha = getDbl("alpha", NULL);

			aBC->pVar[P] = getDbl("pressure", NULL);

			c = sqrt(gamma * aBC->pVar[P] / aBC->pVar[RHO]);
			v = Ma * c;
			aBC->pVar[VX] = v * cos(alpha * pi / 180.0);
			aBC->pVar[VY] = v * sin(alpha * pi / 180.0);
			break;
		case OUTFLOW:
			printf("| BC Type: Outflow\n");
			break;
		case EXACTSOL:
			printf("| BC Type: Exact Function\n");

			aBC->exactFunc = getInt("BCexactFunc", NULL);
			break;
		case PERIODIC:
			printf("| BC Type: Periodic\n");

			aBC->connection = getDblArray("connection", NDIM, NULL);
			break;
		case PRESSURE_OUT:
			printf("| BC Type: Pressure Outlet\n");

			aBC->pVar[P] = getDbl("pressure", NULL);
			break;
		default:
			printf("| ERROR: Illegal boundary condition!\n");
			exit(1);
		}
	}
}

/*
 * set boundary value at x
 */
void boundary(side_t *aSide, double time, double int_pVar[NVAR],
		double ghost_pVar[NVAR], double x[NDIM])
{
	/* extract normal vector */
	double n[NDIM];
	n[X] = aSide->connection->n[X];
	n[Y] = aSide->connection->n[Y];

	/* determine type of boundary condition */
	switch (aSide->BC->BCtype) {
	case SLIPWALL: {
		double VXloc[NDIM], VYloc[NDIM];
		/* rotate into local coordinate system */
		VXloc[X] =   n[X] * int_pVar[VX] + n[Y] * int_pVar[VY];
		VYloc[X] = - n[Y] * int_pVar[VX] + n[X] * int_pVar[VY];

		/* mirror VX and extrapolate VY */
		VXloc[Y] = - VXloc[X];
		VYloc[Y] =   VYloc[X];

		/* backrotate into global coordinate system */
		ghost_pVar[VX] = n[X] * VXloc[Y] - n[Y] * VYloc[Y];
		ghost_pVar[VY] = n[Y] * VXloc[Y] + n[X] * VYloc[Y];

		/* scalar and derived conservative variables */
		ghost_pVar[RHO] = int_pVar[RHO];
		ghost_pVar[P]   = int_pVar[P];

		break;
	}
	case WALL: {
		double VXloc[NDIM], VYloc[NDIM];
		/* rotate into local coordinate system */
		VXloc[X] =   n[X] * int_pVar[VX] + n[Y] * int_pVar[VY];
		VYloc[X] = - n[Y] * int_pVar[VX] + n[X] * int_pVar[VY];

		/* mirror VX and extrapolate VY */
		VXloc[Y] = - VXloc[X];
		if (mu > 0.0) {
			VYloc[Y] = - VYloc[X];
		} else {
			VYloc[Y] =   VYloc[X];
		}

		/* backrotate into global coordinate system */
		ghost_pVar[VX] = n[X] * VXloc[Y] - n[Y] * VYloc[Y];
		ghost_pVar[VY] = n[Y] * VXloc[Y] + n[X] * VYloc[Y];

		/* scalar and derived conservative variables */
		ghost_pVar[RHO] = int_pVar[RHO];
		ghost_pVar[P]   = int_pVar[P];

		break;
	}
	case INFLOW:
		ghost_pVar[RHO] = aSide->BC->pVar[RHO];
		ghost_pVar[VX]  = aSide->BC->pVar[VX];
		ghost_pVar[VY]  = aSide->BC->pVar[VY];
		ghost_pVar[P]   = aSide->BC->pVar[P];

		break;
	case OUTFLOW:
		ghost_pVar[RHO] = int_pVar[RHO];
		ghost_pVar[VX]  = int_pVar[VX];
		ghost_pVar[VY]  = int_pVar[VY];
		ghost_pVar[P]   = int_pVar[P];

		break;
	case CHARACTERISTIC: {
		/* compute Eigenvalues of ghost cell */
		double c = sqrt(gamma * aSide->BC->pVar[P] / aSide->BC->pVar[RHO]);
		double v = n[X] * aSide->BC->pVar[VX] + n[Y] * aSide->BC->pVar[VY];

		/* rotate primitive state into local coordinate system */
		double int_pVarloc[NVAR], ghost_pVarloc[NVAR];
		int_pVarloc[RHO] = int_pVar[RHO];
		int_pVarloc[VX]  =   n[X] * int_pVar[VX] + n[Y] * int_pVar[VY];
		int_pVarloc[VY]  = - n[Y] * int_pVar[VX] + n[X] * int_pVar[VY];
		int_pVarloc[P]   = int_pVar[P];

		ghost_pVarloc[RHO] = aSide->BC->pVar[RHO];
		ghost_pVarloc[VX]  =   n[X] * aSide->BC->pVar[VX] + n[Y] * aSide->BC->pVar[VY];
		ghost_pVarloc[VY]  = - n[Y] * aSide->BC->pVar[VX] + n[X] * aSide->BC->pVar[VY];
		ghost_pVarloc[P]   = aSide->BC->pVar[P];

		/* compute conservative variables of both cells */
		double int_cVar[NVAR], ghost_cVar[NVAR];
		primCons(int_pVarloc, int_cVar);
		primCons(ghost_pVarloc, ghost_cVar);

		/* compute characteristic variables of inner and ghost cell */
		double int_charVar[3], ghost_charVar[3];
		consChar(int_charVar, int_cVar, int_pVar);
		consChar(ghost_charVar, ghost_cVar, int_pVar);

		/* determine characteristic state at boundary */
		if (v + c > 0.0) {
			ghost_charVar[2] = int_charVar[2];
		}
		if (v > 0.0) {
			ghost_charVar[1] = int_charVar[1];
		}
		if (v - c > 0.0) {
			ghost_charVar[0] = int_charVar[0];
		}

		/* determine the conservative state of the ghost cell */
		charCons(ghost_charVar, ghost_cVar, int_pVar);
		if (v > 0.0) {
			ghost_cVar[MY] = int_cVar[MY];
		}

		/* determine the primitive state of the ghost cell */
		consPrim(ghost_cVar, ghost_pVar);

		/* rotate the primitive state into the global coordinate system */
		double VXloc = ghost_pVar[VX];
		double VYloc = ghost_pVar[VY];
		ghost_pVar[VX] = n[X] * VXloc - n[Y] * VYloc;
		ghost_pVar[VY] = n[Y] * VXloc + n[X] * VYloc;

		break;
	}
	case EXACTSOL:
		exactFunc(aSide->BC->exactFunc, x, time, ghost_pVar);
		break;
	case PRESSURE_OUT: {
		double c = sqrt(gamma * int_pVar[P] / int_pVar[RHO]);
		double v = n[X] * int_pVar[VX] + n[Y] * int_pVar[VY];

		double p;
		if (v / c < 1.0) {
			p = aSide->BC->pVar[P];
		} else {
			p = int_pVar[P];
		}

		ghost_pVar[RHO] = int_pVar[RHO] * p / int_pVar[P];
		ghost_pVar[VX]  = int_pVar[VX];
		ghost_pVar[VY]  = int_pVar[VY];
		ghost_pVar[P]   = p;

		break;
	}
	}
}

/*
 * set the ghost values at sides
 */
void setBCatSides(double time)
{
	#pragma omp parallel for
	for (long iSide = 0; iSide < nBCsides; ++iSide) {
		side_t *gSide = BCside[iSide];
		side_t *aSide = gSide->connection;
		elem_t *aElem = aSide->elem;

		double x[NDIM];
		x[X] = aSide->GP[X] + aElem->bary[X];
		x[Y] = aSide->GP[Y] + aElem->bary[Y];

		boundary(gSide, time, aSide->pVar, gSide->pVar, x);
	}
}

/*
 * set the ghost values at elements
 */
void setBCatBarys(double time)
{
	#pragma omp parallel for
	for (long iSide = 0; iSide < nBCsides; ++iSide) {
		side_t *gSide = BCside[iSide];
		elem_t *gElem = gSide->elem;
		side_t *aSide = gSide->connection;
		elem_t *aElem = aSide->elem;

		boundary(gSide, time, aElem->pVar, gElem->pVar, gElem->bary);
	}
}

/*
 * free all memory that was allocated for the boundary conditions
 */
void freeBoundary(void)
{
	boundary_t *aBC = firstBC;
	while (aBC) {
		if (aBC->next) {
			boundary_t *tmp = aBC;
			aBC = aBC->next;
			free(tmp);
		} else {
			free(aBC);
			break;
		}
	}
	free(BCside);
}
