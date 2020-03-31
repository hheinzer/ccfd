/*
 * fluxCalculation.c
 *
 * Created: Tue 31 Mar 2020 05:18:40 PM CEST
 * Author : hhh
 */

#include <math.h>

#include "main.h"
#include "mesh.h"
#include "equation.h"
#include "exactRiemann.h"

/*
 * Godunov flux, which is the exact flux
 */
void flux_god(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{
	double cL = sqrt(gamma * pL / rhoL);
	double cR = sqrt(gamma * pR / rhoR);

	double rho, vx, p;
	exactRiemann(rhoL, rhoR, &rho, vxL, vxR, &vx, pL, pR, &p, cL, cR, 0.0);

	double vy;
	if (vx > 0.0) {
		vy = vyL;
	} else {
		vy = vyR;
	}

	fluxLoc[0] = rho * vx;
	fluxLoc[1] = rho * vx * vy + p;
	fluxLoc[2] = rho * vx * vy;
	fluxLoc[3] = vx * (gamma / gamma1 * p + 0.5 * rho * (vx * vx + vy * vy));
}

/*
 * Roe flux
 */
void flux_roe(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * HLL flux
 */
void flux_hll(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * HLLE flux
 */
void flux_hlle(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * HLLC flux
 */
void flux_hllc(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * Lax-Friedrich flux
 */
void flux_lxf(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * Steger-Warming flux
 */
void flux_stw(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * Central flux
 */
void flux_cen(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * AUSMD flux
 */
void flux_ausmd(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * AUSMDV flux
 */
void flux_ausmdv(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * Van Leer flux
 */
void flux_vanleer(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[NVAR])
{

}

/*
 * select the convective flux
 */
void convectiveFlux(double rhoL, double rhoR,
		    double vxL,  double vxR,
		    double vyL,  double vyR,
		    double pL,   double pR,
		    double fluxLoc[NVAR])
{
	switch (iFlux) {
	case GOD:
		flux_god(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case ROE:
		flux_roe(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case HLL:
		flux_hll(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case HLLE:
		flux_hlle(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case HLLC:
		flux_hllc(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case LXF:
		flux_lxf(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case STW:
		flux_stw(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case CEN:
		flux_cen(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case AUSMD:
		flux_ausmd(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case AUSMDV:
		flux_ausmdv(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	case VANLEER:
		flux_vanleer(rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR, fluxLoc);
		break;
	}
}

/*
 * calculation of left and right state, the velocity vector is transformed into
 * the normal system of the cell interfaces
 * finishes with backrotation
 */
void fluxCalculation(void)
{
	#pragma omp parallel for
	for (long iSide = 0; iSide < nSides; ++iSide) {
		side_t *aSide = side[iSide];

		/* extract left state */
		double pVar[NVAR];
		pVar[RHO] = aSide->pVar[RHO];
		pVar[VX]  = aSide->pVar[VX];
		pVar[VY]  = aSide->pVar[VY];
		pVar[P]   = aSide->pVar[P];

		/* rotate it into normal direction */
		double pVarL[NVAR];
		pVarL[RHO] = pVar[RHO];
		pVarL[VX]  =   aSide->n[X] * pVar[VX] + aSide->n[Y] * pVar[VY];
		pVarL[VY]  = - aSide->n[Y] * pVar[VX] + aSide->n[X] * pVar[VY];
		pVarL[P]   = pVar[P];

		/* extract right state */
		pVar[RHO] = aSide->connection->pVar[RHO];
		pVar[VX]  = aSide->connection->pVar[VX];
		pVar[VY]  = aSide->connection->pVar[VY];
		pVar[P]   = aSide->connection->pVar[P];

		/* rotate it into normal direction */
		double pVarR[NVAR];
		pVarR[RHO] = pVar[RHO];
		pVarR[VX]  =   aSide->n[X] * pVar[VX] + aSide->n[Y] * pVar[VY];
		pVarR[VY]  = - aSide->n[Y] * pVar[VX] + aSide->n[X] * pVar[VY];
		pVarR[P]   = pVar[P];

		/* calculate flux */
		double fluxLoc[NVAR];
		convectiveFlux(pVarL[RHO], pVarR[RHO],
			       pVarL[VX],  pVarR[VX],
			       pVarL[VY],  pVarR[VY],
			       pVarL[P],   pVarR[P],
			       fluxLoc);
	}
}
