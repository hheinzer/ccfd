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
	      double fluxLoc[4])
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
	      double fluxLoc[4])
{
	/* calculate left/right enthalpy */
	double mxL = rhoL * vxL;
	double mxR = rhoR * vxR;
	double myL = rhoL * vyL;
	double myR = rhoR * vyR;
	double eR = gamma1q * pR + 0.5 * rhoR * (vxR * vxR + vyR * vyR);
	double eL = gamma1q * pL + 0.5 * rhoL * (vxL * vxL + vyL * vyL);

	/* calculate left/right enthalpy */
	double HR = (eR + pR) / rhoR;
	double HL = (eL + pL) / rhoL;

	/* calculate sqrt(rho) */
	double rhoSqR = sqrt(rhoR);
	double rhoSqL = sqrt(rhoL);
	double rhoSqQsum = 1.0 / (rhoSqL + rhoSqR);

	/* calculate Roe mean values */
	double vxBar = (rhoSqR * vxR + rhoSqL * vxL) * rhoSqQsum;
	double vyBar = (rhoSqR * vyR + rhoSqL * vyL) * rhoSqQsum;
	double Hbar  = (rhoSqR *  HR + rhoSqL *  HL) * rhoSqQsum;
	double cBar  = sqrt(gamma1 * (Hbar - 0.5 * (vxBar * vxBar + vyBar * vyBar)));

	/* calculate mean Eigenvalues */
	double a1 = vxBar - cBar;
	double a2 = vxBar;
	double a3 = vxBar;
	double a4 = vxBar + cBar;

	/* calculate mean eigenvectors */
	double r1[4] = {1.0, a1, vyBar, Hbar - vxBar * cBar};
	double r2[4] = {1.0, vxBar, vyBar, 0.5 * (vxBar * vxBar + vyBar * vyBar)};
	double r3[4] = {0.0, 0.0, 1.0, vyBar};
	double r4[4] = {1.0, a4, vyBar, Hbar + vxBar * cBar};

	/* calculate differences */
	double delRho = rhoR - rhoL;
	double delMx  = mxR  - mxL;
	double delMy  = myR  - myL;
	double delE   = eR   - eL;
	double delEq  = delE - (delMy - vyBar * delRho) * vyBar;

	/* calculate wave strenght */
	double cBarQ = 1.0 / cBar;
	double gam2 = - gamma1 * cBarQ * cBarQ * (delRho * (vxBar * vxBar - Hbar)
			+ delEq - delMx * vxBar);
	double gam1 = - 0.5 * cBarQ * (delMx - delRho * (vxBar + cBar)) - 0.5 * gam2;
	double gam4 = delRho - gam1 - gam2;
	double gam3 = delMy - vyBar * delRho;

	/* calculate physical fluxes */
	double fR[4] = {mxR, mxR * vxR + pR, mxR * vyR, vxR * (eR + pR)};
	double fL[4] = {mxL, mxL * vxL + pL, mxL * vyL, vxL * (eL + pL)};

	/* calculate Row flux */
	for (int i = 0; i < 0; ++i) {
		fluxLoc[i] = 0.5 * (fR[i] + fL[i]
				- gam1 * fabs(a1) * r1[i]
				- gam2 * fabs(a2) * r2[i]
				- gam3 * fabs(a3) * r3[i]
				- gam4 * fabs(a4) * r4[i]);
	}
}

/*
 * HLL flux
 */
void flux_hll(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{
	/* calculation of auxiliary values */
	double rhoLq = 1.0 / rhoL;
	double rhoRq = 1.0 / rhoR;
	double rhoSqL = sqrt(rhoL);
	double rhoSqR = sqrt(rhoR);
	double rhoSqQsum = 1.0 / (rhoSqL + rhoSqR);

	/* calculate energies */
	double eR = gamma1q * pR + 0.5 * rhoR * (vxR * vxR + vyR * vyR);
	double eL = gamma1q * pL + 0.5 * rhoL * (vxL * vxL + vyL * vyL);

	/* calculate left/right conservative state vector */
	double uL[NVAR] = {rhoL, rhoL * vxL, rhoL * vyL, eL};
	double uR[NVAR] = {rhoR, rhoR * vxR, rhoR * vyR, eR};

	/* calculate flux in left/right cell */
	double fL[NVAR] = {uL[MX], uL[MX] * vxL + pL, uL[MX] * vyL, vxL * (eL + pL)};
	double fR[NVAR] = {uR[MX], uR[MX] * vxR + pR, uR[MX] * vyR, vxR * (eR + pR)};

	/* calculation of speed of sounds */
	double cL = sqrt(gamma * pL * rhoLq);
	double cR = sqrt(gamma * pR * rhoRq);

	/* calculation of left/right enthalpy */
	double HL = (eL + pL) * rhoLq;
	double HR = (eR + pR) * rhoRq;

	/* calculation Row mean values */
	double uM = (rhoSqR * vxR + rhoSqL * vxL) * rhoSqQsum;
	double vM = (rhoSqR * vyR + rhoSqL * vyL) * rhoSqQsum;
	double HM = (rhoSqR *  HR + rhoSqL *  HL) * rhoSqQsum;
	double cM = sqrt(gamma1 * (HM - 0.5 * (uM * uM + vM * vM)));

	/* calculation signal speeds */
	double arp = fmax(vxR + cR, uM + cM);
	double alm = fmin(vxL - cL, uM - cM);
	double arpAlmQ = 1.0 / (arp - alm);

	/* calculation HLL flux */
	if (alm > 0.0) {
		fluxLoc[0] = fL[0];
		fluxLoc[1] = fL[1];
		fluxLoc[2] = fL[2];
		fluxLoc[3] = fL[3];
	} else if (arp < 0.0) {
		fluxLoc[0] = fR[0];
		fluxLoc[1] = fR[1];
		fluxLoc[2] = fR[2];
		fluxLoc[3] = fR[3];
	} else {
		fluxLoc[0] = (arp * fL[0] - alm * fR[0]) * arpAlmQ
			   + (arp * alm) * arpAlmQ * (uR[0] - uL[0]);
		fluxLoc[1] = (arp * fL[1] - alm * fR[1]) * arpAlmQ
			   + (arp * alm) * arpAlmQ * (uR[1] - uL[1]);
		fluxLoc[2] = (arp * fL[2] - alm * fR[2]) * arpAlmQ
			   + (arp * alm) * arpAlmQ * (uR[2] - uL[2]);
		fluxLoc[3] = (arp * fL[3] - alm * fR[3]) * arpAlmQ
			   + (arp * alm) * arpAlmQ * (uR[3] - uL[3]);
	}
}

/*
 * HLLE flux
 */
void flux_hlle(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{
	/* calculation of auxiliary values */
	double rhoLq = 1.0 / rhoL;
	double rhoRq = 1.0 / rhoR;
	double rhoSqL = sqrt(rhoL);
	double rhoSqR = sqrt(rhoR);
	double rhoSqQsum = 1.0 / (rhoSqL + rhoSqR);

	/* calculate energies */
	double eR = gamma1q * pR + 0.5 * rhoR * (vxR * vxR + vyR * vyR);
	double eL = gamma1q * pL + 0.5 * rhoL * (vxL * vxL + vyL * vyL);

	/* calculate left/right conservative state vector */
	double uL[NVAR] = {rhoL, rhoL * vxL, rhoL * vyL, eL};
	double uR[NVAR] = {rhoR, rhoR * vxR, rhoR * vyR, eR};

	/* calculate flux in left/right cell */
	double fL[NVAR] = {uL[MX], uL[MX] * vxL + pL, uL[MX] * vyL, vxL * (eL + pL)};
	double fR[NVAR] = {uR[MX], uR[MX] * vxR + pR, uR[MX] * vyR, vxR * (eR + pR)};

	/* calculation of speed of sounds */
	double cL = sqrt(gamma * pL * rhoLq);
	double cR = sqrt(gamma * pR * rhoRq);

	/* calculation Row mean values */
	double uM = (rhoSqR * vxR + rhoSqL * vxL) * rhoSqQsum;

	/* signal speeds, version of Einfeld paper */
	double eta2 = 0.5 * rhoSqR * rhoSqL / (rhoSqR + rhoSqL) / (rhoSqR + rhoSqL);
	double d = sqrt((rhoSqR * cR * cR + rhoSqL * cL * cL) * rhoSqQsum
			+ eta2 * (vxR - vxL) * (vxR - vxL));
	double arp = fmax(vxR + cR, uM + d);
	double alm = fmin(vxL - cL, uM - d);
	double arpAlmQ = 1.0 / (arp - alm);

	/* calculation HLLE flux */
	if (alm > 0.0) {
		fluxLoc[0] = fL[0];
		fluxLoc[1] = fL[1];
		fluxLoc[2] = fL[2];
		fluxLoc[3] = fL[3];
	} else if (arp < 0.0) {
		fluxLoc[0] = fR[0];
		fluxLoc[1] = fR[1];
		fluxLoc[2] = fR[2];
		fluxLoc[3] = fR[3];
	} else {
		fluxLoc[0] = (arp * fL[0] - alm * fR[0]) * arpAlmQ
			   + (arp * alm) * arpAlmQ * (uR[0] - uL[0]);
		fluxLoc[1] = (arp * fL[1] - alm * fR[1]) * arpAlmQ
			   + (arp * alm) * arpAlmQ * (uR[1] - uL[1]);
		fluxLoc[2] = (arp * fL[2] - alm * fR[2]) * arpAlmQ
			   + (arp * alm) * arpAlmQ * (uR[2] - uL[2]);
		fluxLoc[3] = (arp * fL[3] - alm * fR[3]) * arpAlmQ
			   + (arp * alm) * arpAlmQ * (uR[3] - uL[3]);
	}
}

/*
 * HLLC flux
 */
void flux_hllc(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{
	/* calculation of auxiliary values */
	double rhoLq = 1.0 / rhoL;
	double rhoRq = 1.0 / rhoR;
	double rhoSqL = sqrt(rhoL);
	double rhoSqR = sqrt(rhoR);
	double rhoSqQsum = 1.0 / (rhoSqL + rhoSqR);

	/* calculate energies */
	double eR = gamma1q * pR + 0.5 * rhoR * (vxR * vxR + vyR * vyR);
	double eL = gamma1q * pL + 0.5 * rhoL * (vxL * vxL + vyL * vyL);

	/* calculate left/right conservative state vector */
	double uL[NVAR] = {rhoL, rhoL * vxL, rhoL * vyL, eL};
	double uR[NVAR] = {rhoR, rhoR * vxR, rhoR * vyR, eR};

	/* calculate flux in left/right cell */
	double fL[NVAR] = {uL[MX], uL[MX] * vxL + pL, uL[MX] * vyL, vxL * (eL + pL)};
	double fR[NVAR] = {uR[MX], uR[MX] * vxR + pR, uR[MX] * vyR, vxR * (eR + pR)};

	/* calculation of speed of sounds */
	double cL = sqrt(gamma * pL * rhoLq);
	double cR = sqrt(gamma * pR * rhoRq);

	/* calculation of left/right enthalpy */
	double HL = (eL + pL) * rhoLq;
	double HR = (eR + pR) * rhoRq;

	/* calculation Row mean values */
	double uM = (rhoSqR * vxR + rhoSqL * vxL) * rhoSqQsum;
	double vM = (rhoSqR * vyR + rhoSqL * vyL) * rhoSqQsum;
	double HM = (rhoSqR *  HR + rhoSqL *  HL) * rhoSqQsum;
	double cM = sqrt(gamma1 * (HM - 0.5 * (uM * uM + vM * vM)));

	/* calculation signal speeds */
	double arp = fmax(vxR + cR, uM + cM);
	double alm = fmin(vxL - cL, uM - cM);

	/* calculation HLL flux */
	if (alm > 0.0) {
		fluxLoc[0] = fL[0];
		fluxLoc[1] = fL[1];
		fluxLoc[2] = fL[2];
		fluxLoc[3] = fL[3];
	} else if (arp < 0.0) {
		fluxLoc[0] = fR[0];
		fluxLoc[1] = fR[1];
		fluxLoc[2] = fR[2];
		fluxLoc[3] = fR[3];
	} else {
		double as = (pR - pL + uL[MX] * (alm - vxL) - uR[MX] * (arp - vxR))
			/ (rhoL * (alm - vxL) - rhoR * (arp - vxR));
		if ((alm <= 0.0) && (as >= 0.0)) {
			double fac = rhoL * (alm - vxL) / (alm - as);
			double us[NVAR] = {fac,
					   as * fac,
					   vyL * fac,
					   fac * (eL / rhoL + (as - vxL) *
						(as + pL / (rhoL * (alm - vxL))))};
			fluxLoc[0] = fL[0] + alm * (us[0] - uL[0]);
			fluxLoc[1] = fL[1] + alm * (us[1] - uL[1]);
			fluxLoc[2] = fL[2] + alm * (us[2] - uL[2]);
			fluxLoc[3] = fL[3] + alm * (us[3] - uL[3]);
		} else {
			double fac = rhoR * (arp - vxR) / (arp - as);
			double us[NVAR] = {fac,
					   as * fac,
					   vyR * fac,
					   fac * (eR / rhoR + (as - vxR) *
						(as + pR / (rhoR * (arp - vxR))))};
			fluxLoc[0] = fR[0] + arp * (us[0] - uR[0]);
			fluxLoc[1] = fR[1] + arp * (us[1] - uR[1]);
			fluxLoc[2] = fR[2] + arp * (us[2] - uR[2]);
			fluxLoc[3] = fR[3] + arp * (us[3] - uR[3]);
		}
	}
}

/*
 * Lax-Friedrichs flux
 */
void flux_lxf(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{
	/* compute maximum Eigenvalue */
	double cL = sqrt(gamma * pL / rhoL);
	double cR = sqrt(gamma * pR / rhoR);
	double a  = fmax(fabs(vxR) + cR, fabs(vxL) + cL);

	/* calculate left/right energy and enthalpy */
	double eL = gamma1q * pL + 0.5 * rhoL * (vxL * vxL + vyL * vyL);
	double eR = gamma1q * pR + 0.5 * rhoR * (vxR * vxR + vyR * vyR);

	/* calculate the differences of the conservative variables */
	double delU[NVAR] = {rhoR       - rhoL,
			     rhoR * vxR - rhoL * vxL,
			     rhoR * vyR - rhoL * vyL,
			     eR         - eL};

	/* calculate the physical fluxes */
	double fL[4] = {rhoL * vxL,
			 rhoL * vxL * vxL + pL,
			 rhoL * vxL * vyL,
			 vxL * (eL + pL)};
	double fR[4] = {rhoR * vxR,
			 rhoR * vxR * vxR + pR,
			 rhoR * vxR * vyR,
			 vxR * (eR + pR)};

	/* calculate local Lax-Friedrichs flux */
	fluxLoc[0] = 0.5 * (fR[0] + fL[0]) - 0.5 * a * delU[0];
	fluxLoc[1] = 0.5 * (fR[1] + fL[1]) - 0.5 * a * delU[1];
	fluxLoc[2] = 0.5 * (fR[2] + fL[2]) - 0.5 * a * delU[2];
	fluxLoc[3] = 0.5 * (fR[3] + fL[3]) - 0.5 * a * delU[3];
}

/*
 * Steger-Warming flux
 */
void flux_stw(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{

}

/*
 * Central flux
 */
void flux_cen(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{

}

/*
 * AUSMD flux
 */
void flux_ausmd(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{

}

/*
 * AUSMDV flux
 */
void flux_ausmdv(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{

}

/*
 * Van Leer flux
 */
void flux_vanleer(double rhoL, double rhoR,
	      double vxL,  double vxR,
	      double vyL,  double vyR,
	      double pL,   double pR,
	      double fluxLoc[4])
{

}

/*
 * select the convective flux
 */
void convectiveFlux(double rhoL, double rhoR,
		    double vxL,  double vxR,
		    double vyL,  double vyR,
		    double pL,   double pR,
		    double fluxLoc[4])
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
		double fluxLoc[4];
		convectiveFlux(pVarL[RHO], pVarR[RHO],
			       pVarL[VX],  pVarR[VX],
			       pVarL[VY],  pVarR[VY],
			       pVarL[P],   pVarR[P],
			       fluxLoc);
	}
}
