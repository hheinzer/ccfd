/*
 * equationOfState.c
 *
 * Created: Sat 28 Mar 2020 09:45:30 PM CET
 * Author : hhh
 */

#include <math.h>
#include <float.h>

#include "main.h"
#include "equation.h"

/*
 * convert primitive variables into conservative
 */
void primCons(const double pVar[NVAR], double cVar[NVAR])
{
	cVar[RHO] = pVar[RHO];
	cVar[MX]  = pVar[VX] * pVar[RHO];
	cVar[MY]  = pVar[VY] * pVar[RHO];
	cVar[E]   = gamma1q * pVar[P]
		+ 0.5 * (cVar[MX] * pVar[VX] + cVar[MY] * pVar[VY]);
}

/*
 * convert conservative variables into primitive
 */
void consPrim(const double cVar[NVAR], double pVar[NVAR])
{
	pVar[RHO] = cVar[RHO];
	pVar[VX]  = cVar[MX] / cVar[RHO];
	pVar[VY]  = cVar[MY] / cVar[RHO];
	pVar[P]   = gamma1 *
		(cVar[E] - 0.5 * (cVar[MX] * pVar[VX] + cVar[MY] * pVar[VY]));

	/* check if density or pressure are negative */
	if (pVar[RHO] < DBL_EPSILON) {
		pVar[RHO] = 0.0;
	}
	if (pVar[P] < DBL_EPSILON) {
		pVar[P] = 0.0;
	}
}

/*
 * convert conservative variables to characteristic
 */
void consChar(double charac[3], double cVar[NVAR], double pVarRef[NVAR])
{
	double c = sqrt(gamma * pVarRef[P] / pVarRef[RHO]);
	double u = pVarRef[VX];
	double H = gamma1q * c * c + 0.5 * u * u;
	double phi = u * u - 2.0 * H;
	double a1 = 1.0 / (2.0 * c * phi);
	double a2 = 1.0 / phi;
	double a3 = u * c;

	double cVar1D[3] = {cVar[RHO], cVar[VX], cVar[E]};
	double K[3][3] = {
		{a1 * u * (phi - a3), -a1 * (phi - 2.0 * a3), -a2},
		{a2 * (u * u + phi), -2.0 * u * a2, 2.0 * a2},
		{-a1 * u * (a3 + phi), a1 * (phi + 2.0 * a3), -a2}
	};
	for (int i = 0; i < 3; ++i) {
		charac[i] = 0.0;
		for (int j = 0; j < 3; ++j) {
			charac[i] += K[i][j] * cVar1D[j];
		}
	}
}

/*
 * convert characteristic variables to conservative
 */
void charCons(double charac[3], double cVar[NVAR], double pVarRef[NVAR])
{
	double u = pVarRef[VX];
	double c = sqrt(gamma * pVarRef[P] / pVarRef[RHO]);
	double H = gamma1q * c * c + 0.5 * u * u;
	double K[3][3] = {
		{1.0, 1.0, 1.0},
		{u - c, u, u + c},
		{H - u * c, 0.5 * u * u, H + u * c}
	};

	double cVar1D[3] = {0.0};
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			cVar1D[i] += K[i][j] * charac[j];
		}
	}
	cVar[RHO] = cVar1D[0];
	cVar[MX]  = cVar1D[1];
	cVar[E]   = cVar1D[2];
}
