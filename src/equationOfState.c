/** \file
 *
 * \brief Contains conversion functions between the different variable types
 *
 * \author hhh
 * \date Sat 28 Mar 2020 09:45:30 PM CET
 */

#include <math.h>
#include <float.h>

#include "main.h"
#include "equation.h"

/**
 * \brief Convert primitive variables into conservative variables
 * \param[in] pVar Primitive variable vector
 * \param[out] cVar Conservative variable vector
 */
void primCons(const double pVar[NVAR], double cVar[NVAR])
{
	cVar[RHO] = pVar[RHO];
	cVar[MX]  = pVar[VX] * pVar[RHO];
	cVar[MY]  = pVar[VY] * pVar[RHO];
	cVar[E]   = gam1q * pVar[P]
		+ 0.5 * (cVar[MX] * pVar[VX] + cVar[MY] * pVar[VY]);
}

/** \brief Convert conservative variables into primitive variables
 *
 * This function is used during reconstruction, therefore it has to be checked
 * if the resulting primitive variables are negative. It that is the case, they
 * are set to zero.
 *
 * \param[in] cVar Conservative variable vector
 * \param[out] pVar Primitive variable vector
 */
void consPrim(const double cVar[NVAR], double pVar[NVAR])
{
	pVar[RHO] = cVar[RHO];
	pVar[VX]  = cVar[MX] / cVar[RHO];
	pVar[VY]  = cVar[MY] / cVar[RHO];
	pVar[P]   = gam1 *
		(cVar[E] - 0.5 * (cVar[MX] * pVar[VX] + cVar[MY] * pVar[VY]));

	/* check if density or pressure are negative */
	if (pVar[RHO] < DBL_EPSILON) {
		pVar[RHO] = 0.0;
	}
	if (pVar[P] < DBL_EPSILON) {
		pVar[P] = 0.0;
	}
}

/**
 * \brief Convert conservative variables to characteristic variables
 * \param[in] cVar Conservative variable vector
 * \param[out] charac Characteristic variable vector
 * \param[in] pVarRef Reference primitive variable vector
 */
void consChar(double cVar[NVAR], double charac[3], double pVarRef[NVAR])
{
	double c = sqrt(gam * pVarRef[P] / pVarRef[RHO]);
	double u = pVarRef[VX];
	double H = gam1q * c * c + 0.5 * u * u;
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

/**
 * \brief Convert characteristic variables to conservative variables
 * \param[in] charac Characteristic variable vector
 * \param[out] cVar Conservative variable vector
 * \param[in] pVarRef Reference primitive variable vector
 */
void charCons(double charac[3], double cVar[NVAR], double pVarRef[NVAR])
{
	double u = pVarRef[VX];
	double c = sqrt(gam * pVarRef[P] / pVarRef[RHO]);
	double H = gam1q * c * c + 0.5 * u * u;
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
