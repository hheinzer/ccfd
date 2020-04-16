/*
 * exactFunction.c
 *
 * Created: Sat 28 Mar 2020 05:31:15 PM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "mesh.h"
#include "equation.h"
#include "initialCondition.h"
#include "equationOfState.h"
#include "exactRiemann.h"

void exactFunc(int iExactFunc, double x[NDIM], double time, double pVar[NVAR])
{
	switch (iExactFunc) {
	case 1: {
		/* Richtmyer-Meshkov Instability (independent of time, only
		 * initial condition) */
		pVar[RHO] = 1.0;
		pVar[VX] = 0.0;
		pVar[VY] = 0.0;
		pVar[P] = 1.0;

		double xLen = xMax - xMin;
		double yLen = yMax - yMin;
		if (x[X] >= 0.3 * xLen + 1.0/30.0 * xLen * cos(2.0 * pi *
					3.0 / yLen * x[Y])) {
			pVar[RHO] = 0.25;
		}

		if ((x[X] <= 0.1 * xLen) && (x[X] >= 1.0/30.0 * xLen)) {
			pVar[RHO] = 4.22;
			pVar[P] = 4.9;
		}
		break;
	}
	case 2: {
		/* Gaussian pressure pulse (independent of time, only initial
		 * condition) */
		pVar[RHO] = 1.0;
		pVar[VX] = 0.0;
		pVar[VY] = 0.0;

		double hWidth = fmin((xMax - xMin), (yMax - yMin)) / 100.0 * 6.0;
		double amplitude = 1.0;

		double xPeak = xMin + 0.5 * (xMax - xMin);
		double yPeak = yMin + 0.5 * (yMax - yMin);
		double expLog = exp(- log(2.0));
		double hWidth2 = hWidth * hWidth;
		double dr2 = (x[X] - xPeak) * (x[X] - xPeak)
			   + (x[Y] - yPeak) * (x[Y] - yPeak);

		pVar[P] = 1.0 + amplitude * pow(expLog, dr2 / hWidth2);
		break;
	}
	case 3: {
		/* sinewave: test case for convergence tests */
		double freq = 1.0;
		double amplitude = 0.1;
		double omega = pi * freq;
		double a = 2.0 * pi;
		double resu[NVAR];
		resu[0] = resu[1] = resu[2] =
			2.0 + amplitude * sin(omega * (x[X] + x[Y]) - a * time);
		resu[3] = resu[0] * resu[0];
		consPrim(resu, pVar);
		break;
	}
	case 4: {
		/* double mach reflection: undisturbed shock for initial
		 * condition and upper boundary */
		double factor = pi / (0.2 * dxRef);
		double cSum[NVAR], cDiff[NVAR];
		cSum[0] = 9.4;
		cSum[1] = 57.157676649772950686;
		cSum[2] = - 33.0;
		cSum[3] = 566.0;
		cDiff[0] = 6.6;
		cDiff[1] = 57.157676649772950686;
		cDiff[2] = - 33.0;
		cDiff[3] = 561.0;

		double cVar[NVAR];
		for (int i = 0; i < NVAR; ++i) {
			cVar[i] = 0.5 * (cSum[i] - cDiff[i] * tanh((x[X] - (1.0 / 6.0 + (20.0 * time + x[Y]) * sqrt3q)) * factor));
		}

		pVar[0] = cVar[0];
		pVar[1] = cVar[1] / cVar[0];
		pVar[2] = cVar[2] / cVar[0];
		pVar[3] = gam1 * (cVar[3] - 0.5 * (cVar[1] * cVar[1] + cVar[2] * cVar[2]) / cVar[0]);
		break;
	}
	case 5: {
		/* 1D Riemann problem in x direction */
		double rho_l, u_l, p_l, c_l;
		rho_l = refState[0][RHO];
		u_l = refState[0][VX];
		p_l = refState[0][P];
		c_l = sqrt(gam * p_l / rho_l);

		double rho_r, u_r, p_r, c_r;
		rho_r = refState[1][RHO];
		u_r = refState[1][VX];
		p_r = refState[1][P];
		c_r = sqrt(gam * p_r / rho_r);

		if (time == 0.0) {
			if (x[X] <= rp1Dinterface) {
				pVar[RHO] = rho_l;
				pVar[VX] = u_l;
				pVar[P] = p_l;
			} else {
				pVar[RHO] = rho_r;
				pVar[VX] = u_r;
				pVar[P] = p_r;
			}
		} else {
			double s = (x[X] - rp1Dinterface) / time;
			exactRiemann(rho_l, rho_r, &pVar[RHO],
					u_l, u_r, &pVar[VX],
					p_l, p_r, &pVar[P],
					c_l, c_r, s);
		}
		pVar[VY] = 0.0;
		break;
	}
	case 6: {
		/* 1D sine wave */
		double freq = 1.0;
		double amplitude = 0.00001;
		double omega = pi * freq;
		double cons0[NVAR] = {1.0, 0.1, 0.0, 1.0};
		double prims0[NVAR];
		consPrim(cons0, prims0);
		double c0 = sqrt(gam * prims0[P] / prims0[RHO]);
		double H0 = (cons0[E] + prims0[P]) / prims0[RHO];
		double R[NVAR][NVAR] = {
			{1.0, 1.0, 0.0, 1.0},
			{prims0[VX] - c0, prims0[VX], 0.0, prims0[VX] + c0},
			{prims0[VY], prims0[VY], 1.0, prims0[VY]},
			{H0 - prims0[VX] * c0, (prims0[VX] * prims0[VX] + prims0[VY] * prims0[VY]) / 2.0, prims0[VY], H0 + prims0[VX] * c0}
		};

		double resu[NVAR];
		double vec[NVAR] = {0.0, 0.0, 0.0, amplitude * sin(omega * (x[X] - (c0 + prims0[VX]) * time))};
		for (int i = 0; i < NVAR; ++i) {
			resu[i] = cons0[i];
			for (int j = 0; j < NVAR; ++j) {
				resu[i] += R[i][j] * vec[j];
			}
		}
		consPrim(resu, pVar);
		break;
	}
	default:
		printf("| ERROR: Exact Function Unknown\n");
		exit(1);
	}
}
