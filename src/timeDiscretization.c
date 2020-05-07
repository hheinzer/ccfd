/** \file
 *
 * \brief Contains the functions for performing the time stepping process
 *
 * \author hhh
 * \date Sat 21 Mar 2020 07:52:42 PM CET
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <float.h>

#include "main.h"
#include "timeDiscretization.h"
#include "readInTools.h"
#include "output.h"
#include "mesh.h"
#include "equation.h"
#include "analyze.h"
#include "linearSolver.h"
#include "equationOfState.h"
#include "finiteVolume.h"
#include "memTools.h"

/* extern variables */
double	cfl;				/**< Courant-Friedrichs-Lewy number */
double	dfl;				/**< diffusive Courant-Friedrichs-Lewy number */
double	t;				/**< global calculation time */

double	timeOverall;			/**< overall time */

int	timeOrder;			/**< order of the time integration */
bool	isTimeStep1D;			/**< flag for 1D problem */

bool	isStationary;			/**< flag for stationary problem */
long	maxIter;			/**< maximum number of iterations */
double	stopTime;			/**< simulation end time */
long	iniIterationNumber;		/**< initial iteration number */
double	startTime;			/**< starting time */
double	abortResidual;			/**< residual at which to abort the calculation */
int	abortVariable;			/**< abort variable */
char	abortVariableName[4];		/**< string of the abort variable */
double	clAbortResidual;		/**< CL abort residual */
double  cdAbortResidual;		/**< CD abort residual */
bool	doAbortOnClResidual;		/**< CL abort flag */
bool	doAbortOnCdResidual;		/**< CD abort flag */
bool	isRestart;			/**< restart flag */
double	restartTime;			/**< calculation time for restart */

int	printIter;			/**< iterations after which to output */
double	printTime;			/**< calculation time after which to output */

int	nRKstages;			/**< number of Runge-Kutta stages */
double	RKcoeff[6] = {0.0};		/**< array of Runge-Kutta coefficients */
bool	isImplicit;			/**< implicit calculation flag */

/* local variables */
double **deltaX;			/**< variable used in implicit calculation */
double **Q;				/**< variable used in implicit calculation */
double **F_X0;				/**< variable used in implicit calculation */
double **F_XK;				/**< variable used in implicit calculation */

/**
 * \brief Initialize the time discretization
 */
void initTimeDisc(void)
{
	printf("\nTime Discretization:\n");
	cfl = getDbl("CFL", "0.9");
	dfl = getDbl("DFL", "0.9");
	isTimeStep1D = getBool("timeStep1D", "F");
	timeOrder = getInt("timeOrder", "1");
	if (timeOrder > 3) {
		printf("| ERROR: Temporal Discretization Order must be either 1, 2 or 3\n");
		exit(1);
	}

	isImplicit = getBool("implicit", "F");
	if (isImplicit) {
		deltaX = dyn2DdblArray(NVAR, nElems);
		Q = dyn2DdblArray(NVAR, nElems);
		F_X0 = dyn2DdblArray(NVAR, nElems);
		F_XK = dyn2DdblArray(NVAR, nElems);
	}

	nRKstages = getInt("nRKstages", "1");
	if (nRKstages > 5) {
		printf("| ERROR: No more than 5 RK Stages are Possible\n");
		exit(1);
	}
	if ((nRKstages == 1) && (timeOrder == 1)) {
		printf("| Time Stepping Scheme: Euler Time Stepping\n");
		RKcoeff[1] = 1.0;
	} else {
		printf("| Time Stepping Scheme: RK Time Stepping\n");
		switch (nRKstages) {
		case 3:
			switch (timeOrder) {
			case 1:
				RKcoeff[1] = 0.1418;
				RKcoeff[2] = 0.4;
				RKcoeff[3] = 1.0;
				break;
			case 2:
				RKcoeff[1] = 0.1918;
				RKcoeff[2] = 0.4929;
				RKcoeff[3] = 1.0;
				break;
			case 3:
				RKcoeff[1] = 0.2884;
				RKcoeff[2] = 0.5010;
				RKcoeff[3] = 1.0;
				break;
			}
			break;
		case 4:
			switch (timeOrder) {
			case 1:
				RKcoeff[1] = 0.0833;
				RKcoeff[2] = 0.2069;
				RKcoeff[3] = 0.4265;
				RKcoeff[4] = 1.0;
				break;
			case 2:
				RKcoeff[1] = 0.1084;
				RKcoeff[2] = 0.2602;
				RKcoeff[3] = 0.5052;
				RKcoeff[4] = 1.0;
				break;
			case 3:
				RKcoeff[1] = 0.1666;
				RKcoeff[2] = 0.3027;
				RKcoeff[3] = 0.5275;
				RKcoeff[4] = 1.0;
				break;
			}
			break;
		case 5:
			switch (timeOrder) {
			case 1:
				RKcoeff[1] = 0.0533;
				RKcoeff[2] = 0.1263;
				RKcoeff[3] = 0.2375;
				RKcoeff[4] = 0.4414;
				RKcoeff[5] = 1.0;
				break;
			case 2:
				RKcoeff[1] = 0.0695;
				RKcoeff[2] = 0.1602;
				RKcoeff[3] = 0.2898;
				RKcoeff[4] = 0.5060;
				RKcoeff[5] = 1.0;
				break;
			case 3:
				RKcoeff[1] = 0.1067;
				RKcoeff[2] = 0.1979;
				RKcoeff[3] = 0.3232;
				RKcoeff[4] = 0.5201;
				RKcoeff[5] = 1.0;
				break;
			}
			break;
		default:
			printf("| ERROR: Wrong Number of RK Coefficients\n");
			exit(1);
		}

		char line[STRLEN] = "";
		sprintf(line, "{%g", RKcoeff[1]);
		for (int iRK = 2; iRK <= nRKstages; ++iRK) {
			sprintf(line + strlen(line), ", %g", RKcoeff[iRK]);
		}
		sprintf(line + strlen(line), "}");
		printf("| %16s[%d] = %27s\n", "RKcoeff", nRKstages, line);
	}

	/* stationary computation */
	if (isStationary) {
		doAbortOnClResidual = doAbortOnCdResidual = false;

		abortResidual = getDbl("abortResidual", "1e-6");

		abortVariable = getInt("abortVariable", "1") - 1;
		switch (abortVariable) {
		case RHO:
			strcpy(abortVariableName, "RHO");
			break;
		case MX:
			strcpy(abortVariableName, "MX");
			break;
		case MY:
			strcpy(abortVariableName, "MY");
			break;
		case E:
			strcpy(abortVariableName, "E");
			break;
		}

		clAbortResidual = getDbl("cl_abortResidual", "0.0");
		cdAbortResidual = getDbl("cd_abortResidual", "0.0");

		printf("| Stationary Problem\n");
		if ((clAbortResidual == 0.0) && (cdAbortResidual == 0)) {
			printf("| Abort Residual '%s': %g\n",
					abortVariableName, abortResidual);
		} else if (clAbortResidual > 0.0) {
			printf("| Cl Abort Residual: %g\n", clAbortResidual);
			doAbortOnClResidual = true;
		} else if (cdAbortResidual > 0.0) {
			printf("| Cd Abort Residual: %g\n", cdAbortResidual);
			doAbortOnCdResidual = true;
		} else {
			printf("| ERROR: Wrong Definition of Abort Residual\n");
			exit(1);
		}
	} else {
		printf("| Transient Problem\n");
	}

	maxIter = getInt("maxIter", "100000");
	stopTime = getDbl("tEnd", NULL);

	if (isRestart) {
		t = restartTime;
	} else {
		t = 0.0;
	}

	if (isStationary) {
		printf("| Initial Iteration Number: %ld\n", iniIterationNumber);
	} else {
		printf("| Start Time: %g\n", t);
	}

	printIter = (iniIterationNumber / IOiterInterval + 1) * IOiterInterval;
	printTime = (floor(t / IOtimeInterval) + 1) * IOtimeInterval;
}

/**
 * \brief Compute the time step
 * \param[in] pTime The print time interval
 * \param[out] dt The resulting time step
 * \param[out] viscousTimeStepDominates Flag for if the viscous time step is
 *	dominating
 */
void calcTimeStep(double pTime, double *dt, bool *viscousTimeStepDominates)
{
	/* calculate local timestep for each cell */
	*viscousTimeStepDominates = false;
	if (isTimeStep1D) {
		double dtMax = 1e150;
		#pragma omp parallel for reduction(min:dtMax)
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];
			double a = sqrt(gam * aElem->pVar[P] / aElem->pVar[RHO]);
			double dtConv = cfl * aElem->sy / (fabs(aElem->pVar[VX]) + a);
			if (!isfinite(dtConv)) {
				printf("| Convective Time Step NaN\n");
				exit(1);
			}
			if (mu > 1e-10) {
				printf("| TimeStep1D not implemented for Navier Stokes. Set mu = 0 or turn off timeStep1D\n");
				exit(1);
			}
			dtMax = fmin(dtMax, dtConv);
		}

		*dt = dtMax;
	} else {
		double gamPrMax = fmax(4.0 / 3.0, gam / Pr);
		double dtConvMax = 1e150;
		#pragma omp parallel for reduction(min:dtConvMax)
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];
			/* convective time step */
			double a = sqrt(gam * aElem->pVar[P] / aElem->pVar[RHO]);
			double sumSpectralRadii = (fabs(aElem->pVar[VX]) + a) * aElem->sx
						+ (fabs(aElem->pVar[VY]) + a) * aElem->sy;
			double dtConv = cfl * aElem->area / sumSpectralRadii;
			if (!isfinite(dtConv)) {
				printf("| Convective Time Step NaN\n");
				exit(1);
			}
			dtConvMax = fmin(dtConvMax, dtConv);
		}

		/* viscous time step */
		double dtViscMax = 1e150;
		if (mu > 1e-10) {
			#pragma omp parallel for reduction(min:dtViscMax)
			for (long iElem = 0; iElem < nElems; ++iElem) {
				elem_t *aElem = elem[iElem];
				double sumSpectralRadii
					= gamPrMax * mu * aElem->sx * aElem->sx
					+ gamPrMax * mu * aElem->sy * aElem->sy;
				double dtVisc = dfl * aElem->pVar[RHO] * aElem->pVar[RHO]
					* aElem->area * aElem->area
					/ (4.0 * sumSpectralRadii);
				if (!isfinite(dtVisc)) {
					printf("| Viscous Time Step NaN\n");
					exit(1);
				}
				dtViscMax = fmin(dtViscMax, dtVisc);
			}
		}

		*dt = fmin(dtConvMax, dtViscMax);
		if (dtViscMax < dtConvMax) {
			*viscousTimeStepDominates = true;
		}
	}

	/* special treatment for data output and stoptime */
	if ((t + *dt > pTime) || (t + *dt > stopTime)) {
		*dt = fmin(pTime, stopTime) - t;
	} else if ((t + 1.5 * *dt > pTime) || (t + 1.5 * *dt > stopTime)) {
		*dt = 0.5 * (fmin(pTime, stopTime) - t);
	}

	/* set local time step for each cell to the global time step */
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		aElem->dt = *dt;
	}
}

/**
 * \brief Performs explicit time step using Euler scheme
 * \param[in] time Computation time at calculation
 * \param[in] dt Time step at calculation
 * \param[out] resIter Residual vector for time step
 */
void explicitTimeStepEuler(double time, double dt, double resIter[NVAR + 2])
{
	fvTimeDerivative(time);

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		aElem->cVar[RHO] += dt * aElem->u_t[RHO];
		aElem->cVar[MX]  += dt * aElem->u_t[MX];
		aElem->cVar[MY]  += dt * aElem->u_t[MY];
		aElem->cVar[E]   += dt * aElem->u_t[E];

		consPrim(aElem->cVar, aElem->pVar);
	}

	globalResidual(resIter);
}

/**
 * \brief Performs explicit time step using Runge-Kutta scheme `nRKstages` stages
 * \param[in] time Computation time at calculation
 * \param[in] dt Time step at calculation
 * \param[out] resIter Residual vector for time step
 */
void explicitTimeStepRK(double time, double dt, double resIter[NVAR + 2])
{
	/* save the initial solution as needed for the RK scheme */
	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];
		aElem->cVarStage[RHO] = aElem->cVar[RHO];
		aElem->cVarStage[MX]  = aElem->cVar[MX];
		aElem->cVarStage[MY]  = aElem->cVar[MY];
		aElem->cVarStage[E]   = aElem->cVar[E];
	}

	/* loop over the RK stages */
	for (int iStage = 1; iStage <= nRKstages; ++iStage) {
		double dtStage = RKcoeff[iStage - 1] * dt;
		fvTimeDerivative(time + dtStage);

		/* time update of conservative variables */
		#pragma omp parallel for
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];

			aElem->cVar[RHO] = aElem->cVarStage[RHO]
				+ RKcoeff[iStage] * dt * aElem->u_t[RHO];

			aElem->cVar[MX]  = aElem->cVarStage[MX]
				+ RKcoeff[iStage] * dt * aElem->u_t[MX];

			aElem->cVar[MY]  = aElem->cVarStage[MY]
				+ RKcoeff[iStage] * dt * aElem->u_t[MY];

			aElem->cVar[E]   = aElem->cVarStage[E]
				+ RKcoeff[iStage] * dt * aElem->u_t[E];

			consPrim(aElem->cVar, aElem->pVar);
		}
	}

	globalResidual(resIter);
}

/** \brief Euler implicit time integration
 *
 * The non-linear equations require the use of a Newton method with internal
 * sub-iteration, using a GMRES method.
 *
 * \param[in] time Computation time at calculation
 * \param[in] dt Time step at calculation
 * \param[out] resIter Residual vector for time step
 */
void implicitTimeStep(double time, double dt, double resIter[NVAR + 2])
{
	/* input parameters for Newton */
	double alpha = 1.0;
	double beta = 1.0;

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		Q[RHO][iElem] = aElem->cVar[RHO];
		Q[MX][iElem]  = aElem->cVar[MX];
		Q[MY][iElem]  = aElem->cVar[MY];
		Q[E][iElem]   = aElem->cVar[E];

		consPrim(aElem->cVar, aElem->pVar);
	}

	/* Newton */
	time = t + beta * dt;

	fvTimeDerivative(time);

	#pragma omp parallel for
	for (long iElem = 0; iElem < nElems; ++iElem) {
		elem_t *aElem = elem[iElem];

		F_X0[RHO][iElem] = aElem->cVar[RHO] - Q[RHO][iElem]
			- alpha * dt * aElem->u_t[RHO];

		F_X0[MX][iElem]  = aElem->cVar[MX]  - Q[MX][iElem]
			- alpha * dt * aElem->u_t[MX];

		F_X0[MY][iElem]  = aElem->cVar[MY]  - Q[MY][iElem]
			- alpha * dt * aElem->u_t[MY];

		F_X0[E][iElem]   = aElem->cVar[E]   - Q[E][iElem]
			- alpha * dt * aElem->u_t[E];

		XK[RHO][iElem] = aElem->cVar[RHO];
		XK[MX][iElem]  = aElem->cVar[MX];
		XK[MY][iElem]  = aElem->cVar[MY];
		XK[E][iElem]   = aElem->cVar[E];

		R_XK[RHO][iElem] = aElem->u_t[RHO];
		R_XK[MX][iElem]  = aElem->u_t[MX];
		R_XK[MY][iElem]  = aElem->u_t[MY];
		R_XK[E][iElem]   = aElem->u_t[E];

		F_XK[RHO][iElem] = F_X0[RHO][iElem];
		F_XK[MX][iElem]  = F_X0[MX][iElem];
		F_XK[MY][iElem]  = F_X0[MY][iElem];
		F_XK[E][iElem]   = F_X0[E][iElem];
	}

	/* vector dot product */
	double norm2_F_X0 = vectorDotProduct(F_X0, F_X0);

	/* preparation for matrix vector multiplication */
	double eps2newtonLoc =
		fmax(1e-8 * 1e-8 * nElems * eps2newton / norm2_F_X0,
		     eps2newton);
	double norm2_F_XK = norm2_F_X0;

	nInnerNewton = 0;

	/* Newton iterations */
	double abortCritGMRES, norm2_F_XK_old = norm2_F_XK, gamA, gamB;
	while ((norm2_F_XK > eps2newtonLoc * norm2_F_X0) && (nInnerNewton < nNewtonIter)) {
		if (nInnerNewton == 0) {
			abortCritGMRES = 0.999;
			norm2_F_XK_old = norm2_F_XK;
		} else {
			gamA = gamEW * norm2_F_XK / norm2_F_XK_old;

			if (gamEW * abortCritGMRES * abortCritGMRES < 0.1) {
				gamB = fmin(0.999, gamA);
			} else {
				gamB = fmin(0.999, fmax(gamA,
					gamEW * abortCritGMRES * abortCritGMRES));
			}

			abortCritGMRES = fmin(0.999, fmax(gamB,
					0.5 * eps2newton_sq / sqrt(norm2_F_XK)));

			norm2_F_XK_old = norm2_F_XK;
		}

		nInnerNewton++;

		GMRES_M(time, dt, alpha, F_XK, sqrt(norm2_F_XK),
				&abortCritGMRES, deltaX);

		#pragma omp parallel for
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];

			XK[RHO][iElem] += deltaX[RHO][iElem];
			XK[MX][iElem]  += deltaX[MX][iElem];
			XK[MY][iElem]  += deltaX[MY][iElem];
			XK[E][iElem]   += deltaX[E][iElem];

			aElem->cVar[RHO] = XK[RHO][iElem];
			aElem->cVar[MX]  = XK[MX][iElem];
			aElem->cVar[MY]  = XK[MY][iElem];
			aElem->cVar[E]   = XK[E][iElem];

			consPrim(aElem->cVar, aElem->pVar);
		}

		fvTimeDerivative(time);

		#pragma omp parallel for
		for (long iElem = 0; iElem < nElems; ++iElem) {
			elem_t *aElem = elem[iElem];

			R_XK[RHO][iElem] = aElem->u_t[RHO];
			R_XK[MX][iElem]  = aElem->u_t[MX];
			R_XK[MY][iElem]  = aElem->u_t[MY];
			R_XK[E][iElem]   = aElem->u_t[E];

			F_XK[RHO][iElem] = aElem->cVar[RHO] - Q[RHO][iElem] - alpha * dt * aElem->u_t[RHO];
			F_XK[MX][iElem]  = aElem->cVar[MX]  - Q[MX][iElem]  - alpha * dt * aElem->u_t[MX];
			F_XK[MY][iElem]  = aElem->cVar[MY]  - Q[MY][iElem]  - alpha * dt * aElem->u_t[MY];
			F_XK[E][iElem]   = aElem->cVar[E]   - Q[E][iElem]   - alpha * dt * aElem->u_t[E];
		}

		norm2_F_XK = vectorDotProduct(F_XK, F_XK);
	}

	nNewtonIterGlobal += nInnerNewton;

	if (nInnerNewton == nNewtonIter) {
		printf("| %g\n", eps2newton);
		printf("| Newton NOT converged with %d Newton iteraions\n", nInnerNewton);
		printf("| Norm / Norm_R0 = %g\n", norm2_F_XK / norm2_F_X0);
		exit(1);
	}

	globalResidual(resIter);
}

/** \brief Main time discretization loop
 *
 * Selection of temporal integration method, as well as management of data
 * output and analysis tools.
 */
void timeDisc(void)
{
	bool hasConverged = (isStationary ? false : true);

	/* write initial condition to disk */
	printf("\nWriting Initial Condition to Disk:\n");
	dataOutput(t, iniIterationNumber);
	printf("| done.\n");

	/* main program loop */
	printf("\nStarting Computation:\n");

	#ifdef _OPENMP
		if (omp_get_max_threads() > 1) {
			printf("| OpenMP Enabled: Running on %d threads\n",
					omp_get_max_threads());
		} else {
			printf("| OpenMP Enabled: Running on 1 thread\n");
		}
	#else
		printf("| OpenMP Disabled: Running on 1 thread\n");
	#endif

	long start = iniIterationNumber + 1;
	double tStart = CPU_TIME();
	double tIOstart = tStart;
	bool viscousTimeStepDominates;
	double dt;
	calcTimeStep(printTime, &dt, &viscousTimeStepDominates);
	printf("| Initial Time Step: %.10g\n", dt);
	if (viscousTimeStepDominates) {
		printf("| Viscous Time Step Dominates!\n");
	}

	/* loop over all iterations */
	long iter;
	for (iter = start; iter <= maxIter; ++iter) {
		calcTimeStep(printTime, &dt, &viscousTimeStepDominates);

		/* main computation loop */
		double resIter[NVAR + 2] = {0.0};
		if (!isImplicit) {
			if ((timeOrder == 1) && (nRKstages == 1)) {
				explicitTimeStepEuler(t, dt, resIter);
			} else {
				explicitTimeStepRK(t, dt, resIter);
			}
		} else {
			implicitTimeStep(t, dt, resIter);
		}
		t += dt;

		/* analyze results */
		analyze(t, iter, resIter);

		/* end time abort criterion */
		if (stopTime - t <= 1e-15) {
			printf("\nTime Limit Reached - Computation Complete!\n");
			printf("| Final Time      : %.10g\n", t);
			printf("| Iteration Number: %ld\n", iter);
			printf("| Writing Final State to Disk\n");
			dataOutput(t, iter);
			finalizeDataOutput();

			/* error calculation (for exact function) */
			if (hasExactSolution) {
				calcErrors(t);
			}

			break;
		}

		/* residual abort criterion */
		if (isStationary) {
			if (fabs(resIter[abortVariable]) <= abortResidual) {
				if (doAbortOnClResidual) {
					if (fabs(resIter[4]) <= clAbortResidual) {
						printf("\nConverged in CL - Calculation complete\n");
						printf("| Iteration number: %ld\n", iter);
						printf("| Writing final state to disk\n");

						if (doCalcWing) {
							printf("| Final Lift and Drag Coefficients:\n");
							printf("|   CL: %.10g\n", wing.cl);
							printf("|   CD: %.10g\n", wing.cd);
						}

						hasConverged = true;
					}
				} else if (doAbortOnCdResidual) {
					if (fabs(resIter[5]) <= cdAbortResidual) {
						printf("\nConverged in CD - Calculation complete\n");
						printf("| Iteration number: %ld\n", iter);
						printf("| Writing final state to disk\n");

						if (doCalcWing) {
							printf("| Final Lift and Drag Coefficients:\n");
							printf("|   CL: %.10g\n", wing.cl);
							printf("|   CD: %.10g\n", wing.cd);
						}

						hasConverged = true;
					}
				} else {
					printf("\nConverged in '%s' - Calculation complete\n", abortVariableName);
					printf("| Iteration number: %ld\n", iter);
					printf("| Writing final state to disk\n");

					if (doCalcWing) {
						printf("| Final Lift and Drag Coefficients:\n");
						printf("|   CL: %.10g\n", wing.cl);
						printf("|   CD: %.10g\n", wing.cd);
					}

					hasConverged = true;
				}
			}

			if (hasConverged) {
				dataOutput(t, iter);
				finalizeDataOutput();

				if (hasExactSolution) {
					calcErrors(t);
				}

				break;
			}
		}

		/* data output */
		if ((printTime - t <= 1e-15) || (iter == printIter)) {
			printf("\nData Output at Iteration %ld\n", iter);
			double tIOend = CPU_TIME();
			printf("| Time since last I/O: %.10g s\n", tIOend - tIOstart);
			tIOstart = tIOend;

			if (isStationary) {
				if (doCalcWing) {
					printf("| Residuals: %s: %20.14e\n",
						abortVariableName, resIter[abortVariable]);
					printf("|            CL : %20.14e\n", resIter[4]);
					printf("|            CD : %20.14e\n", resIter[5]);
				} else {
					printf("|            RHO: %20.14e\n", resIter[RHO]);
					printf("|            MX : %20.14e\n", resIter[MX]);
					printf("|            MY : %20.14e\n", resIter[MY]);
					printf("|            E  : %20.14e\n", resIter[E]);
				}
			} else {
				printf("| Time     : %.10g\n", t);
				calcTimeStep(printTime + 1e150, &dt, &viscousTimeStepDominates);
				printf("| Time Step: %.10g\n", dt);
				if (viscousTimeStepDominates) {
					printf("| Viscous Time Step Dominates!\n");
				}
			}

			if (hasExactSolution) {
				calcErrors(t);
			}

			/* data output */
			dataOutput(t, iter);
			finalizeDataOutput();

			if (printTime - t <= 1e-15) {
				printTime += IOtimeInterval;
			}
			if (iter == printIter) {
				printIter += IOiterInterval;
			}
		}
	}

	double tEnd = CPU_TIME();

	/* error handler for the case that the maximum iteration number was reached */
	if (iter > maxIter) {
		printf("\nERROR - ERROR - ERROR\n");
		printf("| Maximum Iteration Number Reached. Calculation Aborted!\n");
		printf("| Final Time %g has not been reached.\n", stopTime);
		printf("| Current Time: %.10g\n", t);
		printf("| Final State will be written to disk\n");

		dataOutput(t, iter - 1);
		finalizeDataOutput();

		if (hasExactSolution) {
			calcErrors(t);
		}
	}

	/* standard output */
	printf("\nComputation Time: %.10g s\n", tEnd - tStart);
	if (isImplicit) {
		printf("| Newton Iterations: %d\n", nNewtonIterGlobal);
		printf("| GMRES Iterations : %d\n", nGMRESiterGlobal);
	}

	/* close all open files */
	if (isStationary) {
		fclose(resFile);
	}

	if (recordPoint.nPoints > 0) {
		for (int iPt = 0; iPt < recordPoint.nPoints; ++iPt) {
			fclose(recordPoint.ioFile[iPt]);
		}
	}

	/* free memory that is allocated for implicit calculation */
	if (isImplicit) {
		free(deltaX);
		free(Q);
		free(F_X0);
		free(F_XK);
	}
}
