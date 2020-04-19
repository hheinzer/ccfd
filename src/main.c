/** \file
 *
 * \brief Contains the main function of `ccfd`
 *
 * \author hhh
 * \date Sat 21 Mar 2020 10:41:51 AM CET
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "readInTools.h"
#include "timeDiscretization.h"
#include "output.h"
#include "equation.h"
#include "boundary.h"
#include "mesh.h"
#include "initialCondition.h"
#include "finiteVolume.h"
#include "linearSolver.h"
#include "analyze.h"

/** \brief Main function
 *
 * This function starts the program by initializing all the necessary global
 * variables and initializing the time discretization loop. The function
 * finishes by deallocating all the allocated memory.
 *
 * \param[in] argc The number of command line arguments passed to `main`
 * \param[in] argv The argument vector containing the command line arguments
 * \return 0 = Success, 1 = Error during execution
 */
int main(int argc, char *argv[])
{
	printf("=============================================================\n");
	printf("                            C C F D                          \n");
	printf("=============================================================\n");
	printf("   Solution of the two-dimensional Euler and Navier-Stokes   \n");
	printf("     equations using an unstructured finite volume solver    \n");
	printf("                                                             \n");
	printf("   This program may be used under the conditions of the GPL  \n");
	printf("                                                             \n");
	printf("  Created at the Institute of Aerodynamics and Gas Dynamics  \n");
	printf("           at the University of Stuttgart, Germany           \n");
	printf("               http://www.iag.uni-stuttgart.de               \n");
	printf("                                                             \n");
	printf("           Recreated in C by Heinz Heinrich Heinzer          \n");
	printf("=============================================================\n");

	/* read parameter file and make the command list */
	fillCmds(argv[1]);
	isStationary = getBool("stationary", "T");

	/* check command line arguments */
	switch (argc) {
	case 2:
		/* no restart required */
		isRestart = false;
		iniIterationNumber = 0;
		startTime = 0.0;
		break;
	case 3:
		/* Restart calculation */
		isRestart = true;
		const char *restartFileName = argv[2];
		strcpy(strIniCondFile, restartFileName);

		printf("\nRestarting from file '%s':\n", restartFileName);
		FILE *restartFile = fopen(restartFileName, "r");
		if (restartFile == NULL) {
			printf("| ERROR: No restart file found\n");
			exit(1);
		}
		fclose(restartFile);

		/* get restart time */
		char *string = strstr(restartFileName, "_") + 1;
		while (strstr(string, "_")) {
			string = strstr(string, "_") + 1;
		}
		double timeStamp = strtod(string, NULL);
		printf("%g\n", timeStamp);

		if (isStationary) {
			iniIterationNumber = (int)timeStamp;
			startTime = 0.0;
		} else {
			iniIterationNumber = 0;
			startTime = timeStamp;
			restartTime = timeStamp;
		}
		break;
	default:
		printf("ERROR: Wrong number of arguments, must be 1 or 2\n");
		exit(1);
	}

	/* initialization routines */
	initOutput();
	initEquation();
	initBoundary();
	initMesh();
	initInitialCondition();
	initFV();
	initTimeDisc();
	initLinearSolver();
	outputTimes = NULL;

	/* setting initial condition */
	setInitialCondition();

	/* initialize c_a, c_w, and c_p calculation as well as record points */
	initAnalyze();

	/* print ignored commands */
	ignoredCmds();

	/* start time stepping routine */
	timeDisc();

	/* clean that memory, like you should */
	freeMesh();
	freeBoundary();
	freeOutputTimes();
	freeInitialCondition();
	freeAnalyze();
	freeLinearSolver();
}
