/*
 * main.c
 *
 * Created: Sat 21 Mar 2020 10:41:51 AM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "readInTools.h"

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
	printf("=============================================================\n\n");

	/* read parameter file and make the command list */
	fillCmds(argv[1]);

	/* check command line arguments */
	switch (argc) {
	case 2:
		//isRestart = false;
		//iniIterationNumber = 0;
		//startTime = 0.0;
		break;
	default:
		printf("ERROR: Wrong number of arguments, must be 1 or 2\n");
		exit(1);
	}
}
