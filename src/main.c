/*
 * main.c
 *
 * Created: Sat 21 Mar 2020 10:41:51 AM CET
 * Author : hhh
 */

#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "readInTools.h"

int main(int argc, char *argv[])
{
	/* greeting */
	printf(" =========================================================== \n");
	printf("                            C C F D                          \n");
	printf(" =========================================================== \n");
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
	printf(" =========================================================== \n\n");


	/* read in parameter file */
	char *iniFileName = "../calc/sod.ini";
	fillCmds(iniFileName);
}

/*
 * test readInTools
 *
char *aStr = getStr("filename", NULL);
int aCount = countKeys("meshbctype", 3);
printf("%i\n", aCount);
int anInt = getInt("meshbctype", "440");
double aDbl = getDbl("tEnd", "0.334");
bool aBool = getBool("exactsolution", "F");
int *anIntArray = getIntArray("nBCsegments", 4, "3,4,2,1");
double *anDblArray = getDblArray("xmax", 2, "3.2,1.9");
ignoredCmds();
*/
