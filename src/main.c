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
	/* read in ini file */
	char *iniFileName = "../calc/sod.ini";
	fillCmds(iniFileName);

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
}
