/*
 * readInTools.h
 *
 * Created: Sat 21 Mar 2020 10:51:13 AM CET
 * Author : hhh
 */

#ifndef READINTOOLS_H
#define READINTOOLS_H

#include "main.h"

void fillCmds(char iniFileName[STRLEN]);
char *getStr(char *key, char *proposal);
int countKeys(const char *key, const int proposal);
int getInt(const char *key, const char *proposal);
double getDbl(const char *key, const char *proposal);
bool getBool(const char *key, const char *proposal);
int *getIntArray(const char *key, const int N, const char *proposal);
double *getDblArray(const char *key, const int N, const char *proposal);
void ignoredCmds(void);

#endif

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
