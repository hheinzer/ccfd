/** \file
 *
 * \author hhh
 * \date Sat 21 Mar 2020 10:51:13 AM CET
 */

#ifndef READINTOOLS_H
#define READINTOOLS_H

#include <stdbool.h>

#include "main.h"

void fillCmds(char iniFileName[STRLEN]);
char *getStr(const char *key, const char *proposal);
int countKeys(const char *key, const int proposal);
int getInt(const char *key, const char *proposal);
double getDbl(const char *key, const char *proposal);
bool getBool(const char *key, const char *proposal);
int *getIntArray(const char *key, const int N, const char *proposal);
double *getDblArray(const char *key, const int N, const char *proposal);
void ignoredCmds(void);

#endif
