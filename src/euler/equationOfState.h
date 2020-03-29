/*
 * equationOfState.h
 *
 * Created: Sat 28 Mar 2020 09:45:50 PM CET
 * Author : hhh
 */

#ifndef EQUATIONOFSTATE_H
#define EQUATIONOFSTATE_H

#include "main.h"

void primCons(double pVar[NVAR], double cVar[NVAR]);
void consPrim(double cVar[NVAR], double pVar[NVAR]);
void consChar(double charac[3], double cVar[NVAR], double pVarRef[NVAR]);
void charCons(double charac[3], double cVar[NVAR], double pVarRef[NVAR]);

#endif
