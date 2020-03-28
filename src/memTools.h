/*
 * memTools.h
 *
 * Created: Fri 27 Mar 2020 02:38:45 PM CET
 * Author : hhh
 */

#ifndef MEMTOOLS_H
#define MEMTOOLS_H

long **dyn2DintArray(long I, long J);
double **dyn2DdblArray(long I, long J);
long ***dyn3DintArray(long I, long J, long K);
double ***dyn3DdblArray(long I, long J, long K);
double ****dyn4DdblArray(long I, long J, long K, long L);

#endif
