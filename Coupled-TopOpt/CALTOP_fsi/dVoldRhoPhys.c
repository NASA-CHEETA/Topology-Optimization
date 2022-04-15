#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "omp.h"

/**
 * @brief dVolume/drhoPhys = EVolume
 * 
 * @param eleVol 
 * @param gradVol 
 * @param ne 
 */
void dVoldRhoPhys(double *eleVol, double *gradVol, ITG ne)
{
    // copy eleVol to gradVol
    int num_cpus = 0;
    char *env;

    // get num of threads declaration, if any
    env = getenv("OMP_NUM_THREADS");
    if (num_cpus == 0)
    {
        if (env)
            num_cpus = atoi(env);
        if (num_cpus < 1)
        {
            num_cpus = 1;
        }
    }

#pragma omp parallel for num_threads(num_cpus)
    for (int i = 0; i < ne; i++)
    {
        gradVol[i] = eleVol[i];
    }
}