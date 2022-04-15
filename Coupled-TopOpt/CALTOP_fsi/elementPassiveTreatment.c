#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "omp.h"

/**
 * @brief Find a set named passive, and assign the values for those element set in a vector to setToValue
 * 
 * @param nset 
 * @param ialset 
 * @param set 
 * @param istartset 
 * @param iendset 
 * @param vector Input vector
 * @param setToValue vector[i] = setToValue
 */
void elementPassiveTreatment(ITG nset, ITG *ialset, char *set, ITG *istartset,
                             ITG *iendset, double *vector, double setToValue)
{
    char passiveSet[] = "PASSIVE";
    char *dummyNameset;
    ITG istart;
    ITG iend;
    ITG elementNumber;

    // For OPENMP
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

    for (int i = 0; i < nset; i++)
    {
        dummyNameset = getSubstring2(set, i * 81 + 1, strlen(passiveSet));
        if (strcmp(dummyNameset, passiveSet) == 0)
        {
            printf("Note: Passive element set found. Setting values to %f.\n", setToValue);
            istart = istartset[i] - 1;
            iend = iendset[i] - 1;

#pragma omp parallel for num_threads(num_cpus) private(elementNumber)
            for (ITG iset = istart; iset <= iend; iset++)
            {
                elementNumber = ialset[iset];
                vector[elementNumber - 1] = setToValue;
            }
        }
    }
    SFREE(dummyNameset);
}