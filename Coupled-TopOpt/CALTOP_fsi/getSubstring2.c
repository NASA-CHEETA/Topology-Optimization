#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "CalculiX.h"

/**
 * @brief Get the ith set name from SET list
 *
 * @param string
 * @param position
 * @param length
 * @return char*
 */
char *getSubstring2(char *string, int position, int length)
{
    char *p;
    int c;

    p = malloc(length + 1);
    if (p == NULL)
    {
        printf("\nERROR: Unable to allocate memory in getSubstring.\n");
        exit(1);
    }

    for (c = 0; c < length; c++)
    {
        *(p + c) = *(string + position - 1);
        string++;
    }
    *(p + c) = '\0';
    return p;
}