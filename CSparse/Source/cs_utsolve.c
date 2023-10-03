// CSparse/Source/cs_utsolve: x=U'\b where x and b are dense
// CSparse, Copyright (c) 2006-2022, Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: LGPL-2.1+
#include "cs.h"
/* solve U'x=b where x and b are dense.  x=b on input, solution on output. */
csi cs_utsolve_singular(const cs *U, double *x, double tolerance)
{
    csi p, j, n, *Up, *Ui ;
    double *Ux ;
    if (!CS_CSC (U) || !x) return (0) ;                     /* check inputs */
    n = U->n ; Up = U->p ; Ui = U->i ; Ux = U->x ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Up [j] ; p < Up [j+1]-1 ; p++)
        {
            x [j] -= Ux [p] * x [Ui [p]] ;
        }
        double pivot = Ux[Up[j + 1] - 1];
        x[j] /= fabs(pivot) > tolerance ? pivot : 1.0;
    }
    return (1) ;
}

csi cs_utsolve (const cs *U, double *x)
{
    return cs_utsolve_singular(U, x, 0.0);
}
