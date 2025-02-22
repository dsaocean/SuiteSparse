//------------------------------------------------------------------------------
// GB_assign_zombie1: delete all entries in C(:,j) for GB_assign
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2023, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

// JIT: not needed.  Only one variant possible.

// C(:,j)<!> = anything: GrB_Row_assign or GrB_Col_assign with an empty
// complemented mask requires all entries in the C(:,j) vector to be deleted.
// C must be sparse or hypersparse.

// C->iso is not affected.

#include "GB_assign.h"
#include "GB_assign_zombie.h"

GrB_Info GB_assign_zombie1
(
    GrB_Matrix C,
    const int64_t j
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    ASSERT (!GB_IS_FULL (C)) ;
    ASSERT (!GB_IS_BITMAP (C)) ;
    ASSERT (GB_ZOMBIES_OK (C)) ;
    ASSERT (GB_JUMBLED_OK (C)) ;
    ASSERT (!GB_PENDING (C)) ;

    //--------------------------------------------------------------------------
    // get C(:,j)
    //--------------------------------------------------------------------------

    int64_t *restrict Ci = C->i ;
    const int64_t *restrict Cp = C->p ;
    int64_t pC_start, pC_end ;

    if (C->h != NULL)
    { 
        // C is hypersparse
        GB_hyper_hash_lookup (Cp, C->Y->p, C->Y->i, C->Y->x, C->Y->vdim-1,
            j, &pC_start, &pC_end) ;
    }
    else
    { 
        // C is sparse
        pC_start = Cp [j] ;
        pC_end   = Cp [j+1] ;
    }

    int64_t cjnz = pC_end - pC_start ;
    int64_t nzombies = C->nzombies ;

    //--------------------------------------------------------------------------
    // determine the number of threads to use
    //--------------------------------------------------------------------------

    int nthreads_max = GB_Context_nthreads_max ( ) ;
    double chunk = GB_Context_chunk ( ) ;
    int nthreads = GB_nthreads (cjnz, chunk, nthreads_max) ;

    //--------------------------------------------------------------------------
    // C(:,j) = empty
    //--------------------------------------------------------------------------

    int64_t pC ;
    #pragma omp parallel for num_threads(nthreads) schedule(static) \
        reduction(+:nzombies)
    for (pC = pC_start ; pC < pC_end ; pC++)
    {
        int64_t i = Ci [pC] ;
        if (!GB_IS_ZOMBIE (i))
        { 
            // delete C(i,j) by marking it as a zombie
            nzombies++ ;
            Ci [pC] = GB_FLIP (i) ;
        }
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    C->nzombies = nzombies ;
    return (GrB_SUCCESS) ;
}

