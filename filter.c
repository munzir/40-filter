/* -*- mode: C; c-basic-offset: 4  -*- */
/*
 * Copyright (c) 2010, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     * Neither the name of the Georgia Tech Research Corporation nor
 *       the names of its contributors may be used to endorse or
 *       promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GEORGIA TECH RESEARCH CORPORATION ''AS
 * IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GEORGIA
 * TECH RESEARCH CORPORATION BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#include <stdlib.h>
#include <assert.h>
#include <somatic/util.h>
#include "filter.h"

/*-------------------------*/
/* Finite Impulse Response */
/*-------------------------*/
double filter_fir( filter_fir_t *f, double x ) {
    f->X[ f->i ] = x;
    double y = 0;
    for( size_t i = 0; i < f->n; i ++ ) {
        size_t j = (i+f->i) % f->n;
        assert( j < f->n );
        y += f->b[i] * f->X[j];
    }
    f->i = (f->i+1) % f->n;
    return y;
}

void filter_fir_n( double *y, filter_fir_t *f, double *x, size_t n ) {
    for( size_t i = 0; i < n; i ++ ) {
        y[i] = filter_fir( &f[i], x[i] );
    }
}



/*--------*/
/* Kalman */
/*--------*/

// fortran prototype
void filter_kalman_predict_( double *x, int *n_x, 
                             double *E, double *A, 
                             double *B, double *u, int *n_u, 
                             double *R );
// fortran prototype
void filter_kalman_correct_( double *x, int *n_x,
                             double *E, double *z, int *n_z,
                             double *C, double *Q );

void filter_kalman_init( filter_kalman_t *kf, size_t n_x, size_t n_u, size_t n_z ) {
    kf->n_x = n_x;
    kf->n_u = n_u;
    kf->n_z = n_z;

    kf->x = SOMATIC_NEW_AR( double, n_x );
    kf->u = SOMATIC_NEW_AR( double, n_u );
    kf->z = SOMATIC_NEW_AR( double, n_z );

    kf->A = SOMATIC_NEW_AR( double, n_x*n_x );
    kf->B = SOMATIC_NEW_AR( double, n_x*n_u );
    kf->C = SOMATIC_NEW_AR( double, n_z*n_x );

    kf->E = SOMATIC_NEW_AR( double, n_x*n_x );
    kf->R = SOMATIC_NEW_AR( double, n_x*n_x );
    kf->Q = SOMATIC_NEW_AR( double, n_z*n_z );
}

void filter_kalman_destroy( filter_kalman_t *kf ) {
    free( kf->x );
    free( kf->u );
    free( kf->z );

    free( kf->A );
    free( kf->B );
    free( kf->C );

    free( kf->E );
    free( kf->R );
    free( kf->Q );
}

void filter_kalman_predict( filter_kalman_t *kf ) {
    int n_x = (int)kf->n_x;
    int n_u = (int)kf->n_u;
    filter_kalman_predict_( kf->x, &n_x, kf->E, kf->A,
                            kf->B, kf->u, &n_u, kf->R );
}
void filter_kalman_correct( filter_kalman_t *kf ) {
    int n_x = (int)kf->n_x;
    int n_z = (int)kf->n_z;
    filter_kalman_correct_( kf->x, &n_x, kf->E, 
                            kf->z, &n_z, kf->C, kf->Q );
}
