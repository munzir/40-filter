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


#ifndef FILTER_H
#define FILTER_H

/*-------------------------*/
/* Finite Impulse Response */
/*-------------------------*/

/// Finite Impulse Response Filter
typedef struct {
    double *X;
    double *b;
    size_t n;
    size_t i;
    size_t order;
} filter_fir_t;

/// Finite Impulse Response filter a scalar value
double filter_fir( filter_fir_t *f, double x );

/// Finite Impulse Response filter n values
void filter_fir_n( double *y, filter_fir_t *f, double *x, size_t n );


/*--------*/
/* Kalman */
/*--------*/

/// kalman filter.
/// matrices are COLUMN-MAJOR
typedef struct {
    double *x;  ///< Mean state
    double *u;  ///< Input
    double *z;  ///< Measurement
    double *A;  ///< Process Model
    double *B;  ///< Input model
    double *C;  ///< Measurement model
    double *E;  ///< Covariance
    double *R;  ///< Process Noise Covariance
    double *Q;  ///< Measurement Noise Covariance
    size_t n_x; ///< size of x
    size_t n_u; ///< size of u
    size_t n_z; ///< size of z
} filter_kalman_t;

/// Allocate all the arrays
void filter_kalman_init( filter_kalman_t *kf, size_t n_x, size_t n_u, size_t n_z );

/// Free all the arrays
void filter_kalman_destroy( filter_kalman_t *kf );

/// Predict step
void filter_kalman_predict( filter_kalman_t *kf );

/// Correct step
void filter_kalman_correct( filter_kalman_t *kf );


#endif
