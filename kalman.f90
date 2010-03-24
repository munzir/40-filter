!! Copyright (c) 2010, Georgia Tech Research Corporation
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!     * Redistributions of source code must retain the above
!!       copyright notice, this list of conditions and the following
!!       disclaimer.
!!     * Redistributions in binary form must reproduce the above
!!       copyright notice, this list of conditions and the following
!!       disclaimer in the documentation and/or other materials
!!       provided with the distribution.
!!     * Neither the name of the Georgia Tech Research Corporation nor
!!       the names of its contributors may be used to endorse or
!!       promote products derived from this software without specific
!!       prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY GEORGIA TECH RESEARCH CORPORATION ''AS
!! IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GEORGIA
!! TECH RESEARCH CORPORATION BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
!! OF THE POSSIBILITY OF SUCH DAMAGE.

!> \file kalman.f90
!! \brief kalman filter
!! \author Neil T. Dantam
  
Module kalman_filter
  Implicit None
Contains
  Function kf_matrix_invert( m, n, A ) result(info)
    Use f77_lapack
    Use f95_lapack
    integer, intent(in) :: m,n
    real(8), intent(inout), dimension(m,n) :: A
    integer :: info

    integer, dimension(m) :: ipiv  ! lu factor pivots
    real(8), dimension(1) :: swork !
        
    !! LU Factor
    Call la_getrf(m, n, A, m, ipiv, info)
    !! invert
    ! first get optimal work array size
    Call la_getri( n, A, m, ipiv, swork, -1, info ) 
    ! now perform the inversion
    info = really_invert( int( swork(1) ) )
  Contains
    Function really_invert( lwork ) result(info)
      integer, intent(in) :: lwork
      integer :: info
      real(8), dimension(lwork) :: work ! work array
      Call la_getri( n, A, m, ipiv, work, lwork, info )
    End Function really_invert
  End Function kf_matrix_invert
  
  Subroutine kf_predict( x, E, A, B, u, R )
    real(8), dimension(:), intent(inout) :: x    ! mean
    real(8), dimension(:,:), intent(inout) :: E  ! covariance
    real(8), dimension(:,:), intent(in) :: A     
    real(8), dimension(:,:), intent(in) :: B
    real(8), dimension(:,:), intent(in) :: R     ! process noise covariance
    real(8), dimension(:), intent(in) :: u       ! input
    ! x = A*x + B*u
    x = matmul(A,x) + matmul(B,u)
    ! E = A * E * A**T + R
    E = matmul( matmul(A,E), transpose(A) ) + R
  End Subroutine kf_predict

  Function kf_correct( x, E, z, C, Q ) result(info)
    real(8), dimension(:), intent(inout) :: x    ! mean
    real(8), dimension(:,:), intent(inout) :: E  ! covariance
    real(8), dimension(:), intent(in) :: z       ! measurement
    real(8), dimension(:,:), intent(in) :: C     ! measurement model
    real(8), dimension(:,:), intent(in) :: Q     ! measurement noise covariance
    integer :: info

    real(8), dimension( size(x), size(z) ) :: K
    real(8), dimension( size(x), size(x) ) :: Kp
    real(8), dimension( size(x), size(x) ) :: Ident
    integer :: i

    ! K = E * C**T * (C * E * C**T + Q)**-1
    Kp = matmul( matmul(C, E), transpose(C) ) + Q
    info = kf_matrix_invert(size(x), size(x), Kp)
    K = matmul( matmul(E, transpose(C)), Kp )

    ! x = x + K * (z - C*x)
    x = x + matmul( K, z - matmul(C,x) )

    ! E = (I - K*C) * E 
    Ident = 0
    Forall ( i = 1:size(x) )
       Ident(i,i) = 1
    End Forall
    E = matmul( Ident - matmul(K,C), E  )
    
  End Function kf_correct
End Module kalman_filter

Subroutine filter_kalman_predict( x, n_x, E, A, B, u, n_u, R )
  Use kalman_filter
  Implicit None
  real(8), intent(inout), dimension(n_x) :: x
  real(8), intent(inout), dimension(n_x,n_x) :: E
  real(8), intent(in), dimension(n_x,n_x) :: A, R
  real(8), intent(in), dimension(n_x,n_u) :: B
  real(8), intent(in), dimension(n_u) :: u
  integer,  intent(in) :: n_x, n_u 
  Call kf_predict( x, E, A, B, u, R )
End Subroutine filter_kalman_predict

Function filter_kalman_correct( x, n_x, E, z, n_z, C, Q ) result(info)
  Use kalman_filter
  Implicit None
  integer, intent(in) :: n_x, n_z
  real(8), intent(inout), dimension(n_x) :: x
  real(8), intent(in), dimension(n_z) :: z
  real(8), intent(inout), dimension(n_x, n_x) :: E
  real(8), intent(in), dimension(n_z, n_x) :: C
  real(8), intent(in), dimension(n_z, n_z) :: Q
  integer :: info
  info = kf_correct( x, E, z, C, Q )
End Function filter_kalman_correct
