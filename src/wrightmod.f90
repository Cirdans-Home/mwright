! BSD 3-Clause License
!
! Copyright (c) 2021, Fabio Durastante, Lidia Aceto
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
!    list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
module wrightmod
  ! Module containg the implementation of the routines for computing the Wright
  ! function on the real line.
  use iso_fortran_env, only : real32, real64, real128, error_unit
  use brentmod, only: fminbnd
  implicit none

  ! Constants
  real(real32),  parameter ::  S_PI = 3.14159265_real32
  real(real64),  parameter ::  D_PI = 3.141592653589793_real64
  real(real128), parameter ::  T_PI = 3.1415926535897932384626433832795_real128

  real(real32),  parameter :: sone = 1.0_real32
  real(real64),  parameter :: done = 1.0_real64
  real(real128), parameter :: tone = 1.0_real128

  complex(real32),  parameter :: sonei = (0.0,1.0_real32)
  complex(real64),  parameter :: donei = (0.0,1.0_real64)
  complex(real128), parameter :: tonei = (0.0,1.0_real128)

  ! Interface
  interface wright
    module procedure dwright
    module procedure twright
    module procedure swright
    ! Vectorized version
    module procedure dwright_vec
    module procedure swright_vec
    module procedure twright_vec
  end interface wright

  ! Only the external interface is available when importing the module,
  ! auxiliary constants and the implementation are private
  private :: dwright, swright, twright, dwright_vec, swright_vec, twright_vec
  private :: S_PI, D_PI, T_PI, sone, done, tone, sonei, donei, tonei

contains

  function dwright(x,t,lambda,mu,N) result(w)
    implicit none

    real(real64) :: x ! The `x` in the argument :math:`-|x|t^\lambda`
    real(real64) :: t ! The `t` in the argument :math:`-|x|t^\lambda` ( :math:`t > 0` )
    real(real64) :: lambda ! First parameter of the Wright function
    real(real64) :: mu ! Second parameter of the Wright function
    real(real64) :: w ! :math:`t^{\mu - 1} W_{\lambda,\mu} (-|x|t^{\lambda})`
    integer, optional :: N ! Optional, number of quadrature points

    ! Internal variables
    integer :: N_,k
    real(real64) :: h,gamma,xi,c,l,ltol
    complex(real64) :: uk,zk,zpk,sk

    l    = -log(epsilon(x))
    ltol = -log(1.0e-15_real64)
    if (present(N)) then
      N_ = N
      h = (4.0_real64*l)/(D_PI*real(N_*N_,real64))
      gamma = (D_PI*D_PI*real(N_*N_,real64))/(16.0_real64*t*l);
    else
      if (mu < 2.0_real64) then
        N_ = floor((sqrt(2.0_real64*l*ltol))/D_PI,kind=kind(N))
        h = (4.0_real64*l)/(D_PI*N_*N_)
        gamma = (D_PI*D_PI*N_*N_)/(16.0_real64*t*l)
      else if (mu == 2.0_real64) then
        N_ = fminbnd(1.0_real64 - 1.0_real64/(l-ltol) + 10.0*epsilon(l), &
          & 0.99_real64, 0.5_real64, 200.0_real64, &
          & epsilon(xi), epsilon(xi), 1.0D-4, fobjective2, c)
        xi = 2.0_real64/(1.0_real64 + log(-log((l-ltol)* &
          & (1.0_real64 - c)*(1.0_real64 - c)))/ltol)
        h = ((2.0_real64+xi*c)*l)/(D_PI*N_*N_)
        gamma = (D_PI*D_PI*N_*N_)/(((2.0_real64+xi*c)**2)*t*l)
      else
        N_ = fminbnd(0.1_real64, 0.99_real64, 0.5_real64, 200.0_real64, &
          & epsilon(xi), epsilon(xi), 1.0D-4, fobjective, c)
        xi = 2.0_real64/(done + ((2.0_real64-mu)/ltol)*log(done-c))
        h = ((2.0_real64+xi*c)*l)/(D_PI*N_*N_)
        gamma = (D_PI*D_PI*N_*N_)/(((2.0_real64+xi*c)**2)*t*l)
      end if
    end if

    sk = (0.0_real64,0.0_real64)
    quadrature: do k=-N_,N_
      uk = real(k,real64)*h
      zk = gamma*(donei*uk + done)**2
      zpk = (2.0_real64)*gamma*donei*(donei*uk + done)
      sk = sk + exp(zk*t)*(zk**(-mu))*exp(-abs(x)*zk**(-lambda))*zpk
    end do quadrature
    sk = h*sk/((2.0_real64)*D_PI*donei)

    w = RealPart(sk)

  contains

    function fobjective(xarg) result(f)
      use iso_fortran_env, only: real64
      implicit none

      real(real64) :: xarg
      real(real64) :: f

      f = (sqrt(l*ltol)/D_PI)* &
            & sqrt( done + (done/xarg)*(done + ((2.0_real64 - mu)/ltol) &
            & *log(done-xarg)))

    end function fobjective

    function fobjective2(xarg) result(f)
      use iso_fortran_env, only: real64
      implicit none

      real(real64) :: xarg
      real(real64) :: f

      f = (sqrt(l*ltol)/D_PI)* &
            & sqrt( done + (done/xarg)*(done + &
            & log(-log( (l-ltol)*(done-xarg)*(done-xarg)))/ltol))

    end function fobjective2

  end function dwright

  function swright(x,t,lambda,mu,N) result(w)
    implicit none

    real(real32) :: x ! The `x` in the argument :math:`-|x|t^\lambda`
    real(real32) :: t ! The `t` in the argument :math:`-|x|t^\lambda` ( :math:`t > 0` )
    real(real32) :: lambda ! First parameter of the Wright function
    real(real32) :: mu ! Second parameter of the Wright function
    real(real32) :: w ! :math:`t^{\mu - 1} W_{\lambda,\mu} (-|x|t^{\lambda})`
    integer, optional :: N ! Optional, number of quadrature points

    ! Internal variables
    integer :: N_,k
    real(real32) :: h,gamma,xi,c,l,ltol
    complex(real32) :: uk,zk,zpk,sk

    l    = -log(epsilon(x))
    ltol = -log(1.0e-6_real32)
    if (present(N)) then
      N_ = N
      h = (4.0_real32*l)/(S_PI*real(N_*N_,real32))
      gamma = (S_PI*S_PI*real(N_*N_,real32))/(16.0_real32*t*l);
    else
      if (mu <= 2.0_real32) then
        N_ = floor((sqrt(2.0_real32*l*ltol))/S_PI,kind=kind(N))
        h = (4.0_real32*l)/(S_PI*N_*N_)
        gamma = (S_PI*S_PI*N_*N_)/(16.0_real32*t*l)
      else if (mu == 2.0_real32) then
        N_ = fminbnd(1.0_real32 - 1.0_real32/(l-ltol) + 10.0*epsilon(l), &
          & 0.99_real32, 0.5_real32, 200.0_real32, &
          & epsilon(xi), epsilon(xi), real(1.0D-4,real32), fobjective2, c)
        xi = 2.0_real32/(1.0_real32 + log(-log((l-ltol)* &
          & (1.0_real32 - c)*(1.0_real32 - c)))/ltol)
        h = ((2.0_real32+xi*c)*l)/(D_PI*N_*N_)
        gamma = (S_PI*S_PI*N_*N_)/(((2.0_real32+xi*c)**2)*t*l)
      else
        N_ = fminbnd(0.1_real32, 0.99_real32, 0.5_real32, 200.0_real32, &
          & epsilon(xi), epsilon(xi), real(1.0D-4,real32), fobjective, c)
        xi = 2.0_real32/(sone + ((2.0_real32-mu)/ltol)*log(sone-c))
        h = ((2.0_real32+xi*c)*l)/(S_PI*N_*N_)
        gamma = (S_PI*S_PI*N_*N_)/(((2.0_real32+xi*c)**2)*t*l)
      end if

    end if

    sk = (0.0_real32,0.0_real32)
    quadrature: do k=-N_,N_
      uk = real(k,real32)*h
      zk = gamma*(sonei*uk + done)**2
      zpk = (2.0_real32)*gamma*sonei*(sonei*uk + sone)
      sk = sk + exp(zk*t)*(zk**(-mu))*exp(-abs(x)*zk**(-lambda))*zpk
    end do quadrature
    sk = h*sk/((2.0_real32)*S_PI*donei)

    w = RealPart(sk)

  contains

    function fobjective(xarg) result(f)
      use iso_fortran_env, only: real32
      implicit none

      real(real32) :: xarg
      real(real32) :: f

      f = (sqrt(l*ltol)/S_PI)* &
            & sqrt( sone + (sone/xarg)*(sone + ((2.0_real32 - mu)/ltol) &
            & *log(sone-xarg)))

    end function fobjective

    function fobjective2(xarg) result(f)
      use iso_fortran_env, only: real32
      implicit none

      real(real32) :: xarg
      real(real32) :: f

      f = (sqrt(l*ltol)/S_PI)* &
            & sqrt( sone + (sone/xarg)*(sone + &
            & log(-log( (l-ltol)*(sone-xarg)*(sone-xarg)))/ltol))

    end function fobjective2

  end function swright

  function twright(x,t,lambda,mu,N) result(w)
    implicit none

    real(real128) :: x ! The `x` in the argument :math:`-|x|t^\lambda`
    real(real128) :: t ! The `t` in the argument :math:`-|x|t^\lambda` ( :math:`t > 0` )
    real(real128) :: lambda ! First parameter of the Wright function
    real(real128) :: mu ! Second parameter of the Wright function
    real(real128) :: w ! :math:`t^{\mu - 1} W_{\lambda,\mu} (-|x|t^{\lambda})`
    integer, optional :: N ! Optional, number of quadrature points

    ! Internal variables
    integer :: N_,k
    real(real128) :: h,gamma,xi,c,l,ltol
    complex(real128) :: uk,zk,zpk,sk

    l    = -log(epsilon(x))
    ltol = -log(1.0e-30_real128)
    if (present(N)) then
      N_ = N
      h = (4.0_real128*l)/(S_PI*real(N_*N_,real128))
      gamma = (T_PI*T_PI*real(N_*N_,real128))/(16.0_real128*t*l);
    else
      if (mu < 2.0_real128) then
        N_ = floor((sqrt(2.0_real128*l*ltol))/T_PI,kind=kind(N))
        h = (4.0_real32*l)/(T_PI*N_*N_)
        gamma = (T_PI*T_PI*N_*N_)/(16.0_real128*t*l)
      else if (mu == 2.0_real128) then
        N_ = fminbnd(1.0_real128 - 1.0_real128/(l-ltol) + 10.0*epsilon(l), &
          & 0.99_real128, 0.5_real128, 200.0_real128, &
          & epsilon(xi), epsilon(xi), real(1.0D-4,real128), fobjective2, c)
        xi = 2.0_real128/(1.0_real128 + log(-log((l-ltol)* &
          & (1.0_real128 - c)*(1.0_real128 - c)))/ltol)
        h = ((2.0_real128+xi*c)*l)/(T_PI*N_*N_)
        gamma = (T_PI*T_PI*N_*N_)/(((2.0_real128+xi*c)**2)*t*l)
      else
        N_ = fminbnd(0.1_real128, 0.99_real128, 0.5_real128, 200.0_real128, &
          & epsilon(xi), epsilon(xi), real(1.0D-4,real128), fobjective, c)
        xi = 2.0_real128/(tone + ((2.0_real128-mu)/ltol)*log(tone-c))
        h = ((2.0_real128+xi*c)*l)/(T_PI*N_*N_)
        gamma = (T_PI*T_PI*N_*N_)/(((2.0_real128+xi*c)**2)*t*l)
      end if

    end if

    sk = (0.0_real128,0.0_real128)
    quadrature: do k=-N_,N_
      uk = real(k,real128)*h
      zk = gamma*(tonei*uk + done)**2
      zpk = (2.0_real128)*gamma*tonei*(tonei*uk + tone)
      sk = sk + exp(zk*t)*(zk**(-mu))*exp(-abs(x)*zk**(-lambda))*zpk
    end do quadrature
    sk = h*sk/((2.0_real128)*T_PI*tonei)

    w = RealPart(sk)

  contains

    function fobjective(xarg) result(f)
      use iso_fortran_env, only: real128
      implicit none

      real(real128) :: xarg
      real(real128) :: f

      f = (sqrt(l*ltol)/T_PI)* &
            & sqrt( tone + (tone/xarg)*(tone + ((2.0_real128 - mu)/ltol) &
            & *log(tone-xarg)))

    end function fobjective

    function fobjective2(xarg) result(f)
      use iso_fortran_env, only: real128
      implicit none

      real(real128) :: xarg
      real(real128) :: f

      f = (sqrt(l*ltol)/T_PI)* &
            & sqrt( tone + (tone/xarg)*(tone + &
            & log(-log( (l-ltol)*(tone-xarg)*(tone-xarg)))/ltol))

    end function fobjective2

  end function twright

  ! Vectorized versions
  function dwright_vec(x,t,lambda,mu,N) result(w)
    implicit none

    real(real64), dimension(:) :: x  ! The `x` in the argument :math:`-|x|t^\lambda`
    real(real64) :: t ! The `t` in the argument :math:`-|x|t^\lambda` ( :math:`t > 0` )
    real(real64) :: lambda ! First parameter of the Wright function
    real(real64) :: mu ! Second parameter of the Wright function
    real(real64), dimension(:), allocatable :: w ! :math:`t^{\mu - 1} W_{\lambda,\mu} (-|x|t^{\lambda})`
    integer, optional :: N ! Optional to arbitrarly set the number of quadrature points

    ! Internal variables
    integer :: N_,k,ll,m,info
    real(real64) :: h,gamma,xi,c,l,ltol
    complex(real64) :: sk
    complex(real64), allocatable, dimension(:) :: uk,zk,zpk

    m = size(x)

    if (.not.allocated(w)) &
      & allocate(w(m))
    if (size(w).ne.m) then
      write(error_unit,'("ERROR: Mismatch in I/O Array Lenght")')
    end if

    l    = -log(epsilon(x))
    ltol = -log(1.0e-15_real64)
    if (present(N)) then
      N_ = N
      h = (4.0_real64*l)/(D_PI*real(N_*N_,real64))
      gamma = (D_PI*D_PI*real(N_*N_,real64))/(16.0_real64*t*l);
    else
      if (mu < 2.0_real64) then
        N_ = floor((sqrt(2.0_real64*l*ltol))/D_PI,kind=kind(N))
        h = (4.0_real64*l)/(D_PI*N_*N_)
        gamma = (D_PI*D_PI*N_*N_)/(16.0_real64*t*l)
      else if (mu == 2.0_real64) then
        N_ = fminbnd(1.0_real64 - 1.0_real64/(l-ltol) + 10.0*epsilon(l), &
          & 0.99_real64, 0.5_real64, 200.0_real64, &
          & epsilon(xi), epsilon(xi), 1.0D-4, fobjective2, c)
        xi = 2.0_real64/(1.0_real64 + log(-log((l-ltol)* &
          & (1.0_real64 - c)*(1.0_real64 - c)))/ltol)
        h = ((2.0_real64+xi*c)*l)/(D_PI*N_*N_)
        gamma = (D_PI*D_PI*N_*N_)/(((2.0_real64+xi*c)**2)*t*l)
      else
        N_ = fminbnd(0.1_real64, 0.99_real64, 0.5_real64, 200.0_real64, &
          & epsilon(xi), epsilon(xi), 1.0D-4, fobjective, c)
        xi = 2.0_real64/(done + ((2.0_real64-mu)/ltol)*log(done-c))
        h = ((2.0_real64+xi*c)*l)/(D_PI*N_*N_)
        gamma = (D_PI*D_PI*N_*N_)/(((2.0_real64+xi*c)**2)*t*l)
      end if

    end if

    allocate(uk(2*N_+1),zk(2*N_+1),zpk(2*N_+1),stat=info)
    if( info.ne.0 ) &
      & write(error_unit,'("ERROR: In allocation temporary vectors")')

    do k=-N_,N_
      uk(k+N_+1) = real(k,real64)*h
      zk(k+N_+1) = gamma*(donei*uk(k+N_+1) + done)*(donei*uk(k+N_+1) + done)
      zpk(k+N_+1) = (2.0_real64)*gamma*donei*(donei*uk(k+N_+1) + done)
    end do

    do ll=1,m
        sk = sum(exp(zk*t - abs(x(ll))*(zk**(-lambda)))*(zk**(-mu))*zpk)
        w(ll) = RealPart( h*sk/((2.0_real64)*D_PI*donei))
    end do

    deallocate(uk,zk,zpk)

  contains

    function fobjective(xarg) result(f)
      use iso_fortran_env, only: real64
      implicit none

      real(real64) :: xarg
      real(real64) :: f

      f = (sqrt(l*ltol)/D_PI)* &
            & sqrt( done + (done/xarg)*(done + ((2.0_real64 - mu)/ltol) &
            & *log(done-xarg)))

    end function fobjective

    function fobjective2(xarg) result(f)
      use iso_fortran_env, only: real64
      implicit none

      real(real64) :: xarg
      real(real64) :: f

      f = (sqrt(l*ltol)/D_PI)* &
            & sqrt( done + (done/xarg)*(done + &
            & log(-log( (l-ltol)*(done-xarg)*(done-xarg)))/ltol))

    end function fobjective2

  end function dwright_vec

  function swright_vec(x,t,lambda,mu,N) result(w)
    implicit none

    real(real32) :: x(:) ! The `x` in the argument :math:`-|x|t^\lambda`
    real(real32) :: t ! The `t` in the argument :math:`-|x|t^\lambda` ( :math:`t > 0` )
    real(real32) :: lambda ! First parameter of the Wright function
    real(real32) :: mu ! Second parameter of the Wright function
    real(real32), dimension(:), allocatable :: w ! :math:`t^{\mu - 1} W_{\lambda,\mu} (-|x|t^{\lambda})`
    integer, optional :: N ! Optional to arbitrarly set the number of quadrature points

    ! Internal variables
    integer :: N_,k,ll,m,info
    real(real32) :: h,gamma,xi,c,l,ltol
    complex(real32) :: sk
    complex(real32), allocatable, dimension(:) :: uk,zk,zpk

    m = size(x)

    if (.not.allocated(w)) &
      & allocate(w(m))
    if (size(w).ne.m) then
      write(error_unit,'("ERROR: Mismatch in I/O Array Lenght")')
    end if

    l    = -log(epsilon(x))
    ltol = -log(1.0e-6_real32)
    if (present(N)) then
      N_ = N
      h = (4.0_real32*l)/(D_PI*real(N_*N_,real32))
      gamma = (D_PI*D_PI*real(N_*N_,real32))/(16.0_real32*t*l)
    else
      if (mu < 2.0_real32) then
        N_ = floor((sqrt(2.0_real32*l*ltol))/S_PI,kind=kind(N))
        h = (4.0_real32*l)/(S_PI*N_*N_)
        gamma = (S_PI*S_PI*N_*N_)/(16.0_real32*t*l)
      else if (mu == 2.0_real32) then
        N_ = fminbnd(1.0_real32 - 1.0_real32/(l-ltol) + 10.0*epsilon(l), &
          & 0.99_real32, 0.5_real32, 200.0_real32, &
          & epsilon(xi), epsilon(xi), real(1.0D-4,real32), fobjective2, c)
        xi = 2.0_real32/(1.0_real32 + log(-log((l-ltol)* &
          & (1.0_real32 - c)*(1.0_real32 - c)))/ltol)
        h = ((2.0_real32+xi*c)*l)/(D_PI*N_*N_)
        gamma = (S_PI*S_PI*N_*N_)/(((2.0_real32+xi*c)**2)*t*l)
      else
        N_ = fminbnd(0.1_real32, 0.99_real32, 0.5_real32, 200.0_real32, &
          & epsilon(xi), epsilon(xi), real(1.0D-4,real32), fobjective, c)
        xi = 2.0_real32/(sone + ((2.0_real32-mu)/ltol)*log(sone-c))
        h = ((2.0_real32+xi*c)*l)/(S_PI*N_*N_)
        gamma = (S_PI*S_PI*N_*N_)/(((2.0_real32+xi*c)**2)*t*l)
      end if

    end if

    allocate(uk(2*N_+1),zk(2*N_+1),zpk(2*N_+1),stat=info)
    if( info.ne.0 ) &
      & write(error_unit,'("ERROR: In allocation temporary vectors")')

    do k=-N_,N_
      uk(k+N_+1) = real(k,real32)*h
      zk(k+N_+1) = gamma*(donei*uk(k+N_+1) + done)*(donei*uk(k+N_+1) + done)
      zpk(k+N_+1) = (2.0_real32)*gamma*donei*(donei*uk(k+N_+1) + done)
    end do

    do ll=1,m
        sk = sum(exp(zk*t - abs(x(ll))*(zk**(-lambda)))*(zk**(-mu))*zpk)
        w(ll) = RealPart(h*sk/((2.0_real32)*S_PI*donei))
    end do

    deallocate(uk,zk,zpk)

  contains

    function fobjective(xarg) result(f)
      use iso_fortran_env, only: real32
      implicit none

      real(real32) :: xarg
      real(real32) :: f

      f = (sqrt(l*ltol)/S_PI)* &
            & sqrt( sone + (sone/xarg)*(sone + ((2.0_real32 - mu)/ltol) &
            & *log(sone-xarg)))

    end function fobjective

    function fobjective2(xarg) result(f)
      use iso_fortran_env, only: real32
      implicit none

      real(real32) :: xarg
      real(real32) :: f

      f = (sqrt(l*ltol)/S_PI)* &
            & sqrt( sone + (sone/xarg)*(sone + &
            & log(-log( (l-ltol)*(sone-xarg)*(sone-xarg)))/ltol))

    end function fobjective2

  end function swright_vec

  function twright_vec(x,t,lambda,mu,N) result(w)
    implicit none

    real(real128) :: x(:) ! The `x` in the argument :math:`-|x|t^\lambda`
    real(real128) :: t ! The `t` in the argument :math:`-|x|t^\lambda` ( :math:`t > 0` )
    real(real128) :: lambda ! First parameter of the Wright function
    real(real128) :: mu ! Second parameter of the Wright function
    real(real128), dimension(:), allocatable :: w ! :math:`t^{\mu - 1} W_{\lambda,\mu} (-|x|t^{\lambda})`
    integer, optional :: N ! Optional to arbitrarly set the number of quadrature points

    ! Internal variables
    integer :: N_,k,ll,m,info
    real(real128) :: h,gamma,xi,c,l,ltol
    complex(real128) :: sk
    complex(real128), allocatable, dimension(:) :: uk,zk,zpk

    m = size(x)

    if (.not.allocated(w)) &
      & allocate(w(m))
    if (size(w).ne.m) then
      write(error_unit,'("ERROR: Mismatch in I/O Array Lenght")')
    end if

    l    = -log(epsilon(x))
    ltol = -log(1.0e-30_real128)
    if (present(N)) then
      N_ = N
      h = (4.0_real128*l)/(D_PI*real(N_*N_,real128))
      gamma = (D_PI*D_PI*real(N_*N_,real128))/(16.0_real128*t*l);
    else
      if (mu < 2.0_real128) then
        N_ = floor((sqrt(2.0_real128*l*ltol))/T_PI,kind=kind(N))
        h = (4.0_real128*l)/(T_PI*N_*N_)
        gamma = (T_PI*T_PI*N_*N_)/(16.0_real128*t*l)
      else if (mu == 2.0_real128) then
        N_ = fminbnd(1.0_real128 - 1.0_real128/(l-ltol) + 10.0*epsilon(l), &
          & 0.99_real128, 0.5_real128, 200.0_real128, &
          & epsilon(xi), epsilon(xi), real(1.0D-4,real128), fobjective2, c)
        xi = 2.0_real128/(1.0_real128 + log(-log((l-ltol)* &
          & (1.0_real128 - c)*(1.0_real128 - c)))/ltol)
        h = ((2.0_real128+xi*c)*l)/(T_PI*N_*N_)
        gamma = (T_PI*T_PI*N_*N_)/(((2.0_real128+xi*c)**2)*t*l)
      else
        N_ = fminbnd(0.1_real128, 0.99_real128, 0.5_real128, 200.0_real128, &
          & epsilon(xi), epsilon(xi), real(1.0D-4,real128), fobjective, c)
        xi = 2.0_real128/(tone + ((2.0_real128-mu)/ltol)*log(tone-c))
        h = ((2.0_real128+xi*c)*l)/(T_PI*N_*N_)
        gamma = (T_PI*T_PI*N_*N_)/(((2.0_real128+xi*c)**2)*t*l)
      end if

    end if

    allocate(uk(2*N_+1),zk(2*N_+1),zpk(2*N_+1),stat=info)
    if( info.ne.0 ) &
      & write(error_unit,'("ERROR: In allocation temporary vectors")')

    do k=-N_,N_
      uk(k+N_+1) = real(k,real128)*h
      zk(k+N_+1) = gamma*(donei*uk(k+N_+1) + done)*(donei*uk(k+N_+1) + done)
      zpk(k+N_+1) = (2.0_real128)*gamma*donei*(donei*uk(k+N_+1) + done)
    end do

    do ll=1,m
        sk = sum(exp(zk*t - abs(x(ll))*(zk**(-lambda)))*(zk**(-mu))*zpk)
        w(ll) = RealPart(h*sk/((2.0_real128)*T_PI*donei))
    end do

    deallocate(uk,zk,zpk)

  contains

    function fobjective(xarg) result(f)
      use iso_fortran_env, only: real128
      implicit none

      real(real128) :: xarg
      real(real128) :: f

      f = (sqrt(l*ltol)/T_PI)* &
            & sqrt( tone + (tone/xarg)*(tone + ((2.0_real128 - mu)/ltol) &
            & *log(tone-xarg)))

    end function fobjective

    function fobjective2(xarg) result(f)
      use iso_fortran_env, only: real128
      implicit none

      real(real128) :: xarg
      real(real128) :: f

      f = (sqrt(l*ltol)/T_PI)* &
            & sqrt( tone + (tone/xarg)*(tone + &
            & log(-log( (l-ltol)*(tone-xarg)*(tone-xarg)))/ltol))

    end function fobjective2

  end function twright_vec

end module wrightmod
