! BSD 3-Clause License
!
!  Copyright (c) 2021, Fabio Durastante, Lidia Aceto
!  All rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:
!
!  1. Redistributions of source code must retain the above copyright notice, this
!     list of conditions and the following disclaimer.
!
!  2. Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!  3. Neither the name of the copyright holder nor the names of its
!     contributors may be used to endorse or promote products derived from
!     this software without specific prior written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
program mainarditest
   ! This test program produces information relative to the error analysis in
   ! double and quadruple precision.
   use iso_fortran_env
   use wrightmod
   implicit none

   ! Parameters
   integer, parameter      :: N = 25
   real(real64)  :: D_PI = 3.141592653589793_real64
   real(real128) :: T_PI = 3.1415926535897932384626433832795_real128

   write(output_unit,'("Test for closed forms of the Mainardi function.")')
   ! Test in double precision
   call mainarditest1()

   ! Test in quadruple precision for λ = -0.5 μ = 0.5
   call mainarditest2(17)
   call mainarditest2(23)
   call mainarditest2(30)
   call mainarditest2(37)

   ! Test in quadruple precision for λ = 0 μ = 0
   call mainarditest3(17)
   call mainarditest3(23)
   call mainarditest3(30)
   call mainarditest3(37)

contains

   subroutine mainarditest1()
      ! This subroutine implements the numerical test of accuracy for the
      ! routine in double precision with respect to a known function.

      ! Values
      real(real64) :: x(N), true(N)
      real(real64), allocatable :: w(:)
      real(real64) lambda, mu, h
      real(real64) :: t = 1.0_real64
      ! Auxiliary variables
      integer i

      h = 10.0_real64/(real(N,real64)-1.0)
      do i=1,N
         x(i) = -5.0_real64 + real((i-1),real64)*h
      end do

      write(output_unit,'("(real64) With λ = -0.5 μ = 0.5 W = exp(-z^2/4)sqrt(pi)")')
      lambda = -0.5_real64
      mu = 0.5_real64
      w = wright(x, t, lambda, mu)
      true = exp(-(abs(x)**2)/4.0_real64)/sqrt(D_PI)
      open (unit = 99, file = "mainardi1.out")
      do i=1,N
         write(99,'(f20.16,",",f20.16,",",f20.16)')x(i),w(i),true(i)
      end do
      close(99)

      deallocate(w)

   end subroutine mainarditest1

   subroutine mainarditest2(qnodes)
      ! This subroutine implements the numerical test of accuracy for the
      ! routine in triple precision with respect to a known function.
      integer, optional :: qnodes

      ! Values
      real(real128) :: x(N), true(N)
      real(real128), allocatable :: w(:)
      real(real128) lambda, mu, h
      real(real128) :: t = 1.0_real128
      ! Auxiliary variables
      integer i
      character(len=10) :: file_id
      character(len=50) :: file_name

      h = 10.0_real128/(real(N,real128)-1.0)
      do i=1,N
         x(i) = -5.0_real128 + real((i-1),real64)*h
      end do

      write(output_unit,'("(real128) With λ = -0.5 μ = 0.5 W = exp(-z^2/4)sqrt(pi)")')
      lambda = -0.5_real128
      mu = 0.5_real128
      w = wright(x, t, lambda, mu, qnodes)
      true = exp(-(abs(x)**2)/4.0_real128)/sqrt(T_PI)

      open (unit = 99, file = "mainardi2.out")
      do i=1,N
         write(99,'(f20.16,",",f40.36,",",f40.36)')x(i),w(i),true(i)
      end do
      close(99)

      if ( present(qnodes) ) then
         write(file_id,'(i2)')qnodes
      else
         write(file_id,'("37")')
      end if
      file_name = "mainardi-m12-quadruple-" // trim(file_id) // ".out"


      open (unit = 99, file = trim(file_name))
      do i=1,N
         write(99,'(f20.16,",",f40.36)')x(i),abs(w(i)-true(i))/true(i)
      end do
      close(99)

      deallocate(w)

   end subroutine mainarditest2

   subroutine mainarditest3(qnodes)
      ! This subroutine implements the numerical test of accuracy for the
      ! routine in triple precision with respect to a known function.
      integer, optional :: qnodes

      ! Values
      real(real128) :: x(N), true(N)
      real(real128), allocatable :: w(:)
      real(real128) lambda, mu, h
      real(real128) :: t = 1.0_real128
      ! Auxiliary variables
      integer i
      character(len=10) :: file_id
      character(len=50) :: file_name

      h = 10.0_real128/(real(N,real128)-1.0)
      do i=1,N
         x(i) = -5.0_real128 + real((i-1),real64)*h
      end do

      write(output_unit,'("(real128) With λ = 0 μ = 1.0 W = exp(-|z|)")')
      lambda = 0.0_real128
      mu = 1.0_real128
      w = wright(x, t, lambda, mu, qnodes)
      true = exp(-abs(x))

      open (unit = 99, file = "mainardi3.out")
      do i=1,N
         write(99,'(f20.16,",",f40.36,",",f40.36)')x(i),w(i),true(i)
      end do
      close(99)

      if ( present(qnodes) ) then
         write(file_id,'(i2)')qnodes
      else
         write(file_id,'("37")')
      end if
      file_name = "mainardi-m0-quadruple-" // trim(file_id) // ".out"


      open (unit = 99, file = trim(file_name))
      do i=1,N
         write(99,'(f20.16,",",f40.36)')x(i),abs(w(i)-true(i))/true(i)
      end do
      close(99)

      deallocate(w)

   end subroutine mainarditest3

end program
