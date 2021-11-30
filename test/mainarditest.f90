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
   use iso_fortran_env
   use wrightmod
   implicit none

   ! Parameters
   integer, parameter      :: N = 200
   real(real64), parameter :: D_PI = 3.1415926535897932384626433832795
   real(real64), parameter :: t = 1.0;
   ! Values
   real(real64) :: x(N), true(N)
   real(real64), allocatable :: w(:)
   real(real64) lambda, mu, h
   ! Auxiliary variables
   integer i

   h = 10.0/(real(N,real64)-1.0)
   do i=1,N
      x(i) = -5.0 + (i-1)*h
   end do

   write(output_unit,'("Test for closed forms of the Mainardi function.")')

   write(output_unit,'("1) λ = -0.5 μ = -0.5 W = exp(-z^2/4)sqrt(pi)")')
   lambda = -0.5
   mu = 0.5
   w = wright(x, t, lambda, mu)
   true = exp(-abs(x)**2/4.0)/sqrt(D_PI)
   open (unit = 99, file = "mainardi1.out")
   do i=1,N
      write(99,'(f20.16,",",f20.16,",",f20.16)')x(i),w(i),true(i)
   end do
   close(99)

   write(output_unit,'("2) λ = -0.5 μ = 0.0 W = abs(z) exp(-abs(z).^2/4)./(2 sqrt(pi))")')
   lambda = -0.5
   mu = 0.0
   w = wright(x, t, lambda, mu)
   true = abs(x)*exp(-abs(x)**(2.0/4.0))/(2.0*sqrt(D_PI));
   open (unit = 99, file = "mainardi2.out")
   do i=1,N
      write(99,'(f20.16,",",f20.16,",",f20.16)')x(i),w(i),true(i)
   end do
   close(99)

end program
