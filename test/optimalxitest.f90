!! BSD 3-Clause License
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
program optimalxitest
  ! This test program uses Brent's method to compute the value of the optimal
  ! :math:`\xi`, and thus :math:`N`, for the computation of the Wright function.

  use iso_fortran_env, only: real32, real64, real128, output_unit
  use brentmod, only: fminbnd
  implicit none

  real(real64) :: c,N
  real(real64) :: mu = 3

  write(output_unit,'("Optimal c/N Test")')

  N = fminbnd(0.1_real64, 0.99_real64, 0.5_real64, 200.0_real64, &
    & epsilon(c), epsilon(c), 1.0D-3, fobjective, c)

  write(output_unit,'("c = ",F20.16," N = ",I3)')c,floor(N)

contains

  function fobjective(carg) result(f)
    use iso_fortran_env, only: real64
    implicit none

    real(real64) :: carg
    real(real64) :: f
    ! Constants and local variables
    real(real64) :: pi = 3.141592653589793
    real(real64) :: l,ltol

    l    = - log(epsilon(c))
    ltol = - log(1.0E-15)

    f = (sqrt(l*ltol)/pi)* &
          & sqrt(1.0_real64 + (1.0_real64/carg)*(1.0_real64 + ((2.0_real64 - mu)/ltol) &
          & *log(1.0_real64-carg)))

    write(output_unit,'("f = ",F20.16," c = ",F20.16)')f,carg

  end function fobjective

end program optimalxitest
