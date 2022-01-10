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
program quadtest
  ! This test program checks the approximation in quadruple precision for values
  ! that do not have a closed form expression.
  ! Usage is:
  !
  !- ``./quadtest N λ μ``
  !
  ! All three parameters are optional. The first one sets the number of
  ! quadrature points, λ μ are the parameters of the Wright function.
  ! The test can be called as:
  !
  !- ``./quadtest N``
  !- ``./quadtest N λ μ``
  !- ``./quadtest λ μ``
  !- ``./quadtest``
  !
  use iso_fortran_env, only : real32, real64, real128, error_unit, output_unit
  use wrightmod
  implicit none

  ! Integer parameter, number of points
  integer, parameter :: N = 100
  ! Integer number of quadrature points
  integer :: Nq
  character(len=25) :: outfilename

  ! Values
  real(real128) :: x(N)
  real(real128), allocatable :: w(:), warb(:)
  real(real128) lambda, mu, h
  real(real128) :: t = 1.0_real128
  character(len = 500) :: commandline
  ! Auxiliary variables
  integer i, exitstat,nlineread
  ! Input of the options from file
  character(len=12), dimension(:), allocatable :: args
  integer :: num_args, ix
  num_args = command_argument_count()
  if (num_args > 0) then
    allocate(args(num_args))
    do ix = 1, num_args
      call get_command_argument(ix,args(ix))
    end do
  end if

  select case (num_args)
  case (1)
    read(args(1),'(I4)')Nq
  case (2)
    read(args(1),'(F20.16)')lambda
    read(args(2),'(F20.16)')mu
  case (3)
    read(args(1),'(I4)')Nq
    read(args(2),'(F20.16)')lambda
    read(args(3),'(F20.16)')mu
  case default
    lambda = -0.5
    mu = 0.5
  end select

  write(output_unit,'("Quadruple precision progam with λ = ",F4.2," μ = ",F4.2)') &
    & lambda, mu

  h = 10.0_real128/(real(N,real128)-1.0)
  do i=1,N
     x(i) = -5.0_real128 + real((i-1),real128)*h
  end do
  ! Write on file
  write(output_unit,'("Write evaluation points")')
  open (unit = 99, file = "quadinput.inp")
  write(99,'(i5)')N
  do i=1,N
     write(99,'(F48.40)')-abs(x(i))
  end do
  close(99)

  ! Compute the value in the Fortran routine
  write(output_unit,'("Compute Wright function")')
  if ( (num_args == 1).or.(num_args == 3)  ) then
    w = wright(x, t, lambda, mu, Nq)
  else
    w = wright(x, t, lambda, mu)
  end if
  ! Compute the value with ARB
  write(commandline,&
    & '("./wrighttest ",F20.4," ",F20.4," 1000 256 quadinput.inp")')&
    & lambda,mu
  write(output_unit,'("Executing ARB: ",A)')trim(commandline)
  call execute_command_line(commandline, exitstat=exitstat)
  if (exitstat /= 0) then
    write(error_unit,'("Error in external command")')
  end if

  ! Read ARB values from file
  write(output_unit,'("Read ARB Wright function values")')
  allocate(warb(N))
  open (unit = 99, file = "wrighttest.out")
  do i=1,N
    read(99, '(F40.39)', iostat=exitstat)warb(i)
  end do
  close(99)

  ! Write the error
  if ( (num_args == 1).or.(num_args == 3) ) then
    write(outfilename,'("quadtest-",I2,".out")')Nq
  else
    write(outfilename,'("quadtest.out")')
  end if
  write(output_unit,'("Write output: ",A)')trim(outfilename)
  open (unit = 99, file = trim(outfilename))
  do i=1,N
    write(99,'(F20.16,",",F40.39,","F50.39,",",F50.39)')&
      & x(i),w(i),warb(i),abs(w(i)-warb(i))
  end do
  close(99)

  ! Clean memory
  deallocate(w,warb)
  if (allocated(args)) then
    deallocate(args)
  end if

end program quadtest
