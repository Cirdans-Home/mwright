module brentmod
! This module contains the implementation of Brent's algorithm for the
! minimization of a function :math:`f(x)` in the interval :math:`[a,b]`.
!
! The implementation is the one from:
!
!.. note::
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
! See: `https://people.math.sc.edu/Burkardt/f_src/brent/brent.html <https://people.math.sc.edu/Burkardt/f_src/brent/brent.html>`_
!
  use iso_fortran_env, only: real32, real64, real128
  implicit none

  interface fminbnd
    module procedure dglomin
    module procedure sglomin
    module procedure tglomin
  end interface

  private :: sglomin,slocal_min, slocal_min_rc, szero_rc, szero, & 
   & dglomin,dlocal_min, dlocal_min_rc, dzero_rc, dzero, &
   & tglomin,tlocal_min, tlocal_min_rc, tzero_rc, tzero

contains

  function dglomin ( a, b, c, m, machep, e, t, f, x )

!*****************************************************************************80
!
!! GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    This function assumes that F(X) is twice continuously differentiable
!    over [A,B] and that F''(X) <= M for all X in [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real64) A, B, the endpoints of the interval.
!    It must be the case that A < B.
!
!    Input, real(real64) C, an initial guess for the global
!    minimizer.  If no good guess is known, C = A or B is acceptable.
!
!    Input, real(real64) M, the bound on the second derivative.
!
!    Input, real(real64) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real(real64) E, a positive tolerance, a bound for the
!    absolute error in the evaluation of F(X) for any X in [A,B].
!
!    Input, real(real64) T, a positive error tolerance.
!
!    Input, external real(real64) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose global minimum is being sought.
!
!    Output, real(real64) X, the estimated value of the abscissa
!    for which F attains its global minimum value in [A,B].
!
!    Output, real(real64) GLOMIN, the value F(X).
!
  implicit none

  real(real64) a
  real(real64) a0
  real(real64) a2
  real(real64) a3
  real(real64) b
  real(real64) c
  real(real64) d0
  real(real64) d1
  real(real64) d2
  real(real64) e
  real(real64) f
  real(real64) dglomin
  real(real64) h
  integer ( kind = 4 ) k
  real(real64) m
  real(real64) m2
  real(real64) machep
  real(real64) p
  real(real64) q
  real(real64) qs
  real(real64) r
  real(real64) s
  real(real64) sc
  real(real64) t
  real(real64) x
  real(real64) y
  real(real64) y0
  real(real64) y1
  real(real64) y2
  real(real64) y3
  real(real64) yb
  real(real64) z0
  real(real64) z1
  real(real64) z2

  a0 = b
  x = a0
  a2 = a
  y0 = f ( b )
  yb = y0
  y2 = f ( a )
  y = y2

  if ( y0 < y ) then
    y = y0
  else
    x = a
  end if

  if ( m <= 0.0D+00 .or. b <= a ) then
    dglomin = y
    return
  end if

  m2 = 0.5D+00 * ( 1.0D+00 + 16.0D+00 * machep ) * m

  if ( c <= a .or. b <= c ) then
    sc = 0.5D+00 * ( a + b )
  else
    sc = c
  end if

  y1 = f ( sc )
  k = 3
  d0 = a2 - sc
  h = 9.0D+00 / 11.0D+00

  if ( y1 < y ) then
    x = sc
    y = y1
  end if

  do

    d1 = a2 - a0
    d2 = sc - a0
    z2 = b - a2
    z0 = y2 - y1
    z1 = y2 - y0
    r = d1 * d1 * z0 - d0 * d0 * z1
    p = r
    qs = 2.0D+00 * ( d0 * z1 - d1 * z0 )
    q = qs

    if ( k < 1000000 .or. y2 <= y ) then

      do

        if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
          z2 * m2 * r * ( z2 * q - r ) ) then
          a3 = a2 + r / q
          y3 = f ( a3 )

          if ( y3 < y ) then
            x = a3
            y = y3
          end if
        end if

        k = mod ( 1611 * k, 1048576 )
        q = 1.0D+00
        r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

        if ( z2 <= r ) then
          exit
        end if

      end do

    else

      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

      do while ( r < z2 )

        if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
          z2 * m2 * r * ( z2 * q - r ) ) then
          a3 = a2 + r / q
          y3 = f ( a3 )

          if ( y3 < y ) then
            x = a3
            y = y3
          end if
        end if

        k = mod ( 1611 * k, 1048576 )
        q = 1.0D+00
        r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

      end do

    end if

    r = m2 * d0 * d1 * d2
    s = sqrt ( ( ( y2 - y ) + t ) / m2 )
    h = 0.5D+00 * ( 1.0D+00 + h )
    p = h * ( p + 2.0D+00 * r * s )
    q = q + 0.5D+00 * qs
    r = - 0.5D+00 * ( d0 + ( z0 + 2.01D+00 * e ) / ( d0 * m2 ) )

    if ( r < s .or. d0 < 0.0D+00 ) then
      r = a2 + s
    else
      r = a2 + r
    end if

    if ( 0.0D+00 < p * q ) then
      a3 = a2 + p / q
    else
      a3 = r
    end if

    do

      a3 = max ( a3, r )

      if ( b <= a3 ) then
        a3 = b
        y3 = yb
      else
        y3 = f ( a3 )
      end if

      if ( y3 < y ) then
        x = a3
        y = y3
      end if

      d0 = a3 - a2

      if ( a3 <= r ) then
        exit
      end if

      p = 2.0D+00 * ( y2 - y3 ) / ( m * d0 )

      if ( ( 1.0D+00 + 9.0D+00 * machep ) * d0 <= abs ( p ) ) then
        exit
      end if

      if ( 0.5D+00 * m2 * ( d0 * d0 + p * p ) <= &
        ( y2 - y ) + ( y3 - y ) + 2.0D+00 * t ) then
        exit
      end if

      a3 = 0.5D+00 * ( a2 + a3 )
      h = 0.9D+00 * h

    end do

    if ( b <= a3 ) then
      exit
    end if

    a0 = sc
    sc = a2
    a2 = a3
    y0 = y1
    y1 = y2
    y2 = y3

  end do

  dglomin = y

  return
end
function dlocal_min ( a, b, eps, t, f, x )

!*****************************************************************************80
!
!! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then convergence is superlinear, and usually of the order of
!    about 1.324....
!
!    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.
!
!    If F is a unimodal function and the computed values of F are always
!    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
!    LOCAL_MIN approximates the abscissa of the global minimum of F on the
!    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
!
!    If F is not unimodal, then LOCAL_MIN may approximate a local, but
!    perhaps non-global, minimum to the same accuracy.
!
!    Thanks to Jonathan Eggleston for pointing out a correction to the
!    golden section step, 01 July 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real64) A, B, the endpoints of the interval.
!
!    Input, real(real64) EPS, a positive relative error tolerance.
!    EPS should be no smaller than twice the relative machine precision,
!    and preferably not much less than the square root of the relative
!    machine precision.
!
!    Input, real(real64) T, a positive absolute error tolerance.
!
!    Input, external real(real64) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Output, real(real64) X, the estimated value of an abscissa
!    for which F attains a local minimum value in [A,B].
!
!    Output, real(real64) LOCAL_MIN, the value F(X).
!
  implicit none

  real(real64) a
  real(real64) b
  real(real64) c
  real(real64) d
  real(real64) e
  real(real64) eps
  real(real64) f
  real(real64) fu
  real(real64) fv
  real(real64) fw
  real(real64) fx
  real(real64) dlocal_min
  real(real64) m
  real(real64) p
  real(real64) q
  real(real64) r
  real(real64) sa
  real(real64) sb
  real(real64) t
  real(real64) t2
  real(real64) tol
  real(real64) u
  real(real64) v
  real(real64) w
  real(real64) x
!
!  C is the square of the inverse of the golden ratio.
!
  c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

  sa = a
  sb = b
  x = sa + c * ( b - a )
  w = x
  v = w
  e = 0.0D+00
  fx = f ( x )
  fw = fx
  fv = fw

  do

    m = 0.5D+00 * ( sa + sb )
    tol = eps * abs ( x ) + t
    t2 = 2.0D+00 * tol
!
!  Check the stopping criterion.
!
    if ( abs ( x - m ) <= t2 - 0.5D+00 * ( sb - sa ) ) then
      exit
    end if
!
!  Fit a parabola.
!
    r = 0.0D+00
    q = r
    p = q

    if ( tol < abs ( e ) ) then

      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0D+00 * ( q - r )

      if ( 0.0D+00 < q ) then
        p = - p
      end if

      q = abs ( q )

      r = e
      e = d

    end if

    if ( abs ( p ) < abs ( 0.5D+00 * q * r ) .and. &
         q * ( sa - x ) < p .and. &
         p < q * ( sb - x ) ) then
!
!  Take the parabolic interpolation step.
!
      d = p / q
      u = x + d
!
!  F must not be evaluated too close to A or B.
!
      if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

        if ( x < m ) then
          d = tol
        else
          d = - tol
        end if

      end if
!
!  A golden-section step.
!
    else

      if ( x < m ) then
        e = sb - x
      else
        e = sa - x
      end if

      d = c * e

    end if
!
!  F must not be evaluated too close to X.
!
    if ( tol <= abs ( d ) ) then
      u = x + d
    else if ( 0.0D+00 < d ) then
      u = x + tol
    else
      u = x - tol
    end if

    fu = f ( u )
!
!  Update A, B, V, W, and X.
!
    if ( fu <= fx ) then

      if ( u < x ) then
        sb = x
      else
        sa = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if ( u < x ) then
        sa = u
      else
        sb = u
      end if

      if ( fu <= fw .or. w == x ) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if ( fu <= fv .or. v == x .or. v == w ) then
        v = u
        fv = fu
      end if

    end if

  end do

  dlocal_min = fx

  return
end
subroutine dlocal_min_rc ( a, b, arg, status, value )

!*****************************************************************************80
!
!! LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    This routine seeks an approximation to the point where a function
!    F attains a minimum on the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324...
!
!    The routine is a revised version of the Brent local minimization
!    algorithm, using reverse communication.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters
!
!    Input/output, real(real64) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for an interval containing the minimizer.  It is required
!    that A < B.
!
!    Output, real(real64) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS dzero, ARG is the routine's
!    estimate for the function minimizer.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to dzero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated minimizer.
!
!    Input, real(real64) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
!  Local parameters:
!
!    C is the squared inverse of the golden ratio.
!
!    EPS is the square root of the relative machine precision.
!
  implicit none

  real(real64) a
  real(real64) arg
  real(real64) b
  real(real64), save :: c
  real(real64), save :: d
  real(real64), save :: e
  real(real64), save :: eps
  real(real64), save :: fu
  real(real64), save :: fv
  real(real64), save :: fw
  real(real64), save :: fx
  real(real64), save :: midpoint
  real(real64), save :: p
  real(real64), save :: q
  real(real64), save :: r
  integer ( kind = 4 ) status
  real(real64), save :: tol
  real(real64), save :: tol1
  real(real64), save :: tol2
  real(real64), save :: u
  real(real64), save :: v
  real(real64) value
  real(real64), save :: w
  real(real64), save :: x
!
!  STATUS (INPUT) = 0, startup.
!
  if ( status == 0 ) then

    if ( b <= a ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LOCAL_MIN_RC - Fatal error!'
      write ( *, '(a)' ) '  A < B is required, but'
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b
      status = -1
      stop 1
    end if

    c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

    eps = sqrt ( epsilon ( eps ) )
    tol = epsilon ( tol )

    v = a + c * ( b - a )
    w = v
    x = v
    e = 0.0D+00

    status = 1
    arg = x

    return
!
!  STATUS (INPUT) = 1, return with initial function value of FX.
!
  else if ( status == 1 ) then

    fx = value
    fv = fx
    fw = fx
!
!  STATUS (INPUT) = 2 or more, update the data.
!
  else if ( 2 <= status ) then

    fu = value

    if ( fu <= fx ) then

      if ( x <= u ) then
        a = x
      else
        b = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if ( u < x ) then
        a = u
      else
        b = u
      end if

      if ( fu <= fw .or. w == x ) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if ( fu <= fv .or. v == x .or. v == w ) then
        v = u
        fv = fu
      end if

    end if

  end if
!
!  Take the next step.
!
  midpoint = 0.5D+00 * ( a + b )
  tol1 = eps * abs ( x ) + tol / 3.0D+00
  tol2 = 2.0D+00 * tol1
!
!  If the stopping criterion is satisfied, we can exit.
!
  if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
    status = 0
    return
  end if
!
!  Is golden-section necessary?
!
  if ( abs ( e ) <= tol1 ) then

    if ( midpoint <= x ) then
      e = a - x
    else
      e = b - x
    end if

    d = c * e
!
!  Consider fitting a parabola.
!
  else

    r = ( x - w ) * ( fx - fv )
    q = ( x - v ) * ( fx - fw )
    p = ( x - v ) * q - ( x - w ) * r
    q = 2.0D+00 * ( q - r )
    if ( 0.0D+00 < q ) then
      p = - p
    end if
    q = abs ( q )
    r = e
    e = d
!
!  Choose a golden-section step if the parabola is not advised.
!
    if ( &
      ( abs ( 0.5D+00 * q * r ) <= abs ( p ) ) .or. &
      ( p <= q * ( a - x ) ) .or. &
      ( q * ( b - x ) <= p ) ) then

      if ( midpoint <= x ) then
        e = a - x
      else
        e = b - x
      end if

      d = c * e
!
!  Choose a parabolic interpolation step.
!
    else

      d = p / q
      u = x + d

      if ( ( u - a ) < tol2 ) then
        d = sign ( tol1, midpoint - x )
      end if

      if ( ( b - u ) < tol2 ) then
        d = sign ( tol1, midpoint - x )
      end if

    end if

  end if
!
!  F must not be evaluated too close to X.
!
  if ( tol1 <= abs ( d ) ) then
    u = x + d
  end if

  if ( abs ( d ) < tol1 ) then
    u = x + sign ( tol1, d )
  end if
!
!  Request value of F(U).
!
  arg = u
  status = status + 1

  return
end

function dzero ( a, b, machep, t, f )

!*****************************************************************************80
!
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the dzero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real64) A, B, the endpoints of the change of
!    sign interval.
!
!    Input, real(real64) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real(real64) T, a positive error tolerance.
!
!    Input, external real(real64) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose dzero is being sought.
!
!    Output, real(real64) ZERO, the estimated value of a dzero of
!    the function F.
!
  implicit none

  real(real64) a
  real(real64) b
  real(real64) c
  real(real64) d
  real(real64) e
  real(real64) f
  real(real64) fa
  real(real64) fb
  real(real64) fc
  real(real64) m
  real(real64) machep
  real(real64) p
  real(real64) q
  real(real64) r
  real(real64) s
  real(real64) sa
  real(real64) sb
  real(real64) t
  real(real64) tol
  real(real64) dzero
!
!  Make local copies of A and B.
!
  sa = a
  sb = b
  fa = f ( sa )
  fb = f ( sb )

  c = sa
  fc = fa
  e = sb - sa
  d = e

  do

    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0D+00 * machep * abs ( sb ) + t
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb / fa

      if ( sa == c ) then

        p = 2.0D+00 * m * s
        q = 1.0D+00 - s

      else

        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if

    fb = f ( sb )

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end do

  dzero = sb

  return
end
subroutine dzero_rc ( a, b, t, arg, status, value )

!*****************************************************************************80
!
!! ZERO_RC seeks the root of a function F(X) using reverse communication.
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the dzero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    The routine is a revised version of the Brent dzero finder
!    algorithm, using reverse communication.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real64) A, B, the endpoints of the change of
!    sign interval.
!
!    Input, real(real64) T, a positive error tolerance.
!
!    Output, real(real64) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS dzero, ARG is the routine's
!    estimate for the function's dzero.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to dzero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated dzero
!
!    Input, real(real64) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
  implicit none

  real(real64) a
  real(real64) arg
  real(real64) b
  real(real64), save :: c
  real(real64), save :: d
  real(real64), save :: e
  real(real64), save :: fa
  real(real64), save :: fb
  real(real64), save :: fc
  real(real64) m
  real(real64), save :: machep
  real(real64) p
  real(real64) q
  real(real64) r
  real(real64) s
  real(real64), save :: sa
  real(real64), save :: sb
  integer ( kind = 4 ) status
  real(real64) t
  real(real64) tol
  real(real64) value
!
!  Input STATUS = 0.
!  Initialize, request F(A).
!
  if ( status == 0 ) then

    machep = epsilon ( a )

    sa = a
    sb = b
    e = sb - sa
    d = e

    status = 1
    arg = a
    return
!
!  Input STATUS = 1.
!  Receive F(A), request F(B).
!
  else if ( status == 1 ) then

    fa = value

    status = 2
    arg = sb
    return
!
!  Input STATUS = 2
!  Receive F(B).
!
  else if ( status == 2 ) then

    fb = value

    if ( 0.0D+00 < fa * fb ) then
      status = -1
      return
    end if

    c = sa
    fc = fa

  else

    fb = value

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end if
!
!  Compute the next point at which a function value is requested.
!
  if ( abs ( fc ) < abs ( fb ) ) then

    sa = sb
    sb = c
    c = sa
    fa = fb
    fb = fc
    fc = fa

  end if

  tol = 2.0D+00 * machep * abs ( sb ) + t
  m = 0.5D+00 * ( c - sb )

  if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
    status = 0
    arg = sb
    return
  end if

  if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

    e = m
    d = e

  else

    s = fb / fa

    if ( sa == c ) then

      p = 2.0D+00 * m * s
      q = 1.0D+00 - s

    else

      q = fa / fc
      r = fb / fc
      p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
      q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

    end if

    if ( 0.0D+00 < p ) then
      q = - q
    else
      p = - p
    end if

    s = e
    e = d

    if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
      p < abs ( 0.5D+00 * s * q ) ) then
      d = p / q
    else
      e = m
      d = e
    end if

  end if

  sa = sb
  fa = fb

  if ( tol < abs ( d ) ) then
    sb = sb + d
  else if ( 0.0D+00 < m ) then
    sb = sb + tol
  else
    sb = sb - tol
  end if

  arg = sb
  status = status + 1

  return
end

function sglomin ( a, b, c, m, machep, e, t, f, x )

!*****************************************************************************80
!
!! GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    This function assumes that F(X) is twice continuously differentiable
!    over [A,B] and that F''(X) <= M for all X in [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real32) A, B, the endpoints of the interval.
!    It must be the case that A < B.
!
!    Input, real(real32) C, an initial guess for the global
!    minimizer.  If no good guess is known, C = A or B is acceptable.
!
!    Input, real(real32) M, the bound on the second derivative.
!
!    Input, real(real32) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real(real32) E, a positive tolerance, a bound for the
!    absolute error in the evaluation of F(X) for any X in [A,B].
!
!    Input, real(real32) T, a positive error tolerance.
!
!    Input, external real(real32) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose global minimum is being sought.
!
!    Output, real(real32) X, the estimated value of the abscissa
!    for which F attains its global minimum value in [A,B].
!
!    Output, real(real32) GLOMIN, the value F(X).
!
implicit none

real(real32) a
real(real32) a0
real(real32) a2
real(real32) a3
real(real32) b
real(real32) c
real(real32) d0
real(real32) d1
real(real32) d2
real(real32) e
real(real32) f
real(real32) sglomin
real(real32) h
integer ( kind = 4 ) k
real(real32) m
real(real32) m2
real(real32) machep
real(real32) p
real(real32) q
real(real32) qs
real(real32) r
real(real32) s
real(real32) sc
real(real32) t
real(real32) x
real(real32) y
real(real32) y0
real(real32) y1
real(real32) y2
real(real32) y3
real(real32) yb
real(real32) z0
real(real32) z1
real(real32) z2

a0 = b
x = a0
a2 = a
y0 = f ( b )
yb = y0
y2 = f ( a )
y = y2

if ( y0 < y ) then
  y = y0
else
  x = a
end if

if ( m <= 0.0D+00 .or. b <= a ) then
  sglomin = y
  return
end if

m2 = 0.5D+00 * ( 1.0D+00 + 16.0D+00 * machep ) * m

if ( c <= a .or. b <= c ) then
  sc = 0.5D+00 * ( a + b )
else
  sc = c
end if

y1 = f ( sc )
k = 3
d0 = a2 - sc
h = 9.0D+00 / 11.0D+00

if ( y1 < y ) then
  x = sc
  y = y1
end if

do

  d1 = a2 - a0
  d2 = sc - a0
  z2 = b - a2
  z0 = y2 - y1
  z1 = y2 - y0
  r = d1 * d1 * z0 - d0 * d0 * z1
  p = r
  qs = 2.0D+00 * ( d0 * z1 - d1 * z0 )
  q = qs

  if ( k < 1000000 .or. y2 <= y ) then

    do

      if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
        z2 * m2 * r * ( z2 * q - r ) ) then
        a3 = a2 + r / q
        y3 = f ( a3 )

        if ( y3 < y ) then
          x = a3
          y = y3
        end if
      end if

      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

      if ( z2 <= r ) then
        exit
      end if

    end do

  else

    k = mod ( 1611 * k, 1048576 )
    q = 1.0D+00
    r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

    do while ( r < z2 )

      if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
        z2 * m2 * r * ( z2 * q - r ) ) then
        a3 = a2 + r / q
        y3 = f ( a3 )

        if ( y3 < y ) then
          x = a3
          y = y3
        end if
      end if

      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

    end do

  end if

  r = m2 * d0 * d1 * d2
  s = sqrt ( ( ( y2 - y ) + t ) / m2 )
  h = 0.5D+00 * ( 1.0D+00 + h )
  p = h * ( p + 2.0D+00 * r * s )
  q = q + 0.5D+00 * qs
  r = - 0.5D+00 * ( d0 + ( z0 + 2.01D+00 * e ) / ( d0 * m2 ) )

  if ( r < s .or. d0 < 0.0D+00 ) then
    r = a2 + s
  else
    r = a2 + r
  end if

  if ( 0.0D+00 < p * q ) then
    a3 = a2 + p / q
  else
    a3 = r
  end if

  do

    a3 = max ( a3, r )

    if ( b <= a3 ) then
      a3 = b
      y3 = yb
    else
      y3 = f ( a3 )
    end if

    if ( y3 < y ) then
      x = a3
      y = y3
    end if

    d0 = a3 - a2

    if ( a3 <= r ) then
      exit
    end if

    p = 2.0D+00 * ( y2 - y3 ) / ( m * d0 )

    if ( ( 1.0D+00 + 9.0D+00 * machep ) * d0 <= abs ( p ) ) then
      exit
    end if

    if ( 0.5D+00 * m2 * ( d0 * d0 + p * p ) <= &
      ( y2 - y ) + ( y3 - y ) + 2.0D+00 * t ) then
      exit
    end if

    a3 = 0.5D+00 * ( a2 + a3 )
    h = 0.9D+00 * h

  end do

  if ( b <= a3 ) then
    exit
  end if

  a0 = sc
  sc = a2
  a2 = a3
  y0 = y1
  y1 = y2
  y2 = y3

end do

sglomin = y

return
end
function slocal_min ( a, b, eps, t, f, x )

!*****************************************************************************80
!
!! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then convergence is superlinear, and usually of the order of
!    about 1.324....
!
!    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.
!
!    If F is a unimodal function and the computed values of F are always
!    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
!    LOCAL_MIN approximates the abscissa of the global minimum of F on the
!    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
!
!    If F is not unimodal, then LOCAL_MIN may approximate a local, but
!    perhaps non-global, minimum to the same accuracy.
!
!    Thanks to Jonathan Eggleston for pointing out a correction to the
!    golden section step, 01 July 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real32) A, B, the endpoints of the interval.
!
!    Input, real(real32) EPS, a positive relative error tolerance.
!    EPS should be no smaller than twice the relative machine precision,
!    and preferably not much less than the square root of the relative
!    machine precision.
!
!    Input, real(real32) T, a positive absolute error tolerance.
!
!    Input, external real(real32) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Output, real(real32) X, the estimated value of an abscissa
!    for which F attains a local minimum value in [A,B].
!
!    Output, real(real32) LOCAL_MIN, the value F(X).
!
implicit none

real(real32) a
real(real32) b
real(real32) c
real(real32) d
real(real32) e
real(real32) eps
real(real32) f
real(real32) fu
real(real32) fv
real(real32) fw
real(real32) fx
real(real32) slocal_min
real(real32) m
real(real32) p
real(real32) q
real(real32) r
real(real32) sa
real(real32) sb
real(real32) t
real(real32) t2
real(real32) tol
real(real32) u
real(real32) v
real(real32) w
real(real32) x
!
!  C is the square of the inverse of the golden ratio.
!
c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

sa = a
sb = b
x = sa + c * ( b - a )
w = x
v = w
e = 0.0D+00
fx = f ( x )
fw = fx
fv = fw

do

  m = 0.5D+00 * ( sa + sb )
  tol = eps * abs ( x ) + t
  t2 = 2.0D+00 * tol
!
!  Check the stopping criterion.
!
  if ( abs ( x - m ) <= t2 - 0.5D+00 * ( sb - sa ) ) then
    exit
  end if
!
!  Fit a parabola.
!
  r = 0.0D+00
  q = r
  p = q

  if ( tol < abs ( e ) ) then

    r = ( x - w ) * ( fx - fv )
    q = ( x - v ) * ( fx - fw )
    p = ( x - v ) * q - ( x - w ) * r
    q = 2.0D+00 * ( q - r )

    if ( 0.0D+00 < q ) then
      p = - p
    end if

    q = abs ( q )

    r = e
    e = d

  end if

  if ( abs ( p ) < abs ( 0.5D+00 * q * r ) .and. &
       q * ( sa - x ) < p .and. &
       p < q * ( sb - x ) ) then
!
!  Take the parabolic interpolation step.
!
    d = p / q
    u = x + d
!
!  F must not be evaluated too close to A or B.
!
    if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

      if ( x < m ) then
        d = tol
      else
        d = - tol
      end if

    end if
!
!  A golden-section step.
!
  else

    if ( x < m ) then
      e = sb - x
    else
      e = sa - x
    end if

    d = c * e

  end if
!
!  F must not be evaluated too close to X.
!
  if ( tol <= abs ( d ) ) then
    u = x + d
  else if ( 0.0D+00 < d ) then
    u = x + tol
  else
    u = x - tol
  end if

  fu = f ( u )
!
!  Update A, B, V, W, and X.
!
  if ( fu <= fx ) then

    if ( u < x ) then
      sb = x
    else
      sa = x
    end if

    v = w
    fv = fw
    w = x
    fw = fx
    x = u
    fx = fu

  else

    if ( u < x ) then
      sa = u
    else
      sb = u
    end if

    if ( fu <= fw .or. w == x ) then
      v = w
      fv = fw
      w = u
      fw = fu
    else if ( fu <= fv .or. v == x .or. v == w ) then
      v = u
      fv = fu
    end if

  end if

end do

slocal_min = fx

return
end
subroutine slocal_min_rc ( a, b, arg, status, value )

!*****************************************************************************80
!
!! LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    This routine seeks an approximation to the point where a function
!    F attains a minimum on the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324...
!
!    The routine is a revised version of the Brent local minimization
!    algorithm, using reverse communication.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters
!
!    Input/output, real(real32) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for an interval containing the minimizer.  It is required
!    that A < B.
!
!    Output, real(real32) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS szero, ARG is the routine's
!    estimate for the function minimizer.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to szero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated minimizer.
!
!    Input, real(real32) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
!  Local parameters:
!
!    C is the squared inverse of the golden ratio.
!
!    EPS is the square root of the relative machine precision.
!
implicit none

real(real32) a
real(real32) arg
real(real32) b
real(real32), save :: c
real(real32), save :: d
real(real32), save :: e
real(real32), save :: eps
real(real32), save :: fu
real(real32), save :: fv
real(real32), save :: fw
real(real32), save :: fx
real(real32), save :: midpoint
real(real32), save :: p
real(real32), save :: q
real(real32), save :: r
integer ( kind = 4 ) status
real(real32), save :: tol
real(real32), save :: tol1
real(real32), save :: tol2
real(real32), save :: u
real(real32), save :: v
real(real32) value
real(real32), save :: w
real(real32), save :: x
!
!  STATUS (INPUT) = 0, startup.
!
if ( status == 0 ) then

  if ( b <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOCAL_MIN_RC - Fatal error!'
    write ( *, '(a)' ) '  A < B is required, but'
    write ( *, '(a,g14.6)' ) '  A = ', a
    write ( *, '(a,g14.6)' ) '  B = ', b
    status = -1
    stop 1
  end if

  c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

  eps = sqrt ( epsilon ( eps ) )
  tol = epsilon ( tol )

  v = a + c * ( b - a )
  w = v
  x = v
  e = 0.0D+00

  status = 1
  arg = x

  return
!
!  STATUS (INPUT) = 1, return with initial function value of FX.
!
else if ( status == 1 ) then

  fx = value
  fv = fx
  fw = fx
!
!  STATUS (INPUT) = 2 or more, update the data.
!
else if ( 2 <= status ) then

  fu = value

  if ( fu <= fx ) then

    if ( x <= u ) then
      a = x
    else
      b = x
    end if

    v = w
    fv = fw
    w = x
    fw = fx
    x = u
    fx = fu

  else

    if ( u < x ) then
      a = u
    else
      b = u
    end if

    if ( fu <= fw .or. w == x ) then
      v = w
      fv = fw
      w = u
      fw = fu
    else if ( fu <= fv .or. v == x .or. v == w ) then
      v = u
      fv = fu
    end if

  end if

end if
!
!  Take the next step.
!
midpoint = 0.5D+00 * ( a + b )
tol1 = eps * abs ( x ) + tol / 3.0D+00
tol2 = 2.0D+00 * tol1
!
!  If the stopping criterion is satisfied, we can exit.
!
if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
  status = 0
  return
end if
!
!  Is golden-section necessary?
!
if ( abs ( e ) <= tol1 ) then

  if ( midpoint <= x ) then
    e = a - x
  else
    e = b - x
  end if

  d = c * e
!
!  Consider fitting a parabola.
!
else

  r = ( x - w ) * ( fx - fv )
  q = ( x - v ) * ( fx - fw )
  p = ( x - v ) * q - ( x - w ) * r
  q = 2.0D+00 * ( q - r )
  if ( 0.0D+00 < q ) then
    p = - p
  end if
  q = abs ( q )
  r = e
  e = d
!
!  Choose a golden-section step if the parabola is not advised.
!
  if ( &
    ( abs ( 0.5D+00 * q * r ) <= abs ( p ) ) .or. &
    ( p <= q * ( a - x ) ) .or. &
    ( q * ( b - x ) <= p ) ) then

    if ( midpoint <= x ) then
      e = a - x
    else
      e = b - x
    end if

    d = c * e
!
!  Choose a parabolic interpolation step.
!
  else

    d = p / q
    u = x + d

    if ( ( u - a ) < tol2 ) then
      d = sign ( tol1, midpoint - x )
    end if

    if ( ( b - u ) < tol2 ) then
      d = sign ( tol1, midpoint - x )
    end if

  end if

end if
!
!  F must not be evaluated too close to X.
!
if ( tol1 <= abs ( d ) ) then
  u = x + d
end if

if ( abs ( d ) < tol1 ) then
  u = x + sign ( tol1, d )
end if
!
!  Request value of F(U).
!
arg = u
status = status + 1

return
end

function szero ( a, b, machep, t, f )

!*****************************************************************************80
!
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the szero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real32) A, B, the endpoints of the change of
!    sign interval.
!
!    Input, real(real32) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real(real32) T, a positive error tolerance.
!
!    Input, external real(real32) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose szero is being sought.
!
!    Output, real(real32) ZERO, the estimated value of a szero of
!    the function F.
!
implicit none

real(real32) a
real(real32) b
real(real32) c
real(real32) d
real(real32) e
real(real32) f
real(real32) fa
real(real32) fb
real(real32) fc
real(real32) m
real(real32) machep
real(real32) p
real(real32) q
real(real32) r
real(real32) s
real(real32) sa
real(real32) sb
real(real32) t
real(real32) tol
real(real32) szero
!
!  Make local copies of A and B.
!
sa = a
sb = b
fa = f ( sa )
fb = f ( sb )

c = sa
fc = fa
e = sb - sa
d = e

do

  if ( abs ( fc ) < abs ( fb ) ) then

    sa = sb
    sb = c
    c = sa
    fa = fb
    fb = fc
    fc = fa

  end if

  tol = 2.0D+00 * machep * abs ( sb ) + t
  m = 0.5D+00 * ( c - sb )

  if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
    exit
  end if

  if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

    e = m
    d = e

  else

    s = fb / fa

    if ( sa == c ) then

      p = 2.0D+00 * m * s
      q = 1.0D+00 - s

    else

      q = fa / fc
      r = fb / fc
      p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
      q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

    end if

    if ( 0.0D+00 < p ) then
      q = - q
    else
      p = - p
    end if

    s = e
    e = d

    if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
      p < abs ( 0.5D+00 * s * q ) ) then
      d = p / q
    else
      e = m
      d = e
    end if

  end if

  sa = sb
  fa = fb

  if ( tol < abs ( d ) ) then
    sb = sb + d
  else if ( 0.0D+00 < m ) then
    sb = sb + tol
  else
    sb = sb - tol
  end if

  fb = f ( sb )

  if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
       ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
    c = sa
    fc = fa
    e = sb - sa
    d = e
  end if

end do

szero = sb

return
end
subroutine szero_rc ( a, b, t, arg, status, value )

!*****************************************************************************80
!
!! ZERO_RC seeks the root of a function F(X) using reverse communication.
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the szero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    The routine is a revised version of the Brent szero finder
!    algorithm, using reverse communication.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real32) A, B, the endpoints of the change of
!    sign interval.
!
!    Input, real(real32) T, a positive error tolerance.
!
!    Output, real(real32) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS szero, ARG is the routine's
!    estimate for the function's szero.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to szero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated szero
!
!    Input, real(real32) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
implicit none

real(real32) a
real(real32) arg
real(real32) b
real(real32), save :: c
real(real32), save :: d
real(real32), save :: e
real(real32), save :: fa
real(real32), save :: fb
real(real32), save :: fc
real(real32) m
real(real32), save :: machep
real(real32) p
real(real32) q
real(real32) r
real(real32) s
real(real32), save :: sa
real(real32), save :: sb
integer ( kind = 4 ) status
real(real32) t
real(real32) tol
real(real32) value
!
!  Input STATUS = 0.
!  Initialize, request F(A).
!
if ( status == 0 ) then

  machep = epsilon ( a )

  sa = a
  sb = b
  e = sb - sa
  d = e

  status = 1
  arg = a
  return
!
!  Input STATUS = 1.
!  Receive F(A), request F(B).
!
else if ( status == 1 ) then

  fa = value

  status = 2
  arg = sb
  return
!
!  Input STATUS = 2
!  Receive F(B).
!
else if ( status == 2 ) then

  fb = value

  if ( 0.0D+00 < fa * fb ) then
    status = -1
    return
  end if

  c = sa
  fc = fa

else

  fb = value

  if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
       ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
    c = sa
    fc = fa
    e = sb - sa
    d = e
  end if

end if
!
!  Compute the next point at which a function value is requested.
!
if ( abs ( fc ) < abs ( fb ) ) then

  sa = sb
  sb = c
  c = sa
  fa = fb
  fb = fc
  fc = fa

end if

tol = 2.0D+00 * machep * abs ( sb ) + t
m = 0.5D+00 * ( c - sb )

if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
  status = 0
  arg = sb
  return
end if

if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

  e = m
  d = e

else

  s = fb / fa

  if ( sa == c ) then

    p = 2.0D+00 * m * s
    q = 1.0D+00 - s

  else

    q = fa / fc
    r = fb / fc
    p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
    q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

  end if

  if ( 0.0D+00 < p ) then
    q = - q
  else
    p = - p
  end if

  s = e
  e = d

  if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
    p < abs ( 0.5D+00 * s * q ) ) then
    d = p / q
  else
    e = m
    d = e
  end if

end if

sa = sb
fa = fb

if ( tol < abs ( d ) ) then
  sb = sb + d
else if ( 0.0D+00 < m ) then
  sb = sb + tol
else
  sb = sb - tol
end if

arg = sb
status = status + 1

return
end

function tglomin ( a, b, c, m, machep, e, t, f, x )

!*****************************************************************************80
!
!! GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    This function assumes that F(X) is twice continuously differentiable
!    over [A,B] and that F''(X) <= M for all X in [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real128) A, B, the endpoints of the interval.
!    It must be the case that A < B.
!
!    Input, real(real128) C, an initial guess for the global
!    minimizer.  If no good guess is known, C = A or B is acceptable.
!
!    Input, real(real128) M, the bound on the second derivative.
!
!    Input, real(real128) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real(real128) E, a positive tolerance, a bound for the
!    absolute error in the evaluation of F(X) for any X in [A,B].
!
!    Input, real(real128) T, a positive error tolerance.
!
!    Input, external real(real128) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose global minimum is being sought.
!
!    Output, real(real128) X, the estimated value of the abscissa
!    for which F attains its global minimum value in [A,B].
!
!    Output, real(real128) GLOMIN, the value F(X).
!
implicit none

real(real128) a
real(real128) a0
real(real128) a2
real(real128) a3
real(real128) b
real(real128) c
real(real128) d0
real(real128) d1
real(real128) d2
real(real128) e
real(real128) f
real(real128) tglomin
real(real128) h
integer ( kind = 4 ) k
real(real128) m
real(real128) m2
real(real128) machep
real(real128) p
real(real128) q
real(real128) qs
real(real128) r
real(real128) s
real(real128) sc
real(real128) t
real(real128) x
real(real128) y
real(real128) y0
real(real128) y1
real(real128) y2
real(real128) y3
real(real128) yb
real(real128) z0
real(real128) z1
real(real128) z2

a0 = b
x = a0
a2 = a
y0 = f ( b )
yb = y0
y2 = f ( a )
y = y2

if ( y0 < y ) then
  y = y0
else
  x = a
end if

if ( m <= 0.0D+00 .or. b <= a ) then
  tglomin = y
  return
end if

m2 = 0.5D+00 * ( 1.0D+00 + 16.0D+00 * machep ) * m

if ( c <= a .or. b <= c ) then
  sc = 0.5D+00 * ( a + b )
else
  sc = c
end if

y1 = f ( sc )
k = 3
d0 = a2 - sc
h = 9.0D+00 / 11.0D+00

if ( y1 < y ) then
  x = sc
  y = y1
end if

do

  d1 = a2 - a0
  d2 = sc - a0
  z2 = b - a2
  z0 = y2 - y1
  z1 = y2 - y0
  r = d1 * d1 * z0 - d0 * d0 * z1
  p = r
  qs = 2.0D+00 * ( d0 * z1 - d1 * z0 )
  q = qs

  if ( k < 1000000 .or. y2 <= y ) then

    do

      if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
        z2 * m2 * r * ( z2 * q - r ) ) then
        a3 = a2 + r / q
        y3 = f ( a3 )

        if ( y3 < y ) then
          x = a3
          y = y3
        end if
      end if

      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

      if ( z2 <= r ) then
        exit
      end if

    end do

  else

    k = mod ( 1611 * k, 1048576 )
    q = 1.0D+00
    r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

    do while ( r < z2 )

      if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < &
        z2 * m2 * r * ( z2 * q - r ) ) then
        a3 = a2 + r / q
        y3 = f ( a3 )

        if ( y3 < y ) then
          x = a3
          y = y3
        end if
      end if

      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * real ( k, kind = 8 )

    end do

  end if

  r = m2 * d0 * d1 * d2
  s = sqrt ( ( ( y2 - y ) + t ) / m2 )
  h = 0.5D+00 * ( 1.0D+00 + h )
  p = h * ( p + 2.0D+00 * r * s )
  q = q + 0.5D+00 * qs
  r = - 0.5D+00 * ( d0 + ( z0 + 2.01D+00 * e ) / ( d0 * m2 ) )

  if ( r < s .or. d0 < 0.0D+00 ) then
    r = a2 + s
  else
    r = a2 + r
  end if

  if ( 0.0D+00 < p * q ) then
    a3 = a2 + p / q
  else
    a3 = r
  end if

  do

    a3 = max ( a3, r )

    if ( b <= a3 ) then
      a3 = b
      y3 = yb
    else
      y3 = f ( a3 )
    end if

    if ( y3 < y ) then
      x = a3
      y = y3
    end if

    d0 = a3 - a2

    if ( a3 <= r ) then
      exit
    end if

    p = 2.0D+00 * ( y2 - y3 ) / ( m * d0 )

    if ( ( 1.0D+00 + 9.0D+00 * machep ) * d0 <= abs ( p ) ) then
      exit
    end if

    if ( 0.5D+00 * m2 * ( d0 * d0 + p * p ) <= &
      ( y2 - y ) + ( y3 - y ) + 2.0D+00 * t ) then
      exit
    end if

    a3 = 0.5D+00 * ( a2 + a3 )
    h = 0.9D+00 * h

  end do

  if ( b <= a3 ) then
    exit
  end if

  a0 = sc
  sc = a2
  a2 = a3
  y0 = y1
  y1 = y2
  y2 = y3

end do

tglomin = y

return
end
function tlocal_min ( a, b, eps, t, f, x )

!*****************************************************************************80
!
!! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then convergence is superlinear, and usually of the order of
!    about 1.324....
!
!    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.
!
!    If F is a unimodal function and the computed values of F are always
!    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
!    LOCAL_MIN approximates the abscissa of the global minimum of F on the
!    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
!
!    If F is not unimodal, then LOCAL_MIN may approximate a local, but
!    perhaps non-global, minimum to the same accuracy.
!
!    Thanks to Jonathan Eggleston for pointing out a correction to the
!    golden section step, 01 July 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real128) A, B, the endpoints of the interval.
!
!    Input, real(real128) EPS, a positive relative error tolerance.
!    EPS should be no smaller than twice the relative machine precision,
!    and preferably not much less than the square root of the relative
!    machine precision.
!
!    Input, real(real128) T, a positive absolute error tolerance.
!
!    Input, external real(real128) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Output, real(real128) X, the estimated value of an abscissa
!    for which F attains a local minimum value in [A,B].
!
!    Output, real(real128) LOCAL_MIN, the value F(X).
!
implicit none

real(real128) a
real(real128) b
real(real128) c
real(real128) d
real(real128) e
real(real128) eps
real(real128) f
real(real128) fu
real(real128) fv
real(real128) fw
real(real128) fx
real(real128) tlocal_min
real(real128) m
real(real128) p
real(real128) q
real(real128) r
real(real128) sa
real(real128) sb
real(real128) t
real(real128) t2
real(real128) tol
real(real128) u
real(real128) v
real(real128) w
real(real128) x
!
!  C is the square of the inverse of the golden ratio.
!
c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

sa = a
sb = b
x = sa + c * ( b - a )
w = x
v = w
e = 0.0D+00
fx = f ( x )
fw = fx
fv = fw

do

  m = 0.5D+00 * ( sa + sb )
  tol = eps * abs ( x ) + t
  t2 = 2.0D+00 * tol
!
!  Check the stopping criterion.
!
  if ( abs ( x - m ) <= t2 - 0.5D+00 * ( sb - sa ) ) then
    exit
  end if
!
!  Fit a parabola.
!
  r = 0.0D+00
  q = r
  p = q

  if ( tol < abs ( e ) ) then

    r = ( x - w ) * ( fx - fv )
    q = ( x - v ) * ( fx - fw )
    p = ( x - v ) * q - ( x - w ) * r
    q = 2.0D+00 * ( q - r )

    if ( 0.0D+00 < q ) then
      p = - p
    end if

    q = abs ( q )

    r = e
    e = d

  end if

  if ( abs ( p ) < abs ( 0.5D+00 * q * r ) .and. &
       q * ( sa - x ) < p .and. &
       p < q * ( sb - x ) ) then
!
!  Take the parabolic interpolation step.
!
    d = p / q
    u = x + d
!
!  F must not be evaluated too close to A or B.
!
    if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

      if ( x < m ) then
        d = tol
      else
        d = - tol
      end if

    end if
!
!  A golden-section step.
!
  else

    if ( x < m ) then
      e = sb - x
    else
      e = sa - x
    end if

    d = c * e

  end if
!
!  F must not be evaluated too close to X.
!
  if ( tol <= abs ( d ) ) then
    u = x + d
  else if ( 0.0D+00 < d ) then
    u = x + tol
  else
    u = x - tol
  end if

  fu = f ( u )
!
!  Update A, B, V, W, and X.
!
  if ( fu <= fx ) then

    if ( u < x ) then
      sb = x
    else
      sa = x
    end if

    v = w
    fv = fw
    w = x
    fw = fx
    x = u
    fx = fu

  else

    if ( u < x ) then
      sa = u
    else
      sb = u
    end if

    if ( fu <= fw .or. w == x ) then
      v = w
      fv = fw
      w = u
      fw = fu
    else if ( fu <= fv .or. v == x .or. v == w ) then
      v = u
      fv = fu
    end if

  end if

end do

tlocal_min = fx

return
end
subroutine tlocal_min_rc ( a, b, arg, status, value )

!*****************************************************************************80
!
!! LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
!
!  Discussion:
!
!    This routine seeks an approximation to the point where a function
!    F attains a minimum on the interval (A,B).
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much
!    slower than that for a Fibonacci search.  If F has a continuous
!    second derivative which is positive at the minimum (which is not
!    at A or B), then convergence is superlinear, and usually of the
!    order of about 1.324...
!
!    The routine is a revised version of the Brent local minimization
!    algorithm, using reverse communication.
!
!    It is worth stating explicitly that this routine will NOT be
!    able to detect a minimizer that occurs at either initial endpoint
!    A or B.  If this is a concern to the user, then the user must
!    either ensure that the initial interval is larger, or to check
!    the function value at the returned minimizer against the values
!    at either endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters
!
!    Input/output, real(real128) A, B.  On input, the left and right
!    endpoints of the initial interval.  On output, the lower and upper
!    bounds for an interval containing the minimizer.  It is required
!    that A < B.
!
!    Output, real(real128) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS tzero, ARG is the routine's
!    estimate for the function minimizer.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to tzero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated minimizer.
!
!    Input, real(real128) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
!  Local parameters:
!
!    C is the squared inverse of the golden ratio.
!
!    EPS is the square root of the relative machine precision.
!
implicit none

real(real128) a
real(real128) arg
real(real128) b
real(real128), save :: c
real(real128), save :: d
real(real128), save :: e
real(real128), save :: eps
real(real128), save :: fu
real(real128), save :: fv
real(real128), save :: fw
real(real128), save :: fx
real(real128), save :: midpoint
real(real128), save :: p
real(real128), save :: q
real(real128), save :: r
integer ( kind = 4 ) status
real(real128), save :: tol
real(real128), save :: tol1
real(real128), save :: tol2
real(real128), save :: u
real(real128), save :: v
real(real128) value
real(real128), save :: w
real(real128), save :: x
!
!  STATUS (INPUT) = 0, startup.
!
if ( status == 0 ) then

  if ( b <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOCAL_MIN_RC - Fatal error!'
    write ( *, '(a)' ) '  A < B is required, but'
    write ( *, '(a,g14.6)' ) '  A = ', a
    write ( *, '(a,g14.6)' ) '  B = ', b
    status = -1
    stop 1
  end if

  c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

  eps = sqrt ( epsilon ( eps ) )
  tol = epsilon ( tol )

  v = a + c * ( b - a )
  w = v
  x = v
  e = 0.0D+00

  status = 1
  arg = x

  return
!
!  STATUS (INPUT) = 1, return with initial function value of FX.
!
else if ( status == 1 ) then

  fx = value
  fv = fx
  fw = fx
!
!  STATUS (INPUT) = 2 or more, update the data.
!
else if ( 2 <= status ) then

  fu = value

  if ( fu <= fx ) then

    if ( x <= u ) then
      a = x
    else
      b = x
    end if

    v = w
    fv = fw
    w = x
    fw = fx
    x = u
    fx = fu

  else

    if ( u < x ) then
      a = u
    else
      b = u
    end if

    if ( fu <= fw .or. w == x ) then
      v = w
      fv = fw
      w = u
      fw = fu
    else if ( fu <= fv .or. v == x .or. v == w ) then
      v = u
      fv = fu
    end if

  end if

end if
!
!  Take the next step.
!
midpoint = 0.5D+00 * ( a + b )
tol1 = eps * abs ( x ) + tol / 3.0D+00
tol2 = 2.0D+00 * tol1
!
!  If the stopping criterion is satisfied, we can exit.
!
if ( abs ( x - midpoint ) <= ( tol2 - 0.5D+00 * ( b - a ) ) ) then
  status = 0
  return
end if
!
!  Is golden-section necessary?
!
if ( abs ( e ) <= tol1 ) then

  if ( midpoint <= x ) then
    e = a - x
  else
    e = b - x
  end if

  d = c * e
!
!  Consider fitting a parabola.
!
else

  r = ( x - w ) * ( fx - fv )
  q = ( x - v ) * ( fx - fw )
  p = ( x - v ) * q - ( x - w ) * r
  q = 2.0D+00 * ( q - r )
  if ( 0.0D+00 < q ) then
    p = - p
  end if
  q = abs ( q )
  r = e
  e = d
!
!  Choose a golden-section step if the parabola is not advised.
!
  if ( &
    ( abs ( 0.5D+00 * q * r ) <= abs ( p ) ) .or. &
    ( p <= q * ( a - x ) ) .or. &
    ( q * ( b - x ) <= p ) ) then

    if ( midpoint <= x ) then
      e = a - x
    else
      e = b - x
    end if

    d = c * e
!
!  Choose a parabolic interpolation step.
!
  else

    d = p / q
    u = x + d

    if ( ( u - a ) < tol2 ) then
      d = sign ( tol1, midpoint - x )
    end if

    if ( ( b - u ) < tol2 ) then
      d = sign ( tol1, midpoint - x )
    end if

  end if

end if
!
!  F must not be evaluated too close to X.
!
if ( tol1 <= abs ( d ) ) then
  u = x + d
end if

if ( abs ( d ) < tol1 ) then
  u = x + sign ( tol1, d )
end if
!
!  Request value of F(U).
!
arg = u
status = status + 1

return
end

function tzero ( a, b, machep, t, f )

!*****************************************************************************80
!
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the tzero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real128) A, B, the endpoints of the change of
!    sign interval.
!
!    Input, real(real128) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real(real128) T, a positive error tolerance.
!
!    Input, external real(real128) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose tzero is being sought.
!
!    Output, real(real128) ZERO, the estimated value of a tzero of
!    the function F.
!
implicit none

real(real128) a
real(real128) b
real(real128) c
real(real128) d
real(real128) e
real(real128) f
real(real128) fa
real(real128) fb
real(real128) fc
real(real128) m
real(real128) machep
real(real128) p
real(real128) q
real(real128) r
real(real128) s
real(real128) sa
real(real128) sb
real(real128) t
real(real128) tol
real(real128) tzero
!
!  Make local copies of A and B.
!
sa = a
sb = b
fa = f ( sa )
fb = f ( sb )

c = sa
fc = fa
e = sb - sa
d = e

do

  if ( abs ( fc ) < abs ( fb ) ) then

    sa = sb
    sb = c
    c = sa
    fa = fb
    fb = fc
    fc = fa

  end if

  tol = 2.0D+00 * machep * abs ( sb ) + t
  m = 0.5D+00 * ( c - sb )

  if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
    exit
  end if

  if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

    e = m
    d = e

  else

    s = fb / fa

    if ( sa == c ) then

      p = 2.0D+00 * m * s
      q = 1.0D+00 - s

    else

      q = fa / fc
      r = fb / fc
      p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
      q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

    end if

    if ( 0.0D+00 < p ) then
      q = - q
    else
      p = - p
    end if

    s = e
    e = d

    if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
      p < abs ( 0.5D+00 * s * q ) ) then
      d = p / q
    else
      e = m
      d = e
    end if

  end if

  sa = sb
  fa = fb

  if ( tol < abs ( d ) ) then
    sb = sb + d
  else if ( 0.0D+00 < m ) then
    sb = sb + tol
  else
    sb = sb - tol
  end if

  fb = f ( sb )

  if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
       ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
    c = sa
    fc = fa
    e = sb - sa
    d = e
  end if

end do

tzero = sb

return
end
subroutine tzero_rc ( a, b, t, arg, status, value )

!*****************************************************************************80
!
!! ZERO_RC seeks the root of a function F(X) using reverse communication.
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the tzero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    The routine is a revised version of the Brent tzero finder
!    algorithm, using reverse communication.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real(real128) A, B, the endpoints of the change of
!    sign interval.
!
!    Input, real(real128) T, a positive error tolerance.
!
!    Output, real(real128) ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS tzero, ARG is the routine's
!    estimate for the function's tzero.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to tzero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated tzero
!
!    Input, real(real128) VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
implicit none

real(real128) a
real(real128) arg
real(real128) b
real(real128), save :: c
real(real128), save :: d
real(real128), save :: e
real(real128), save :: fa
real(real128), save :: fb
real(real128), save :: fc
real(real128) m
real(real128), save :: machep
real(real128) p
real(real128) q
real(real128) r
real(real128) s
real(real128), save :: sa
real(real128), save :: sb
integer ( kind = 4 ) status
real(real128) t
real(real128) tol
real(real128) value
!
!  Input STATUS = 0.
!  Initialize, request F(A).
!
if ( status == 0 ) then

  machep = epsilon ( a )

  sa = a
  sb = b
  e = sb - sa
  d = e

  status = 1
  arg = a
  return
!
!  Input STATUS = 1.
!  Receive F(A), request F(B).
!
else if ( status == 1 ) then

  fa = value

  status = 2
  arg = sb
  return
!
!  Input STATUS = 2
!  Receive F(B).
!
else if ( status == 2 ) then

  fb = value

  if ( 0.0D+00 < fa * fb ) then
    status = -1
    return
  end if

  c = sa
  fc = fa

else

  fb = value

  if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
       ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
    c = sa
    fc = fa
    e = sb - sa
    d = e
  end if

end if
!
!  Compute the next point at which a function value is requested.
!
if ( abs ( fc ) < abs ( fb ) ) then

  sa = sb
  sb = c
  c = sa
  fa = fb
  fb = fc
  fc = fa

end if

tol = 2.0D+00 * machep * abs ( sb ) + t
m = 0.5D+00 * ( c - sb )

if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
  status = 0
  arg = sb
  return
end if

if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

  e = m
  d = e

else

  s = fb / fa

  if ( sa == c ) then

    p = 2.0D+00 * m * s
    q = 1.0D+00 - s

  else

    q = fa / fc
    r = fb / fc
    p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
    q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

  end if

  if ( 0.0D+00 < p ) then
    q = - q
  else
    p = - p
  end if

  s = e
  e = d

  if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
    p < abs ( 0.5D+00 * s * q ) ) then
    d = p / q
  else
    e = m
    d = e
  end if

end if

sa = sb
fa = fb

if ( tol < abs ( d ) ) then
  sb = sb + d
else if ( 0.0D+00 < m ) then
  sb = sb + tol
else
  sb = sb - tol
end if

arg = sb
status = status + 1

return
end


end module brentmod
