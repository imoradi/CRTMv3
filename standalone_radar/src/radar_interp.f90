!
! radar_interp
!
! Four-point interpolation used by the CRTM cloud optics lookup tables
! (two overlapping quadratic Lagrange polynomials blended linearly between
! the two central points; CRTM_Interpolation, ORDER=2). Reproducing it keeps
! the bulk-LUT path identical to CRTM. The routines return the weights of the
! four grid points and their derivatives with respect to the interpolation
! point, so that value = SUM(w*z) and d(value)/dx = SUM(dw*z).
! The derivative is obtained with a complex step, which is exact to machine
! precision because the weights are rational functions of x.
! Points outside the grid are clamped to the end points, where the
! derivative is zero.
!
MODULE radar_interp
  USE radar_kinds, ONLY: fp, ZERO, ONE, CSTEP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interp_weights
  PUBLIC :: interp_weights_c

CONTAINS

  ! Index of the first of the four points (CRTM Find_Random_Index)
  PURE INTEGER FUNCTION first_index(x, xi) RESULT(i1)
    REAL(fp), INTENT(IN) :: x(:), xi
    INTEGER :: k, n
    n = SIZE(x)
    DO k = 1, n
      IF (xi <= x(k)) EXIT
    END DO
    i1 = MIN(MAX(1, k-2), n-3)
  END FUNCTION first_index

  ! Weights in complex arithmetic (x4: the four grid points, xi complex)
  PURE SUBROUTINE weights_c(x4, xi, w)
    REAL(fp),    INTENT(IN)  :: x4(4)
    COMPLEX(fp), INTENT(IN)  :: xi
    COMPLEX(fp), INTENT(OUT) :: w(4)
    COMPLEX(fp) :: dxl(3), dxr(3), lpl(3), lpr(3), wl, wr
    REAL(fp)    :: dl(3), dr(3)
    dxl = xi - x4(1:3)
    dxr = xi - x4(2:4)
    dl = [x4(1)-x4(2), x4(1)-x4(3), x4(2)-x4(3)]
    dr = [x4(2)-x4(3), x4(2)-x4(4), x4(3)-x4(4)]
    lpl(1) = dxl(2)*dxl(3)/( dl(1)*dl(2))
    lpl(2) = dxl(1)*dxl(3)/(-dl(1)*dl(3))
    lpl(3) = dxl(1)*dxl(2)/( dl(2)*dl(3))
    lpr(1) = dxr(2)*dxr(3)/( dr(1)*dr(2))
    lpr(2) = dxr(1)*dxr(3)/(-dr(1)*dr(3))
    lpr(3) = dxr(1)*dxr(2)/( dr(2)*dr(3))
    IF (REAL(xi) < x4(2)) THEN
      wr = (0.0_fp, 0.0_fp)
    ELSE IF (REAL(xi) > x4(3)) THEN
      wr = (1.0_fp, 0.0_fp)
    ELSE
      wr = dxl(2)/(-dl(3))
    END IF
    wl = 1.0_fp - wr
    w(1) = wl*lpl(1)
    w(2) = wl*lpl(2) + wr*lpr(1)
    w(3) = wl*lpl(3) + wr*lpr(2)
    w(4) = wr*lpr(3)
  END SUBROUTINE weights_c

  !
  ! interp_weights(x, xi, i1, w, dw)
  !   x  : grid (monotonically increasing, at least 4 points)
  !   xi : interpolation point
  !   i1 : first grid index used (points i1:i1+3)
  !   w  : weights, dw: d(w)/d(xi)
  !
  PURE SUBROUTINE interp_weights(x, xi, i1, w, dw)
    REAL(fp), INTENT(IN)  :: x(:), xi
    INTEGER,  INTENT(OUT) :: i1
    REAL(fp), INTENT(OUT) :: w(4), dw(4)
    REAL(fp)    :: xc, h
    LOGICAL     :: clamped
    COMPLEX(fp) :: wc(4)
    INTEGER :: n
    n  = SIZE(x)
    xc = MIN(MAX(xi, x(1)), x(n))
    clamped = (xi < x(1)) .OR. (xi > x(n))
    i1 = first_index(x, xc)
    h  = CSTEP*MAX(ABS(xc), 1.0e-300_fp)
    CALL weights_c(x(i1:i1+3), CMPLX(xc, h, fp), wc)
    w  = REAL(wc)
    IF (clamped) THEN
      dw = ZERO
    ELSE
      dw = AIMAG(wc)/h
    END IF
  END SUBROUTINE interp_weights

  ! Complex version for use inside complex-step code (weights only)
  PURE SUBROUTINE interp_weights_c(x, xi, i1, w)
    REAL(fp),    INTENT(IN)  :: x(:)
    COMPLEX(fp), INTENT(IN)  :: xi
    INTEGER,     INTENT(OUT) :: i1
    COMPLEX(fp), INTENT(OUT) :: w(4)
    COMPLEX(fp) :: xc
    INTEGER :: n
    n = SIZE(x)
    IF (REAL(xi) < x(1)) THEN
      xc = CMPLX(x(1), 0.0_fp, fp)
    ELSE IF (REAL(xi) > x(n)) THEN
      xc = CMPLX(x(n), 0.0_fp, fp)
    ELSE
      xc = xi
    END IF
    i1 = first_index(x, REAL(xc))
    CALL weights_c(x(i1:i1+3), xc, w)
  END SUBROUTINE interp_weights_c

END MODULE radar_interp
