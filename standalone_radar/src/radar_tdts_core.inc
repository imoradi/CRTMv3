!
! CRTM_Active_MS_TDTS_core.inc
!
! Generic body of the time-dependent two-stream (TDTS) wide-angle multiple
! scattering calculation for a spaceborne radar with a Gaussian antenna
! pattern. This file is included twice by CRTM_Active_MS_TDTS.F90:
!   TDTS_RK = REAL(fp)    -> forward model
!   TDTS_RK = COMPLEX(fp) -> complex-step tangent-linear model
! All branch decisions use the real part of the variables so that both
! versions follow the same path.
!
! Fortran translation of tdts.c, wide_angle.c and wide_angle_regrid.c from
! the "multiscatter" package version 1.2.11:
!   Copyright (C) 2006-2011 University of Reading
!   Author: Robin Hogan <r.j.hogan@reading.ac.uk>
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!       http://www.apache.org/licenses/LICENSE-2.0
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
! Changes from the original: radar case only (no forward lobe, Gaussian
! receiver, one field of view); regridding works on layer edges rather than
! gate mid-points; written in Fortran with generic real/complex arithmetic.
!
! Reference: Hogan, R. J., and A. Battaglia, 2008: Fast lidar and radar
!   multiple-scattering models. Part II: Wide-angle scattering using the
!   time-dependent two-stream approximation. J. Atmos. Sci., 65, 3636-3651.
!

  !------------------------------------------------------------------------------
  ! Transfer coefficients Delta0..Delta5 (Hogan and Battaglia 2008, appendix),
  ! layer-mean one-way transmittance and transport mean free path.
  ! Arrays are indexed -1:m (one padding point each side).
  !------------------------------------------------------------------------------
  SUBROUTINE TDTS_DELTAS_NAME(m, drange, ext, ssa, g, &
                              delta0, delta1, delta2, delta3, delta4, delta5, &
                              transmittance, transport_mfp)
    INTEGER,  INTENT(IN)  :: m
    REAL(fp), INTENT(IN)  :: drange
    TDTS_RK,  INTENT(IN)  :: ext(0:m-1), ssa(0:m-1), g(0:m-1)
    TDTS_RK,  INTENT(OUT) :: delta0(-1:m), delta1(-1:m), delta2(-1:m)
    TDTS_RK,  INTENT(OUT) :: delta3(-1:m), delta4(-1:m), delta5(-1:m)
    TDTS_RK,  INTENT(OUT) :: transmittance(-1:m), transport_mfp(-1:m)
    INTEGER :: i
    TDTS_RK :: optical_depth, mfp_trans, mfp_abs, f_trans, f_abs, f
    TDTS_RK :: diffusion, diffusion1, diffusion2

    delta0 = ZERO; delta1 = ZERO; delta2 = ZERO
    delta3 = ZERO; delta4 = ZERO; delta5 = ZERO
    transmittance = ZERO; transport_mfp = ZERO
    optical_depth = ZERO
    transmittance(-1) = ONE
    DO i = 0, m-1
      IF (REAL(ext(i),fp) > ZERO) THEN
        ! Mean free paths for transport and absorption
        mfp_trans = ONE/(ext(i)*(ONE-ssa(i)*g(i)))
        mfp_abs   = ONE/(ext(i)*(ONE-ssa(i)))
        ! Factors
        f_trans = EXP(-drange/mfp_trans)
        f_abs   = EXP(-drange/mfp_abs)
        f = MS_MU1*(mfp_trans/drange - f_trans/(ONE-f_trans))
        ! Diffusion terms
        diffusion  = MS_MU1*ssa(i)*f_abs*SQRT(mfp_trans/(3.0_fp*drange))
        diffusion1 = diffusion*EXP(-3.7_fp*(mfp_trans/drange)**0.75_fp)
        diffusion2 = diffusion*EXP(-3.7_fp*mfp_trans/drange)
        ! Transfer terms
        delta0(i) = f_trans*(ONE-MS_MU1) + (f_abs-f_trans)*(POINT_5-f) - diffusion1
        delta1(i) = (f_abs-f_trans)*(ONE-f)*POINT_5 - diffusion2
        delta2(i) = MS_MU1*f_trans + (f_abs-f_trans)*f + diffusion1*POINT_5
        delta3(i) = (f_abs-f_trans)*f*0.25_fp + diffusion2*POINT_5
        delta4(i) = diffusion1*POINT_5
        delta5(i) = delta3(i)
        ! Remove small negative terms while conserving energy
        IF (REAL(delta0(i),fp) < ZERO) THEN
          delta2(i) = delta2(i) - POINT_5*delta0(i)
          delta4(i) = delta4(i) - POINT_5*delta0(i)
          delta0(i) = ZERO
        END IF
        IF (REAL(delta1(i),fp) < ZERO) THEN
          delta3(i) = delta3(i) - POINT_5*delta1(i)
          delta5(i) = delta5(i) - POINT_5*delta1(i)
          delta1(i) = ZERO
        END IF
        ! Transmittance and diffusivity
        transmittance(i) = EXP(-optical_depth)*(ONE-EXP(-ext(i)*drange))/(ext(i)*drange)
        optical_depth = optical_depth + ext(i)*drange
        transport_mfp(i) = ONE/((ONE-ssa(i)*g(i))*ext(i))
      ELSE
        delta0(i) = ONE - MS_MU1
        delta1(i) = ZERO
        delta2(i) = MS_MU1
        delta3(i) = ZERO
        delta4(i) = ZERO
        delta5(i) = ZERO
        transmittance(i) = EXP(-optical_depth)
        transport_mfp(i) = 1.0e6_fp*drange
      END IF
    END DO
    delta1(-1) = MS_MU1; delta1(m) = MS_MU1
    delta3(-1) = ZERO;   delta3(m) = ZERO
    delta4(-1) = ZERO;   delta5(m) = ZERO
    transmittance(m) = transmittance(m-1)
    transport_mfp(-1) = 1.0e6_fp*drange
    transport_mfp(m)  = 1.0e6_fp*drange
  END SUBROUTINE TDTS_DELTAS_NAME


  !------------------------------------------------------------------------------
  ! Lateral expansion of the photon distribution in one time step from the
  ! Ornstein-Furth formula (Hogan and Battaglia 2008, Eqs. 19, 22, 23); the
  ! fast version that is accurate to 1.5%.
  !------------------------------------------------------------------------------
  FUNCTION TDTS_EXPANSION_NAME(drange, lt, varw0, I, V) RESULT(expansion)
    REAL(fp), INTENT(IN) :: drange
    TDTS_RK,  INTENT(IN) :: lt, varw0, I, V
    TDTS_RK :: expansion
    TDTS_RK :: n, norm_varw, sqrt_var

    IF (REAL(I,fp) <= ZERO) THEN
      expansion = ZERO
      RETURN
    ELSE IF (REAL(V,fp) <= REAL(I*varw0,fp)) THEN
      n = ZERO
    ELSE
      norm_varw = (V/I-varw0)/(lt*lt)
      IF (REAL(norm_varw,fp) < ONE) THEN
        sqrt_var = SQRT(norm_varw*TWO*MS_ONE_OVER_D_FACTOR)
        n = sqrt_var*(ONE+0.2_fp*sqrt_var)
      ELSE
        n = norm_varw*MS_ONE_OVER_D_FACTOR + ONE - 0.24_fp/(norm_varw*norm_varw)
      END IF
    END IF
    expansion = MS_D_FACTOR*lt*(drange + lt*(EXP(-n-drange/lt)-EXP(-n)))
  END FUNCTION TDTS_EXPANSION_NAME


  !------------------------------------------------------------------------------
  ! Time-dependent two-stream integration on a regular grid of m gates of
  ! thickness drange (port of ms_tdts, radar Gaussian receiver). The result
  ! is added to bscat (apparent backscatter, m-1 sr-1).
  !   fp2: receiver footprint radius squared at each gate (m2)
  !------------------------------------------------------------------------------
  SUBROUTINE TDTS_CORE_NAME(m, drange, ext, ssa, g, src_power_in, src_power_out, &
                            src_width2, fp2, bscat)
    INTEGER,  INTENT(IN)     :: m
    REAL(fp), INTENT(IN)     :: drange
    TDTS_RK,  INTENT(IN)     :: ext(0:m-1), ssa(0:m-1), g(0:m-1)
    TDTS_RK,  INTENT(IN)     :: src_power_in(0:m-1), src_power_out(0:m-1), src_width2(0:m-1)
    REAL(fp), INTENT(IN)     :: fp2(0:m-1)
    TDTS_RK,  INTENT(IN OUT) :: bscat(0:m-1)
    TDTS_RK, DIMENSION(-1:m) :: Iin, Iout, Iin_next, Iout_next
    TDTS_RK, DIMENSION(-1:m) :: Vin, Vout, Vin_next, Vout_next
    TDTS_RK, DIMENSION(-1:m) :: delta0, delta1, delta2, delta3, delta4, delta5
    TDTS_RK, DIMENSION(-1:m) :: transmittance, transport_mfp
    TDTS_RK :: expansion, bscat_inc
    INTEGER :: nt, it, i, max_i, ibscat

    Iin = ZERO; Iout = ZERO; Iin_next = ZERO; Iout_next = ZERO
    Vin = ZERO; Vout = ZERO; Vin_next = ZERO; Vout_next = ZERO

    CALL TDTS_DELTAS_NAME(m, drange, ext, ssa, g, &
                          delta0, delta1, delta2, delta3, delta4, delta5, &
                          transmittance, transport_mfp)
    nt = 2*m
    Time_Loop: DO it = 0, nt-1
      ! Skip regions not yet reached by the beam, or that cannot scatter
      ! back to the receiver within the time period of interest
      max_i = MIN(it, m-1)
      IF (it > nt-m) max_i = nt-it

      Gate_Loop: DO i = 0, max_i
        ! Step the photon energy forward in time
        Iin_next(i)  = Iin(i)*delta0(i) + Iout(i)*delta1(i) + Iin(i+1)*delta2(i+1) &
                     + Iout(i-1)*delta3(i-1) + Iin(i-1)*delta4(i-1) + Iout(i+1)*delta5(i+1)
        Iout_next(i) = Iout(i)*delta0(i) + Iin(i)*delta1(i) + Iout(i-1)*delta2(i-1) &
                     + Iin(i+1)*delta3(i+1) + Iout(i+1)*delta4(i+1) + Iin(i-1)*delta5(i-1)
        ! Step the variance forward in time
        expansion = TDTS_EXPANSION_NAME(drange, transport_mfp(i), src_width2(i), Iin(i), Vin(i))
        Vin_next(i)  = Vin(i)*delta0(i) + Vout(i)*delta1(i) + Vin(i+1)*delta2(i+1) &
                     + Vout(i-1)*delta3(i-1) + Vin(i-1)*delta4(i-1) + Vout(i+1)*delta5(i+1) &
                     + Iin(i)*expansion
        expansion = TDTS_EXPANSION_NAME(drange, transport_mfp(i), src_width2(i), Iout(i), Vout(i))
        Vout_next(i) = Vout(i)*delta0(i) + Vin(i)*delta1(i) + Vout(i-1)*delta2(i-1) &
                     + Vin(i+1)*delta3(i+1) + Vout(i+1)*delta4(i+1) + Vin(i-1)*delta5(i-1) &
                     + Iout(i)*expansion
        ! Increment the apparent backscatter
        ibscat = (it+i)/2
        IF (ibscat >= m) CYCLE Gate_Loop
        IF (REAL(Vin_next(i),fp) > ZERO .AND. REAL(Vout_next(i),fp) > ZERO) THEN
          IF (REAL(g(i),fp) >= ONE/(3.0_fp*MS_MU1)) THEN
            bscat_inc = transmittance(i)*ONE_OVER_FOUR_PI &
                      * (Iin_next(i)*ssa(i)*ext(i)*TWO &
                         / (ONE + Vin_next(i)/(Iin_next(i)*fp2(i))))
          ELSE
            bscat_inc = transmittance(i)*ONE_OVER_FOUR_PI &
                      * (Iin_next(i)*ssa(i)*ext(i)*(ONE+3.0_fp*g(i)*MS_MU1) &
                         / (ONE + Vin_next(i)/(Iin_next(i)*fp2(i))) &
                       + Iout_next(i)*ssa(i)*ext(i)*(ONE-3.0_fp*g(i)*MS_MU1) &
                         / (ONE + Vout_next(i)/(Iout_next(i)*fp2(i))))
          END IF
          bscat(ibscat) = bscat(ibscat) + bscat_inc
        END IF
      END DO Gate_Loop

      ! Add the source power in the first half of the time simulated
      IF (it < m) THEN
        Iin_next(it)  = Iin_next(it)  + src_power_in(it)
        Iout_next(it) = Iout_next(it) + src_power_out(it)
        Vin_next(it)  = Vin_next(it)  + src_power_in(it)*src_width2(it)
        Vout_next(it) = Vout_next(it) + src_power_out(it)*src_width2(it)
      END IF

      ! Swap current and next time levels
      CALL TDTS_SWAP_NAME(Iin,  Iin_next)
      CALL TDTS_SWAP_NAME(Iout, Iout_next)
      CALL TDTS_SWAP_NAME(Vin,  Vin_next)
      CALL TDTS_SWAP_NAME(Vout, Vout_next)
    END DO Time_Loop
  END SUBROUTINE TDTS_CORE_NAME

  SUBROUTINE TDTS_SWAP_NAME(a, b)
    TDTS_RK, INTENT(IN OUT) :: a(:), b(:)
    TDTS_RK :: t(SIZE(a))
    t = a; a = b; b = t
  END SUBROUTINE TDTS_SWAP_NAME


  !------------------------------------------------------------------------------
  ! Wide-angle multiple scattering on an irregular layer grid.
  !
  ! Port of ms_wide_angle_regrid + ms_wide_angle (MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE)
  ! + ms_tdts for a radar with a Gaussian antenna of the same 1/e half-width
  ! for transmitter and receiver unless rho_fov differs.
  !
  ! Inputs (layers ordered away from the radar, k = 1..n):
  !   edge(0:n)  distance from the radar to the layer edges along the beam (m)
  !   ext(n)     extinction coefficient (m-1), gas + hydrometeors, unscaled
  !   ssa(n)     single-scatter albedo, unscaled
  !   g(n)       asymmetry factor
  !   rho_tr     transmitter 1/e half-width (rad)
  !   rho_fov    receiver 1/e half-width (rad)
  !   cbe        coherent backscatter enhancement factor (1 to 2)
  ! Output:
  !   bscat(n)   layer-mean apparent backscatter from wide-angle multiple
  !              scattering (m-1 sr-1); add to the single-scattering part.
  !------------------------------------------------------------------------------
  SUBROUTINE TDTS_PIPE_NAME(n, edge, ext, ssa, g, rho_tr, rho_fov, cbe, bscat)
    INTEGER,  INTENT(IN)  :: n
    REAL(fp), INTENT(IN)  :: edge(0:n)
    TDTS_RK,  INTENT(IN)  :: ext(n), ssa(n), g(n)
    REAL(fp), INTENT(IN)  :: rho_tr, rho_fov, cbe
    TDTS_RK,  INTENT(OUT) :: bscat(n)
    INTEGER  :: first_gate, i, j, m
    REAL(fp) :: footprint, drange_new, dri, scat_od, a, b, ov, dist
    TDTS_RK  :: mfp, od2way, ss_multiplier, od, optical_depth, src_power
    TDTS_RK, ALLOCATABLE :: ext_new(:), ssa_new(:), g_new(:), bscat_new(:)
    TDTS_RK, ALLOCATABLE :: src_in(:), src_out(:), src_w2(:)
    REAL(fp), ALLOCATABLE :: fp2(:)

    bscat = ZERO
    IF (n < 3) RETURN

    ! Find the first gate where the mean free path is small compared with
    ! the footprint, and the two-way optical depth down to it
    first_gate = -1
    od2way = ZERO
    DO i = 1, n-1
      IF (REAL(ext(i),fp) > ZERO .AND. REAL(ssa(i),fp) > ZERO) THEN
        footprint = ABS(edge(i-1)+edge(i))*POINT_5*rho_fov
        mfp = ONE/(ext(i)*ssa(i)*(ONE-g(i)))
        IF (REAL(mfp,fp) < footprint*MS_MFP_FOOTPRINT_THRESHOLD) THEN
          first_gate = i
          EXIT
        END IF
      END IF
      od2way = od2way + TWO*ext(i)*(edge(i)-edge(i-1))
    END DO
    IF (first_gate < 0 .OR. first_gate >= n-1) RETURN
    ss_multiplier = (ONE + (rho_tr*rho_tr)/(rho_fov*rho_fov))*EXP(-od2way)

    ! Regular grid spacing: the thinnest scattering layer within the first
    ! MS_OD_DRANGE_THRESHOLD scattering optical depths of the cloud
    drange_new = edge(first_gate) - edge(first_gate-1)
    scat_od = drange_new*REAL(ext(first_gate)*ssa(first_gate),fp)
    DO i = first_gate+1, n
      dri = edge(i) - edge(i-1)
      IF (drange_new > dri .AND. REAL(ext(i),fp) > ZERO .AND. REAL(ssa(i),fp) > ZERO) &
        drange_new = dri
      scat_od = scat_od + drange_new*REAL(ext(i)*ssa(i),fp)
      IF (scat_od > MS_OD_DRANGE_THRESHOLD) EXIT
    END DO
    m = NINT((edge(n)-edge(first_gate-1))/drange_new)
    IF (m > TDTS_MAX_GATES) THEN
      m = TDTS_MAX_GATES
      drange_new = (edge(n)-edge(first_gate-1))/REAL(m,fp)
    END IF
    IF (m < 3) RETURN

    ALLOCATE(ext_new(0:m-1), ssa_new(0:m-1), g_new(0:m-1), bscat_new(0:m-1), &
             src_in(0:m-1), src_out(0:m-1), src_w2(0:m-1), fp2(0:m-1))
    ext_new = ZERO; ssa_new = ZERO; g_new = ZERO; bscat_new = ZERO

    ! Interpolate the optical properties on to the regular grid, conserving
    ! optical depth, scattering optical depth and its g-weighted integral
    DO j = 0, m-1
      a = edge(first_gate-1) + drange_new*REAL(j,fp)
      b = a + drange_new
      DO i = first_gate, n
        ov = MIN(b, edge(i)) - MAX(a, edge(i-1))
        IF (ov <= ZERO) CYCLE
        ext_new(j) = ext_new(j) + ov*ext(i)
        ssa_new(j) = ssa_new(j) + ov*ext(i)*ssa(i)
        g_new(j)   = g_new(j)   + ov*ext(i)*ssa(i)*g(i)
      END DO
      IF (REAL(ssa_new(j),fp) > ZERO) THEN
        g_new(j) = g_new(j)/ssa_new(j)
      ELSE
        g_new(j) = 0.65_fp
      END IF
      IF (REAL(ext_new(j),fp) > ZERO) THEN
        ssa_new(j) = ssa_new(j)/ext_new(j)
      ELSE
        ssa_new(j) = ONE
      END IF
      ext_new(j) = ext_new(j)/drange_new
    END DO

    ! Source terms for the two-stream calculation (no forward lobe)
    optical_depth = ZERO
    DO j = 0, m-1
      dist = edge(first_gate-1) + drange_new*(REAL(j,fp)+POINT_5)
      od = ext_new(j)*drange_new
      src_power = cbe*ss_multiplier*ssa_new(j)*EXP(-optical_depth)*(ONE-EXP(-od))
      IF (REAL(od,fp) > ZERO) optical_depth = optical_depth + od
      IF (REAL(g_new(j),fp) >= ONE/(3.0_fp*MS_MU1)) THEN
        src_in(j)  = ZERO
        src_out(j) = src_power
      ELSE IF (REAL(g_new(j),fp) <= -ONE/(3.0_fp*MS_MU1)) THEN
        src_in(j)  = src_power
        src_out(j) = ZERO
      ELSE
        src_in(j)  = src_power*POINT_5*(ONE-3.0_fp*g_new(j)*MS_MU1)
        src_out(j) = src_power*POINT_5*(ONE+3.0_fp*g_new(j)*MS_MU1)
      END IF
      src_w2(j) = rho_tr*rho_tr*dist*dist
      fp2(j)    = rho_fov*rho_fov*dist*dist
    END DO

    CALL TDTS_CORE_NAME(m, drange_new, ext_new, ssa_new, g_new, src_in, src_out, &
                        src_w2, fp2, bscat_new)

    ! Map the backscatter back on to the layers, conserving its integral
    DO i = first_gate, n
      DO j = 0, m-1
        a = edge(first_gate-1) + drange_new*REAL(j,fp)
        b = a + drange_new
        ov = MIN(b, edge(i)) - MAX(a, edge(i-1))
        IF (ov <= ZERO) CYCLE
        bscat(i) = bscat(i) + ov*bscat_new(j)
      END DO
      bscat(i) = bscat(i)/(edge(i)-edge(i-1))
    END DO

    DEALLOCATE(ext_new, ssa_new, g_new, bscat_new, src_in, src_out, src_w2, fp2)
  END SUBROUTINE TDTS_PIPE_NAME
