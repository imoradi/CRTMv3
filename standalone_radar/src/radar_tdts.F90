!
! radar_tdts (copied from CRTM_Active_MS_TDTS, CRTMv3 feature/multiple-scattering)
!
! Wide-angle multiple scattering for spaceborne radar using the
! time-dependent two-stream (TDTS) method of Hogan and Battaglia (2008).
!
! The forward model and its complex-step tangent-linear share the generic
! code in radar_tdts_core.inc. The adjoint is hand coded; the
! adjoint of the two-stream time stepping follows tdts_AD.c of the
! "multiscatter" package.
!
! Fortran translation of parts of the "multiscatter" package version 1.2.11:
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
!
! Reference: Hogan, R. J., and A. Battaglia, 2008: Fast lidar and radar
!   multiple-scattering models. Part II: Wide-angle scattering using the
!   time-dependent two-stream approximation. J. Atmos. Sci., 65, 3636-3651.
!
! CREATION HISTORY:
!       Written by:     Isaac Moradi - NASA GSFC, November 2026
!

MODULE radar_tdts

  USE radar_kinds, ONLY: fp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: TDTS_Wide_Angle
  PUBLIC :: TDTS_Wide_Angle_TL
  PUBLIC :: TDTS_Wide_Angle_AD

  REAL(fp), PARAMETER :: ZERO = 0.0_fp, ONE = 1.0_fp, TWO = 2.0_fp, POINT_5 = 0.5_fp
  REAL(fp), PARAMETER :: ONE_OVER_FOUR_PI = 0.0795774715459477_fp
  ! Cosine of the two-stream angle
  REAL(fp), PARAMETER :: MS_MU1 = 0.5_fp
  ! Factor in the Ornstein-Furth formula (4/3) and its reciprocal
  REAL(fp), PARAMETER :: MS_D_FACTOR = 1.33333333_fp
  REAL(fp), PARAMETER :: MS_ONE_OVER_D_FACTOR = 0.75_fp
  ! Wide-angle scattering starts where the scattering mean free path is less
  ! than this many times the footprint radius
  REAL(fp), PARAMETER :: MS_MFP_FOOTPRINT_THRESHOLD = 50.0_fp
  ! Scattering optical depth over which the finest layer sets the grid
  REAL(fp), PARAMETER :: MS_OD_DRANGE_THRESHOLD = 5.0_fp
  ! Upper limit on the number of points of the regular grid
  INTEGER,  PARAMETER :: TDTS_MAX_GATES = 500
  ! Complex step for the tangent-linear
  REAL(fp), PARAMETER :: CSTEP = 1.0e-20_fp

CONTAINS

! ---------------------------------------------------------------------------
! Generic code, real version (forward model)
! ---------------------------------------------------------------------------
#define TDTS_RK REAL(fp)
#define TDTS_DELTAS_NAME tdts_deltas_r
#define TDTS_EXPANSION_NAME tdts_expansion_r
#define TDTS_CORE_NAME tdts_core_r
#define TDTS_SWAP_NAME tdts_swap_r
#define TDTS_PIPE_NAME tdts_pipe_r
#include "radar_tdts_core.inc"
#undef TDTS_RK
#undef TDTS_DELTAS_NAME
#undef TDTS_EXPANSION_NAME
#undef TDTS_CORE_NAME
#undef TDTS_SWAP_NAME
#undef TDTS_PIPE_NAME

! ---------------------------------------------------------------------------
! Generic code, complex version (complex-step tangent-linear)
! ---------------------------------------------------------------------------
#define TDTS_RK COMPLEX(fp)
#define TDTS_DELTAS_NAME tdts_deltas_c
#define TDTS_EXPANSION_NAME tdts_expansion_c
#define TDTS_CORE_NAME tdts_core_c
#define TDTS_SWAP_NAME tdts_swap_c
#define TDTS_PIPE_NAME tdts_pipe_c
#include "radar_tdts_core.inc"
#undef TDTS_RK
#undef TDTS_DELTAS_NAME
#undef TDTS_EXPANSION_NAME
#undef TDTS_CORE_NAME
#undef TDTS_SWAP_NAME
#undef TDTS_PIPE_NAME


!--------------------------------------------------------------------------------
!
! NAME:
!       TDTS_Wide_Angle
!
! PURPOSE:
!       Layer-mean apparent backscatter (m-1 sr-1) from wide-angle multiple
!       scattering, to be added to the attenuated single-scattering return.
!       See TDTS_PIPE_NAME in radar_tdts_core.inc for the arguments.
!
!--------------------------------------------------------------------------------
  SUBROUTINE TDTS_Wide_Angle(n, edge, ext, ssa, g, rho_tr, rho_fov, cbe, bscat)
    INTEGER,  INTENT(IN)  :: n
    REAL(fp), INTENT(IN)  :: edge(0:n), ext(n), ssa(n), g(n)
    REAL(fp), INTENT(IN)  :: rho_tr, rho_fov, cbe
    REAL(fp), INTENT(OUT) :: bscat(n)
    CALL tdts_pipe_r(n, edge, ext, ssa, g, rho_tr, rho_fov, cbe, bscat)
  END SUBROUTINE TDTS_Wide_Angle


!--------------------------------------------------------------------------------
!
! NAME:
!       TDTS_Wide_Angle_TL
!
! PURPOSE:
!       Tangent-linear of TDTS_Wide_Angle by the complex-step method, which
!       gives the derivative to machine precision with the same code path
!       as the forward model.
!
!--------------------------------------------------------------------------------
  SUBROUTINE TDTS_Wide_Angle_TL(n, edge, ext, ssa, g, rho_tr, rho_fov, cbe, &
                                ext_TL, ssa_TL, g_TL, bscat, bscat_TL)
    INTEGER,  INTENT(IN)  :: n
    REAL(fp), INTENT(IN)  :: edge(0:n), ext(n), ssa(n), g(n)
    REAL(fp), INTENT(IN)  :: rho_tr, rho_fov, cbe
    REAL(fp), INTENT(IN)  :: ext_TL(n), ssa_TL(n), g_TL(n)
    REAL(fp), INTENT(OUT) :: bscat(n), bscat_TL(n)
    COMPLEX(fp) :: ext_c(n), ssa_c(n), g_c(n), bscat_c(n)

    ext_c = CMPLX(ext, CSTEP*ext_TL, fp)
    ssa_c = CMPLX(ssa, CSTEP*ssa_TL, fp)
    g_c   = CMPLX(g,   CSTEP*g_TL,   fp)
    CALL tdts_pipe_c(n, edge, ext_c, ssa_c, g_c, rho_tr, rho_fov, cbe, bscat_c)
    bscat    = REAL(bscat_c, fp)
    bscat_TL = AIMAG(bscat_c)/CSTEP
  END SUBROUTINE TDTS_Wide_Angle_TL


!--------------------------------------------------------------------------------
!
! NAME:
!       TDTS_Wide_Angle_AD
!
! PURPOSE:
!       Adjoint of TDTS_Wide_Angle. The adjoints of ext, ssa and g are
!       accumulated (not initialised) in ext_AD, ssa_AD and g_AD.
!
!--------------------------------------------------------------------------------
  SUBROUTINE TDTS_Wide_Angle_AD(n, edge, ext, ssa, g, rho_tr, rho_fov, cbe, &
                                bscat_AD, ext_AD, ssa_AD, g_AD)
    INTEGER,  INTENT(IN)     :: n
    REAL(fp), INTENT(IN)     :: edge(0:n), ext(n), ssa(n), g(n)
    REAL(fp), INTENT(IN)     :: rho_tr, rho_fov, cbe
    REAL(fp), INTENT(IN)     :: bscat_AD(n)
    REAL(fp), INTENT(IN OUT) :: ext_AD(n), ssa_AD(n), g_AD(n)

    INTEGER  :: first_gate, i, j, m
    REAL(fp) :: footprint, drange_new, dri, scat_od, a, b, ov, dist
    REAL(fp) :: mfp, od2way, ss_multiplier, od2way_AD, ss_multiplier_AD
    REAL(fp) :: c1, e_od, e_OD_cum, sp_AD, OD_AD_sum
    REAL(fp), ALLOCATABLE :: Ej(:), Sj(:), Gj(:), ext_new(:), ssa_new(:), g_new(:)
    REAL(fp), ALLOCATABLE :: od(:), OD_cum(:), src_power(:)
    REAL(fp), ALLOCATABLE :: src_in(:), src_out(:), src_w2(:), fp2(:)
    REAL(fp), ALLOCATABLE :: Ej_AD(:), Sj_AD(:), Gj_AD(:)
    REAL(fp), ALLOCATABLE :: ext_new_AD(:), ssa_new_AD(:), g_new_AD(:), bscat_new_AD(:)
    REAL(fp), ALLOCATABLE :: src_in_AD(:), src_out_AD(:), od_AD(:), OD_cum_AD(:)

    IF (n < 3) RETURN

    ! ---------------- Forward recomputation ----------------
    first_gate = -1
    od2way = ZERO
    DO i = 1, n-1
      IF (ext(i) > ZERO .AND. ssa(i) > ZERO) THEN
        footprint = ABS(edge(i-1)+edge(i))*POINT_5*rho_fov
        mfp = ONE/(ext(i)*ssa(i)*(ONE-g(i)))
        IF (mfp < footprint*MS_MFP_FOOTPRINT_THRESHOLD) THEN
          first_gate = i
          EXIT
        END IF
      END IF
      od2way = od2way + TWO*ext(i)*(edge(i)-edge(i-1))
    END DO
    IF (first_gate < 0 .OR. first_gate >= n-1) RETURN
    ss_multiplier = (ONE + (rho_tr*rho_tr)/(rho_fov*rho_fov))*EXP(-od2way)

    drange_new = edge(first_gate) - edge(first_gate-1)
    scat_od = drange_new*ext(first_gate)*ssa(first_gate)
    DO i = first_gate+1, n
      dri = edge(i) - edge(i-1)
      IF (drange_new > dri .AND. ext(i) > ZERO .AND. ssa(i) > ZERO) drange_new = dri
      scat_od = scat_od + drange_new*ext(i)*ssa(i)
      IF (scat_od > MS_OD_DRANGE_THRESHOLD) EXIT
    END DO
    m = NINT((edge(n)-edge(first_gate-1))/drange_new)
    IF (m > TDTS_MAX_GATES) THEN
      m = TDTS_MAX_GATES
      drange_new = (edge(n)-edge(first_gate-1))/REAL(m,fp)
    END IF
    IF (m < 3) RETURN

    ALLOCATE(Ej(0:m-1), Sj(0:m-1), Gj(0:m-1), ext_new(0:m-1), ssa_new(0:m-1), g_new(0:m-1), &
             od(0:m-1), OD_cum(0:m-1), src_power(0:m-1), src_in(0:m-1), src_out(0:m-1), &
             src_w2(0:m-1), fp2(0:m-1), Ej_AD(0:m-1), Sj_AD(0:m-1), Gj_AD(0:m-1), &
             ext_new_AD(0:m-1), ssa_new_AD(0:m-1), g_new_AD(0:m-1), bscat_new_AD(0:m-1), &
             src_in_AD(0:m-1), src_out_AD(0:m-1), od_AD(0:m-1), OD_cum_AD(0:m-1))
    Ej = ZERO; Sj = ZERO; Gj = ZERO
    DO j = 0, m-1
      a = edge(first_gate-1) + drange_new*REAL(j,fp)
      b = a + drange_new
      DO i = first_gate, n
        ov = MIN(b, edge(i)) - MAX(a, edge(i-1))
        IF (ov <= ZERO) CYCLE
        Ej(j)  = Ej(j)  + ov*ext(i)
        Sj(j)  = Sj(j)  + ov*ext(i)*ssa(i)
        Gj(j) = Gj(j) + ov*ext(i)*ssa(i)*g(i)
      END DO
      IF (Sj(j) > ZERO) THEN
        g_new(j) = Gj(j)/Sj(j)
      ELSE
        g_new(j) = 0.65_fp
      END IF
      IF (Ej(j) > ZERO) THEN
        ssa_new(j) = Sj(j)/Ej(j)
      ELSE
        ssa_new(j) = ONE
      END IF
      ext_new(j) = Ej(j)/drange_new
    END DO

    c1 = cbe*ss_multiplier
    e_OD_cum = ZERO
    DO j = 0, m-1
      dist = edge(first_gate-1) + drange_new*(REAL(j,fp)+POINT_5)
      od(j) = ext_new(j)*drange_new
      OD_cum(j) = e_OD_cum
      src_power(j) = c1*ssa_new(j)*EXP(-OD_cum(j))*(ONE-EXP(-od(j)))
      IF (od(j) > ZERO) e_OD_cum = e_OD_cum + od(j)
      IF (g_new(j) >= ONE/(3.0_fp*MS_MU1)) THEN
        src_in(j) = ZERO; src_out(j) = src_power(j)
      ELSE IF (g_new(j) <= -ONE/(3.0_fp*MS_MU1)) THEN
        src_in(j) = src_power(j); src_out(j) = ZERO
      ELSE
        src_in(j)  = src_power(j)*POINT_5*(ONE-3.0_fp*g_new(j)*MS_MU1)
        src_out(j) = src_power(j)*POINT_5*(ONE+3.0_fp*g_new(j)*MS_MU1)
      END IF
      src_w2(j) = rho_tr*rho_tr*dist*dist
      fp2(j)    = rho_fov*rho_fov*dist*dist
    END DO

    ! ---------------- Adjoint ----------------
    ! Back-mapping on to the layers
    bscat_new_AD = ZERO
    DO i = first_gate, n
      DO j = 0, m-1
        a = edge(first_gate-1) + drange_new*REAL(j,fp)
        b = a + drange_new
        ov = MIN(b, edge(i)) - MAX(a, edge(i-1))
        IF (ov <= ZERO) CYCLE
        bscat_new_AD(j) = bscat_new_AD(j) + ov*bscat_AD(i)/(edge(i)-edge(i-1))
      END DO
    END DO

    ! Two-stream time stepping
    ext_new_AD = ZERO; ssa_new_AD = ZERO; g_new_AD = ZERO
    src_in_AD = ZERO; src_out_AD = ZERO
    CALL tdts_core_AD(m, drange_new, ext_new, ssa_new, g_new, src_in, src_out, src_w2, fp2, &
                      bscat_new_AD, ext_new_AD, ssa_new_AD, g_new_AD, src_in_AD, src_out_AD)

    ! Source terms
    od_AD = ZERO; OD_cum_AD = ZERO; ss_multiplier_AD = ZERO
    DO j = m-1, 0, -1
      IF (g_new(j) >= ONE/(3.0_fp*MS_MU1)) THEN
        sp_AD = src_out_AD(j)
      ELSE IF (g_new(j) <= -ONE/(3.0_fp*MS_MU1)) THEN
        sp_AD = src_in_AD(j)
      ELSE
        sp_AD = src_in_AD(j)*POINT_5*(ONE-3.0_fp*g_new(j)*MS_MU1) &
              + src_out_AD(j)*POINT_5*(ONE+3.0_fp*g_new(j)*MS_MU1)
        g_new_AD(j) = g_new_AD(j) + src_power(j)*POINT_5*3.0_fp*MS_MU1*(src_out_AD(j)-src_in_AD(j))
      END IF
      e_od = EXP(-od(j))
      ssa_new_AD(j)    = ssa_new_AD(j)    + sp_AD*c1*EXP(-OD_cum(j))*(ONE-e_od)
      ss_multiplier_AD = ss_multiplier_AD + sp_AD*cbe*ssa_new(j)*EXP(-OD_cum(j))*(ONE-e_od)
      od_AD(j)         = od_AD(j)         + sp_AD*c1*ssa_new(j)*EXP(-OD_cum(j))*e_od
      OD_cum_AD(j)     = OD_cum_AD(j)     - sp_AD*src_power(j)
    END DO
    ! OD_cum(j) = sum of od(l) for l < j (od >= 0)
    OD_AD_sum = ZERO
    DO j = m-1, 0, -1
      IF (od(j) > ZERO) od_AD(j) = od_AD(j) + OD_AD_sum
      OD_AD_sum = OD_AD_sum + OD_cum_AD(j)
    END DO
    ext_new_AD = ext_new_AD + od_AD*drange_new

    ! Regridding
    DO j = 0, m-1
      Ej_AD(j) = ext_new_AD(j)/drange_new
      Sj_AD(j) = ZERO
      Gj_AD(j) = ZERO
      IF (Ej(j) > ZERO) THEN
        Sj_AD(j) = ssa_new_AD(j)/Ej(j)
        Ej_AD(j) = Ej_AD(j) - ssa_new_AD(j)*Sj(j)/(Ej(j)*Ej(j))
      END IF
      IF (Sj(j) > ZERO) THEN
        Gj_AD(j) = g_new_AD(j)/Sj(j)
        Sj_AD(j) = Sj_AD(j) - g_new_AD(j)*Gj(j)/(Sj(j)*Sj(j))
      END IF
      a = edge(first_gate-1) + drange_new*REAL(j,fp)
      b = a + drange_new
      DO i = first_gate, n
        ov = MIN(b, edge(i)) - MAX(a, edge(i-1))
        IF (ov <= ZERO) CYCLE
        ext_AD(i) = ext_AD(i) + ov*(Ej_AD(j) + ssa(i)*Sj_AD(j) + ssa(i)*g(i)*Gj_AD(j))
        ssa_AD(i) = ssa_AD(i) + ov*ext(i)*(Sj_AD(j) + g(i)*Gj_AD(j))
        g_AD(i)   = g_AD(i)   + ov*ext(i)*ssa(i)*Gj_AD(j)
      END DO
    END DO

    ! Two-way optical depth down to the first wide-angle gate
    od2way_AD = -ss_multiplier_AD*ss_multiplier
    DO i = 1, first_gate-1
      ext_AD(i) = ext_AD(i) + od2way_AD*TWO*(edge(i)-edge(i-1))
    END DO

    DEALLOCATE(Ej, Sj, Gj, ext_new, ssa_new, g_new, od, OD_cum, src_power, src_in, src_out, &
               src_w2, fp2, Ej_AD, Sj_AD, Gj_AD, ext_new_AD, ssa_new_AD, g_new_AD, bscat_new_AD, &
               src_in_AD, src_out_AD, od_AD, OD_cum_AD)
  END SUBROUTINE TDTS_Wide_Angle_AD


  !------------------------------------------------------------------------------
  ! Adjoint of the transfer coefficients (tdts_deltas_r), including the
  ! correction for negative Delta0/Delta1.
  !------------------------------------------------------------------------------
  SUBROUTINE tdts_deltas_AD(m, drange, ext, ssa, g, transmittance, &
                            delta0_AD, delta1_AD, delta2_AD, delta3_AD, delta4_AD, delta5_AD, &
                            transmittance_AD, transport_mfp_AD, ext_AD, ssa_AD, g_AD)
    INTEGER,  INTENT(IN)     :: m
    REAL(fp), INTENT(IN)     :: drange
    REAL(fp), INTENT(IN)     :: ext(0:m-1), ssa(0:m-1), g(0:m-1), transmittance(-1:m)
    REAL(fp), INTENT(IN OUT) :: delta0_AD(-1:m), delta1_AD(-1:m), delta2_AD(-1:m)
    REAL(fp), INTENT(IN OUT) :: delta3_AD(-1:m), delta4_AD(-1:m), delta5_AD(-1:m)
    REAL(fp), INTENT(IN OUT) :: transmittance_AD(-1:m), transport_mfp_AD(-1:m)
    REAL(fp), INTENT(IN OUT) :: ext_AD(0:m-1), ssa_AD(0:m-1), g_AD(0:m-1)
    INTEGER  :: i
    REAL(fp) :: mt, ft, fa, f, sq, dif, d1e, d2e, dif1, dif2, d0r, d1r, u, eu
    REAL(fp) :: ft_AD, fa_AD, f_AD, dif_AD, dif1_AD, dif2_AD, mt_AD, d3r_AD, d0r_AD, d1r_AD
    REAL(fp) :: OD_AD_sum

    transmittance_AD(m-1) = transmittance_AD(m-1) + transmittance_AD(m)
    transmittance_AD(m) = ZERO
    OD_AD_sum = ZERO
    DO i = m-1, 0, -1
      IF (ext(i) > ZERO) THEN
        ! Forward quantities
        mt  = ONE/(ext(i)*(ONE-ssa(i)*g(i)))
        ft  = EXP(-drange/mt)
        fa  = EXP(-drange*ext(i)*(ONE-ssa(i)))
        f   = MS_MU1*(mt/drange - ft/(ONE-ft))
        sq  = SQRT(mt/(3.0_fp*drange))
        dif = MS_MU1*ssa(i)*fa*sq
        d1e = EXP(-3.7_fp*(mt/drange)**0.75_fp)
        d2e = EXP(-3.7_fp*mt/drange)
        dif1 = dif*d1e
        dif2 = dif*d2e
        d0r = ft*(ONE-MS_MU1) + (fa-ft)*(POINT_5-f) - dif1
        d1r = (fa-ft)*(ONE-f)*POINT_5 - dif2

        ! Transmittance: T(i) = exp(-OD(i))*(1-exp(-u))/u, u = ext*drange
        u  = ext(i)*drange
        eu = EXP(-u)
        ext_AD(i) = ext_AD(i) + transmittance_AD(i)*(transmittance(i)/(ONE-eu))*(eu - (ONE-eu)/u)*drange
        ! (The above uses exp(-OD)=T*u/(1-eu))
        ext_AD(i) = ext_AD(i) + OD_AD_sum*drange
        OD_AD_sum = OD_AD_sum - transmittance(i)*transmittance_AD(i)
        transmittance_AD(i) = ZERO

        ! Transport mean free path (same expression as mt)
        mt_AD = transport_mfp_AD(i)
        transport_mfp_AD(i) = ZERO

        ! Corrections for negative terms
        d0r_AD = delta0_AD(i)
        d1r_AD = delta1_AD(i)
        IF (d1r < ZERO) d1r_AD = -POINT_5*(delta3_AD(i) + delta5_AD(i))
        IF (d0r < ZERO) d0r_AD = -POINT_5*(delta2_AD(i) + delta4_AD(i))
        ! delta5 = delta3 before correction
        d3r_AD = delta3_AD(i) + delta5_AD(i)

        ft_AD = (POINT_5+f-MS_MU1)*d0r_AD + POINT_5*(f-ONE)*d1r_AD &
              + (MS_MU1-f)*delta2_AD(i) - 0.25_fp*f*d3r_AD
        fa_AD = (POINT_5-f)*d0r_AD + POINT_5*(ONE-f)*d1r_AD &
              + f*delta2_AD(i) + 0.25_fp*f*d3r_AD
        f_AD  = (ft-fa)*d0r_AD + POINT_5*(ft-fa)*d1r_AD &
              + (fa-ft)*delta2_AD(i) + 0.25_fp*(fa-ft)*d3r_AD
        dif1_AD = -d0r_AD + POINT_5*delta2_AD(i) + POINT_5*delta4_AD(i)
        dif2_AD = -d1r_AD + POINT_5*d3r_AD
        delta0_AD(i) = ZERO; delta1_AD(i) = ZERO; delta2_AD(i) = ZERO
        delta3_AD(i) = ZERO; delta4_AD(i) = ZERO; delta5_AD(i) = ZERO

        ! Diffusion terms
        dif_AD = d1e*dif1_AD + d2e*dif2_AD
        mt_AD = mt_AD - 3.7_fp*dif2*dif2_AD/drange &
              - 3.7_fp*0.75_fp*dif1*mt**(-0.25_fp)*drange**(-0.75_fp)*dif1_AD
        ssa_AD(i) = ssa_AD(i) + MS_MU1*sq*fa*dif_AD
        fa_AD = fa_AD + MS_MU1*sq*ssa(i)*dif_AD
        mt_AD = mt_AD + POINT_5*dif*dif_AD/mt

        ! Factors
        mt_AD = mt_AD + MS_MU1*f_AD/drange
        ft_AD = ft_AD - MS_MU1*f_AD/((ONE-ft)*(ONE-ft))
        ext_AD(i) = ext_AD(i) - fa*drange*(ONE-ssa(i))*fa_AD
        ssa_AD(i) = ssa_AD(i) + fa*drange*ext(i)*fa_AD
        mt_AD = mt_AD + ft*drange*ft_AD/(mt*mt)

        ! Mean free path for transport
        ext_AD(i) = ext_AD(i) - mt*mt_AD/ext(i)
        ssa_AD(i) = ssa_AD(i) + mt*mt_AD*g(i)/(ONE-ssa(i)*g(i))
        g_AD(i)   = g_AD(i)   + mt*mt_AD*ssa(i)/(ONE-ssa(i)*g(i))
      ELSE
        delta0_AD(i) = ZERO; delta1_AD(i) = ZERO; delta2_AD(i) = ZERO
        delta3_AD(i) = ZERO; delta4_AD(i) = ZERO; delta5_AD(i) = ZERO
        ! T(i) = exp(-OD(i))
        OD_AD_sum = OD_AD_sum - transmittance(i)*transmittance_AD(i)
        transmittance_AD(i) = ZERO
        transport_mfp_AD(i) = ZERO
      END IF
    END DO
  END SUBROUTINE tdts_deltas_AD


  !------------------------------------------------------------------------------
  ! Adjoint of tdts_expansion_r with respect to lt, I and V
  !------------------------------------------------------------------------------
  SUBROUTINE tdts_expansion_AD(drange, lt, varw0, I, V, expansion_AD, lt_AD, I_AD, V_AD)
    REAL(fp), INTENT(IN)     :: drange, lt, varw0, I, V, expansion_AD
    REAL(fp), INTENT(IN OUT) :: lt_AD, I_AD, V_AD
    REAL(fp) :: n, norm_varw, sqrt_var, exp1, exp2, n_AD, norm_varw_AD

    IF (I <= ZERO) RETURN
    norm_varw = ZERO
    sqrt_var  = ZERO
    IF (V <= I*varw0) THEN
      n = ZERO
    ELSE
      norm_varw = (V/I-varw0)/(lt*lt)
      IF (norm_varw < ONE) THEN
        sqrt_var = SQRT(norm_varw*TWO*MS_ONE_OVER_D_FACTOR)
        n = sqrt_var*(ONE+0.2_fp*sqrt_var)
      ELSE
        n = norm_varw*MS_ONE_OVER_D_FACTOR + ONE - 0.24_fp/(norm_varw*norm_varw)
      END IF
    END IF
    exp1 = EXP(-n-drange/lt)
    exp2 = EXP(-n)
    lt_AD = lt_AD + MS_D_FACTOR*(drange + TWO*lt*(exp1-exp2) + drange*exp1)*expansion_AD
    n_AD = MS_D_FACTOR*lt*lt*(exp2-exp1)*expansion_AD
    IF (V <= I*varw0) RETURN
    IF (norm_varw < ONE) THEN
      norm_varw_AD = SQRT(MS_ONE_OVER_D_FACTOR*POINT_5/norm_varw)*(ONE+0.4_fp*sqrt_var)*n_AD
    ELSE
      norm_varw_AD = n_AD*(MS_ONE_OVER_D_FACTOR + (TWO*0.24_fp)/(norm_varw*norm_varw*norm_varw))
    END IF
    V_AD  = V_AD  + norm_varw_AD/(I*lt*lt)
    I_AD  = I_AD  - norm_varw_AD*V/(I*I*lt*lt)
    lt_AD = lt_AD - norm_varw_AD*TWO*norm_varw/lt
  END SUBROUTINE tdts_expansion_AD


  !------------------------------------------------------------------------------
  ! Adjoint of tdts_core_r (radar, Gaussian receiver). Follows ms_tdts_AD.
  ! The adjoint of src_width2 is not needed (it depends on geometry only).
  !------------------------------------------------------------------------------
  SUBROUTINE tdts_core_AD(m, drange, ext, ssa, g, src_power_in, src_power_out, src_width2, fp2, &
                          bscat_AD, ext_AD, ssa_AD, g_AD, src_power_in_AD, src_power_out_AD)
    INTEGER,  INTENT(IN)     :: m
    REAL(fp), INTENT(IN)     :: drange
    REAL(fp), INTENT(IN)     :: ext(0:m-1), ssa(0:m-1), g(0:m-1)
    REAL(fp), INTENT(IN)     :: src_power_in(0:m-1), src_power_out(0:m-1), src_width2(0:m-1)
    REAL(fp), INTENT(IN)     :: fp2(0:m-1), bscat_AD(0:m-1)
    REAL(fp), INTENT(IN OUT) :: ext_AD(0:m-1), ssa_AD(0:m-1), g_AD(0:m-1)
    REAL(fp), INTENT(IN OUT) :: src_power_in_AD(0:m-1), src_power_out_AD(0:m-1)

    REAL(fp), ALLOCATABLE :: Iin(:,:), Iout(:,:), Vin(:,:), Vout(:,:)
    REAL(fp), DIMENSION(-1:m) :: delta0, delta1, delta2, delta3, delta4, delta5
    REAL(fp), DIMENSION(-1:m) :: transmittance, transport_mfp
    REAL(fp), DIMENSION(-1:m) :: delta0_AD, delta1_AD, delta2_AD, delta3_AD, delta4_AD, delta5_AD
    REAL(fp), DIMENSION(-1:m) :: transmittance_AD, transport_mfp_AD
    REAL(fp), DIMENSION(-1:m) :: Iin_AD, Iout_AD, Vin_AD, Vout_AD
    REAL(fp), DIMENSION(-1:m) :: Iin_next_AD, Iout_next_AD, Vin_next_AD, Vout_next_AD
    REAL(fp) :: expansion, expansion_AD, inc_AD, fin, fout, main_in, main_out, cin, cout, tas
    INTEGER  :: nt, it, i, max_i, ibscat

    nt = 2*m
    ALLOCATE(Iin(-1:m,0:nt), Iout(-1:m,0:nt), Vin(-1:m,0:nt), Vout(-1:m,0:nt))
    Iin = ZERO; Iout = ZERO; Vin = ZERO; Vout = ZERO

    CALL tdts_deltas_r(m, drange, ext, ssa, g, delta0, delta1, delta2, delta3, delta4, delta5, &
                       transmittance, transport_mfp)

    ! Forward pass storing the full history
    DO it = 0, nt-1
      max_i = MIN(it, m-1)
      IF (it > nt-m) max_i = nt-it
      DO i = 0, max_i
        Iin(i,it+1)  = Iin(i,it)*delta0(i) + Iout(i,it)*delta1(i) + Iin(i+1,it)*delta2(i+1) &
                     + Iout(i-1,it)*delta3(i-1) + Iin(i-1,it)*delta4(i-1) + Iout(i+1,it)*delta5(i+1)
        Iout(i,it+1) = Iout(i,it)*delta0(i) + Iin(i,it)*delta1(i) + Iout(i-1,it)*delta2(i-1) &
                     + Iin(i+1,it)*delta3(i+1) + Iout(i+1,it)*delta4(i+1) + Iin(i-1,it)*delta5(i-1)
        expansion = tdts_expansion_r(drange, transport_mfp(i), src_width2(i), Iin(i,it), Vin(i,it))
        Vin(i,it+1)  = Vin(i,it)*delta0(i) + Vout(i,it)*delta1(i) + Vin(i+1,it)*delta2(i+1) &
                     + Vout(i-1,it)*delta3(i-1) + Vin(i-1,it)*delta4(i-1) + Vout(i+1,it)*delta5(i+1) &
                     + Iin(i,it)*expansion
        expansion = tdts_expansion_r(drange, transport_mfp(i), src_width2(i), Iout(i,it), Vout(i,it))
        Vout(i,it+1) = Vout(i,it)*delta0(i) + Vin(i,it)*delta1(i) + Vout(i-1,it)*delta2(i-1) &
                     + Vin(i+1,it)*delta3(i+1) + Vout(i+1,it)*delta4(i+1) + Vin(i-1,it)*delta5(i-1) &
                     + Iout(i,it)*expansion
      END DO
      IF (it < m) THEN
        Iin(it,it+1)  = Iin(it,it+1)  + src_power_in(it)
        Iout(it,it+1) = Iout(it,it+1) + src_power_out(it)
        Vin(it,it+1)  = Vin(it,it+1)  + src_power_in(it)*src_width2(it)
        Vout(it,it+1) = Vout(it,it+1) + src_power_out(it)*src_width2(it)
      END IF
    END DO

    ! Adjoint pass
    delta0_AD = ZERO; delta1_AD = ZERO; delta2_AD = ZERO
    delta3_AD = ZERO; delta4_AD = ZERO; delta5_AD = ZERO
    transmittance_AD = ZERO; transport_mfp_AD = ZERO
    Iin_AD = ZERO; Iout_AD = ZERO; Vin_AD = ZERO; Vout_AD = ZERO
    Iin_next_AD = ZERO; Iout_next_AD = ZERO; Vin_next_AD = ZERO; Vout_next_AD = ZERO

    DO it = nt-1, 0, -1
      max_i = MIN(it, m-1)
      IF (it > nt-m) max_i = nt-it
      IF (it < m) THEN
        ! Restore the values before the source was added
        Iin(it,it+1)  = Iin(it,it+1)  - src_power_in(it)
        Iout(it,it+1) = Iout(it,it+1) - src_power_out(it)
        Vin(it,it+1)  = Vin(it,it+1)  - src_power_in(it)*src_width2(it)
        Vout(it,it+1) = Vout(it,it+1) - src_power_out(it)*src_width2(it)
        src_power_in_AD(it)  = src_power_in_AD(it)  + Iin_next_AD(it)  + src_width2(it)*Vin_next_AD(it)
        src_power_out_AD(it) = src_power_out_AD(it) + Iout_next_AD(it) + src_width2(it)*Vout_next_AD(it)
      END IF

      DO i = max_i, 0, -1
        ! Apparent backscatter increment
        ibscat = (it+i)/2
        IF (ibscat < m) THEN
          IF (Vin(i,it+1) > ZERO .AND. Vout(i,it+1) > ZERO) THEN
            inc_AD = bscat_AD(ibscat)
            tas = transmittance(i)*ext(i)*ssa(i)
            IF (g(i) >= ONE/(3.0_fp*MS_MU1)) THEN
              fin = ONE/(ONE + Vin(i,it+1)/(Iin(i,it+1)*fp2(i)))
              main_in = ONE_OVER_FOUR_PI*TWO*Iin(i,it+1)*fin
              transmittance_AD(i) = transmittance_AD(i) + main_in*ssa(i)*ext(i)*inc_AD
              ext_AD(i) = ext_AD(i) + main_in*ssa(i)*transmittance(i)*inc_AD
              ssa_AD(i) = ssa_AD(i) + main_in*ext(i)*transmittance(i)*inc_AD
              Iin_next_AD(i) = Iin_next_AD(i) + ONE_OVER_FOUR_PI*TWO*tas*inc_AD &
                             * (fin + Vin(i,it+1)*fin*fin/(Iin(i,it+1)*fp2(i)))
              Vin_next_AD(i) = Vin_next_AD(i) - ONE_OVER_FOUR_PI*TWO*tas*inc_AD*fin*fin/fp2(i)
            ELSE
              fin  = ONE/(ONE + Vin(i,it+1)/(Iin(i,it+1)*fp2(i)))
              fout = ONE/(ONE + Vout(i,it+1)/(Iout(i,it+1)*fp2(i)))
              main_in  = ONE_OVER_FOUR_PI*Iin(i,it+1)*fin
              main_out = ONE_OVER_FOUR_PI*Iout(i,it+1)*fout
              cin  = ONE + 3.0_fp*g(i)*MS_MU1
              cout = ONE - 3.0_fp*g(i)*MS_MU1
              transmittance_AD(i) = transmittance_AD(i) + ssa(i)*ext(i)*inc_AD*(main_in*cin + main_out*cout)
              ext_AD(i) = ext_AD(i) + ssa(i)*transmittance(i)*inc_AD*(main_in*cin + main_out*cout)
              ssa_AD(i) = ssa_AD(i) + ext(i)*transmittance(i)*inc_AD*(main_in*cin + main_out*cout)
              g_AD(i)   = g_AD(i)   + 3.0_fp*MS_MU1*inc_AD*tas*(main_in - main_out)
              Iin_next_AD(i)  = Iin_next_AD(i)  + ONE_OVER_FOUR_PI*cin*tas*inc_AD &
                              * (fin + Vin(i,it+1)*fin*fin/(Iin(i,it+1)*fp2(i)))
              Vin_next_AD(i)  = Vin_next_AD(i)  - ONE_OVER_FOUR_PI*cin*tas*inc_AD*fin*fin/fp2(i)
              Iout_next_AD(i) = Iout_next_AD(i) + ONE_OVER_FOUR_PI*cout*tas*inc_AD &
                              * (fout + Vout(i,it+1)*fout*fout/(Iout(i,it+1)*fp2(i)))
              Vout_next_AD(i) = Vout_next_AD(i) - ONE_OVER_FOUR_PI*cout*tas*inc_AD*fout*fout/fp2(i)
            END IF
          END IF
        END IF

        ! Vout_next
        expansion = tdts_expansion_r(drange, transport_mfp(i), src_width2(i), Iout(i,it), Vout(i,it))
        Vout_AD(i)     = Vout_AD(i)     + delta0(i)*Vout_next_AD(i)
        delta0_AD(i)   = delta0_AD(i)   + Vout(i,it)*Vout_next_AD(i)
        Vin_AD(i)      = Vin_AD(i)      + delta1(i)*Vout_next_AD(i)
        delta1_AD(i)   = delta1_AD(i)   + Vin(i,it)*Vout_next_AD(i)
        Vout_AD(i-1)   = Vout_AD(i-1)   + delta2(i-1)*Vout_next_AD(i)
        delta2_AD(i-1) = delta2_AD(i-1) + Vout(i-1,it)*Vout_next_AD(i)
        Vin_AD(i+1)    = Vin_AD(i+1)    + delta3(i+1)*Vout_next_AD(i)
        delta3_AD(i+1) = delta3_AD(i+1) + Vin(i+1,it)*Vout_next_AD(i)
        Vout_AD(i+1)   = Vout_AD(i+1)   + delta4(i+1)*Vout_next_AD(i)
        delta4_AD(i+1) = delta4_AD(i+1) + Vout(i+1,it)*Vout_next_AD(i)
        Vin_AD(i-1)    = Vin_AD(i-1)    + delta5(i-1)*Vout_next_AD(i)
        delta5_AD(i-1) = delta5_AD(i-1) + Vin(i-1,it)*Vout_next_AD(i)
        Iout_AD(i)     = Iout_AD(i)     + expansion*Vout_next_AD(i)
        expansion_AD   = Iout(i,it)*Vout_next_AD(i)
        Vout_next_AD(i) = ZERO
        CALL tdts_expansion_AD(drange, transport_mfp(i), src_width2(i), Iout(i,it), Vout(i,it), &
                               expansion_AD, transport_mfp_AD(i), Iout_AD(i), Vout_AD(i))

        ! Vin_next
        expansion = tdts_expansion_r(drange, transport_mfp(i), src_width2(i), Iin(i,it), Vin(i,it))
        Vin_AD(i)      = Vin_AD(i)      + delta0(i)*Vin_next_AD(i)
        delta0_AD(i)   = delta0_AD(i)   + Vin(i,it)*Vin_next_AD(i)
        Vout_AD(i)     = Vout_AD(i)     + delta1(i)*Vin_next_AD(i)
        delta1_AD(i)   = delta1_AD(i)   + Vout(i,it)*Vin_next_AD(i)
        Vin_AD(i+1)    = Vin_AD(i+1)    + delta2(i+1)*Vin_next_AD(i)
        delta2_AD(i+1) = delta2_AD(i+1) + Vin(i+1,it)*Vin_next_AD(i)
        Vout_AD(i-1)   = Vout_AD(i-1)   + delta3(i-1)*Vin_next_AD(i)
        delta3_AD(i-1) = delta3_AD(i-1) + Vout(i-1,it)*Vin_next_AD(i)
        Vin_AD(i-1)    = Vin_AD(i-1)    + delta4(i-1)*Vin_next_AD(i)
        delta4_AD(i-1) = delta4_AD(i-1) + Vin(i-1,it)*Vin_next_AD(i)
        Vout_AD(i+1)   = Vout_AD(i+1)   + delta5(i+1)*Vin_next_AD(i)
        delta5_AD(i+1) = delta5_AD(i+1) + Vout(i+1,it)*Vin_next_AD(i)
        Iin_AD(i)      = Iin_AD(i)      + expansion*Vin_next_AD(i)
        expansion_AD   = Iin(i,it)*Vin_next_AD(i)
        Vin_next_AD(i) = ZERO
        CALL tdts_expansion_AD(drange, transport_mfp(i), src_width2(i), Iin(i,it), Vin(i,it), &
                               expansion_AD, transport_mfp_AD(i), Iin_AD(i), Vin_AD(i))

        ! Iout_next
        Iout_AD(i)     = Iout_AD(i)     + delta0(i)*Iout_next_AD(i)
        delta0_AD(i)   = delta0_AD(i)   + Iout(i,it)*Iout_next_AD(i)
        Iin_AD(i)      = Iin_AD(i)      + delta1(i)*Iout_next_AD(i)
        delta1_AD(i)   = delta1_AD(i)   + Iin(i,it)*Iout_next_AD(i)
        Iout_AD(i-1)   = Iout_AD(i-1)   + delta2(i-1)*Iout_next_AD(i)
        delta2_AD(i-1) = delta2_AD(i-1) + Iout(i-1,it)*Iout_next_AD(i)
        Iin_AD(i+1)    = Iin_AD(i+1)    + delta3(i+1)*Iout_next_AD(i)
        delta3_AD(i+1) = delta3_AD(i+1) + Iin(i+1,it)*Iout_next_AD(i)
        Iout_AD(i+1)   = Iout_AD(i+1)   + delta4(i+1)*Iout_next_AD(i)
        delta4_AD(i+1) = delta4_AD(i+1) + Iout(i+1,it)*Iout_next_AD(i)
        Iin_AD(i-1)    = Iin_AD(i-1)    + delta5(i-1)*Iout_next_AD(i)
        delta5_AD(i-1) = delta5_AD(i-1) + Iin(i-1,it)*Iout_next_AD(i)
        Iout_next_AD(i) = ZERO

        ! Iin_next
        Iin_AD(i)      = Iin_AD(i)      + delta0(i)*Iin_next_AD(i)
        delta0_AD(i)   = delta0_AD(i)   + Iin(i,it)*Iin_next_AD(i)
        Iout_AD(i)     = Iout_AD(i)     + delta1(i)*Iin_next_AD(i)
        delta1_AD(i)   = delta1_AD(i)   + Iout(i,it)*Iin_next_AD(i)
        Iin_AD(i+1)    = Iin_AD(i+1)    + delta2(i+1)*Iin_next_AD(i)
        delta2_AD(i+1) = delta2_AD(i+1) + Iin(i+1,it)*Iin_next_AD(i)
        Iout_AD(i-1)   = Iout_AD(i-1)   + delta3(i-1)*Iin_next_AD(i)
        delta3_AD(i-1) = delta3_AD(i-1) + Iout(i-1,it)*Iin_next_AD(i)
        Iin_AD(i-1)    = Iin_AD(i-1)    + delta4(i-1)*Iin_next_AD(i)
        delta4_AD(i-1) = delta4_AD(i-1) + Iin(i-1,it)*Iin_next_AD(i)
        Iout_AD(i+1)   = Iout_AD(i+1)   + delta5(i+1)*Iin_next_AD(i)
        delta5_AD(i+1) = delta5_AD(i+1) + Iout(i+1,it)*Iin_next_AD(i)
        Iin_next_AD(i) = ZERO
      END DO

      ! Move to the previous time step. Gates that were not updated in this
      ! step keep their adjoint (their value was carried over unchanged only
      ! in the sense that it was never read again; see tdts_core_r).
      Iin_next_AD  = Iin_AD;  Iin_AD  = ZERO
      Iout_next_AD = Iout_AD; Iout_AD = ZERO
      Vin_next_AD  = Vin_AD;  Vin_AD  = ZERO
      Vout_next_AD = Vout_AD; Vout_AD = ZERO
    END DO

    CALL tdts_deltas_AD(m, drange, ext, ssa, g, transmittance, &
                        delta0_AD, delta1_AD, delta2_AD, delta3_AD, delta4_AD, delta5_AD, &
                        transmittance_AD, transport_mfp_AD, ext_AD, ssa_AD, g_AD)
    DEALLOCATE(Iin, Iout, Vin, Vout)
  END SUBROUTINE tdts_core_AD

END MODULE radar_tdts
