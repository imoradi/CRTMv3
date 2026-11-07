!
! radar_operator
!
! Standalone radar reflectivity operator (forward, tangent-linear, adjoint
! and K-matrix), independent of CRTM. It reproduces the CRTM active-sensor
! calculation (CRTM_Active_Sensor) when the CRTM-compatible settings are used
! (the defaults), and adds:
!   - spaceborne (looking down) and ground-based (looking up) radars
!   - reflectivity with and without attenuation, linear and dBZ
!   - vertical (zenith-angle independent) layer and cumulative optical
!     depths, so attenuation can be applied later for any zenith angle
!   - full or approximate gas and hydrometeor attenuation
!   - bulk lookup tables, single particles with a PSD, or bulk optical
!     properties supplied directly per layer
!
! Conventions
!   Layer 1 is the top of the profile; z_int(0) is the top interface and
!   z_int(n) the bottom (m). Inputs per layer: pressure p (hPa),
!   temperature t (K), specific humidity q (kg/kg) and hydrometeor water
!   content rho(:,s) (kg m-3) for each species s.
!   Reflectivity: mm^6 m^-3 (linear) and dBZ.
!   Optical depths are one-way, vertical (not divided by cos(zenith)).
!   Two-way transmittance to a gate:  trans2 = exp(-2*tau_cum/cos(zenith)).
!
MODULE radar_operator
  USE radar_kinds
  USE radar_gas,     ONLY: GAS_MODEL_NONE, GAS_MODEL_ITU_LBL, GAS_MODEL_ITU_APPROX, GAS_MODEL_USER, &
                           gas_absorption_deriv, gas_absorption
  USE radar_species, ONLY: radar_species_type, N_OPT, I_EXT, I_SCA, I_BSC, I_FSCA, I_GSCA, I_ATT, &
                           species_optics_deriv, species_optics, species_destroy, SOURCE_NONE
  USE radar_tdts,    ONLY: TDTS_Wide_Angle, TDTS_Wide_Angle_TL, TDTS_Wide_Angle_AD
  IMPLICIT NONE
  PRIVATE

  ! Platforms
  INTEGER, PARAMETER, PUBLIC :: PLATFORM_SPACEBORNE = 1   ! looking down from above the profile
  INTEGER, PARAMETER, PUBLIC :: PLATFORM_GROUND     = 2   ! looking up from the bottom of the profile
  ! Hydrometeor attenuation
  INTEGER, PARAMETER, PUBLIC :: ATTEN_NONE   = 0          ! gas only
  INTEGER, PARAMETER, PUBLIC :: ATTEN_FULL   = 1          ! extinction from the scattering properties
  INTEGER, PARAMETER, PUBLIC :: ATTEN_APPROX = 2          ! ITU-R P.838/P.840 type formulas
  ! Treatment of scattering in the attenuation (ATTEN_FULL only)
  INTEGER, PARAMETER, PUBLIC :: SCALING_NONE      = 0     ! full extinction (single scattering)
  INTEGER, PARAMETER, PUBLIC :: SCALING_DELTA     = 1     ! delta-truncated extinction (CRTM default)
  INTEGER, PARAMETER, PUBLIC :: SCALING_EFFECTIVE = 2     ! CRTM ACTIVE_MS_EFFECTIVE_ATTENUATION
  ! |K|^2 of water
  INTEGER, PARAMETER, PUBLIC :: KW_TURNER_273 = 1         ! Turner et al. (2016) at 273.15 K (CRTM)
  INTEGER, PARAMETER, PUBLIC :: KW_FIXED      = 2         ! config%kw2
  ! Outputs for the K-matrix
  INTEGER, PARAMETER, PUBLIC :: OUT_Z_LIN = 1, OUT_Z_DBZ = 2, OUT_ZATT_LIN = 3, OUT_ZATT_DBZ = 4
  INTEGER, PARAMETER, PUBLIC :: OUT_TAU_CUM = 5, OUT_PIA = 6
  REAL(fp), PARAMETER, PUBLIC :: MISSING_DBZ = -9999.0_fp

  PUBLIC :: radar_config_type, radar_operator_type, radar_profile_type, radar_output_type
  PUBLIC :: radar_k_type
  PUBLIC :: radar_init, radar_add_species, radar_destroy
  PUBLIC :: radar_profile_create, radar_profile_zero, radar_output_create, radar_output_zero
  PUBLIC :: radar_k_create
  PUBLIC :: radar_forward, radar_tl, radar_ad, radar_k
  PUBLIC :: water_kw2_turner

  ! ---------------------------------------------------------------------------
  TYPE :: radar_config_type
    REAL(fp) :: frequency          = 13.6_fp          ! GHz
    INTEGER  :: platform           = PLATFORM_SPACEBORNE
    INTEGER  :: gas_model          = GAS_MODEL_ITU_LBL
    INTEGER  :: hydro_attenuation  = ATTEN_FULL
    INTEGER  :: scaling            = SCALING_DELTA
    INTEGER  :: n_streams          = 4                ! delta truncation order (CRTM n_Streams)
    REAL(fp) :: eta                = 0.5_fp           ! SCALING_EFFECTIVE forward fraction
    ! Wide-angle multiple scattering (Hogan and Battaglia 2008 TDTS), spaceborne only
    LOGICAL  :: multiple_scattering = .FALSE.
    REAL(fp) :: altitude           = -1.0_fp          ! m, <=0: 407 km (<50 GHz) or 705 km
    REAL(fp) :: beamwidth          = -1.0_fp          ! deg, <=0: 0.71 (<50 GHz) or 0.108
    REAL(fp) :: cbe                = 1.0_fp           ! coherent backscatter enhancement
    INTEGER  :: kw_model           = KW_TURNER_273
    REAL(fp) :: kw2                = 0.93_fp
    ! Fraction of the gate's own layer included in the attenuation to the gate
    ! (1 = whole layer, as CRTM; 0.5 = to the layer centre)
    REAL(fp) :: gate_fraction      = 1.0_fp
    ! Species with zero single-scatter albedo give no backscatter (CRTM)
    LOGICAL  :: crtm_compat        = .TRUE.
    ! Layers whose water path (kg m-2) is below this are skipped (CRTM)
    REAL(fp) :: wc_threshold       = 1.0e-6_fp
    ! dBZ TL/AD: zero below tlad_min_dbz; Z held >= tlad_cap_dbz in the gain
    REAL(fp) :: tlad_min_dbz       = 0.0_fp
    REAL(fp) :: tlad_cap_dbz       = 0.0_fp
  END TYPE radar_config_type

  TYPE :: radar_operator_type
    LOGICAL  :: initialized = .FALSE.
    TYPE(radar_config_type) :: config
    REAL(fp) :: wavelength = ZERO        ! m
    REAL(fp) :: kw2 = ZERO
    REAL(fp) :: radar_constant = ZERO    ! Z (mm6 m-3) = radar_constant * backscatter (m-1)
    INTEGER  :: n_species = 0
    TYPE(radar_species_type), ALLOCATABLE :: species(:)
  END TYPE radar_operator_type

  ! Profile (also used for TL and AD perturbations: t, q, rho, tau_gas, bulk_*)
  TYPE :: radar_profile_type
    INTEGER  :: n_layers = 0, n_species = 0
    REAL(fp) :: zenith = ZERO            ! beam zenith angle (deg); for a ground radar 90 - elevation
    REAL(fp), ALLOCATABLE :: p(:), t(:), q(:)
    REAL(fp), ALLOCATABLE :: z_int(:)    ! (0:n) m
    REAL(fp), ALLOCATABLE :: rho(:,:)    ! (n, n_species) kg m-3
    REAL(fp), ALLOCATABLE :: tau_gas(:)  ! GAS_MODEL_USER: vertical one-way gas optical depth
    ! Bulk optical properties supplied directly (added to the species)
    LOGICAL  :: has_bulk = .FALSE.
    REAL(fp), ALLOCATABLE :: bulk_ext(:), bulk_sca(:), bulk_bsc(:), bulk_g(:)   ! m-1, m-1, m-1, -
  END TYPE radar_profile_type

  TYPE :: radar_output_type
    INTEGER :: n_layers = 0
    REAL(fp), ALLOCATABLE :: z_lin(:), z_dbz(:)          ! unattenuated
    REAL(fp), ALLOCATABLE :: zatt_lin(:), zatt_dbz(:)    ! attenuated along the beam (zenith angle)
    REAL(fp), ALLOCATABLE :: tau(:)                      ! vertical layer optical depth used for attenuation
    REAL(fp), ALLOCATABLE :: tau_gas(:), tau_hydro(:)    ! its gas and hydrometeor parts
    REAL(fp), ALLOCATABLE :: tau_cum(:)                  ! vertical one-way optical depth radar->gate
    REAL(fp), ALLOCATABLE :: trans2(:)                   ! two-way transmittance along the beam
    REAL(fp), ALLOCATABLE :: pia(:)                      ! two-way path-integrated attenuation (dB)
    REAL(fp), ALLOCATABLE :: z_ms(:)                     ! multiple-scattering contribution (mm6 m-3)
  END TYPE radar_output_type

  TYPE :: radar_k_type
    INTEGER :: n_layers = 0, n_species = 0
    REAL(fp), ALLOCATABLE :: t(:,:), q(:,:)        ! (gate, layer)
    REAL(fp), ALLOCATABLE :: rho(:,:,:)            ! (gate, layer, species)
  END TYPE radar_k_type

  ! Internal: base state of one profile
  TYPE :: work_type
    INTEGER  :: n = 0, ns = 0
    REAL(fp) :: mu = ONE
    REAL(fp), ALLOCATABLE :: dz(:)
    REAL(fp), ALLOCATABLE :: ag(:), dag_dt(:), dag_dq(:)              ! gas (m-1)
    LOGICAL,  ALLOCATABLE :: active(:,:)
    REAL(fp), ALLOCATABLE :: opt(:,:,:), dopt_drho(:,:,:), dopt_dt(:,:,:)   ! (N_OPT, n, ns)
    REAL(fp), ALLOCATABLE :: e(:), s(:), f(:), g(:), b(:), a(:)       ! layer totals (m-1)
    LOGICAL,  ALLOCATABLE :: delta_on(:)
    REAL(fp), ALLOCATABLE :: kappa(:), tau(:), scum(:), trans2(:), z(:), zatt(:), zms(:)
    ! TDTS
    REAL(fp) :: rho_fov = ZERO, cbe = ONE
    REAL(fp), ALLOCATABLE :: edge(:), ms_w(:), ms_g(:)
    LOGICAL,  ALLOCATABLE :: ms_scat(:)
  END TYPE work_type

CONTAINS

  ! ===========================================================================
  ! Setup
  ! ===========================================================================

  ! |K|^2 from the Turner et al. (2016) liquid water permittivity (as CRTM)
  FUNCTION water_kw2_turner(freq_ghz, temp_k) RESULT(kw2)
    REAL(fp), INTENT(IN) :: freq_ghz, temp_k
    REAL(fp) :: kw2
    REAL(fp) :: tc, eps_s, d1, d2, t1, t2, w, den1, den2, ere, eim
    COMPLEX(fp) :: eps, k
    REAL(fp), PARAMETER :: a1 = 8.111e+1_fp, a2 = 2.025e+0_fp, b1 = 4.434e-3_fp, b2 = 1.073e-2_fp
    REAL(fp), PARAMETER :: c1 = 1.302e-13_fp, c2 = 1.012e-14_fp, dd1 = 6.627e+2_fp, dd2 = 6.089e+2_fp
    REAL(fp), PARAMETER :: tcc = 133.1383_fp
    REAL(fp), PARAMETER :: e0 = 87.9144_fp, e1 = 0.404399_fp, e2 = 9.58726e-4_fp, e3 = -1.32802e-6_fp
    tc = temp_k - 273.15_fp
    eps_s = e0 + e1*tc + e2*tc**2 + e3*tc**3
    d1 = a1*EXP(-b1*tc)
    d2 = a2*EXP(-b2*tc)
    t1 = c1*EXP(dd1/(tc + tcc))
    t2 = c2*EXP(dd2/(tc + tcc))
    w  = TWO*PI*freq_ghz*1.0e9_fp
    den1 = ONE + (w*t1)**2
    den2 = ONE + (w*t2)**2
    ere = eps_s - w**2*(t1**2*d1/den1 + t2**2*d2/den2)
    eim = w*(t1*d1/den1 + t2*d2/den2)
    eps = CMPLX(ere, eim, fp)
    k = (eps - ONE)/(eps + TWO)
    kw2 = ABS(k)**2
  END FUNCTION water_kw2_turner

  !
  ! radar_init: set up the operator for one frequency. Species are added
  ! afterwards with radar_add_species (in the order of the rho columns).
  !
  SUBROUTINE radar_init(op, config, n_species)
    TYPE(radar_operator_type), INTENT(OUT) :: op
    TYPE(radar_config_type),   INTENT(IN)  :: config
    INTEGER,                   INTENT(IN)  :: n_species
    op%config = config
    op%wavelength = C_LIGHT/(config%frequency*1.0e9_fp)
    IF (config%kw_model == KW_FIXED) THEN
      op%kw2 = config%kw2
    ELSE
      op%kw2 = water_kw2_turner(config%frequency, 273.15_fp)
    END IF
    op%radar_constant = M6_TO_MM6*op%wavelength**4/(PI**5*op%kw2)
    op%n_species = n_species
    ALLOCATE(op%species(n_species))
    IF (op%config%multiple_scattering .AND. op%config%platform /= PLATFORM_SPACEBORNE) THEN
      WRITE(*,*) 'radar_init: multiple scattering is only available for spaceborne radars; disabled'
      op%config%multiple_scattering = .FALSE.
    END IF
    op%initialized = .TRUE.
  END SUBROUTINE radar_init

  SUBROUTINE radar_add_species(op, index, sp)
    TYPE(radar_operator_type), INTENT(INOUT) :: op
    INTEGER,                   INTENT(IN)    :: index
    TYPE(radar_species_type),  INTENT(IN)    :: sp
    op%species(index) = sp
    op%species(index)%zero_bsc_if_no_scattering = op%config%crtm_compat
  END SUBROUTINE radar_add_species

  SUBROUTINE radar_destroy(op)
    TYPE(radar_operator_type), INTENT(INOUT) :: op
    INTEGER :: i
    IF (ALLOCATED(op%species)) THEN
      DO i = 1, SIZE(op%species)
        CALL species_destroy(op%species(i))
      END DO
      DEALLOCATE(op%species)
    END IF
    op%initialized = .FALSE.
  END SUBROUTINE radar_destroy

  SUBROUTINE radar_profile_create(prof, n_layers, n_species, with_bulk)
    TYPE(radar_profile_type), INTENT(OUT) :: prof
    INTEGER,                  INTENT(IN)  :: n_layers, n_species
    LOGICAL, OPTIONAL,        INTENT(IN)  :: with_bulk
    prof%n_layers = n_layers
    prof%n_species = n_species
    ALLOCATE(prof%p(n_layers), prof%t(n_layers), prof%q(n_layers), prof%z_int(0:n_layers), &
             prof%rho(n_layers, n_species), prof%tau_gas(n_layers))
    prof%has_bulk = .FALSE.
    IF (PRESENT(with_bulk)) prof%has_bulk = with_bulk
    ALLOCATE(prof%bulk_ext(n_layers), prof%bulk_sca(n_layers), prof%bulk_bsc(n_layers), prof%bulk_g(n_layers))
    CALL radar_profile_zero(prof)
  END SUBROUTINE radar_profile_create

  SUBROUTINE radar_profile_zero(prof)
    TYPE(radar_profile_type), INTENT(INOUT) :: prof
    prof%p = ZERO; prof%t = ZERO; prof%q = ZERO; prof%z_int = ZERO; prof%rho = ZERO
    prof%tau_gas = ZERO
    prof%bulk_ext = ZERO; prof%bulk_sca = ZERO; prof%bulk_bsc = ZERO; prof%bulk_g = ZERO
  END SUBROUTINE radar_profile_zero

  SUBROUTINE radar_output_create(out, n_layers)
    TYPE(radar_output_type), INTENT(OUT) :: out
    INTEGER,                 INTENT(IN)  :: n_layers
    out%n_layers = n_layers
    ALLOCATE(out%z_lin(n_layers), out%z_dbz(n_layers), out%zatt_lin(n_layers), out%zatt_dbz(n_layers), &
             out%tau(n_layers), out%tau_gas(n_layers), out%tau_hydro(n_layers), out%tau_cum(n_layers), &
             out%trans2(n_layers), out%pia(n_layers), out%z_ms(n_layers))
    CALL radar_output_zero(out)
  END SUBROUTINE radar_output_create

  SUBROUTINE radar_output_zero(out)
    TYPE(radar_output_type), INTENT(INOUT) :: out
    out%z_lin = ZERO; out%z_dbz = ZERO; out%zatt_lin = ZERO; out%zatt_dbz = ZERO
    out%tau = ZERO; out%tau_gas = ZERO; out%tau_hydro = ZERO; out%tau_cum = ZERO
    out%trans2 = ZERO; out%pia = ZERO; out%z_ms = ZERO
  END SUBROUTINE radar_output_zero

  SUBROUTINE radar_k_create(k, n_layers, n_species)
    TYPE(radar_k_type), INTENT(OUT) :: k
    INTEGER,            INTENT(IN)  :: n_layers, n_species
    k%n_layers = n_layers; k%n_species = n_species
    ALLOCATE(k%t(n_layers, n_layers), k%q(n_layers, n_layers), k%rho(n_layers, n_layers, n_species))
    k%t = ZERO; k%q = ZERO; k%rho = ZERO
  END SUBROUTINE radar_k_create

  ! ===========================================================================
  ! Base state
  ! ===========================================================================
  SUBROUTINE compute_base(op, prof, wk, need_deriv)
    TYPE(radar_operator_type), INTENT(IN)    :: op
    TYPE(radar_profile_type),  INTENT(IN)    :: prof
    TYPE(work_type),           INTENT(INOUT) :: wk
    LOGICAL,                   INTENT(IN)    :: need_deriv
    INTEGER  :: n, ns, k, s
    REAL(fp) :: alt, bw, gam
    TYPE(radar_config_type) :: c

    c  = op%config
    n  = prof%n_layers
    ns = op%n_species
    wk%n = n; wk%ns = ns
    ALLOCATE(wk%dz(n), wk%ag(n), wk%dag_dt(n), wk%dag_dq(n), wk%active(n,ns), &
             wk%opt(N_OPT,n,ns), wk%dopt_drho(N_OPT,n,ns), wk%dopt_dt(N_OPT,n,ns), &
             wk%e(n), wk%s(n), wk%f(n), wk%g(n), wk%b(n), wk%a(n), wk%delta_on(n), &
             wk%kappa(n), wk%tau(n), wk%scum(n), wk%trans2(n), wk%z(n), wk%zatt(n), wk%zms(n), &
             wk%edge(0:n), wk%ms_w(n), wk%ms_g(n), wk%ms_scat(n))

    wk%dz = prof%z_int(0:n-1) - prof%z_int(1:n)
    IF (ANY(wk%dz <= ZERO)) THEN
      WRITE(*,*) 'radar_operator: z_int must decrease with layer index (layer 1 at the top)'
      ERROR STOP 1
    END IF
    wk%mu = COS(prof%zenith*PI/180.0_fp)

    ! Gas absorption
    DO k = 1, n
      IF (c%gas_model == GAS_MODEL_ITU_LBL .OR. c%gas_model == GAS_MODEL_ITU_APPROX) THEN
        IF (need_deriv) THEN
          CALL gas_absorption_deriv(c%gas_model, c%frequency, prof%p(k), prof%t(k), prof%q(k), &
                                    wk%ag(k), wk%dag_dt(k), wk%dag_dq(k))
        ELSE
          wk%ag(k) = gas_absorption(c%gas_model, c%frequency, prof%p(k), prof%t(k), prof%q(k))
          wk%dag_dt(k) = ZERO; wk%dag_dq(k) = ZERO
        END IF
      ELSE
        wk%ag(k) = ZERO; wk%dag_dt(k) = ZERO; wk%dag_dq(k) = ZERO
      END IF
    END DO

    ! Hydrometeor optics
    wk%opt = ZERO; wk%dopt_drho = ZERO; wk%dopt_dt = ZERO
    DO s = 1, ns
      DO k = 1, n
        wk%active(k,s) = prof%rho(k,s)*wk%dz(k) > c%wc_threshold .AND. op%species(s)%source /= SOURCE_NONE
        IF (.NOT. wk%active(k,s)) CYCLE
        IF (need_deriv) THEN
          CALL species_optics_deriv(op%species(s), prof%rho(k,s), prof%t(k), wk%opt(:,k,s), &
                                    wk%dopt_drho(:,k,s), wk%dopt_dt(:,k,s))
        ELSE
          CALL species_optics(op%species(s), prof%rho(k,s), prof%t(k), wk%opt(:,k,s))
        END IF
      END DO
    END DO

    ! Layer totals
    DO k = 1, n
      wk%e(k) = wk%ag(k) + SUM(wk%opt(I_EXT,k,:))
      wk%s(k) = SUM(wk%opt(I_SCA,k,:))
      wk%f(k) = SUM(wk%opt(I_FSCA,k,:))
      wk%g(k) = SUM(wk%opt(I_GSCA,k,:))
      wk%b(k) = SUM(wk%opt(I_BSC,k,:))
      wk%a(k) = SUM(wk%opt(I_ATT,k,:))
      IF (prof%has_bulk) THEN
        wk%e(k) = wk%e(k) + prof%bulk_ext(k)
        wk%s(k) = wk%s(k) + prof%bulk_sca(k)
        wk%g(k) = wk%g(k) + prof%bulk_g(k)*prof%bulk_sca(k)
        wk%b(k) = wk%b(k) + prof%bulk_bsc(k)
      END IF
      ! CRTM applies the delta truncation only to layers with significant scattering
      wk%delta_on(k) = wk%s(k)*wk%dz(k) > 1.0e-10_fp
    END DO

    ! Attenuation coefficient
    DO k = 1, n
      SELECT CASE (c%hydro_attenuation)
      CASE (ATTEN_NONE)
        wk%kappa(k) = wk%ag(k)
      CASE (ATTEN_APPROX)
        wk%kappa(k) = wk%ag(k) + wk%a(k)
        IF (prof%has_bulk) wk%kappa(k) = wk%kappa(k) + prof%bulk_ext(k)
      CASE DEFAULT
        IF (c%multiple_scattering) THEN
          wk%kappa(k) = wk%e(k)
        ELSE
          SELECT CASE (c%scaling)
          CASE (SCALING_DELTA)
            wk%kappa(k) = wk%e(k)
            IF (wk%delta_on(k)) wk%kappa(k) = wk%e(k) - wk%f(k)
          CASE (SCALING_EFFECTIVE)
            IF (wk%delta_on(k)) THEN
              wk%kappa(k) = wk%e(k) - wk%f(k) - c%eta*(wk%s(k) - wk%f(k))
            ELSE
              wk%kappa(k) = wk%e(k) - c%eta*wk%s(k)
            END IF
          CASE DEFAULT
            wk%kappa(k) = wk%e(k)
          END SELECT
        END IF
      END SELECT
    END DO
    wk%tau = wk%kappa*wk%dz
    IF (c%gas_model == GAS_MODEL_USER) wk%tau = wk%tau + prof%tau_gas

    ! Cumulative optical depth radar -> gate
    gam = c%gate_fraction
    IF (c%platform == PLATFORM_GROUND) THEN
      DO k = 1, n
        wk%scum(k) = SUM(wk%tau(k+1:n)) + gam*wk%tau(k)
      END DO
    ELSE
      DO k = 1, n
        wk%scum(k) = SUM(wk%tau(1:k-1)) + gam*wk%tau(k)
      END DO
    END IF
    wk%trans2 = EXP(-TWO*wk%scum/wk%mu)

    wk%z = op%radar_constant*wk%b
    wk%zatt = wk%z*wk%trans2

    ! Wide-angle multiple scattering (unscaled optics, as in CRTM)
    wk%zms = ZERO
    IF (c%multiple_scattering .AND. c%hydro_attenuation == ATTEN_FULL) THEN
      alt = c%altitude
      bw  = c%beamwidth
      IF (alt <= ZERO) alt = MERGE(705000.0_fp, 407000.0_fp, c%frequency > 50.0_fp)
      IF (bw  <= ZERO) bw  = MERGE(0.108_fp, 0.71_fp, c%frequency > 50.0_fp)
      wk%rho_fov = (bw*PI/180.0_fp)/(TWO*SQRT(LOG(TWO)))
      wk%cbe = c%cbe
      wk%edge = (alt - prof%z_int)/wk%mu
      DO k = 1, n
        wk%ms_scat(k) = wk%s(k) > 1.0e-10_fp*wk%e(k) .AND. wk%s(k) > ZERO
        IF (wk%ms_scat(k)) THEN
          wk%ms_w(k) = wk%s(k)/wk%e(k)
          wk%ms_g(k) = wk%g(k)/wk%s(k)
        ELSE
          wk%ms_w(k) = ZERO; wk%ms_g(k) = ZERO
        END IF
      END DO
      CALL TDTS_Wide_Angle(n, wk%edge, wk%e, wk%ms_w, wk%ms_g, wk%rho_fov, wk%rho_fov, wk%cbe, wk%zms)
      wk%zms = op%radar_constant*FOUR_PI*wk%zms
      wk%zatt = wk%zatt + wk%zms
    END IF
  END SUBROUTINE compute_base

  SUBROUTINE tlad_limits(c, zmin, zcap)
    TYPE(radar_config_type), INTENT(IN)  :: c
    REAL(fp),                INTENT(OUT) :: zmin, zcap
    zmin = TEN**(c%tlad_min_dbz/TEN)
    zcap = MAX(TEN**(c%tlad_cap_dbz/TEN), zmin)
  END SUBROUTINE tlad_limits

  ! ===========================================================================
  ! Forward
  ! ===========================================================================
  SUBROUTINE radar_forward(op, prof, out)
    TYPE(radar_operator_type), INTENT(IN)    :: op
    TYPE(radar_profile_type),  INTENT(IN)    :: prof
    TYPE(radar_output_type),   INTENT(INOUT) :: out
    TYPE(work_type) :: wk
    CALL compute_base(op, prof, wk, .FALSE.)
    out%z_lin    = wk%z
    out%zatt_lin = wk%zatt
    WHERE (wk%z > TINY(ONE))
      out%z_dbz = TEN*LOG10(wk%z)
    ELSEWHERE
      out%z_dbz = MISSING_DBZ
    END WHERE
    WHERE (wk%zatt > TINY(ONE))
      out%zatt_dbz = TEN*LOG10(wk%zatt)
    ELSEWHERE
      out%zatt_dbz = MISSING_DBZ
    END WHERE
    out%tau = wk%tau
    IF (op%config%gas_model == GAS_MODEL_USER) THEN
      out%tau_gas = prof%tau_gas
    ELSE
      out%tau_gas = wk%ag*wk%dz
    END IF
    out%tau_hydro = wk%tau - out%tau_gas
    out%tau_cum = wk%scum
    out%trans2 = wk%trans2
    out%pia = TWO*TEN/LOG(TEN)*wk%scum/wk%mu
    out%z_ms = wk%zms
  END SUBROUTINE radar_forward

  ! ===========================================================================
  ! Tangent linear
  ! ===========================================================================
  SUBROUTINE radar_tl(op, prof, prof_tl, out_tl)
    TYPE(radar_operator_type), INTENT(IN)    :: op
    TYPE(radar_profile_type),  INTENT(IN)    :: prof, prof_tl
    TYPE(radar_output_type),   INTENT(INOUT) :: out_tl
    TYPE(work_type) :: wk
    TYPE(radar_config_type) :: c
    INTEGER  :: n, k, s
    REAL(fp) :: zmin, zcap, gam, dx(N_OPT)
    REAL(fp), DIMENSION(prof%n_layers) :: ag_tl, e_tl, s_tl, f_tl, g_tl, b_tl, a_tl, kappa_tl, tau_tl, &
                                          scum_tl, t2_tl, z_tl, zatt_tl, zms_tl, w_tl, gg_tl, bs_tl

    c = op%config
    CALL compute_base(op, prof, wk, .TRUE.)
    n = wk%n
    CALL tlad_limits(c, zmin, zcap)

    ag_tl = wk%dag_dt*prof_tl%t + wk%dag_dq*prof_tl%q
    e_tl = ag_tl; s_tl = ZERO; f_tl = ZERO; g_tl = ZERO; b_tl = ZERO; a_tl = ZERO
    DO s = 1, wk%ns
      DO k = 1, n
        IF (.NOT. wk%active(k,s)) CYCLE
        dx = wk%dopt_drho(:,k,s)*prof_tl%rho(k,s) + wk%dopt_dt(:,k,s)*prof_tl%t(k)
        e_tl(k) = e_tl(k) + dx(I_EXT)
        s_tl(k) = s_tl(k) + dx(I_SCA)
        f_tl(k) = f_tl(k) + dx(I_FSCA)
        g_tl(k) = g_tl(k) + dx(I_GSCA)
        b_tl(k) = b_tl(k) + dx(I_BSC)
        a_tl(k) = a_tl(k) + dx(I_ATT)
      END DO
    END DO
    IF (prof%has_bulk) THEN
      e_tl = e_tl + prof_tl%bulk_ext
      s_tl = s_tl + prof_tl%bulk_sca
      g_tl = g_tl + prof%bulk_g*prof_tl%bulk_sca + prof_tl%bulk_g*prof%bulk_sca
      b_tl = b_tl + prof_tl%bulk_bsc
    END IF

    DO k = 1, n
      SELECT CASE (c%hydro_attenuation)
      CASE (ATTEN_NONE)
        kappa_tl(k) = ag_tl(k)
      CASE (ATTEN_APPROX)
        kappa_tl(k) = ag_tl(k) + a_tl(k)
        IF (prof%has_bulk) kappa_tl(k) = kappa_tl(k) + prof_tl%bulk_ext(k)
      CASE DEFAULT
        IF (c%multiple_scattering) THEN
          kappa_tl(k) = e_tl(k)
        ELSE
          SELECT CASE (c%scaling)
          CASE (SCALING_DELTA)
            kappa_tl(k) = e_tl(k)
            IF (wk%delta_on(k)) kappa_tl(k) = e_tl(k) - f_tl(k)
          CASE (SCALING_EFFECTIVE)
            IF (wk%delta_on(k)) THEN
              kappa_tl(k) = e_tl(k) - f_tl(k) - c%eta*(s_tl(k) - f_tl(k))
            ELSE
              kappa_tl(k) = e_tl(k) - c%eta*s_tl(k)
            END IF
          CASE DEFAULT
            kappa_tl(k) = e_tl(k)
          END SELECT
        END IF
      END SELECT
    END DO
    tau_tl = kappa_tl*wk%dz
    IF (c%gas_model == GAS_MODEL_USER) tau_tl = tau_tl + prof_tl%tau_gas

    gam = c%gate_fraction
    IF (c%platform == PLATFORM_GROUND) THEN
      DO k = 1, n
        scum_tl(k) = SUM(tau_tl(k+1:n)) + gam*tau_tl(k)
      END DO
    ELSE
      DO k = 1, n
        scum_tl(k) = SUM(tau_tl(1:k-1)) + gam*tau_tl(k)
      END DO
    END IF
    t2_tl = -TWO/wk%mu*wk%trans2*scum_tl
    z_tl = op%radar_constant*b_tl
    zatt_tl = z_tl*wk%trans2 + wk%z*t2_tl

    zms_tl = ZERO
    IF (c%multiple_scattering .AND. c%hydro_attenuation == ATTEN_FULL) THEN
      DO k = 1, n
        IF (wk%ms_scat(k)) THEN
          w_tl(k)  = (s_tl(k) - wk%ms_w(k)*e_tl(k))/wk%e(k)
          gg_tl(k) = (g_tl(k) - wk%ms_g(k)*s_tl(k))/wk%s(k)
        ELSE
          w_tl(k) = ZERO; gg_tl(k) = ZERO
        END IF
      END DO
      CALL TDTS_Wide_Angle_TL(n, wk%edge, wk%e, wk%ms_w, wk%ms_g, wk%rho_fov, wk%rho_fov, wk%cbe, &
                              e_tl, w_tl, gg_tl, bs_tl, zms_tl)
      zms_tl = op%radar_constant*FOUR_PI*zms_tl
      zatt_tl = zatt_tl + zms_tl
    END IF

    out_tl%z_lin = z_tl
    out_tl%zatt_lin = zatt_tl
    WHERE (wk%z > zmin)
      out_tl%z_dbz = TEN*z_tl/(MAX(wk%z, zcap)*LOG(TEN))
    ELSEWHERE
      out_tl%z_dbz = ZERO
    END WHERE
    WHERE (wk%zatt > zmin)
      out_tl%zatt_dbz = TEN*zatt_tl/(MAX(wk%zatt, zcap)*LOG(TEN))
    ELSEWHERE
      out_tl%zatt_dbz = ZERO
    END WHERE
    out_tl%tau = tau_tl
    IF (c%gas_model == GAS_MODEL_USER) THEN
      out_tl%tau_gas = prof_tl%tau_gas
    ELSE
      out_tl%tau_gas = ag_tl*wk%dz
    END IF
    out_tl%tau_hydro = tau_tl - out_tl%tau_gas
    out_tl%tau_cum = scum_tl
    out_tl%trans2 = t2_tl
    out_tl%pia = TWO*TEN/LOG(TEN)*scum_tl/wk%mu
    out_tl%z_ms = zms_tl
  END SUBROUTINE radar_tl

  ! ===========================================================================
  ! Adjoint. Adjoint inputs in out_ad: z_lin, z_dbz, zatt_lin, zatt_dbz, tau,
  ! tau_cum, trans2, pia (accumulated, so any combination can be seeded).
  ! Results are accumulated into prof_ad (t, q, rho, tau_gas, bulk_*).
  ! ===========================================================================
  SUBROUTINE radar_ad(op, prof, out_ad, prof_ad)
    TYPE(radar_operator_type), INTENT(IN)    :: op
    TYPE(radar_profile_type),  INTENT(IN)    :: prof
    TYPE(radar_output_type),   INTENT(IN)    :: out_ad
    TYPE(radar_profile_type),  INTENT(INOUT) :: prof_ad
    TYPE(work_type) :: wk
    TYPE(radar_config_type) :: c
    INTEGER  :: n, k, s, j
    REAL(fp) :: zmin, zcap, gam, acc, ox(N_OPT)
    REAL(fp), DIMENSION(prof%n_layers) :: z_ad, zatt_ad, t2_ad, scum_ad, tau_ad, kappa_ad, &
                                          ag_ad, e_ad, s_ad, f_ad, g_ad, b_ad, a_ad, &
                                          bs_ad, w_ad, gg_ad, ext_ad

    c = op%config
    CALL compute_base(op, prof, wk, .TRUE.)
    n = wk%n
    CALL tlad_limits(c, zmin, zcap)

    z_ad = out_ad%z_lin
    zatt_ad = out_ad%zatt_lin
    WHERE (wk%z > zmin) z_ad = z_ad + TEN*out_ad%z_dbz/(MAX(wk%z, zcap)*LOG(TEN))
    WHERE (wk%zatt > zmin) zatt_ad = zatt_ad + TEN*out_ad%zatt_dbz/(MAX(wk%zatt, zcap)*LOG(TEN))

    e_ad = ZERO; s_ad = ZERO; f_ad = ZERO; g_ad = ZERO; b_ad = ZERO; a_ad = ZERO; ag_ad = ZERO

    ! Multiple scattering
    IF (c%multiple_scattering .AND. c%hydro_attenuation == ATTEN_FULL) THEN
      bs_ad = op%radar_constant*FOUR_PI*zatt_ad
      ext_ad = ZERO; w_ad = ZERO; gg_ad = ZERO
      CALL TDTS_Wide_Angle_AD(n, wk%edge, wk%e, wk%ms_w, wk%ms_g, wk%rho_fov, wk%rho_fov, wk%cbe, &
                              bs_ad, ext_ad, w_ad, gg_ad)
      e_ad = e_ad + ext_ad
      DO k = 1, n
        IF (wk%ms_scat(k)) THEN
          ! gg = g/s
          g_ad(k) = g_ad(k) + gg_ad(k)/wk%s(k)
          s_ad(k) = s_ad(k) - gg_ad(k)*wk%ms_g(k)/wk%s(k)
          ! w = s/e
          s_ad(k) = s_ad(k) + w_ad(k)/wk%e(k)
          e_ad(k) = e_ad(k) - w_ad(k)*wk%ms_w(k)/wk%e(k)
        END IF
      END DO
    END IF

    ! zatt = z*trans2 (+ zms)
    z_ad  = z_ad + wk%trans2*zatt_ad
    t2_ad = out_ad%trans2 + wk%z*zatt_ad
    scum_ad = out_ad%tau_cum + TWO*TEN/LOG(TEN)*out_ad%pia/wk%mu - TWO/wk%mu*wk%trans2*t2_ad
    b_ad = b_ad + op%radar_constant*z_ad

    ! Cumulative sums
    gam = c%gate_fraction
    tau_ad = out_ad%tau
    IF (c%platform == PLATFORM_GROUND) THEN
      acc = ZERO
      DO j = 1, n
        ! scum(k) = sum_{j>k} tau(j) + gam*tau(k)
        tau_ad(j) = tau_ad(j) + acc + gam*scum_ad(j)
        acc = acc + scum_ad(j)
      END DO
    ELSE
      acc = ZERO
      DO j = n, 1, -1
        tau_ad(j) = tau_ad(j) + acc + gam*scum_ad(j)
        acc = acc + scum_ad(j)
      END DO
    END IF

    IF (c%gas_model == GAS_MODEL_USER) prof_ad%tau_gas = prof_ad%tau_gas + tau_ad
    kappa_ad = tau_ad*wk%dz

    DO k = 1, n
      SELECT CASE (c%hydro_attenuation)
      CASE (ATTEN_NONE)
        ag_ad(k) = ag_ad(k) + kappa_ad(k)
      CASE (ATTEN_APPROX)
        ag_ad(k) = ag_ad(k) + kappa_ad(k)
        a_ad(k)  = a_ad(k) + kappa_ad(k)
        IF (prof%has_bulk) prof_ad%bulk_ext(k) = prof_ad%bulk_ext(k) + kappa_ad(k)
      CASE DEFAULT
        IF (c%multiple_scattering) THEN
          e_ad(k) = e_ad(k) + kappa_ad(k)
        ELSE
          SELECT CASE (c%scaling)
          CASE (SCALING_DELTA)
            e_ad(k) = e_ad(k) + kappa_ad(k)
            IF (wk%delta_on(k)) f_ad(k) = f_ad(k) - kappa_ad(k)
          CASE (SCALING_EFFECTIVE)
            e_ad(k) = e_ad(k) + kappa_ad(k)
            s_ad(k) = s_ad(k) - c%eta*kappa_ad(k)
            IF (wk%delta_on(k)) f_ad(k) = f_ad(k) - (ONE - c%eta)*kappa_ad(k)
          CASE DEFAULT
            e_ad(k) = e_ad(k) + kappa_ad(k)
          END SELECT
        END IF
      END SELECT
    END DO

    ! Layer totals -> species, gas and bulk
    IF (c%hydro_attenuation == ATTEN_FULL) ag_ad = ag_ad + e_ad
    IF (prof%has_bulk) THEN
      IF (c%hydro_attenuation == ATTEN_FULL) prof_ad%bulk_ext = prof_ad%bulk_ext + e_ad
      prof_ad%bulk_sca = prof_ad%bulk_sca + s_ad + prof%bulk_g*g_ad
      prof_ad%bulk_g   = prof_ad%bulk_g + prof%bulk_sca*g_ad
      prof_ad%bulk_bsc = prof_ad%bulk_bsc + b_ad
    END IF
    DO s = 1, wk%ns
      DO k = 1, n
        IF (.NOT. wk%active(k,s)) CYCLE
        ox(I_EXT) = e_ad(k); ox(I_SCA) = s_ad(k); ox(I_FSCA) = f_ad(k)
        ox(I_GSCA) = g_ad(k); ox(I_BSC) = b_ad(k); ox(I_ATT) = a_ad(k)
        prof_ad%rho(k,s) = prof_ad%rho(k,s) + SUM(wk%dopt_drho(:,k,s)*ox)
        prof_ad%t(k)     = prof_ad%t(k)     + SUM(wk%dopt_dt(:,k,s)*ox)
      END DO
    END DO
    prof_ad%t = prof_ad%t + wk%dag_dt*ag_ad
    prof_ad%q = prof_ad%q + wk%dag_dq*ag_ad
  END SUBROUTINE radar_ad

  ! ===========================================================================
  ! K-matrix (Jacobians) of one output quantity, by repeated adjoint calls.
  ! k%t(i,j) = d out(i) / d t(j), etc.
  ! ===========================================================================
  SUBROUTINE radar_k(op, prof, which, k)
    TYPE(radar_operator_type), INTENT(IN)    :: op
    TYPE(radar_profile_type),  INTENT(IN)    :: prof
    INTEGER,                   INTENT(IN)    :: which
    TYPE(radar_k_type),        INTENT(INOUT) :: k
    TYPE(radar_output_type)  :: out_ad
    TYPE(radar_profile_type) :: prof_ad
    INTEGER :: i, n
    n = prof%n_layers
    CALL radar_output_create(out_ad, n)
    CALL radar_profile_create(prof_ad, n, prof%n_species, prof%has_bulk)
    DO i = 1, n
      CALL radar_output_zero(out_ad)
      CALL radar_profile_zero(prof_ad)
      SELECT CASE (which)
      CASE (OUT_Z_LIN);    out_ad%z_lin(i) = ONE
      CASE (OUT_Z_DBZ);    out_ad%z_dbz(i) = ONE
      CASE (OUT_ZATT_LIN); out_ad%zatt_lin(i) = ONE
      CASE (OUT_TAU_CUM);  out_ad%tau_cum(i) = ONE
      CASE (OUT_PIA);      out_ad%pia(i) = ONE
      CASE DEFAULT;        out_ad%zatt_dbz(i) = ONE
      END SELECT
      CALL radar_ad(op, prof, out_ad, prof_ad)
      k%t(i,:) = prof_ad%t
      k%q(i,:) = prof_ad%q
      k%rho(i,:,:) = prof_ad%rho
    END DO
  END SUBROUTINE radar_k

END MODULE radar_operator
