!
! test_radar_operator
!
! Consistency tests of the standalone radar operator:
!   1. tangent linear vs. finite differences of the forward model
!   2. adjoint vs. tangent linear (dot-product test)
!   3. K-matrix vs. tangent linear (columns)
! for several configurations (spaceborne/ground, attenuation options,
! multiple scattering, bulk LUT and single-particle species, user gas
! optical depth, directly supplied bulk properties).
!
! Usage: test_radar_operator CloudCoeff.nc4
!
PROGRAM test_radar_operator
  USE radar_kinds
  USE radar_gas,     ONLY: GAS_MODEL_ITU_LBL, GAS_MODEL_ITU_APPROX, GAS_MODEL_USER
  USE radar_species
  USE radar_operator
  IMPLICIT NONE
  CHARACTER(512) :: lut
  INTEGER, PARAMETER :: NL = 50, NS = 5
  TYPE(radar_operator_type) :: op
  TYPE(radar_config_type)   :: cfg
  TYPE(radar_profile_type)  :: prof
  INTEGER :: nfail
  REAL(fp) :: freqs(3) = [13.6_fp, 35.5_fp, 94.05_fp]
  INTEGER :: ifreq

  CALL GET_COMMAND_ARGUMENT(1, lut)
  IF (LEN_TRIM(lut) == 0) THEN
    WRITE(*,*) 'usage: test_radar_operator CloudCoeff.nc4'
    STOP 1
  END IF
  nfail = 0
  CALL make_profile(prof)

  DO ifreq = 1, 3
    cfg = radar_config_type()
    cfg%frequency = freqs(ifreq)
    CALL run_case('spaceborne, delta scaling (CRTM default)', cfg, .FALSE.)
    cfg%scaling = SCALING_NONE; cfg%gate_fraction = 0.5_fp; cfg%crtm_compat = .FALSE.
    CALL run_case('spaceborne, full extinction, gate centre', cfg, .FALSE.)
    cfg%platform = PLATFORM_GROUND; prof%zenith = 30.0_fp
    CALL run_case('ground-based, 30 deg zenith', cfg, .FALSE.)
    cfg%hydro_attenuation = ATTEN_APPROX; cfg%gas_model = GAS_MODEL_ITU_APPROX
    CALL run_case('ground-based, approximate attenuation', cfg, .FALSE.)
    cfg = radar_config_type(); cfg%frequency = freqs(ifreq); prof%zenith = 0.0_fp
    cfg%scaling = SCALING_EFFECTIVE
    CALL run_case('spaceborne, effective attenuation', cfg, .FALSE.)
    cfg%scaling = SCALING_DELTA; cfg%multiple_scattering = .TRUE.
    CALL run_case('spaceborne, TDTS multiple scattering', cfg, .FALSE.)
    cfg = radar_config_type(); cfg%frequency = freqs(ifreq); cfg%gas_model = GAS_MODEL_USER
    CALL run_case('user gas optical depth + bulk properties', cfg, .TRUE.)
    cfg = radar_config_type(); cfg%frequency = freqs(ifreq)
    CALL run_case('single-particle rain species (PSD)', cfg, .FALSE., particle_rain=.TRUE.)
  END DO

  IF (nfail > 0) THEN
    WRITE(*,'(a,i0,a)') 'FAILED: ', nfail, ' test(s)'
    STOP 1
  END IF
  WRITE(*,'(a)') 'All tests passed'

CONTAINS

  SUBROUTINE make_profile(pr)
    TYPE(radar_profile_type), INTENT(OUT) :: pr
    INTEGER :: k
    REAL(fp) :: zc, h
    CALL radar_profile_create(pr, NL, NS, with_bulk=.FALSE.)
    h = 18000.0_fp/NL
    DO k = 0, NL
      pr%z_int(k) = 18000.0_fp - k*h
    END DO
    DO k = 1, NL
      zc = 0.5_fp*(pr%z_int(k-1) + pr%z_int(k))
      pr%t(k) = 295.0_fp - 6.5e-3_fp*MIN(zc, 14000.0_fp)
      pr%p(k) = 1010.0_fp*EXP(-zc/8000.0_fp)
      pr%q(k) = 0.016_fp*EXP(-zc/2200.0_fp)
      ! cloud liquid 1-4 km, rain 0-5 km, cloud ice 7-12 km, snow 5-10 km, graupel 4-7 km
      IF (zc > 1000.0_fp .AND. zc < 4000.0_fp) pr%rho(k,1) = 4.0e-4_fp*SIN(PI*(zc-1000.0_fp)/3000.0_fp)
      IF (zc < 5000.0_fp) pr%rho(k,2) = 1.5e-3_fp*(1.0_fp - zc/5000.0_fp) + 1.0e-4_fp
      IF (zc > 7000.0_fp .AND. zc < 12000.0_fp) pr%rho(k,3) = 5.0e-5_fp*SIN(PI*(zc-7000.0_fp)/5000.0_fp)
      IF (zc > 5000.0_fp .AND. zc < 10000.0_fp) pr%rho(k,4) = 6.0e-4_fp*SIN(PI*(zc-5000.0_fp)/5000.0_fp)
      IF (zc > 4000.0_fp .AND. zc < 7000.0_fp) pr%rho(k,5) = 8.0e-4_fp*SIN(PI*(zc-4000.0_fp)/3000.0_fp)
    END DO
    pr%tau_gas = 2.0e-3_fp
    pr%bulk_ext = 0.0_fp; pr%bulk_sca = 0.0_fp; pr%bulk_bsc = 0.0_fp; pr%bulk_g = 0.0_fp
  END SUBROUTINE make_profile

  ! Rayleigh-type liquid spheres for the single-particle path
  SUBROUTINE rain_particles(f, sp)
    REAL(fp), INTENT(IN) :: f
    TYPE(radar_species_type), INTENT(OUT) :: sp
    INTEGER, PARAMETER :: nd = 80
    REAL(fp) :: d(nd), m(nd), ce(nd), cs(nd), cb(nd), g(nd), lam, kw2, x
    INTEGER :: i
    lam = C_LIGHT/(f*1.0e9_fp)
    kw2 = 0.9_fp
    DO i = 1, nd
      d(i) = 1.0e-4_fp*(60.0_fp)**((i-1.0_fp)/(nd-1.0_fp))     ! 0.1 - 6 mm
      m(i) = 1000.0_fp*PI/6.0_fp*d(i)**3
      x = PI*d(i)/lam
      cb(i) = PI**5*kw2*d(i)**6/lam**4/(ONE + x**2)
      cs(i) = TWO/3.0_fp*PI**5*kw2*d(i)**6/lam**4/(ONE + x**2)
      ce(i) = cs(i) + PI**2*d(i)**3/lam*0.3_fp
      g(i)  = 0.6_fp*x**2/(ONE + x**2)
    END DO
    CALL species_from_particles(sp, 'rain_particles', PHASE_LIQUID, f, d, m, ce, cs, cb, g, 8.0e6_fp, 0.0_fp)
  END SUBROUTINE rain_particles

  SUBROUTINE setup(c, particle_rain)
    TYPE(radar_config_type), INTENT(IN) :: c
    LOGICAL, INTENT(IN) :: particle_rain
    TYPE(radar_species_type) :: sp
    CHARACTER(16) :: names(NS) = [CHARACTER(16) :: 'cloud_liquid', 'rain', 'cloud_ice', 'snow', 'graupel']
    INTEGER :: i
    CALL radar_init(op, c, NS)
    DO i = 1, NS
      IF (particle_rain .AND. i == 2) THEN
        CALL rain_particles(c%frequency, sp)
      ELSE
        CALL species_from_lut(sp, TRIM(names(i)), TRIM(lut), c%frequency, c%n_streams)
      END IF
      CALL radar_add_species(op, i, sp)
    END DO
  END SUBROUTINE setup

  SUBROUTINE random_perturbation(pr, d, bulk)
    TYPE(radar_profile_type), INTENT(IN)    :: pr
    TYPE(radar_profile_type), INTENT(INOUT) :: d
    LOGICAL, INTENT(IN) :: bulk
    REAL(fp) :: r(NL)
    INTEGER :: s
    CALL radar_profile_zero(d)
    CALL RANDOM_NUMBER(r); d%t = (r - 0.5_fp)*0.5_fp
    CALL RANDOM_NUMBER(r); d%q = (r - 0.5_fp)*0.05_fp*pr%q
    DO s = 1, NS
      CALL RANDOM_NUMBER(r); d%rho(:,s) = (r - 0.5_fp)*0.05_fp*pr%rho(:,s)
    END DO
    CALL RANDOM_NUMBER(r); d%tau_gas = (r - 0.5_fp)*1.0e-4_fp
    IF (bulk) THEN
      CALL RANDOM_NUMBER(r); d%bulk_ext = (r - 0.5_fp)*0.05_fp*pr%bulk_ext
      CALL RANDOM_NUMBER(r); d%bulk_sca = (r - 0.5_fp)*0.05_fp*pr%bulk_sca
      CALL RANDOM_NUMBER(r); d%bulk_bsc = (r - 0.5_fp)*0.05_fp*pr%bulk_bsc
      CALL RANDOM_NUMBER(r); d%bulk_g   = (r - 0.5_fp)*0.05_fp*pr%bulk_g
    END IF
  END SUBROUTINE random_perturbation

  SUBROUTINE axpy(a, d, pr, res)
    REAL(fp), INTENT(IN) :: a
    TYPE(radar_profile_type), INTENT(IN)    :: d, pr
    TYPE(radar_profile_type), INTENT(INOUT) :: res
    res%t = pr%t + a*d%t; res%q = pr%q + a*d%q; res%rho = pr%rho + a*d%rho
    res%tau_gas = pr%tau_gas + a*d%tau_gas
    res%bulk_ext = pr%bulk_ext + a*d%bulk_ext; res%bulk_sca = pr%bulk_sca + a*d%bulk_sca
    res%bulk_bsc = pr%bulk_bsc + a*d%bulk_bsc; res%bulk_g = pr%bulk_g + a*d%bulk_g
  END SUBROUTINE axpy

  REAL(fp) FUNCTION dot_out(a, b)
    TYPE(radar_output_type), INTENT(IN) :: a, b
    dot_out = SUM(a%z_lin*b%z_lin) + SUM(a%z_dbz*b%z_dbz) + SUM(a%zatt_lin*b%zatt_lin) + &
              SUM(a%zatt_dbz*b%zatt_dbz) + SUM(a%tau*b%tau) + SUM(a%tau_cum*b%tau_cum) + &
              SUM(a%trans2*b%trans2) + SUM(a%pia*b%pia)
  END FUNCTION dot_out

  REAL(fp) FUNCTION dot_prof(a, b)
    TYPE(radar_profile_type), INTENT(IN) :: a, b
    dot_prof = SUM(a%t*b%t) + SUM(a%q*b%q) + SUM(a%rho*b%rho) + SUM(a%tau_gas*b%tau_gas) + &
               SUM(a%bulk_ext*b%bulk_ext) + SUM(a%bulk_sca*b%bulk_sca) + &
               SUM(a%bulk_bsc*b%bulk_bsc) + SUM(a%bulk_g*b%bulk_g)
  END FUNCTION dot_prof

  SUBROUTINE run_case(label, c, bulk, particle_rain)
    CHARACTER(*), INTENT(IN) :: label
    TYPE(radar_config_type), INTENT(IN) :: c
    LOGICAL, INTENT(IN) :: bulk
    LOGICAL, OPTIONAL, INTENT(IN) :: particle_rain
    TYPE(radar_profile_type) :: pr, d, pp, pm, pad
    TYPE(radar_output_type)  :: o, op_, om, otl, oad
    TYPE(radar_k_type) :: kk
    REAL(fp) :: h, err_tl, err_ad, err_k, lhs, rhs, r(NL), num, den, col(NL)
    LOGICAL :: prain
    INTEGER :: j, it
    prain = .FALSE.
    IF (PRESENT(particle_rain)) prain = particle_rain
    CALL setup(c, prain)
    pr = prof
    IF (bulk) THEN
      pr%has_bulk = .TRUE.
      DO j = 1, NL
        IF (pr%rho(j,2) > 0.0_fp) THEN
          pr%bulk_ext(j) = 2.0e-4_fp; pr%bulk_sca(j) = 1.0e-4_fp
          pr%bulk_bsc(j) = 1.0e-7_fp; pr%bulk_g(j) = 0.3_fp
        END IF
      END DO
    END IF
    CALL radar_output_create(o, NL); CALL radar_output_create(op_, NL); CALL radar_output_create(om, NL)
    CALL radar_output_create(otl, NL); CALL radar_output_create(oad, NL)
    CALL radar_forward(op, pr, o)

    ! 1. TL vs central finite difference (linear outputs, tau, transmittance)
    CALL radar_profile_create(d, NL, NS, bulk); d%has_bulk = bulk
    CALL random_perturbation(pr, d, bulk)
    CALL radar_tl(op, pr, d, otl)
    err_tl = 0.0_fp
    pp = pr; pm = pr
    h = 1.0e-4_fp
    CALL axpy( h, d, pr, pp); CALL axpy(-h, d, pr, pm)
    CALL radar_forward(op, pp, op_); CALL radar_forward(op, pm, om)
    err_tl = MAX(rel(otl%zatt_lin, (op_%zatt_lin - om%zatt_lin)/(2*h)), &
                 rel(otl%z_lin,    (op_%z_lin    - om%z_lin   )/(2*h)), &
                 rel(otl%tau_cum,  (op_%tau_cum  - om%tau_cum )/(2*h)), &
                 rel(otl%trans2,   (op_%trans2   - om%trans2  )/(2*h)))

    ! 2. Dot-product test
    CALL radar_output_zero(oad)
    CALL RANDOM_NUMBER(r); oad%zatt_dbz = r
    CALL RANDOM_NUMBER(r); oad%z_dbz = r
    CALL RANDOM_NUMBER(r); oad%zatt_lin = r
    CALL RANDOM_NUMBER(r); oad%z_lin = r
    CALL RANDOM_NUMBER(r); oad%tau_cum = r
    CALL RANDOM_NUMBER(r); oad%tau = r
    CALL RANDOM_NUMBER(r); oad%trans2 = r
    CALL RANDOM_NUMBER(r); oad%pia = r
    CALL radar_profile_create(pad, NL, NS, bulk); pad%has_bulk = bulk
    CALL radar_ad(op, pr, oad, pad)
    lhs = dot_out(otl, oad)
    rhs = dot_prof(d, pad)
    err_ad = ABS(lhs - rhs)/MAX(ABS(lhs), TINY(lhs))

    ! 3. K-matrix (dBZ attenuated) vs TL with unit perturbations
    CALL radar_k_create(kk, NL, NS)
    CALL radar_k(op, pr, OUT_ZATT_DBZ, kk)
    err_k = 0.0_fp
    DO it = 1, 3
      j = 10*it + 5
      CALL radar_profile_zero(d)
      d%rho(j,2) = 1.0e-5_fp + 0.0_fp*it
      d%t(j) = 0.3_fp
      CALL radar_tl(op, pr, d, otl)
      col = kk%rho(:,j,2)*d%rho(j,2) + kk%t(:,j)*d%t(j)
      num = MAXVAL(ABS(col - otl%zatt_dbz)); den = MAX(MAXVAL(ABS(otl%zatt_dbz)), TINY(den))
      err_k = MAX(err_k, num/den)
    END DO

    WRITE(*,'(f7.2,1x,a45,3(a,es9.2),a)') c%frequency, label, '  TL-FD:', err_tl, '  AD-TL:', err_ad, &
          '  K-TL:', err_k, merge('  ok  ', '  FAIL', err_tl < 1.0e-4_fp .AND. err_ad < 1.0e-10_fp .AND. err_k < 1.0e-10_fp)
    IF (.NOT. (err_tl < 1.0e-4_fp .AND. err_ad < 1.0e-10_fp .AND. err_k < 1.0e-10_fp)) nfail = nfail + 1
    CALL radar_destroy(op)
  END SUBROUTINE run_case

  REAL(fp) FUNCTION rel(a, b)
    REAL(fp), INTENT(IN) :: a(:), b(:)
    rel = MAXVAL(ABS(a - b))/MAX(MAXVAL(ABS(b)), TINY(1.0_fp))
  END FUNCTION rel

END PROGRAM test_radar_operator
