!
! radar_simple
!
! Simple interface for a host model (e.g. GEOS). Set up one handle per
! radar frequency once, then call radar_simple_run for each column.
!
!   CALL radar_simple_init(handle, frequency_ghz, cloudcoeff_file, species, &
!                          platform=PLATFORM_SPACEBORNE)
!   CALL radar_simple_run(handle, zenith_deg, p_hpa, t_k, q_kgkg, z_int_m, wc_kgm3, &
!                         z_dbz, zatt_dbz, z_lin=..., zatt_lin=..., tau=..., &
!                         tau_cum=..., pia=...)
!
!   species  : names of the hydrometeor columns of wc_kgm3, e.g.
!              [character(16) :: 'cloud_liquid', 'cloud_ice', 'rain', 'snow', 'graupel']
!   wc_kgm3  : hydrometeor mass content (kg m-3), (n_layers, n_species);
!              radar_mixing_ratio_to_content converts mixing ratios.
!   z_dbz    : reflectivity without attenuation (dBZ)
!   zatt_dbz : reflectivity attenuated along the beam at zenith_deg (dBZ)
!   tau      : vertical one-way optical depth of each layer (gas + hydrometeors)
!   tau_cum  : vertical one-way optical depth from the radar to each gate;
!              the attenuated reflectivity for any zenith angle theta is
!              z_lin*exp(-2*tau_cum/cos(theta))
!   pia      : two-way path-integrated attenuation along the beam (dB)
!
! Layer 1 is the top of the column (z_int(0) is the top interface).
!
MODULE radar_simple
  USE radar_kinds
  USE radar_species,  ONLY: radar_species_type, species_from_lut
  USE radar_operator
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: radar_simple_init, radar_simple_run, radar_simple_destroy
  PUBLIC :: radar_simple_operator
  PUBLIC :: radar_mixing_ratio_to_content
  PUBLIC :: PLATFORM_SPACEBORNE, PLATFORM_GROUND

  INTEGER, PARAMETER :: MAX_HANDLES = 32
  TYPE(radar_operator_type), SAVE, TARGET :: ops(MAX_HANDLES)

CONTAINS

  SUBROUTINE radar_simple_init(handle, frequency, cloudcoeff_file, species, platform, config, habits)
    INTEGER,                 INTENT(OUT) :: handle
    REAL(fp),                INTENT(IN)  :: frequency
    CHARACTER(*),            INTENT(IN)  :: cloudcoeff_file
    CHARACTER(*),            INTENT(IN)  :: species(:)
    INTEGER,       OPTIONAL, INTENT(IN)  :: platform
    TYPE(radar_config_type), OPTIONAL, INTENT(IN) :: config
    CHARACTER(*),  OPTIONAL, INTENT(IN)  :: habits(:)
    TYPE(radar_config_type)  :: c
    TYPE(radar_species_type) :: sp
    INTEGER :: i
    handle = 0
    DO i = 1, MAX_HANDLES
      IF (.NOT. ops(i)%initialized) THEN
        handle = i
        EXIT
      END IF
    END DO
    IF (handle == 0) THEN
      WRITE(*,*) 'radar_simple_init: no free handle'
      ERROR STOP 1
    END IF
    IF (PRESENT(config)) c = config
    c%frequency = frequency
    IF (PRESENT(platform)) c%platform = platform
    CALL radar_init(ops(handle), c, SIZE(species))
    DO i = 1, SIZE(species)
      IF (PRESENT(habits)) THEN
        CALL species_from_lut(sp, TRIM(species(i)), cloudcoeff_file, frequency, c%n_streams, habit=TRIM(habits(i)))
      ELSE
        CALL species_from_lut(sp, TRIM(species(i)), cloudcoeff_file, frequency, c%n_streams)
      END IF
      CALL radar_add_species(ops(handle), i, sp)
    END DO
  END SUBROUTINE radar_simple_init

  ! Access the operator behind a handle (for TL/AD/K or to add particle species)
  FUNCTION radar_simple_operator(handle) RESULT(op)
    INTEGER, INTENT(IN) :: handle
    TYPE(radar_operator_type), POINTER :: op
    op => ops(handle)
  END FUNCTION radar_simple_operator

  SUBROUTINE radar_simple_run(handle, zenith, p, t, q, z_int, wc, z_dbz, zatt_dbz, &
                              z_lin, zatt_lin, tau, tau_cum, pia)
    INTEGER,  INTENT(IN)  :: handle
    REAL(fp), INTENT(IN)  :: zenith, p(:), t(:), q(:), z_int(0:), wc(:,:)
    REAL(fp), INTENT(OUT) :: z_dbz(:), zatt_dbz(:)
    REAL(fp), OPTIONAL, INTENT(OUT) :: z_lin(:), zatt_lin(:), tau(:), tau_cum(:), pia(:)
    TYPE(radar_profile_type) :: prof
    TYPE(radar_output_type)  :: out
    INTEGER :: n
    n = SIZE(p)
    CALL radar_profile_create(prof, n, ops(handle)%n_species)
    prof%zenith = zenith
    prof%p = p; prof%t = t; prof%q = q; prof%z_int = z_int(0:n)
    prof%rho = wc(1:n, 1:ops(handle)%n_species)
    CALL radar_output_create(out, n)
    CALL radar_forward(ops(handle), prof, out)
    z_dbz = out%z_dbz
    zatt_dbz = out%zatt_dbz
    IF (PRESENT(z_lin))    z_lin = out%z_lin
    IF (PRESENT(zatt_lin)) zatt_lin = out%zatt_lin
    IF (PRESENT(tau))      tau = out%tau
    IF (PRESENT(tau_cum))  tau_cum = out%tau_cum
    IF (PRESENT(pia))      pia = out%pia
  END SUBROUTINE radar_simple_run

  SUBROUTINE radar_simple_destroy(handle)
    INTEGER, INTENT(IN) :: handle
    CALL radar_destroy(ops(handle))
  END SUBROUTINE radar_simple_destroy

  ! Mass content (kg m-3) from a mixing ratio (kg/kg of moist air)
  ELEMENTAL FUNCTION radar_mixing_ratio_to_content(qx, p_hpa, t_k, qv) RESULT(wc)
    REAL(fp), INTENT(IN) :: qx, p_hpa, t_k, qv
    REAL(fp) :: wc
    wc = qx*p_hpa*100.0_fp/(R_DRY*t_k*(ONE + (ONE/EPS_WV - ONE)*qv))
  END FUNCTION radar_mixing_ratio_to_content

END MODULE radar_simple
