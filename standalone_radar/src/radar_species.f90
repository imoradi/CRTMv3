!
! radar_species
!
! Optical properties of one hydrometeor species at the radar frequency.
! A species gets its properties from one of three sources:
!
!   SOURCE_LUT      Bulk properties tabulated against water content and
!                   temperature (liquid) or water content and habit (frozen),
!                   read from a CRTM CloudCoeff netCDF file, e.g. the DDA
!                   tables CloudCoeff_DDA_Moradi_*.nc4. The interpolation is
!                   the same as in CRTM_CloudScatter so results match CRTM.
!   SOURCE_PARTICLE Single-particle cross sections (extinction, scattering,
!                   radar backscattering, asymmetry, mass) on a size grid,
!                   integrated over a gamma particle size distribution
!                   N(D) = N0 D^mu exp(-lambda D). The slope lambda is solved
!                   from the water content, so the water content stays the
!                   control variable for the TL/AD.
!   (per-layer bulk properties given directly by the caller are handled in
!    radar_operator, not here.)
!
! For each layer the species returns the volume coefficients (m-1)
!   ext  extinction             sca  scattering
!   bsc  radar backscattering (sigma_b = 4 pi dsigma/dOmega at 180 deg)
!   fsca f*sca, f = delta-truncation fraction (CRTM convention)
!   gsca g*sca, g = asymmetry parameter
!   att  attenuation coefficient from the approximate formulas
! and their derivatives with respect to water content and temperature.
!
MODULE radar_species
  USE radar_kinds,        ONLY: fp, ZERO, ONE, CSTEP
  USE radar_interp,       ONLY: interp_weights_c
  USE radar_hydro_approx, ONLY: APPROX_NONE, APPROX_CLOUD_LIQUID, APPROX_ICE, APPROX_RAIN, &
                                hydro_approx_attenuation_c
  USE netcdf
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: radar_species_type
  PUBLIC :: SOURCE_NONE, SOURCE_LUT, SOURCE_PARTICLE
  PUBLIC :: PHASE_LIQUID, PHASE_FROZEN
  PUBLIC :: PSD_GAMMA_FIXED_N0, PSD_GAMMA, PSD_T08_RAIN, PSD_T08_SNOW, PSD_T08_SNOW_FIELD05, PSD_T08_GRAUPEL, PSD_T08_ICE, &
            PSD_F07_MIDLAT, PSD_F07_TROPICAL
  PUBLIC :: species_from_ssdb, species_default_ssdb
  PUBLIC :: N_OPT, I_EXT, I_SCA, I_BSC, I_FSCA, I_GSCA, I_ATT
  PUBLIC :: species_from_lut
  PUBLIC :: species_from_particles
  PUBLIC :: species_optics
  PUBLIC :: species_optics_deriv
  PUBLIC :: species_destroy
  PUBLIC :: lut_habit_names

  INTEGER, PARAMETER :: SOURCE_NONE = 0, SOURCE_LUT = 1, SOURCE_PARTICLE = 2
  INTEGER, PARAMETER :: PHASE_LIQUID = 0, PHASE_FROZEN = 1
  ! Particle size distributions
  !   PSD_GAMMA_FIXED_N0  N0 D^mu exp(-lambda D), lambda solved so the PSD holds the
  !                       water content with the given particle masses
  !   PSD_GAMMA           as NuRTM: lambda from a power-law mass-size fit to the
  !                       database, then renormalised to the database masses
  !   PSD_T08_*           Thompson et al. (2008) rain, snow (exponential with the
  !                       T04 temperature-dependent intercept, or the two-mode
  !                       Field et al. 2005 form) and graupel, as NuRTM
  INTEGER, PARAMETER :: PSD_GAMMA_FIXED_N0 = 1, PSD_GAMMA = 2
  INTEGER, PARAMETER :: PSD_T08_RAIN = 3, PSD_T08_SNOW = 4, PSD_T08_SNOW_FIELD05 = 5, PSD_T08_GRAUPEL = 6, PSD_T08_ICE = 7
  INTEGER, PARAMETER :: PSD_F07_MIDLAT = 8, PSD_F07_TROPICAL = 9
  ! Field et al. (2005) moment relation coefficients (WRF module_mp_thompson sa, sb)
  REAL(fp), PARAMETER :: F05_SA(10) = [5.065339_fp, -0.062659_fp, -3.032362_fp, 0.029469_fp, -0.000285_fp, &
                                       0.31255_fp, 0.000204_fp, 0.003199_fp, 0.0_fp, -0.015952_fp]
  REAL(fp), PARAMETER :: F05_SB(10) = [0.476221_fp, -0.015896_fp, 0.165977_fp, 0.007468_fp, -0.000141_fp, &
                                       0.060366_fp, 0.000079_fp, 0.000594_fp, 0.0_fp, -0.003577_fp]

  ! Indices of the optical quantities
  INTEGER, PARAMETER :: N_OPT = 6
  INTEGER, PARAMETER :: I_EXT = 1, I_SCA = 2, I_BSC = 3, I_FSCA = 4, I_GSCA = 5, I_ATT = 6

  ! Habits in the CRTM DDA CloudCoeff files (index in the *_S_MW arrays)
  INTEGER, PARAMETER :: N_DDA_HABITS = 18
  CHARACTER(*), PARAMETER :: DDA_HABITS(N_DDA_HABITS) = [ CHARACTER(32) :: &
       'PlateType1', 'ColumnType1', '6BulletRosette', 'Perpendicular4BulletRosette', &
       'Flat3BulletRosette', 'IconCloudIce', 'SectorSnowflake', 'EvansSnowAggregate', &
       '8ColumnAggregate', 'LargePlateAggregate', 'LargeColumnAggregate', &
       'LargeBlockAggregate', 'IconSnow', 'IconHail', 'GemGraupel', 'GemSnow', &
       'GemHail', 'IceSphere' ]

  TYPE :: radar_species_type
    CHARACTER(64) :: name = ''
    INTEGER  :: source = SOURCE_NONE
    INTEGER  :: phase  = PHASE_LIQUID
    INTEGER  :: approx_kind = APPROX_NONE
    ! Non-scattering species (CRTM WATER_CLOUD/ICE_CLOUD): properties taken at
    ! the first water-content grid point and single-scatter albedo zero
    LOGICAL  :: rayleigh_only = .FALSE.
    ! CRTM sets the backscatter to zero when the single-scatter albedo is zero
    LOGICAL  :: zero_bsc_if_no_scattering = .TRUE.
    REAL(fp) :: frequency = ZERO
    ! ---- LUT source (already interpolated to the radar frequency)
    INTEGER  :: nx = 0, nt = 0
    REAL(fp), ALLOCATABLE :: x(:)          ! water content grid (kg m-3)
    REAL(fp), ALLOCATABLE :: tgrid(:)      ! temperature grid (K), liquid only
    REAL(fp), ALLOCATABLE :: ke(:,:)       ! mass extinction (m2 kg-1)   (nx, nt)
    REAL(fp), ALLOCATABLE :: kb(:,:)       ! mass backscatter (m2 kg-1)
    REAL(fp), ALLOCATABLE :: w(:,:)        ! single-scatter albedo
    REAL(fp), ALLOCATABLE :: pc1(:,:)      ! Legendre coefficient l=1 (g = pc1/1.5)
    REAL(fp), ALLOCATABLE :: pcn(:,:)      ! Legendre coefficient l=n_streams
    ! ---- PARTICLE source
    INTEGER  :: nd = 0
    INTEGER  :: psd_type = PSD_GAMMA_FIXED_N0
    REAL(fp) :: n0 = 8.0e6_fp              ! m^-(4+mu)
    REAL(fp) :: mu = ZERO
    REAL(fp) :: a_ms = ZERO, b_ms = 3.0_fp ! mass-size relation m = a D^b (PSD_GAMMA)
    REAL(fp) :: rho_air = 1.225_fp         ! graupel intercept (PSD_T08_GRAUPEL)
    REAL(fp) :: ni_fixed = ZERO            ! PSD_T08_ICE: >0 fixed number concentration (m-3),
                                           ! otherwise Cooper (1986) from temperature
    REAL(fp) :: t_fix = ZERO               ! >0: fixed temperature for the T08 snow PSD (K)
    INTEGER  :: n_streams = 4
    INTEGER  :: ntp = 1                    ! temperatures of the cross sections
    REAL(fp), ALLOCATABLE :: tgp(:)        ! (ntp) K
    REAL(fp), ALLOCATABLE :: d(:)          ! particle size (m) on which the PSD is defined
    REAL(fp), ALLOCATABLE :: dd(:)         ! quadrature weights (m)
    REAL(fp), ALLOCATABLE :: mass(:)       ! particle mass (kg)
    REAL(fp), ALLOCATABLE :: cext(:,:), csca(:,:), cbsc(:,:), asym(:,:)   ! (nd, ntp) m2, m2, m2, -
  END TYPE radar_species_type

CONTAINS

  FUNCTION lut_habit_names() RESULT(names)
    CHARACTER(32) :: names(N_DDA_HABITS)
    names = DDA_HABITS
  END FUNCTION lut_habit_names

  SUBROUTINE nc_check(status, what)
    INTEGER,      INTENT(IN) :: status
    CHARACTER(*), INTENT(IN) :: what
    IF (status /= NF90_NOERR) THEN
      WRITE(*,'(a)') 'radar_species: netCDF error ('//TRIM(what)//'): '//TRIM(nf90_strerror(status))
      ERROR STOP 1
    END IF
  END SUBROUTINE nc_check

  !
  ! species_from_lut
  !   Build a species from a CRTM CloudCoeff netCDF file.
  !   name : 'cloud_liquid' (CRTM WATER_CLOUD), 'rain', 'cloud_ice' (CRTM
  !          ICE_CLOUD), 'snow', 'graupel', 'hail', or any DDA habit name
  !          (e.g. 'LargePlateAggregate') for a frozen species.
  !   frequency (GHz), n_streams sets which Legendre coefficient is used for
  !   the delta truncation (CRTM: 4 -> l=4 of the 4-stream set).
  !   habit (optional) overrides the default habit of snow/graupel/hail/ice.
  !
  SUBROUTINE species_from_lut(sp, name, filename, frequency, n_streams, habit)
    TYPE(radar_species_type), INTENT(OUT) :: sp
    CHARACTER(*),             INTENT(IN)  :: name, filename
    REAL(fp),                 INTENT(IN)  :: frequency
    INTEGER,                  INTENT(IN)  :: n_streams
    CHARACTER(*), OPTIONAL,   INTENT(IN)  :: habit
    INTEGER :: ncid, varid, dimid, nf, nr, ntt, nh, nl, i1, ih, j, k
    INTEGER :: loff
    REAL(fp), ALLOCATABLE :: freq(:), a3(:,:,:), p5(:,:,:,:,:)
    REAL(fp) :: wf(4)
    COMPLEX(fp) :: wfc(4)
    CHARACTER(32) :: hname

    sp%name = name
    sp%source = SOURCE_LUT
    sp%frequency = frequency
    sp%n_streams = n_streams
    SELECT CASE (n_streams)
    CASE (4);  loff = 0
    CASE (6);  loff = 5
    CASE (8);  loff = 12
    CASE (16); loff = 21
    CASE DEFAULT
      WRITE(*,*) 'radar_species: n_streams must be 4, 6, 8 or 16'
      ERROR STOP 1
    END SELECT

    ! Species defaults (match CRTM CLOUD_INDEX_DDA_ARTS)
    hname = ''
    SELECT CASE (TRIM(name))
    CASE ('cloud_liquid', 'cloud_water', 'water')
      sp%phase = PHASE_LIQUID; sp%rayleigh_only = .TRUE.;  sp%approx_kind = APPROX_CLOUD_LIQUID
    CASE ('rain')
      sp%phase = PHASE_LIQUID; sp%rayleigh_only = .FALSE.; sp%approx_kind = APPROX_RAIN
    CASE ('cloud_ice', 'ice')
      sp%phase = PHASE_FROZEN; sp%rayleigh_only = .TRUE.;  sp%approx_kind = APPROX_ICE; hname = 'IceSphere'
    CASE ('snow')
      sp%phase = PHASE_FROZEN; sp%rayleigh_only = .FALSE.; sp%approx_kind = APPROX_ICE; hname = 'SectorSnowflake'
    CASE ('graupel')
      sp%phase = PHASE_FROZEN; sp%rayleigh_only = .FALSE.; sp%approx_kind = APPROX_ICE; hname = 'GemGraupel'
    CASE ('hail')
      sp%phase = PHASE_FROZEN; sp%rayleigh_only = .FALSE.; sp%approx_kind = APPROX_ICE; hname = 'GemHail'
    CASE DEFAULT
      sp%phase = PHASE_FROZEN; sp%rayleigh_only = .FALSE.; sp%approx_kind = APPROX_ICE; hname = name
    END SELECT
    IF (PRESENT(habit)) THEN
      IF (LEN_TRIM(habit) > 0) hname = habit
    END IF
    ih = 0
    IF (sp%phase == PHASE_FROZEN) THEN
      DO j = 1, N_DDA_HABITS
        IF (TRIM(DDA_HABITS(j)) == TRIM(hname)) ih = j
      END DO
      IF (ih == 0) THEN
        WRITE(*,*) 'radar_species: unknown habit '//TRIM(hname)
        ERROR STOP 1
      END IF
    END IF

    CALL nc_check(nf90_open(filename, NF90_NOWRITE, ncid), 'open '//TRIM(filename))
    CALL nc_check(nf90_inq_dimid(ncid, 'n_MW_Frequencies', dimid), 'dim f')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=nf), 'dim f')
    CALL nc_check(nf90_inq_dimid(ncid, 'n_MW_Radii', dimid), 'dim r')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=nr), 'dim r')
    CALL nc_check(nf90_inq_dimid(ncid, 'n_Temperatures', dimid), 'dim t')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=ntt), 'dim t')
    CALL nc_check(nf90_inq_dimid(ncid, 'n_MW_Densities', dimid), 'dim h')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=nh), 'dim h')
    CALL nc_check(nf90_inq_dimid(ncid, 'n_Legendre_Terms', dimid), 'dim l')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=nl), 'dim l')

    ALLOCATE(freq(nf))
    CALL nc_check(nf90_inq_varid(ncid, 'Frequency_MW', varid), 'Frequency_MW')
    CALL nc_check(nf90_get_var(ncid, varid, freq), 'Frequency_MW')
    ! Frequency weights (CRTM clamps to the table range)
    CALL interp_weights_c(freq, CMPLX(frequency, ZERO, fp), i1, wfc)
    wf = REAL(wfc)

    sp%nx = nr
    ALLOCATE(sp%x(nr))
    CALL nc_check(nf90_inq_varid(ncid, 'Water_Density_MW', varid), 'Water_Density_MW')
    CALL nc_check(nf90_get_var(ncid, varid, sp%x), 'Water_Density_MW')

    IF (sp%phase == PHASE_LIQUID) THEN
      sp%nt = ntt
      ALLOCATE(sp%tgrid(ntt))
      CALL nc_check(nf90_inq_varid(ncid, 'Temperature', varid), 'Temperature')
      CALL nc_check(nf90_get_var(ncid, varid, sp%tgrid), 'Temperature')
      ALLOCATE(a3(4, nr, ntt), p5(4, nr, ntt, nl, 1))
      ALLOCATE(sp%ke(nr,ntt), sp%kb(nr,ntt), sp%w(nr,ntt), sp%pc1(nr,ntt), sp%pcn(nr,ntt))
      CALL read3('ke_L_MW', sp%ke)
      CALL read3('kb_L_MW', sp%kb)
      CALL read3('w_L_MW',  sp%w)
      CALL nc_check(nf90_inq_varid(ncid, 'pcoeff_L_MW', varid), 'pcoeff_L_MW')
      CALL nc_check(nf90_get_var(ncid, varid, p5, start=[i1,1,1,1,1], count=[4,nr,ntt,nl,1]), 'pcoeff_L_MW')
      DO k = 1, ntt
        DO j = 1, nr
          sp%pc1(j,k) = SUM(wf*p5(:,j,k,1+loff+1,1))
          sp%pcn(j,k) = SUM(wf*p5(:,j,k,n_streams+loff+1,1))
        END DO
      END DO
    ELSE
      sp%nt = 1
      ALLOCATE(sp%tgrid(1)); sp%tgrid = 273.15_fp
      ALLOCATE(a3(4, nr, nh), p5(4, nr, nh, nl, 1))
      ALLOCATE(sp%ke(nr,1), sp%kb(nr,1), sp%w(nr,1), sp%pc1(nr,1), sp%pcn(nr,1))
      CALL read3h('ke_S_MW', sp%ke(:,1))
      CALL read3h('kb_S_MW', sp%kb(:,1))
      CALL read3h('w_S_MW',  sp%w(:,1))
      CALL nc_check(nf90_inq_varid(ncid, 'pcoeff_S_MW', varid), 'pcoeff_S_MW')
      CALL nc_check(nf90_get_var(ncid, varid, p5, start=[i1,1,1,1,1], count=[4,nr,nh,nl,1]), 'pcoeff_S_MW')
      DO j = 1, nr
        sp%pc1(j,1) = SUM(wf*p5(:,j,ih,1+loff+1,1))
        sp%pcn(j,1) = SUM(wf*p5(:,j,ih,n_streams+loff+1,1))
      END DO
    END IF
    CALL nc_check(nf90_close(ncid), 'close')
    IF (sp%rayleigh_only) THEN
      ! CRTM WATER_CLOUD / ICE_CLOUD: no scattering
      sp%w = ZERO; sp%pc1 = ZERO; sp%pcn = ZERO
    END IF
  CONTAINS
    SUBROUTINE read3(vname, out)
      CHARACTER(*), INTENT(IN)  :: vname
      REAL(fp),     INTENT(OUT) :: out(:,:)
      INTEGER :: jj, kk
      CALL nc_check(nf90_inq_varid(ncid, vname, varid), vname)
      CALL nc_check(nf90_get_var(ncid, varid, a3, start=[i1,1,1], count=[4,nr,ntt]), vname)
      DO kk = 1, ntt
        DO jj = 1, nr
          out(jj,kk) = SUM(wf*a3(:,jj,kk))
        END DO
      END DO
    END SUBROUTINE read3
    SUBROUTINE read3h(vname, out)
      CHARACTER(*), INTENT(IN)  :: vname
      REAL(fp),     INTENT(OUT) :: out(:)
      INTEGER :: jj
      CALL nc_check(nf90_inq_varid(ncid, vname, varid), vname)
      CALL nc_check(nf90_get_var(ncid, varid, a3, start=[i1,1,1], count=[4,nr,nh]), vname)
      DO jj = 1, nr
        out(jj) = SUM(wf*a3(:,jj,ih))
      END DO
    END SUBROUTINE read3h
  END SUBROUTINE species_from_lut

  !
  ! species_from_particles
  !   d     particle size (m) on which the PSD is defined (melted or maximum
  !         dimension, consistent with n0/mu)
  !   mass  particle mass (kg)
  !   cext, csca, cbsc  extinction, scattering and radar backscattering
  !         cross sections (m2) at the radar frequency;
  !         cbsc = 4 pi * dsigma/dOmega(180 deg) = 4 pi Z11(180)
  !   asym  asymmetry parameter
  !   n0, mu  gamma PSD N(D) = n0 D^mu exp(-lambda D)   (m^-4-mu)
  !
  SUBROUTINE species_from_particles(sp, name, phase, frequency, d, mass, cext, csca, cbsc, asym, &
                                    n0, mu, n_streams, approx_kind)
    TYPE(radar_species_type), INTENT(OUT) :: sp
    CHARACTER(*),             INTENT(IN)  :: name
    INTEGER,                  INTENT(IN)  :: phase
    REAL(fp),                 INTENT(IN)  :: frequency, d(:), mass(:), cext(:), csca(:), cbsc(:), asym(:)
    REAL(fp),                 INTENT(IN)  :: n0, mu
    INTEGER,  OPTIONAL,       INTENT(IN)  :: n_streams, approx_kind
    INTEGER :: n
    n = SIZE(d)
    sp%name = name
    sp%source = SOURCE_PARTICLE
    sp%phase = phase
    sp%frequency = frequency
    sp%nd = n
    sp%n0 = n0
    sp%mu = mu
    sp%rayleigh_only = .FALSE.
    IF (PRESENT(n_streams)) sp%n_streams = n_streams
    IF (PRESENT(approx_kind)) THEN
      sp%approx_kind = approx_kind
    ELSE IF (phase == PHASE_LIQUID) THEN
      sp%approx_kind = APPROX_RAIN
    ELSE
      sp%approx_kind = APPROX_ICE
    END IF
    sp%psd_type = PSD_GAMMA_FIXED_N0
    sp%ntp = 1
    ALLOCATE(sp%tgp(1)); sp%tgp = 273.15_fp
    ALLOCATE(sp%d(n), sp%dd(n), sp%mass(n), sp%cext(n,1), sp%csca(n,1), sp%cbsc(n,1), sp%asym(n,1))
    sp%d = d; sp%mass = mass
    sp%cext(:,1) = cext; sp%csca(:,1) = csca; sp%cbsc(:,1) = cbsc; sp%asym(:,1) = asym
    CALL bin_widths(d, sp%dd)
  END SUBROUTINE species_from_particles

  !
  ! species_default_ssdb
  !   The recommended species set built from the NuRTM single-scattering
  !   databases in ssdb_dir, all temperature dependent:
  !     cloud_liquid  LiquidSphere    gamma N0 = 3e35, mu = 4 (d_veq)
  !     rain          LiquidSphere    Thompson (2008) rain (d_veq)
  !     cloud_ice     IceSphere       Thompson (2008) ice, Ni from Cooper (1986)
  !     snow          6BulletRosette  Field et al. (2007) mid-latitude, layer T
  !     graupel       GemGraupel      Field et al. (2007) mid-latitude, layer T
  !     hail          GemHail         Field et al. (2007) mid-latitude, layer T
  !   The cross sections are interpolated in temperature for every species
  !   (database grids: 230-310 K liquid, 190-270 K ice); the frozen PSDs use
  !   the layer temperature (Field 2007 held within -70 to -0.1 C).
  !
  SUBROUTINE species_default_ssdb(sp, name, ssdb_dir, frequency, n_streams)
    TYPE(radar_species_type), INTENT(OUT) :: sp
    CHARACTER(*),             INTENT(IN)  :: name, ssdb_dir
    REAL(fp),                 INTENT(IN)  :: frequency
    INTEGER,        OPTIONAL, INTENT(IN)  :: n_streams
    INTEGER :: ns
    ns = 4
    IF (PRESENT(n_streams)) ns = n_streams
    SELECT CASE (TRIM(name))
    CASE ('cloud_liquid')
      CALL species_from_ssdb(sp, name, TRIM(ssdb_dir)//'/LiquidSphere.nc', frequency, PSD_GAMMA, &
                             n0=3.0e35_fp, mu=4.0_fp, size_param='d_veq', n_streams=ns, &
                             approx_kind=APPROX_CLOUD_LIQUID)
    CASE ('rain')
      CALL species_from_ssdb(sp, name, TRIM(ssdb_dir)//'/LiquidSphere.nc', frequency, PSD_T08_RAIN, &
                             size_param='d_veq', n_streams=ns, approx_kind=APPROX_RAIN)
    CASE ('cloud_ice')
      CALL species_from_ssdb(sp, name, TRIM(ssdb_dir)//'/IceSphere.nc', frequency, PSD_T08_ICE, &
                             n_streams=ns)
    CASE ('snow')
      CALL species_from_ssdb(sp, name, TRIM(ssdb_dir)//'/6BulletRosette.nc', frequency, PSD_F07_MIDLAT, &
                             n_streams=ns)
    CASE ('graupel')
      CALL species_from_ssdb(sp, name, TRIM(ssdb_dir)//'/GemGraupel.nc', frequency, PSD_F07_MIDLAT, &
                             n_streams=ns)
    CASE ('hail')
      CALL species_from_ssdb(sp, name, TRIM(ssdb_dir)//'/GemHail.nc', frequency, PSD_F07_MIDLAT, &
                             n_streams=ns)
    CASE DEFAULT
      WRITE(*,'(a)') 'species_default_ssdb: unknown species '//TRIM(name)
      ERROR STOP 1
    END SELECT
  END SUBROUTINE species_default_ssdb

  ! Trapezoid bin widths for sorted sizes
  SUBROUTINE bin_widths(d, dd)
    REAL(fp), INTENT(IN)  :: d(:)
    REAL(fp), INTENT(OUT) :: dd(:)
    INTEGER :: i, n
    n = SIZE(d)
    dd = ZERO
    IF (n == 1) THEN
      dd = ONE
      RETURN
    END IF
    DO i = 1, n-1
      dd(i)   = dd(i)   + 0.5_fp*(d(i+1) - d(i))
      dd(i+1) = dd(i+1) + 0.5_fp*(d(i+1) - d(i))
    END DO
  END SUBROUTINE bin_widths

  !
  ! species_from_ssdb
  !   Build a species from a NuRTM single-scattering database file (one habit,
  !   totally random orientation; variables frequency (Hz), temperature (K),
  !   d_veq, d_max (m), mass (kg), scat_angle (deg), ext, abs, bsc (m2) and
  !   pha (m2 sr-1, F11 is element 1)). The cross sections are interpolated
  !   linearly to the radar frequency, as NuRTM does, and kept on the
  !   temperature grid; the layer temperature interpolates between them.
  !   The asymmetry parameter is computed from F11.
  !
  !   psd_type   PSD_GAMMA (n0, mu), PSD_T08_RAIN, PSD_T08_SNOW,
  !              PSD_T08_SNOW_FIELD05, PSD_T08_GRAUPEL, PSD_T08_ICE (ni: fixed number
  !              concentration in m-3, default Cooper 1986), PSD_F07_MIDLAT,
  !              PSD_F07_TROPICAL (Field et al. 2007; t_fix = 235 reproduces the
  !              CRTM DDA tables), or PSD_GAMMA_FIXED_N0
  !   size_param 'd_max' (default) or 'd_veq': the size the PSD is defined on
  !   phase      default: liquid when the temperature grid goes above 273.15 K
  !
  SUBROUTINE species_from_ssdb(sp, name, filename, frequency, psd_type, n0, mu, size_param, &
                               n_streams, phase, approx_kind, rho_air, t_fix, ni)
    TYPE(radar_species_type), INTENT(OUT) :: sp
    CHARACTER(*),             INTENT(IN)  :: name, filename
    REAL(fp),                 INTENT(IN)  :: frequency
    INTEGER,                  INTENT(IN)  :: psd_type
    REAL(fp),       OPTIONAL, INTENT(IN)  :: n0, mu, rho_air, t_fix, ni
    CHARACTER(*),   OPTIONAL, INTENT(IN)  :: size_param
    INTEGER,        OPTIONAL, INTENT(IN)  :: n_streams, phase, approx_kind
    INTEGER :: ncid, varid, dimid, nf, nt, nd, na, i1, i2, k, j, ia, st
    REAL(fp), ALLOCATABLE :: freq(:), dveq(:), dmax(:), ang(:), e(:,:), a(:,:), b(:,:), f11(:,:,:)
    REAL(fp) :: fhz, wf, mu_, x, num, den, dth, xm, ym, sxy, sxx
    REAL(fp), ALLOCATABLE :: c1(:,:,:), lx(:), ly(:)
    LOGICAL :: have_bsc

    sp%name = name
    sp%source = SOURCE_PARTICLE
    sp%frequency = frequency
    sp%psd_type = psd_type
    IF (PRESENT(n0)) sp%n0 = n0
    mu_ = ZERO
    IF (PRESENT(mu)) mu_ = mu
    sp%mu = mu_
    IF (PRESENT(rho_air)) sp%rho_air = rho_air
    IF (PRESENT(t_fix)) sp%t_fix = t_fix
    IF (PRESENT(ni)) sp%ni_fixed = ni
    IF (PRESENT(n_streams)) sp%n_streams = n_streams
    sp%rayleigh_only = .FALSE.

    CALL nc_check(nf90_open(filename, NF90_NOWRITE, ncid), 'open '//TRIM(filename))
    CALL nc_check(nf90_inq_dimid(ncid, 'frequency', dimid), 'frequency')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=nf), 'frequency')
    CALL nc_check(nf90_inq_dimid(ncid, 'temperature', dimid), 'temperature')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=nt), 'temperature')
    CALL nc_check(nf90_inq_dimid(ncid, 'size', dimid), 'size')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=nd), 'size')
    CALL nc_check(nf90_inq_dimid(ncid, 'scat_angle', dimid), 'scat_angle')
    CALL nc_check(nf90_inquire_dimension(ncid, dimid, len=na), 'scat_angle')
    ALLOCATE(freq(nf), dveq(nd), dmax(nd), ang(na))
    ALLOCATE(sp%tgp(nt), sp%d(nd), sp%dd(nd), sp%mass(nd))
    ALLOCATE(sp%cext(nd,nt), sp%csca(nd,nt), sp%cbsc(nd,nt), sp%asym(nd,nt))
    CALL nc_check(nf90_inq_varid(ncid, 'frequency', varid), 'frequency')
    CALL nc_check(nf90_get_var(ncid, varid, freq), 'frequency')
    CALL nc_check(nf90_inq_varid(ncid, 'temperature', varid), 'temperature')
    CALL nc_check(nf90_get_var(ncid, varid, sp%tgp), 'temperature')
    CALL nc_check(nf90_inq_varid(ncid, 'd_veq', varid), 'd_veq')
    CALL nc_check(nf90_get_var(ncid, varid, dveq), 'd_veq')
    CALL nc_check(nf90_inq_varid(ncid, 'd_max', varid), 'd_max')
    CALL nc_check(nf90_get_var(ncid, varid, dmax), 'd_max')
    CALL nc_check(nf90_inq_varid(ncid, 'mass', varid), 'mass')
    CALL nc_check(nf90_get_var(ncid, varid, sp%mass), 'mass')
    CALL nc_check(nf90_inq_varid(ncid, 'scat_angle', varid), 'scat_angle')
    CALL nc_check(nf90_get_var(ncid, varid, ang), 'scat_angle')
    sp%ntp = nt
    sp%nd = nd

    ! Frequency bracket (linear interpolation, clamped)
    fhz = frequency*1.0e9_fp
    IF (fhz <= freq(1)) THEN
      i1 = 1; i2 = 1; wf = ZERO
    ELSE IF (fhz >= freq(nf)) THEN
      i1 = nf; i2 = nf; wf = ZERO
    ELSE
      DO i1 = 1, nf-1
        IF (freq(i1+1) >= fhz) EXIT
      END DO
      i2 = i1 + 1
      wf = (fhz - freq(i1))/(freq(i2) - freq(i1))
    END IF

    ALLOCATE(e(nt,nd), a(nt,nd), b(nt,nd), f11(nt,nd,na), c1(nt,nd,na))
    sp%cext = ZERO; sp%csca = ZERO; sp%cbsc = ZERO; sp%asym = ZERO
    st = nf90_inq_varid(ncid, 'bsc', varid)
    have_bsc = (st == NF90_NOERR)
    DO k = 1, 2
      IF (k == 1) THEN
        j = i1; x = ONE - wf
      ELSE
        j = i2; x = wf
      END IF
      CALL nc_check(nf90_inq_varid(ncid, 'ext', varid), 'ext')
      CALL nc_check(nf90_get_var(ncid, varid, e, start=[j,1,1], count=[1,nt,nd]), 'ext')
      CALL nc_check(nf90_inq_varid(ncid, 'abs', varid), 'abs')
      CALL nc_check(nf90_get_var(ncid, varid, a, start=[j,1,1], count=[1,nt,nd]), 'abs')
      CALL nc_check(nf90_inq_varid(ncid, 'pha', varid), 'pha')
      CALL nc_check(nf90_get_var(ncid, varid, c1, start=[j,1,1,1,1], count=[1,nt,nd,na,1]), 'pha')
      IF (have_bsc) THEN
        CALL nc_check(nf90_inq_varid(ncid, 'bsc', varid), 'bsc')
        CALL nc_check(nf90_get_var(ncid, varid, b, start=[j,1,1], count=[1,nt,nd]), 'bsc')
      ELSE
        b = 4.0_fp*ACOS(-ONE)*c1(:,:,na)
      END IF
      IF (k == 1) THEN
        f11 = x*c1
      ELSE
        f11 = f11 + x*c1
      END IF
      sp%cext = sp%cext + x*TRANSPOSE(e)
      sp%csca = sp%csca + x*TRANSPOSE(e - a)
      sp%cbsc = sp%cbsc + x*TRANSPOSE(b)
    END DO
    CALL nc_check(nf90_close(ncid), 'close')
    sp%csca = MAX(sp%csca, ZERO)
    ! Asymmetry parameter from F11 (trapezoid in scattering angle)
    DO k = 1, nt
      DO j = 1, nd
        num = ZERO; den = ZERO
        DO ia = 1, na-1
          dth = (ang(ia+1) - ang(ia))*ACOS(-ONE)/180.0_fp
          x = 0.5_fp*dth
          num = num + x*(f11(k,j,ia)*COSD(ang(ia))*SIND(ang(ia)) + f11(k,j,ia+1)*COSD(ang(ia+1))*SIND(ang(ia+1)))
          den = den + x*(f11(k,j,ia)*SIND(ang(ia)) + f11(k,j,ia+1)*SIND(ang(ia+1)))
        END DO
        IF (den > ZERO) sp%asym(j,k) = num/den
      END DO
    END DO

    ! PSD size parameter
    sp%d = dmax
    IF (PRESENT(size_param)) THEN
      IF (TRIM(size_param) == 'd_veq') sp%d = dveq
    END IF
    CALL bin_widths(sp%d, sp%dd)
    ! Mass-size fit m = a D^b (least squares in log space), used by PSD_GAMMA
    ALLOCATE(lx(nd), ly(nd))
    lx = LOG(sp%d); ly = LOG(sp%mass)
    xm = SUM(lx)/nd; ym = SUM(ly)/nd
    sxy = SUM((lx - xm)*(ly - ym)); sxx = SUM((lx - xm)**2)
    IF (nd > 1 .AND. sxx > ZERO) THEN
      sp%b_ms = sxy/sxx
    ELSE
      sp%b_ms = 3.0_fp
    END IF
    sp%a_ms = EXP(ym - sp%b_ms*xm)

    sp%phase = PHASE_FROZEN
    IF (MAXVAL(sp%tgp) > 273.15_fp) sp%phase = PHASE_LIQUID
    IF (PRESENT(phase)) sp%phase = phase
    IF (PRESENT(approx_kind)) THEN
      sp%approx_kind = approx_kind
    ELSE IF (sp%phase == PHASE_LIQUID .AND. psd_type /= PSD_T08_RAIN .AND. &
             (psd_type == PSD_GAMMA .OR. psd_type == PSD_GAMMA_FIXED_N0) .AND. sp%mu > ZERO) THEN
      sp%approx_kind = APPROX_CLOUD_LIQUID
    ELSE IF (sp%phase == PHASE_LIQUID) THEN
      sp%approx_kind = APPROX_RAIN
    ELSE
      sp%approx_kind = APPROX_ICE
    END IF
  END SUBROUTINE species_from_ssdb

  SUBROUTINE species_destroy(sp)
    TYPE(radar_species_type), INTENT(INOUT) :: sp
    IF (ALLOCATED(sp%x))     DEALLOCATE(sp%x)
    IF (ALLOCATED(sp%tgrid)) DEALLOCATE(sp%tgrid)
    IF (ALLOCATED(sp%ke))    DEALLOCATE(sp%ke, sp%kb, sp%w, sp%pc1, sp%pcn)
    IF (ALLOCATED(sp%d))     DEALLOCATE(sp%d, sp%dd, sp%mass, sp%cext, sp%csca, sp%cbsc, sp%asym)
    IF (ALLOCATED(sp%tgp))   DEALLOCATE(sp%tgp)
    sp%source = SOURCE_NONE
  END SUBROUTINE species_destroy

  !
  ! Optical properties in complex arithmetic.
  !   rho  water content (kg m-3), t temperature (K)
  !   opt(N_OPT) volume coefficients (m-1), see I_* indices
  !
  SUBROUTINE species_optics_c(sp, rho, t, opt)
    TYPE(radar_species_type), INTENT(IN)  :: sp
    COMPLEX(fp),              INTENT(IN)  :: rho, t
    COMPLEX(fp),              INTENT(OUT) :: opt(N_OPT)
    COMPLEX(fp) :: wx(4), wt(4), ke, kb, w, pc1, pcn
    INTEGER :: ix, it, j, k
    opt = (0.0_fp, 0.0_fp)
    IF (REAL(rho) <= 0.0_fp) RETURN
    SELECT CASE (sp%source)
    CASE (SOURCE_LUT)
      IF (sp%rayleigh_only) THEN
        ix = 1
        wx = (0.0_fp, 0.0_fp); wx(1) = (1.0_fp, 0.0_fp)
      ELSE
        CALL interp_weights_c(sp%x, rho, ix, wx)
      END IF
      IF (sp%nt > 1) THEN
        CALL interp_weights_c(sp%tgrid, t, it, wt)
      ELSE
        it = 1
        wt = (0.0_fp, 0.0_fp); wt(1) = (1.0_fp, 0.0_fp)
      END IF
      ke = 0; kb = 0; w = 0; pc1 = 0; pcn = 0
      DO k = 1, MIN(4, sp%nt)
        DO j = 1, 4
          ke  = ke  + wx(j)*wt(k)*sp%ke (ix+j-1, it+k-1)
          kb  = kb  + wx(j)*wt(k)*sp%kb (ix+j-1, it+k-1)
          w   = w   + wx(j)*wt(k)*sp%w  (ix+j-1, it+k-1)
          pc1 = pc1 + wx(j)*wt(k)*sp%pc1(ix+j-1, it+k-1)
          pcn = pcn + wx(j)*wt(k)*sp%pcn(ix+j-1, it+k-1)
        END DO
      END DO
      ! CRTM_CloudScatter quality control
      IF (REAL(ke) <= 0.0_fp) THEN
        ke = 0; kb = 0; w = 0
      END IF
      IF (REAL(w) <= 0.0_fp) THEN
        w = 0; pc1 = 0; pcn = 0
        IF (sp%zero_bsc_if_no_scattering) kb = 0
      END IF
      IF (REAL(w) >= 1.0_fp) w = 1
      opt(I_EXT)  = ke*rho
      opt(I_BSC)  = kb*rho
      opt(I_SCA)  = w*ke*rho
      opt(I_FSCA) = pcn*opt(I_SCA)
      opt(I_GSCA) = pc1/1.5_fp*opt(I_SCA)
    CASE (SOURCE_PARTICLE)
      CALL particle_optics_c(sp, rho, t, opt)
    END SELECT
    opt(I_ATT) = hydro_approx_attenuation_c(sp%approx_kind, sp%frequency, rho, t)
  END SUBROUTINE species_optics_c

  ! PSD integration in complex arithmetic (complex-step derivatives).
  ! Cross sections are interpolated linearly in temperature (clamped).
  SUBROUTINE particle_optics_c(sp, rho, t, opt)
    TYPE(radar_species_type), INTENT(IN)  :: sp
    COMPLEX(fp),              INTENT(IN)  :: rho, t
    COMPLEX(fp),              INTENT(OUT) :: opt(N_OPT)
    COMPLEX(fp) :: nd(sp%nd), sca, gs, wt
    COMPLEX(fp) :: ce(sp%nd), cs(sp%nd), cb(sp%nd), cg(sp%nd)
    INTEGER :: k1, k2
    CALL psd_numbers_c(sp, rho, t, nd)
    ! temperature weights
    IF (sp%ntp == 1) THEN
      k1 = 1; k2 = 1; wt = 0
    ELSE IF (REAL(t) <= sp%tgp(1)) THEN
      k1 = 1; k2 = 1; wt = 0
    ELSE IF (REAL(t) >= sp%tgp(sp%ntp)) THEN
      k1 = sp%ntp; k2 = sp%ntp; wt = 0
    ELSE
      DO k1 = 1, sp%ntp-1
        IF (sp%tgp(k1+1) >= REAL(t)) EXIT
      END DO
      k2 = k1 + 1
      wt = (t - sp%tgp(k1))/(sp%tgp(k2) - sp%tgp(k1))
    END IF
    ce = (1 - wt)*sp%cext(:,k1) + wt*sp%cext(:,k2)
    cs = (1 - wt)*sp%csca(:,k1) + wt*sp%csca(:,k2)
    cb = (1 - wt)*sp%cbsc(:,k1) + wt*sp%cbsc(:,k2)
    cg = (1 - wt)*sp%asym(:,k1)*sp%csca(:,k1) + wt*sp%asym(:,k2)*sp%csca(:,k2)
    opt(I_EXT) = SUM(ce*nd)
    sca        = SUM(cs*nd)
    opt(I_SCA) = sca
    opt(I_BSC) = SUM(cb*nd)
    gs         = SUM(cg*nd)
    opt(I_GSCA) = gs
    ! Delta truncation with a Henyey-Greenstein phase function, written in
    ! the CRTM convention for the Legendre coefficients ((2l+1)/2 g^l)
    IF (REAL(sca) > 0.0_fp) THEN
      opt(I_FSCA) = (2*sp%n_streams + 1)/2.0_fp*(gs/sca)**sp%n_streams*sca
    ELSE
      opt(I_FSCA) = 0
    END IF
  END SUBROUTINE particle_optics_c

  ! Complex min/max that pick a branch by the real part
  ELEMENTAL FUNCTION cmin(a, b) RESULT(c)
    COMPLEX(fp), INTENT(IN) :: a, b
    COMPLEX(fp) :: c
    IF (REAL(a) <= REAL(b)) THEN
      c = a
    ELSE
      c = b
    END IF
  END FUNCTION cmin
  ELEMENTAL FUNCTION cmax(a, b) RESULT(c)
    COMPLEX(fp), INTENT(IN) :: a, b
    COMPLEX(fp) :: c
    IF (REAL(a) >= REAL(b)) THEN
      c = a
    ELSE
      c = b
    END IF
  END FUNCTION cmax

  ! Field et al. (2005) moment relation M_n = a(n,Tc) M2^b(n,Tc)
  FUNCTION field05_moment(n, tc, m2) RESULT(mn)
    REAL(fp),    INTENT(IN) :: n
    COMPLEX(fp), INTENT(IN) :: tc, m2
    COMPLEX(fp) :: mn
    mn = 10.0_fp**poly(F05_SA)*m2**poly(F05_SB)
  CONTAINS
    COMPLEX(fp) FUNCTION poly(c)
      REAL(fp), INTENT(IN) :: c(10)
      poly = c(1) + c(2)*tc + c(3)*n + c(4)*tc*n + c(5)*tc**2 + c(6)*n**2 + &
             c(7)*tc**2*n + c(8)*tc*n**2 + c(9)*tc**3 + c(10)*n**3
    END FUNCTION poly
  END FUNCTION field05_moment

  ! Field et al. (2007) moment relation M_n = A(n) exp(B(n) Tc) M2^C(n):
  ! which = 1 returns A exp(B Tc), which = 3 returns C
  FUNCTION f07_coef(n, tc, which) RESULT(v)
    REAL(fp),    INTENT(IN) :: n
    COMPLEX(fp), INTENT(IN) :: tc
    INTEGER,     INTENT(IN) :: which
    COMPLEX(fp) :: v
    REAL(fp) :: a, b, c
    a = EXP(13.6078_fp - 7.76062_fp*n + 0.478694_fp*n*n)
    b = -0.0360722_fp + 0.0150830_fp*n + 0.00149453_fp*n*n
    c = 0.806856_fp + 0.00581022_fp*n + 0.0456723_fp*n*n
    IF (which == 1) THEN
      v = a*EXP(b*tc)
    ELSE
      v = CMPLX(c, ZERO, fp)
    END IF
  END FUNCTION f07_coef

  ! Number concentration per size bin (m-3) holding the water content rho
  SUBROUTINE psd_numbers_c(sp, rho, t, n)
    TYPE(radar_species_type), INTENT(IN)  :: sp
    COMPLEX(fp),              INTENT(IN)  :: rho, t
    COMPLEX(fp),              INTENT(OUT) :: n(sp%nd)
    REAL(fp), PARAMETER :: PI = 3.14159265358979323846_fp
    REAL(fp) :: lam, dmdl, peak
    COMPLEX(fp) :: lamc, n0, am, tk, tc, m2, m3, r, mtot
    REAL(fp) :: bm
    n = 0
    IF (REAL(rho) <= 0.0_fp) RETURN
    tk = t
    IF (sp%t_fix > ZERO) tk = CMPLX(sp%t_fix, ZERO, fp)
    SELECT CASE (sp%psd_type)
    CASE (PSD_GAMMA_FIXED_N0)
      ! lambda is solved in real arithmetic, and the complex part of rho is
      ! carried through d(lambda)/d(rho) = 1/(dM/dlambda)
      CALL solve_lambda(sp, REAL(rho), lam, dmdl)
      lamc = CMPLX(lam, AIMAG(rho)/dmdl, fp)
      n = sp%n0*sp%d**sp%mu*EXP(-lamc*sp%d)*sp%dd
      RETURN
    CASE (PSD_GAMMA)
      lamc = (sp%n0*sp%a_ms*GAMMA(sp%b_ms + sp%mu + 1)/rho)**(1/(sp%b_ms + sp%mu + 1))
      n = sp%n0*sp%d**sp%mu*EXP(-lamc*sp%d)
    CASE (PSD_T08_RAIN, PSD_T08_SNOW, PSD_T08_GRAUPEL)
      bm = 3.0_fp
      IF (sp%psd_type == PSD_T08_RAIN) THEN
        n0 = 8.0e6_fp; am = PI/6*1000.0_fp
      ELSE IF (sp%psd_type == PSD_T08_SNOW) THEN
        n0 = cmin(CMPLX(2.0e8_fp, ZERO, fp), 2.0e6_fp*EXP(-0.12_fp*cmin(CMPLX(-0.001_fp, ZERO, fp), tk - 273.15_fp)))
        am = PI/6*100.0_fp
      ELSE
        n0 = cmax(CMPLX(1.0e4_fp, ZERO, fp), cmin(200.0_fp/(rho/sp%rho_air), CMPLX(5.0e6_fp, ZERO, fp)))
        am = PI/6*400.0_fp
      END IF
      lamc = (n0*am*GAMMA(bm + 1)/rho)**(1/(bm + 1))
      n = n0*EXP(-lamc*sp%d)
    CASE (PSD_T08_ICE)
      ! Thompson et al. (2008) cloud ice: exponential in D (mu = 0) with
      ! m(D) = pi/6 * 890 D^3. Ni is prognostic in Thompson; here it is fixed or
      ! diagnosed with Cooper (1986), Ni = 5 exp(0.304 (273.15 - T)) m-3,
      ! limited to 250 L-1 as in the Thompson ice nucleation. The slope is held
      ! within the Thompson limits and the PSD renormalised to the water content.
      IF (sp%ni_fixed > ZERO) THEN
        n0 = sp%ni_fixed
      ELSE
        n0 = cmin(CMPLX(250.0e3_fp, ZERO, fp), 5.0_fp*EXP(0.304_fp*(273.15_fp - cmin(tk, CMPLX(273.15_fp, ZERO, fp)))))
      END IF
      am = PI/6*890.0_fp
      lamc = (am*n0*GAMMA(4.0_fp)/rho)**(1.0_fp/3.0_fp)
      ! Thompson limits on the ice slope: 1/(2 D0i + 100 um) <= lambda <= 1/1 um, D0i = 300 um
      IF (REAL(lamc) < 1.0_fp/(2*300.0e-6_fp + 100.0e-6_fp)) lamc = 1.0_fp/(2*300.0e-6_fp + 100.0e-6_fp)
      IF (REAL(lamc) > 1.0e6_fp) lamc = 1.0e6_fp
      n = EXP(-lamc*sp%d)
    CASE (PSD_F07_MIDLAT, PSD_F07_TROPICAL)
      ! Field et al. (2007) moment-based ice/snow PSD with the database mass-size
      ! fit m = a D^b (as used to build the CRTM DDA tables, pyWeather
      ! psd_field.py; those tables use a fixed T = 235 K: t_fix).
      ! Valid for -70 < Tc < 0 C; the temperature is held in that range.
      tc = tk - 273.0_fp
      tc = cmax(cmin(tc, CMPLX(-0.1_fp, ZERO, fp)), CMPLX(-70.0_fp, ZERO, fp))
      m3 = rho/sp%a_ms                                   ! M_beta
      IF (ABS(sp%b_ms - 2.0_fp) < 1.0e-12_fp) THEN
        m2 = m3
      ELSE
        m2 = (m3/f07_coef(sp%b_ms, tc, 1))**(1.0_fp/f07_coef(sp%b_ms, tc, 3))
      END IF
      m3 = f07_coef(3.0_fp, tc, 1)*m2**f07_coef(3.0_fp, tc, 3)
      r = m2/m3
      IF (sp%psd_type == PSD_F07_MIDLAT) THEN
        n = 141.0_fp*EXP(-16.8_fp*r*sp%d) + 102.0_fp*(r*sp%d)**2.07_fp*EXP(-4.82_fp*r*sp%d)
      ELSE
        n = 152.0_fp*EXP(-12.4_fp*r*sp%d) + 3.28_fp*(r*sp%d)**(-0.78_fp)*EXP(-1.94_fp*r*sp%d)
      END IF
      WHERE (REAL(r*sp%d) >= 20.0_fp) n = 0
      n = n*m2**4/m3**3
    CASE (PSD_T08_SNOW_FIELD05)
      tc = cmin(CMPLX(-0.1_fp, ZERO, fp), tk - 273.15_fp)
      m2 = rho/0.069_fp
      m3 = field05_moment(3.0_fp, tc, m2)
      r = m2/m3
      n = m2**4/m3**3*(0.490_fp*EXP(-r*20.78_fp*sp%d) + 17.46_fp*r**0.6357_fp*sp%d**0.6357_fp*EXP(-r*3.29_fp*sp%d))
    END SELECT
    ! Discretise and renormalise to the database masses (NuRTM psd_numbers)
    peak = MAXVAL(REAL(n))
    IF (peak > ZERO .AND. peak < HUGE(peak)) THEN
      n = n*sp%dd/peak
    ELSE
      n = n*sp%dd
    END IF
    mtot = SUM(n*sp%mass)
    IF (.NOT. (REAL(mtot) > ZERO) .OR. .NOT. (ABS(REAL(rho/mtot)) < HUGE(peak))) THEN
      n = 0
      n(1) = rho/sp%mass(1)
      RETURN
    END IF
    n = n*(rho/mtot)
  END SUBROUTINE psd_numbers_c

  ! Solve M(lambda) = rho for the gamma PSD slope (Newton on log(lambda))
  SUBROUTINE solve_lambda(sp, rho, lam, dmdl)
    TYPE(radar_species_type), INTENT(IN)  :: sp
    REAL(fp),                 INTENT(IN)  :: rho
    REAL(fp),                 INTENT(OUT) :: lam, dmdl
    REAL(fp) :: x, xlo, xhi, m, nd(sp%nd), f, fx
    INTEGER :: it
    ! Bracket in log(lambda) (lambda between 1 and 1e8 m-1); M decreases with lambda
    xlo = 0.0_fp; xhi = LOG(1.0e8_fp)
    x = 0.5_fp*(xlo + xhi)
    DO it = 1, 200
      lam = EXP(x)
      nd = sp%n0*sp%d**sp%mu*EXP(-lam*sp%d)*sp%dd
      m = SUM(sp%mass*nd)
      dmdl = -SUM(sp%mass*sp%d*nd)
      f = LOG(MAX(m, TINY(m))) - LOG(rho)
      IF (ABS(f) < 1.0e-13_fp) EXIT
      IF (f > 0.0_fp) THEN
        xlo = x
      ELSE
        xhi = x
      END IF
      ! Newton step in log space: d(log M)/d(log lambda) = lambda*dM/dlambda/M
      fx = lam*dmdl/MAX(m, TINY(m))
      IF (fx < 0.0_fp) THEN
        x = x - f/fx
      ELSE
        x = 0.5_fp*(xlo + xhi)
      END IF
      IF (x <= xlo .OR. x >= xhi) x = 0.5_fp*(xlo + xhi)
    END DO
    lam = EXP(x)
    nd = sp%n0*sp%d**sp%mu*EXP(-lam*sp%d)*sp%dd
    dmdl = -SUM(sp%mass*sp%d*nd)
  END SUBROUTINE solve_lambda

  ! Real-valued optics
  SUBROUTINE species_optics(sp, rho, t, opt)
    TYPE(radar_species_type), INTENT(IN)  :: sp
    REAL(fp),                 INTENT(IN)  :: rho, t
    REAL(fp),                 INTENT(OUT) :: opt(N_OPT)
    COMPLEX(fp) :: oc(N_OPT)
    CALL species_optics_c(sp, CMPLX(rho, ZERO, fp), CMPLX(t, ZERO, fp), oc)
    opt = REAL(oc)
  END SUBROUTINE species_optics

  ! Optics and derivatives with respect to rho and t
  SUBROUTINE species_optics_deriv(sp, rho, t, opt, dopt_drho, dopt_dt)
    TYPE(radar_species_type), INTENT(IN)  :: sp
    REAL(fp),                 INTENT(IN)  :: rho, t
    REAL(fp),                 INTENT(OUT) :: opt(N_OPT), dopt_drho(N_OPT), dopt_dt(N_OPT)
    COMPLEX(fp) :: oc(N_OPT)
    REAL(fp) :: hr, ht
    IF (rho <= ZERO) THEN
      opt = ZERO; dopt_drho = ZERO; dopt_dt = ZERO
      RETURN
    END IF
    hr = CSTEP*rho
    ht = CSTEP*t
    CALL species_optics_c(sp, CMPLX(rho, hr, fp), CMPLX(t, ZERO, fp), oc)
    opt = REAL(oc)
    dopt_drho = AIMAG(oc)/hr
    CALL species_optics_c(sp, CMPLX(rho, ZERO, fp), CMPLX(t, ht, fp), oc)
    dopt_dt = AIMAG(oc)/ht
  END SUBROUTINE species_optics_deriv

END MODULE radar_species
