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
  PUBLIC :: PSD_GAMMA_FIXED_N0
  PUBLIC :: N_OPT, I_EXT, I_SCA, I_BSC, I_FSCA, I_GSCA, I_ATT
  PUBLIC :: species_from_lut
  PUBLIC :: species_from_particles
  PUBLIC :: species_optics
  PUBLIC :: species_optics_deriv
  PUBLIC :: species_destroy
  PUBLIC :: lut_habit_names

  INTEGER, PARAMETER :: SOURCE_NONE = 0, SOURCE_LUT = 1, SOURCE_PARTICLE = 2
  INTEGER, PARAMETER :: PHASE_LIQUID = 0, PHASE_FROZEN = 1
  INTEGER, PARAMETER :: PSD_GAMMA_FIXED_N0 = 1

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
    INTEGER  :: n_streams = 4
    REAL(fp), ALLOCATABLE :: d(:)          ! particle size (m)
    REAL(fp), ALLOCATABLE :: dd(:)         ! quadrature weights (m)
    REAL(fp), ALLOCATABLE :: mass(:)       ! particle mass (kg)
    REAL(fp), ALLOCATABLE :: cext(:), csca(:), cbsc(:), asym(:)   ! m2, m2, m2, -
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
    INTEGER :: n, i
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
    ALLOCATE(sp%d(n), sp%dd(n), sp%mass(n), sp%cext(n), sp%csca(n), sp%cbsc(n), sp%asym(n))
    sp%d = d; sp%mass = mass; sp%cext = cext; sp%csca = csca; sp%cbsc = cbsc; sp%asym = asym
    ! Trapezoidal weights
    sp%dd = ZERO
    DO i = 1, n-1
      sp%dd(i)   = sp%dd(i)   + 0.5_fp*(d(i+1) - d(i))
      sp%dd(i+1) = sp%dd(i+1) + 0.5_fp*(d(i+1) - d(i))
    END DO
  END SUBROUTINE species_from_particles

  SUBROUTINE species_destroy(sp)
    TYPE(radar_species_type), INTENT(INOUT) :: sp
    IF (ALLOCATED(sp%x))     DEALLOCATE(sp%x)
    IF (ALLOCATED(sp%tgrid)) DEALLOCATE(sp%tgrid)
    IF (ALLOCATED(sp%ke))    DEALLOCATE(sp%ke, sp%kb, sp%w, sp%pc1, sp%pcn)
    IF (ALLOCATED(sp%d))     DEALLOCATE(sp%d, sp%dd, sp%mass, sp%cext, sp%csca, sp%cbsc, sp%asym)
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
      CALL particle_optics_c(sp, rho, opt)
    END SELECT
    opt(I_ATT) = hydro_approx_attenuation_c(sp%approx_kind, sp%frequency, rho, t)
  END SUBROUTINE species_optics_c

  ! PSD integration. lambda is solved in real arithmetic, and the complex
  ! part of rho is carried through d(lambda)/d(rho) = 1/(dM/dlambda).
  SUBROUTINE particle_optics_c(sp, rho, opt)
    TYPE(radar_species_type), INTENT(IN)  :: sp
    COMPLEX(fp),              INTENT(IN)  :: rho
    COMPLEX(fp),              INTENT(OUT) :: opt(N_OPT)
    REAL(fp) :: lam, dmdl
    COMPLEX(fp) :: lamc, nd(sp%nd), sca, gs
    CALL solve_lambda(sp, REAL(rho), lam, dmdl)
    lamc = CMPLX(lam, AIMAG(rho)/dmdl, fp)
    nd = sp%n0*sp%d**sp%mu*EXP(-lamc*sp%d)*sp%dd
    opt(I_EXT) = SUM(sp%cext*nd)
    sca        = SUM(sp%csca*nd)
    opt(I_SCA) = sca
    opt(I_BSC) = SUM(sp%cbsc*nd)
    gs         = SUM(sp%asym*sp%csca*nd)
    opt(I_GSCA) = gs
    ! Delta truncation with a Henyey-Greenstein phase function, written in
    ! the CRTM convention for the Legendre coefficients ((2l+1)/2 g^l)
    IF (REAL(sca) > 0.0_fp) THEN
      opt(I_FSCA) = (2*sp%n_streams + 1)/2.0_fp*(gs/sca)**sp%n_streams*sca
    ELSE
      opt(I_FSCA) = 0
    END IF
  END SUBROUTINE particle_optics_c

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
