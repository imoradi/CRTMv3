!
! compare_crtm
!
! Runs the CRTM active-sensor operator and the standalone radar operator on
! the same IFS profiles and writes both to a netCDF file.
!
! Usage: compare_crtm ifs.nc stride zenith_deg out.nc coef_dir/ cloudcoeff_dir/ cloudcoeff_file n_jac nurtm_ssdb_dir [snow_habit]
!   snow_habit: e.g. 6BulletRosette (CRTM and LUT variants); default CRTM SNOW_CLOUD (SectorSnowflake)
!
! Standalone variants (per channel):
!   A  CRTM-compatible settings, gas optical depth taken from CRTM (ODPS)
!   B  CRTM-compatible settings, gas from ITU-R P.676 line-by-line
!   C  approximate attenuation (ITU-R P.676 approx gas, P.838/P.840 hydrometeors)
!   D  physical: full (unscaled) extinction, cloud backscatter included, ITU LBL gas
!   E  NuRTM DDA single particles + Thompson (2008) PSDs, delta scaling, ITU LBL gas
!   F  as E with full (unscaled) extinction
!
PROGRAM compare_crtm
  USE CRTM_Module
  USE CRTM_SpcCoeff,           ONLY: SC
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, CRTM_GeometryInfo_SetValue
  USE CRTM_GeometryInfo,       ONLY: CRTM_GeometryInfo_Compute
  USE CRTM_Predictor_Define,   ONLY: CRTM_Predictor_type, CRTM_Predictor_Create, CRTM_Predictor_Destroy
  USE CRTM_Predictor,          ONLY: CRTM_PVar_type => iVar_type, CRTM_Compute_Predictors
  USE CRTM_AtmAbsorption,      ONLY: CRTM_AAvar_type => iVar_type, CRTM_Compute_AtmAbsorption
  USE CRTM_AtmOptics_Define,   ONLY: CRTM_AtmOptics_type, CRTM_AtmOptics_Create, CRTM_AtmOptics_Destroy, &
                                     CRTM_AtmOptics_Zero
  USE CRTM_AncillaryInput_Define, ONLY: CRTM_AncillaryInput_type
  USE radar_kinds,    ONLY: rfp => fp
  USE radar_gas,      ONLY: GAS_MODEL_ITU_LBL, GAS_MODEL_ITU_APPROX, GAS_MODEL_USER
  USE radar_species,  ONLY: radar_species_type, species_from_lut, species_from_ssdb, &
                            PSD_GAMMA, PSD_T08_RAIN, PSD_T08_SNOW_FIELD05, PSD_T08_ICE, PSD_T08_SNOW
  USE radar_operator
  USE netcdf
  IMPLICIT NONE

  INTEGER, PARAMETER :: NV = 6      ! standalone variants
  INTEGER, PARAMETER :: NSP = 4     ! species: cloud liquid, cloud ice, rain, snow
  CHARACTER(16), PARAMETER :: SPNAME(NSP) = [CHARACTER(16) :: 'cloud_liquid', 'cloud_ice', 'rain', 'snow']
  INTEGER :: CTYPE(NSP) = [WATER_CLOUD, ICE_CLOUD, RAIN_CLOUD, SNOW_CLOUD]
  CHARACTER(32) :: snow_habit = ''
  REAL(fp), PARAMETER :: G0 = 9.80665_fp
  INTEGER, PARAMETER :: CHUNK = 200

  CHARACTER(512) :: ifs_file, out_file, coef_dir, cc_dir, cc_file, arg, ssdb_dir
  INTEGER :: stride, nobs, nlev, nlay, nsel, ncid, varid, dimid, i, j, k, m, ic, iv, s, n_jac
  REAL(fp) :: zenith
  REAL(fp), ALLOCATABLE :: lp(:,:), gp(:,:), tt(:,:), qq(:,:), o3(:,:), hyd(:,:,:)
  INTEGER,  ALLOCATABLE :: sel(:)
  TYPE(CRTM_ChannelInfo_type) :: chinfo(2)
  INTEGER :: nch_tot, err
  REAL(fp) :: freq(3)
  INTEGER :: sensor_of(3), chan_of(3)
  ! Outputs (layer, obs, channel)
  REAL(4), ALLOCATABLE :: c_z(:,:,:), c_zatt(:,:,:), c_tgas(:,:,:)
  REAL(4), ALLOCATABLE :: s_z(:,:,:,:), s_zatt(:,:,:,:), s_tcum(:,:,:,:), s_thyd(:,:,:,:)
  REAL(4), ALLOCATABLE :: o_hgt(:,:), o_t(:,:), o_wc(:,:,:)
  TYPE(radar_operator_type) :: ops(3, NV)
  INTEGER :: id, v(20), nv_

  CALL GET_COMMAND_ARGUMENT(1, ifs_file)
  CALL GET_COMMAND_ARGUMENT(2, arg); READ(arg,*) stride
  CALL GET_COMMAND_ARGUMENT(3, arg); READ(arg,*) zenith
  CALL GET_COMMAND_ARGUMENT(4, out_file)
  CALL GET_COMMAND_ARGUMENT(5, coef_dir)
  CALL GET_COMMAND_ARGUMENT(6, cc_dir)
  CALL GET_COMMAND_ARGUMENT(7, cc_file)
  n_jac = 0
  CALL GET_COMMAND_ARGUMENT(8, arg)
  IF (LEN_TRIM(arg) > 0) READ(arg,*) n_jac
  CALL GET_COMMAND_ARGUMENT(9, ssdb_dir)
  ! Optional snow habit for CRTM and the LUT variants (default: CRTM SNOW_CLOUD = SectorSnowflake)
  CALL GET_COMMAND_ARGUMENT(10, snow_habit)
  IF (TRIM(snow_habit) == '6BulletRosette') CTYPE(4) = SixBulletRosette

  ! ---------------------------------------------------------------- read IFS
  CALL chk(nf90_open(ifs_file, NF90_NOWRITE, ncid))
  CALL chk(nf90_inq_dimid(ncid, 'obs_id', dimid)); CALL chk(nf90_inquire_dimension(ncid, dimid, len=nobs))
  CALL chk(nf90_inq_dimid(ncid, 'level', dimid));  CALL chk(nf90_inquire_dimension(ncid, dimid, len=nlev))
  nlay = nlev - 1
  nsel = (nobs - 1)/stride + 1
  ALLOCATE(sel(nsel))
  sel = [(1 + (i-1)*stride, i = 1, nsel)]
  ALLOCATE(lp(nobs,nlev), gp(nobs,nlev), tt(nobs,nlev), qq(nobs,nlev), o3(nobs,nlev), hyd(nobs,nlev,NSP))
  CALL getv('level_pressure', lp); CALL getv('GEOPOTENTIAL', gp); CALL getv('t', tt)
  CALL getv('q', qq); CALL getv('o3', o3)
  CALL getv('clwc', hyd(:,:,1)); CALL getv('ciwc', hyd(:,:,2))
  CALL getv('crwc', hyd(:,:,3)); CALL getv('cswc', hyd(:,:,4))
  CALL chk(nf90_close(ncid))
  WRITE(*,'(a,i0,a,i0,a,i0)') 'profiles: ', nobs, '  selected: ', nsel, '  layers: ', nlay

  ! ---------------------------------------------------------------- CRTM init
  err = CRTM_Init([CHARACTER(16) :: 'dpr_gpm', 'cpr_cloudsat'], chinfo, &
                  Aerosol_Model='CRTM', AerosolCoeff_Format='Binary', AerosolCoeff_File='AerosolCoeff.bin', &
                  Cloud_Model='CRTM', CloudCoeff_Format='netCDF', CloudCoeff_File=TRIM(cc_file), &
                  File_Path=TRIM(coef_dir), NC_File_Path=TRIM(cc_dir), Quiet=.TRUE.)
  IF (err /= SUCCESS) STOP 'CRTM_Init failed'
  nch_tot = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
  IF (nch_tot /= 3) STOP 'expected 3 channels'
  sensor_of = [1, 1, 2]; chan_of = [1, 2, 1]
  DO ic = 1, 3
    freq(ic) = SC(chinfo(sensor_of(ic))%Sensor_Index)%Frequency(chan_of(ic))
    WRITE(*,'(a,i0,a,a,a,f8.3,a,l1)') 'channel ', ic, ' ', TRIM(chinfo(sensor_of(ic))%Sensor_Id), &
          ' frequency ', freq(ic), ' active=', SC(chinfo(sensor_of(ic))%Sensor_Index)%Is_Active_Sensor
  END DO

  ! ---------------------------------------------------------------- standalone init
  DO ic = 1, 3
    DO iv = 1, NV
      CALL setup_variant(ops(ic,iv), iv, freq(ic))
    END DO
  END DO

  ALLOCATE(c_z(nlay,nsel,3), c_zatt(nlay,nsel,3), c_tgas(nlay,nsel,3))
  ALLOCATE(s_z(nlay,nsel,3,NV), s_zatt(nlay,nsel,3,NV), s_tcum(nlay,nsel,3,NV), s_thyd(nlay,nsel,3,NV))
  ALLOCATE(o_hgt(nlay,nsel), o_t(nlay,nsel), o_wc(nlay,nsel,NSP))

  CALL run_all()
  CALL write_output()
  IF (n_jac > 0) CALL jacobians(n_jac)
  err = CRTM_Destroy(chinfo)

CONTAINS

  SUBROUTINE chk(st)
    INTEGER, INTENT(IN) :: st
    IF (st /= NF90_NOERR) THEN
      WRITE(*,*) TRIM(nf90_strerror(st)); STOP 1
    END IF
  END SUBROUTINE chk

  SUBROUTINE getv(name, a)
    CHARACTER(*), INTENT(IN) :: name
    REAL(fp), INTENT(OUT) :: a(:,:)
    CALL chk(nf90_inq_varid(ncid, name, varid)); CALL chk(nf90_get_var(ncid, varid, a))
  END SUBROUTINE getv

  SUBROUTINE setup_variant(op, iv, f)
    TYPE(radar_operator_type), INTENT(OUT) :: op
    INTEGER, INTENT(IN) :: iv
    REAL(fp), INTENT(IN) :: f
    TYPE(radar_config_type) :: c
    TYPE(radar_species_type) :: sp
    INTEGER :: s
    c%frequency = f
    c%platform = PLATFORM_SPACEBORNE
    SELECT CASE (iv)
    CASE (1); c%gas_model = GAS_MODEL_USER
    CASE (2); c%gas_model = GAS_MODEL_ITU_LBL
    CASE (3); c%gas_model = GAS_MODEL_ITU_APPROX; c%hydro_attenuation = ATTEN_APPROX
    CASE (4); c%gas_model = GAS_MODEL_ITU_LBL; c%scaling = SCALING_NONE; c%crtm_compat = .FALSE.
    CASE (5); c%gas_model = GAS_MODEL_ITU_LBL
    CASE (6); c%gas_model = GAS_MODEL_ITU_LBL; c%scaling = SCALING_NONE
    END SELECT
    CALL radar_init(op, c, NSP)
    DO s = 1, NSP
      IF (iv >= 5) THEN
        ! NuRTM single-particle databases with the PSDs of the NuRTM validation runs
        SELECT CASE (s)
        CASE (1)
          CALL species_from_ssdb(sp, 'cloud_liquid', TRIM(ssdb_dir)//'/LiquidSphere.nc', f, &
                                 PSD_GAMMA, n0=3.0e35_fp, mu=4.0_fp, size_param='d_veq')
        CASE (2)
          CALL species_from_ssdb(sp, 'cloud_ice', TRIM(ssdb_dir)//'/IceSphere.nc', f, &
                                 PSD_T08_ICE)
        CASE (3)
          CALL species_from_ssdb(sp, 'rain', TRIM(ssdb_dir)//'/LiquidSphere.nc', f, &
                                 PSD_T08_RAIN, size_param='d_veq')
        CASE (4)
          ! RADAR_SNOW_PSD=exponential selects the Thompson (2004) exponential snow PSD
          CALL GET_ENVIRONMENT_VARIABLE('RADAR_SNOW_PSD', arg)
          IF (TRIM(arg) == 'exponential') THEN
            CALL species_from_ssdb(sp, 'snow', TRIM(ssdb_dir)//'/6BulletRosette.nc', f, PSD_T08_SNOW)
          ELSE
            CALL species_from_ssdb(sp, 'snow', TRIM(ssdb_dir)//'/6BulletRosette.nc', f, &
                                   PSD_T08_SNOW_FIELD05)
          END IF
        END SELECT
      ELSE
        IF (s == 4 .AND. LEN_TRIM(snow_habit) > 0) THEN
          CALL species_from_lut(sp, TRIM(SPNAME(s)), TRIM(cc_dir)//TRIM(cc_file), f, c%n_streams, &
                                habit=TRIM(snow_habit))
        ELSE
          CALL species_from_lut(sp, TRIM(SPNAME(s)), TRIM(cc_dir)//TRIM(cc_file), f, c%n_streams)
        END IF
      END IF
      CALL radar_add_species(op, s, sp)
    END DO
  END SUBROUTINE setup_variant

  ! Build CRTM and standalone inputs for observation index io
  SUBROUTINE build_inputs(io, atm, prof)
    INTEGER, INTENT(IN) :: io
    TYPE(CRTM_Atmosphere_type), INTENT(INOUT) :: atm
    TYPE(radar_profile_type),   INTENT(INOUT) :: prof
    INTEGER :: k, s
    REAL(fp) :: dp, z(0:nlay)
    atm%Climatology = US_STANDARD_ATMOSPHERE
    atm%Absorber_Id = [H2O_ID, O3_ID]
    atm%Absorber_Units = [MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS]
    atm%Level_Pressure(0:nlay) = lp(io, 1:nlev)
    atm%Level_Pressure(0) = MIN(atm%Level_Pressure(0), 0.005_fp)
    z = gp(io, 1:nlev)/G0
    atm%Height(0:nlay) = z
    ! CRTM_Forward zeroes clouds where Cloud_Fraction < 1e-6: run overcast
    atm%Cloud_Fraction = 1.0_fp
    DO k = 1, nlay
      atm%Pressure(k) = 0.5_fp*(lp(io,k) + lp(io,k+1))
      atm%Temperature(k) = 0.5_fp*(tt(io,k) + tt(io,k+1))
      prof%q(k) = 0.5_fp*(qq(io,k) + qq(io,k+1))
      atm%Absorber(k,1) = 1000.0_fp*prof%q(k)/(1.0_fp - prof%q(k))      ! g/kg mixing ratio
      atm%Absorber(k,2) = MAX(0.5_fp*(o3(io,k) + o3(io,k+1))*28.9644_fp/47.9982_fp*1.0e6_fp, 1.0e-3_fp)  ! ppmv
      dp = (atm%Level_Pressure(k) - atm%Level_Pressure(k-1))*100.0_fp
      DO s = 1, NSP
        atm%Cloud(s)%Type = CTYPE(s)
        atm%Cloud(s)%Effective_Radius(k) = 0.0_fp
        atm%Cloud(s)%Water_Content(k) = MAX(0.5_fp*(hyd(io,k,s) + hyd(io,k+1,s)), 0.0_fp)*dp/G0   ! kg m-2
        prof%rho(k,s) = atm%Cloud(s)%Water_Content(k)/(z(k-1) - z(k))                         ! kg m-3
      END DO
    END DO
    prof%p = atm%Pressure(1:nlay)
    prof%t = atm%Temperature(1:nlay)
    prof%z_int = z
    prof%zenith = zenith
  END SUBROUTINE build_inputs

  ! CRTM gas optical depth (vertical, per layer) for every channel
  SUBROUTINE crtm_gas_od(atm, geo, tgas)
    TYPE(CRTM_Atmosphere_type), INTENT(IN) :: atm
    TYPE(CRTM_Geometry_type),   INTENT(IN) :: geo
    REAL(fp), INTENT(OUT) :: tgas(nlay, 3)
    TYPE(CRTM_GeometryInfo_type) :: gi
    TYPE(CRTM_Predictor_type) :: pred
    TYPE(CRTM_PVar_type) :: pvar
    TYPE(CRTM_AAvar_type) :: aavar
    TYPE(CRTM_AtmOptics_type) :: ao
    TYPE(CRTM_AncillaryInput_type) :: anc
    INTEGER :: ic, sidx
    CALL CRTM_GeometryInfo_SetValue(gi, Geometry=geo)
    CALL CRTM_GeometryInfo_Compute(gi)
    DO ic = 1, 3
      sidx = chinfo(sensor_of(ic))%Sensor_Index
      CALL CRTM_Predictor_Create(pred, atm%n_Layers, sidx)
      CALL CRTM_Compute_Predictors(sidx, atm, gi, anc, pred, pvar)
      CALL CRTM_AtmOptics_Create(ao, atm%n_Layers, 16, 1)
      CALL CRTM_AtmOptics_Zero(ao)
      CALL CRTM_Compute_AtmAbsorption(sidx, chan_of(ic), anc, pred, ao, aavar)
      tgas(:,ic) = ao%Optical_Depth(1:nlay)
      CALL CRTM_AtmOptics_Destroy(ao)
      CALL CRTM_Predictor_Destroy(pred)
    END DO
  END SUBROUTINE crtm_gas_od

  SUBROUTINE run_all()
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm(:)
    TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfc(:)
    TYPE(CRTM_Geometry_type),   ALLOCATABLE :: geo(:)
    TYPE(CRTM_Options_type),    ALLOCATABLE :: opt(:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
    TYPE(radar_profile_type),   ALLOCATABLE :: prof(:)
    TYPE(radar_output_type) :: out
    REAL(fp) :: tgas(nlay,3)
    INTEGER :: i0, i1, np, ip, n, ic, iv, sidx
    CALL radar_output_create(out, nlay)
    DO i0 = 1, nsel, CHUNK
      i1 = MIN(nsel, i0 + CHUNK - 1)
      np = i1 - i0 + 1
      ALLOCATE(atm(np), sfc(np), geo(np), opt(np), prof(np))
      CALL CRTM_Atmosphere_Create(atm, nlay, 2, NSP, 0)
      CALL CRTM_Surface_Create(sfc, 0)
      DO ip = 1, np
        atm(ip)%Add_Extra_Layers = .FALSE.
        CALL radar_profile_create(prof(ip), nlay, NSP)
        CALL build_inputs(sel(i0+ip-1), atm(ip), prof(ip))
        sfc(ip)%Water_Coverage = 1.0_fp; sfc(ip)%Water_Temperature = 290.0_fp
        sfc(ip)%Wind_Speed = 5.0_fp; sfc(ip)%Salinity = 33.0_fp
        CALL CRTM_Geometry_SetValue(geo(ip), Sensor_Zenith_Angle=zenith, Sensor_Scan_Angle=zenith)
        opt(ip)%Use_n_Streams = .TRUE.
        opt(ip)%n_Streams = 4
        o_hgt(:,i0+ip-1) = REAL(0.5_fp*(prof(ip)%z_int(0:nlay-1) + prof(ip)%z_int(1:nlay)), 4)
        o_t(:,i0+ip-1) = REAL(prof(ip)%t, 4)
        o_wc(:,i0+ip-1,:) = REAL(prof(ip)%rho, 4)
      END DO
      DO n = 1, 2
        ALLOCATE(rts(CRTM_ChannelInfo_n_Channels(chinfo(n)), np))
        CALL CRTM_RTSolution_Create(rts, nlay)
        err = CRTM_Forward(atm, sfc, geo, chinfo(n:n), rts, Options=opt)
        IF (err /= SUCCESS) STOP "CRTM_Forward failed"
        DO ic = 1, 3
          IF (sensor_of(ic) /= n) CYCLE
          DO ip = 1, np
            c_z(:,i0+ip-1,ic)    = REAL(rts(chan_of(ic),ip)%Reflectivity(1:nlay), 4)
            c_zatt(:,i0+ip-1,ic) = REAL(rts(chan_of(ic),ip)%Reflectivity_Attenuated(1:nlay), 4)
          END DO
        END DO
        CALL CRTM_RTSolution_Destroy(rts)
        DEALLOCATE(rts)
      END DO
      ! gas OD and standalone runs
      !$OMP PARALLEL DO PRIVATE(ip, tgas, ic, iv, out) SCHEDULE(dynamic)
      DO ip = 1, np
        CALL crtm_gas_od(atm(ip), geo(ip), tgas)
        CALL radar_output_create(out, nlay)
        DO ic = 1, 3
          c_tgas(:,i0+ip-1,ic) = REAL(tgas(:,ic), 4)
          DO iv = 1, NV
            prof(ip)%tau_gas = tgas(:,ic)
            CALL radar_forward(ops(ic,iv), prof(ip), out)
            s_z(:,i0+ip-1,ic,iv) = REAL(out%z_dbz, 4)
            s_zatt(:,i0+ip-1,ic,iv) = REAL(out%zatt_dbz, 4)
            s_tcum(:,i0+ip-1,ic,iv) = REAL(out%tau_cum, 4)
            s_thyd(:,i0+ip-1,ic,iv) = REAL(out%tau_hydro, 4)
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
      CALL CRTM_Atmosphere_Destroy(atm)
      DEALLOCATE(atm, sfc, geo, opt, prof)
      IF (MOD(i0/CHUNK, 10) == 0) WRITE(*,'(a,i0,a,i0)') '  done ', i1, ' of ', nsel
    END DO
    sidx = 0
  END SUBROUTINE run_all

  SUBROUTINE write_output()
    INTEGER :: dl, do_, dc, dv, ds
    CHARACTER(16) :: vn(NV) = [CHARACTER(16) :: 'A_crtmgas', 'B_itulbl', 'C_approx', 'D_physical', &
                                                'E_nurtm_delta', 'F_nurtm_full']
    CALL chk(nf90_create(out_file, IOR(NF90_NETCDF4, NF90_CLOBBER), id))
    CALL chk(nf90_def_dim(id, 'layer', nlay, dl))
    CALL chk(nf90_def_dim(id, 'obs', nsel, do_))
    CALL chk(nf90_def_dim(id, 'channel', 3, dc))
    CALL chk(nf90_def_dim(id, 'variant', NV, dv))
    CALL chk(nf90_def_dim(id, 'species', NSP, ds))
    nv_ = 0
    CALL defv('frequency', [dc], 'GHz'); CALL defv('obs_index', [do_], '1-based index in the IFS file')
    CALL defv('height', [dl,do_], 'm'); CALL defv('temperature', [dl,do_], 'K')
    CALL defv('water_content', [dl,do_,ds], 'kg m-3 (cloud_liquid, cloud_ice, rain, snow)')
    CALL defv('crtm_z', [dl,do_,dc], 'dBZ'); CALL defv('crtm_zatt', [dl,do_,dc], 'dBZ')
    CALL defv('crtm_tau_gas', [dl,do_,dc], 'vertical one-way gas optical depth (ODPS)')
    CALL defv('sa_z', [dl,do_,dc,dv], 'dBZ'); CALL defv('sa_zatt', [dl,do_,dc,dv], 'dBZ')
    CALL defv('sa_tau_cum', [dl,do_,dc,dv], 'vertical one-way optical depth TOA->gate')
    CALL defv('sa_tau_hydro', [dl,do_,dc,dv], 'vertical one-way hydrometeor optical depth of the layer')
    CALL chk(nf90_put_att(id, NF90_GLOBAL, 'variants', &
         'A: CRTM-compatible, CRTM gas OD; B: CRTM-compatible, ITU-R P.676 LBL gas; '// &
         'C: approximate gas (P.676 Annex 2) and hydrometeor (P.838/P.840) attenuation; '// &
         'D: unscaled extinction, cloud backscatter on, LBL gas; '// &
         'E: NuRTM DDA particles + Thompson 2008 PSD, delta scaling, LBL gas; F: as E, unscaled extinction'))
    CALL chk(nf90_put_att(id, NF90_GLOBAL, 'zenith_deg', zenith))
    CALL chk(nf90_enddef(id))
    CALL chk(nf90_put_var(id, v(1), freq)); CALL chk(nf90_put_var(id, v(2), sel))
    CALL chk(nf90_put_var(id, v(3), o_hgt)); CALL chk(nf90_put_var(id, v(4), o_t))
    CALL chk(nf90_put_var(id, v(5), o_wc))
    CALL chk(nf90_put_var(id, v(6), c_z)); CALL chk(nf90_put_var(id, v(7), c_zatt))
    CALL chk(nf90_put_var(id, v(8), c_tgas))
    CALL chk(nf90_put_var(id, v(9), s_z)); CALL chk(nf90_put_var(id, v(10), s_zatt))
    CALL chk(nf90_put_var(id, v(11), s_tcum)); CALL chk(nf90_put_var(id, v(12), s_thyd))
    CALL chk(nf90_close(id))
    WRITE(*,*) 'wrote ', TRIM(out_file), ' (variants: ', (TRIM(vn(i))//' ', i=1,NV), ')'
  END SUBROUTINE write_output

  SUBROUTINE defv(name, dims, units)
    CHARACTER(*), INTENT(IN) :: name, units
    INTEGER, INTENT(IN) :: dims(:)
    nv_ = nv_ + 1
    IF (name == 'frequency') THEN
      CALL chk(nf90_def_var(id, name, NF90_DOUBLE, dims, v(nv_)))
    ELSE IF (name == 'obs_index') THEN
      CALL chk(nf90_def_var(id, name, NF90_INT, dims, v(nv_)))
    ELSE
      CALL chk(nf90_def_var(id, name, NF90_FLOAT, dims, v(nv_), deflate_level=1))
    END IF
    CALL chk(nf90_put_att(id, v(nv_), 'units', units))
  END SUBROUTINE defv

  ! Jacobians of the attenuated dBZ with respect to the water contents (kg m-2)
  ! from the CRTM K-matrix and from the standalone K-matrix (variant A),
  ! for the n most precipitating selected profiles.
  SUBROUTINE jacobians(nj)
    INTEGER, INTENT(IN) :: nj
    TYPE(CRTM_Atmosphere_type) :: atm(1), atm_k(1,1)
    TYPE(CRTM_Surface_type)    :: sfc(1), sfc_k(1,1)
    TYPE(CRTM_Geometry_type)   :: geo(1)
    TYPE(CRTM_Options_type)    :: opt(1)
    TYPE(CRTM_RTSolution_type) :: rts(1,1), rts_k(1,1)
    TYPE(radar_profile_type) :: prof
    TYPE(radar_k_type) :: kk
    TYPE(CRTM_ChannelInfo_type) :: ci(1)
    REAL(fp) :: tgas(nlay,3), dz(nlay)
    REAL(4), ALLOCATABLE :: jc(:,:,:,:,:), js(:,:,:,:,:)
    INTEGER, ALLOCATABLE :: pick(:)
    REAL(fp), ALLOCATABLE :: score(:)
    INTEGER :: ij, ip, ic, g, s, d1, d2, d3, d4, d5, vjc, vjs, vp
    ALLOCATE(score(nsel)); score = SUM(SUM(o_wc(:,:,3:4), DIM=3), DIM=1)
    ALLOCATE(pick(nj))
    DO ij = 1, nj
      pick(ij) = MAXLOC(score, DIM=1); score(pick(ij)) = -1.0_fp
    END DO
    ALLOCATE(jc(nlay,nlay,NSP,3,nj), js(nlay,nlay,NSP,3,nj))
    jc = 0; js = 0
    CALL radar_k_create(kk, nlay, NSP)
    DO ij = 1, nj
      CALL CRTM_Atmosphere_Create(atm, nlay, 2, NSP, 0)
      CALL CRTM_Surface_Create(sfc, 0)
      atm(1)%Add_Extra_Layers = .FALSE.
      CALL radar_profile_create(prof, nlay, NSP)
      CALL build_inputs(sel(pick(ij)), atm(1), prof)
      sfc(1)%Water_Coverage = 1.0_fp; sfc(1)%Water_Temperature = 290.0_fp
      sfc(1)%Wind_Speed = 5.0_fp; sfc(1)%Salinity = 33.0_fp
      CALL CRTM_Geometry_SetValue(geo(1), Sensor_Zenith_Angle=zenith, Sensor_Scan_Angle=zenith)
      opt(1)%Use_n_Streams = .TRUE.; opt(1)%n_Streams = 4
      CALL crtm_gas_od(atm(1), geo(1), tgas)
      dz = prof%z_int(0:nlay-1) - prof%z_int(1:nlay)
      DO ic = 1, 3
        ci(1) = chinfo(sensor_of(ic))
        err = CRTM_ChannelInfo_Subset(ci(1), Channel_Subset=[chan_of(ic)])
        ! standalone
        prof%tau_gas = tgas(:,ic)
        CALL radar_k(ops(ic,1), prof, OUT_ZATT_DBZ, kk)
        DO s = 1, NSP
          DO g = 1, nlay
            js(g,:,s,ic,ij) = REAL(kk%rho(g,:,s)/dz, 4)      ! per kg m-2
          END DO
        END DO
        ! CRTM, one gate at a time
        DO g = 1, nlay
          CALL CRTM_Atmosphere_Create(atm_k(1,1), nlay, 2, NSP, 0)
          CALL CRTM_Surface_Create(sfc_k(1,1), 0)
          CALL CRTM_RTSolution_Create(rts, nlay)
          CALL CRTM_RTSolution_Create(rts_k, nlay)
          CALL CRTM_Atmosphere_Zero(atm_k(1,1))
          CALL CRTM_Surface_Zero(sfc_k(1,1))
          CALL CRTM_RTSolution_Zero(rts_k(1,1))
          rts_k(1,1)%Radiance = 0.0_fp; rts_k(1,1)%Brightness_Temperature = 0.0_fp
          rts_k(1,1)%Reflectivity_Attenuated(g) = 1.0_fp
          err = CRTM_K_Matrix(atm, sfc, rts_k, geo, ci, atm_k, sfc_k, rts, Options=opt)
          IF (err /= SUCCESS) STOP 'CRTM_K_Matrix failed'
          DO s = 1, NSP
            jc(g,:,s,ic,ij) = REAL(atm_k(1,1)%Cloud(s)%Water_Content(1:nlay), 4)
          END DO
          CALL CRTM_RTSolution_Destroy(rts); CALL CRTM_RTSolution_Destroy(rts_k)
          CALL CRTM_Atmosphere_Destroy(atm_k(1,1)); CALL CRTM_Surface_Destroy(sfc_k(1,1))
        END DO
      END DO
      CALL CRTM_Atmosphere_Destroy(atm)
      WRITE(*,'(a,i0,a,i0)') '  jacobian profile ', ij, ' obs ', sel(pick(ij))
    END DO
    CALL chk(nf90_open(out_file, NF90_WRITE, id))
    CALL chk(nf90_redef(id))
    CALL chk(nf90_def_dim(id, 'gate', nlay, d1)); CALL chk(nf90_def_dim(id, 'jlayer', nlay, d2))
    CALL chk(nf90_inq_dimid(id, 'species', d3)); CALL chk(nf90_inq_dimid(id, 'channel', d4))
    CALL chk(nf90_def_dim(id, 'jprofile', nj, d5))
    CALL chk(nf90_def_var(id, 'jac_crtm', NF90_FLOAT, [d1,d2,d3,d4,d5], vjc, deflate_level=1))
    CALL chk(nf90_def_var(id, 'jac_sa', NF90_FLOAT, [d1,d2,d3,d4,d5], vjs, deflate_level=1))
    CALL chk(nf90_def_var(id, 'jac_obs_pos', NF90_INT, [d5], vp))
    CALL chk(nf90_put_att(id, vjc, 'units', 'dBZ (attenuated) per kg m-2 of water content'))
    CALL chk(nf90_enddef(id))
    CALL chk(nf90_put_var(id, vjc, jc)); CALL chk(nf90_put_var(id, vjs, js)); CALL chk(nf90_put_var(id, vp, pick))
    CALL chk(nf90_close(id))
  END SUBROUTINE jacobians

END PROGRAM compare_crtm
