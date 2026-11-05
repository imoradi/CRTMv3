!
! CRTM_Active_Sensor
!
! Module containing the active sensor routines.
!
!
! CREATION HISTORY:
!       Written by:     Isaac Moradi - NASA GSFC March 10, 2022
!                       Isaac.Moradi@NASA.GOV
!

MODULE CRTM_Active_Sensor

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use statements
  USE Type_Kinds,      ONLY: FP
  USE CRTM_Parameters, ONLY: ZERO, ONE, TWO, TEN, ONE_THOUSAND, MISSING_REFL, EPSILON_FP
  USE CRTM_SpcCoeff,            ONLY: SC, &
                                      SpcCoeff_IsMicrowaveSensor , &
                                      SpcCoeff_IsInfraredSensor  , &
                                      SpcCoeff_IsVisibleSensor   , &
                                      SpcCoeff_IsUltravioletSensor
  USE CRTM_Atmosphere_Define,   ONLY: CRTM_Atmosphere_type, &
                                      H2O_ID
  USE CRTM_AtmOptics_Define,     ONLY: CRTM_AtmOptics_type
  USE CRTM_RTSolution_Define,    ONLY: CRTM_RTSolution_type
  USE Spectral_Units_Conversion, ONLY: GHz_to_inverse_cm
  USE Fundamental_Constants,     ONLY: PI
  USE ODPS_CoordinateMapping,     ONLY: Geopotential_Height
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type
  USE Message_Handler         , ONLY: Display_Message, WARNING
  USE CRTM_Options_Define,      ONLY: CRTM_Options_type, &
                                      ACTIVE_MS_NONE, &
                                      ACTIVE_MS_EFFECTIVE_ATTENUATION, &
                                      ACTIVE_MS_TDTS
  USE CRTM_Active_MS_TDTS,      ONLY: TDTS_Wide_Angle, &
                                      TDTS_Wide_Angle_TL, &
                                      TDTS_Wide_Angle_AD
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC  :: CRTM_Compute_Reflectivity
  PUBLIC  :: CRTM_Compute_Reflectivity_TL
  PUBLIC  :: CRTM_Compute_Reflectivity_AD
  PUBLIC  :: Calculate_Height
  PUBLIC  :: Calculate_Cloud_Water_Density

  ! ----------
  ! Parameters
  ! ----------
  REAL(fp), PARAMETER :: POINT_01 = 0.01_fp
  REAL(fp), PARAMETER :: Kw_2_fixed = 0.93_fp
  ! 1.0d818 converts from m^3 to mm^6/m^3 (standard radar reflectivity units)
  REAL(fp), PARAMETER :: M6_MM6 = 1.0d18
  REAL(fp), PARAMETER :: REFLECTIVITY_THRESHOLD = TINY(REAL(fp))
  ! Limits on the dBZ tangent-linear/adjoint. d(dBZ)/dZ = 10/(Z ln10) is
  ! unbounded as Z -> 0, which gives huge Jacobians for near-clear layers.
  !   REFLECTIVITY_TLAD_MIN: below this Z (mm^6 m^-3) the dBZ TL/AD is zero.
  !   REFLECTIVITY_TLAD_CAP: Z in the denominator is held at no less than this,
  !                          so the gain is at most 10/(CAP ln10).
  ! Both depend on the base state only, so TL stays linear and AD exact.
  ! Defaults are 0 dBZ; override with the environment variables
  ! CRTM_ACTIVE_TLAD_MIN_DBZ and CRTM_ACTIVE_TLAD_CAP_DBZ (in dBZ).
  REAL(fp), SAVE :: REFLECTIVITY_TLAD_MIN = ONE
  REAL(fp), SAVE :: REFLECTIVITY_TLAD_CAP = ONE
  LOGICAL,  SAVE :: TLAD_LIMITS_SET = .FALSE.
  REAL(fp), PARAMETER :: MM6_TO_CM6 = 1.0e-6_fp
  REAL(fp), PARAMETER :: KW_TEMPERATURE_K = 273.15_fp
  REAL(fp), PARAMETER :: FOUR_PI = 4.0_fp*3.14159265358979323846_fp
  ! Layers with a single-scatter albedo below this are treated as non-scattering
  ! (matches BS_THRESHOLD in CRTM_AtmOptics_Combine)
  REAL(fp), PARAMETER :: BS_THRESHOLD_MS = 1.0e-10_fp
CONTAINS

!--------------------------------------------------------------------------------
!
! NAME:
!       Calculate_Height
!
! PURPOSE:
!       Subroutine to calculate the height usig atmospheric input profiles
!
! CALLING SEQUENCE:
!       Height = Calculate_Height(Atm)
!
! INPUT ARGUMENTS:
!
!       Atm:            Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
! OUTPUT ARGUMENTS:
!       Height:         Array containing height
!                       UNITS:      m
!                       TYPE:       REAL
!                       DIMENSION:  Array
!
!--------------------------------------------------------------------------------

Function Calculate_Height(Atm) RESULT (Height)

   TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
   REAL(fp) :: Height(0:Atm%n_Layers)

   INTEGER  :: j, H2O_idx

   !** locate the WV absorber ID
   H2O_idx = 1
   DO j = 1,Atm%n_Absorbers
        IF (Atm%Absorber_ID(j) == H2O_ID) H2O_idx = j
   END DO
   ! Calculate the geometric heights of the pressure levels in km
   CALL Geopotential_Height(Atm%Level_Pressure      , & ! Input
                            Atm%Temperature         , & ! Input
                            Atm%Absorber(:, H2O_idx), & ! Input
                            ZERO                    , & ! Input - surface height
                            Height                ) ! Output in km
   
   Height = Height * ONE_THOUSAND ! convert to meters

END FUNCTION Calculate_Height

!--------------------------------------------------------------------------------
!
! NAME:
!       Calculate_Cloud_Water_Density
!
! PURPOSE:
!       Subroutine to calculate cloud water density (mass/volume) and also height if not set yet
!
! CALLING SEQUENCE:
!       CALL Calculate_Cloud_Water_Density(Atm, GeometryInfo)
!
! INPUT ARGUMENTS:
!
!       Atm:            Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN/OUT)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!--------------------------------------------------------------------------------

SUBROUTINE Calculate_Cloud_Water_Density(Atm, &
                                         GeometryInfo)
   TYPE(CRTM_Atmosphere_type), INTENT(IN OUT)     :: Atm
   TYPE(CRTM_GeometryInfo_type), OPTIONAL, INTENT(IN)     :: GeometryInfo
   INTEGER :: n_Layers
   REAL(fp) :: Height(0:Atm%n_Layers), dZ_m(Atm%n_Layers)
   Integer :: n

   n_Layers = Atm%n_Layers

   ! Calculate heights if hasn't been set already
   IF (ALL(Atm%Height .LT. EPSILON_FP)) THEN
       Atm%Height = Calculate_Height(Atm)
   END IF

   dZ_m = (Atm%Height(0:n_Layers-1) - Atm%Height(1:n_Layers)) 

   DO n = 1, Atm%n_Clouds
      Atm%Cloud(n)%Water_Density = Atm%Cloud(n)%Water_Content / dZ_m
   END DO

END SUBROUTINE Calculate_Cloud_Water_Density

!--------------------------------------------------------------------------------
!
! NAME:
!       Water_Permittivity_Turner_2016
!
! PURPOSE:
!       Function to calculate water permittivity using Turnet et al. 2016 Method
!
! CALLING SEQUENCE:
!       CALL Water_Permittivity_Turner_2016(freq_Hz, & ! Input
!                                           temp_K)  ! Input
!                                           RESULT(perm) ! Output
!
! INPUT ARGUMENTS:
!
!       freq_Hz:        Scalar representing frequency in Hz.
!                       UNITS:      Hz
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       temp_K:         Array representing temperature in Kelvin.
!                       UNITS:      Kelvin
!                       TYPE:       REAL(fp)
!                       DIMENSION:  1-D Array
!                       ATTRIBUTES: INTENT(IN)
!
! OOUTUT ARGUMENTS:
!
!       perm:        Array representing water permittivity in a complex form of (real, imaginary).
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  1-D Complex Array
!
!--------------------------------------------------------------------------------

FUNCTION Water_Permittivity_Turner_2016(freq_Hz, temp_K) RESULT(perm)

  REAL(fp), INTENT(IN) :: freq_Hz, temp_K(:)
  REAL(fp), DIMENSION(SIZE(temp_K)) :: perm_re, perm_im
  COMPLEX :: perm(SIZE(temp_K))

  REAL(fp) :: a1, a2, b1, b2, c1, c2, d1, d2
  REAL(fp) :: e0, e1, e2, e3, tc
  REAL(fp), DIMENSION(SIZE(temp_K)) :: delta1, delta2, tau1, tau2
  REAL(fp), DIMENSION(SIZE(temp_K)) :: term11, term12, denom1, denom2
  REAL(fp), DIMENSION(SIZE(temp_K)) :: term21, term22, temp_c, eps_s

  temp_c = temp_k - 273.15_fp

  IF (freq_Hz .GT. 500.d9) THEN
     STOP
  END IF

  a1 = 8.111E+1
  a2 = 2.025E+0

  b1 = 4.434E-3
  b2 = 1.073E-2

  c1 = 1.302E-13
  c2 = 1.012E-14

  d1 = 6.627E+2
  d2 = 6.089E+2

  tc = 133.1383_fp

  e0 = 87.9144_fp
  e1 = 0.404399_fp
  e2 = 9.58726d-4
  e3 = -1.32802d-6

  ! static dielectric constant
  eps_s = e0 + e1 * temp_c + e2 * (temp_c**TWO) + e3 * (temp_c**3.0_fp)

  delta1 = a1 * EXP(-b1 * temp_c)
  delta2 = a2 * EXP(-b2 * temp_c)

  tau1   = c1 * EXP(d1 / (temp_c + tc))
  tau2   = c2 * EXP(d2 / (temp_c + tc))

  denom1 = ONE + (TWO * PI * freq_Hz * tau1)**TWO
  denom2 = ONE + (TWO * PI * freq_Hz * tau2)**TWO

  term11 = ((tau1**TWO) * delta1) / denom1
  term12 = ((tau2**TWO) * delta2) / denom2

  term21 = (tau1 * delta1) / denom1
  term22 = (tau2 * delta2) / denom2

  perm_re = eps_s - ((TWO * PI * freq_Hz)**TWO) * (term11 + term12)
  perm_im = (TWO * PI * freq_Hz) * (term21 + term22)

  perm = cmplx(perm_re, perm_im, FP)

END FUNCTION Water_Permittivity_Turner_2016

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_Reflectivity
!
! PURPOSE:
!       Subroutine to calculate the reflectivity for active sensors
!
! CALLING SEQUENCE:
!       CALL CRTM_Compute_Reflectivity(Atm         , &  ! Input
!                                      AtmOptics   , &  ! Input
!                                      GeometryInfo , &  ! Input
!                                      SensorIndex , &  ! Input
!                                      ChannelIndex, &  ! Input
!                                      RTSolution )     ! Input/Output
!
! INPUT ARGUMENTS:
!
!       Atm:            Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:      Structure containing the combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       RTSolution:     Structure containing the soluition to the RT equation
!                       for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Compute_Reflectivity(Atm         , &  ! Input
                                       AtmOptics   , &  ! Input
                                       GeometryInfo , &  ! Input
                                       SensorIndex , &  ! Input
                                       ChannelIndex, &  ! Input
                                       RTSolution  , &  ! Input/Output
                                       Options     )    ! Optional input
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER                   , INTENT(IN)     :: SensorIndex
    INTEGER                   , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_AtmOptics_type) , INTENT(IN OUT)     :: AtmOptics
    TYPE(CRTM_RTSolution_type), INTENT(IN OUT) :: RTSolution
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN) :: Options

    REAL(fp) :: Frequency, Wavenumber, Wavelength_m
    REAL(fp) :: Reflectivity(AtmOptics%n_Layers)
    REAL(fp) :: Reflectivity_Attenuated(AtmOptics%n_Layers)
    REAL(fp) :: Transmittance(AtmOptics%n_Layers)
    REAL(fp) :: Optical_Depth(AtmOptics%n_Layers)
    REAL(fp) :: P1(AtmOptics%n_Layers)
    REAL(fp) :: Height(0:AtmOptics%n_Layers)
    REAL(fp) :: dZ_m(AtmOptics%n_Layers)
    REAL(fp) :: Temp_K(AtmOptics%n_Layers)
    REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Kw_2, perm_re, perm_im
    COMPLEX :: perm(AtmOptics%n_Layers), kw(AtmOptics%n_Layers)
    INTEGER :: k, n_Layers
    REAL(fp) :: temp_sum
    REAL(fp) :: Z_ms(AtmOptics%n_Layers)

    n_Layers = Atm%n_Layers

    ! Calculate heights if hasn't been set already
    IF (ALL(Atm%Height .LT. EPSILON_FP)) THEN
        Height = Calculate_Height(Atm)
    ELSE
        Height = Atm%Height
    ENDIF
    dZ_m = (Height(0:n_Layers-1) - Height(1:n_Layers))

    IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) THEN
       Frequency = SC(SensorIndex)%Frequency(ChannelIndex) ! GHz
       Wavelength_m = POINT_01 / GHz_to_Inverse_cm( Frequency )
    ELSE IF( SpcCoeff_IsInfraredSensor(SC(SensorIndex)) ) THEN
       Wavenumber = SC(SensorIndex)%Wavenumber(ChannelIndex) ! 1/cm
       Wavelength_m = POINT_01 / Wavenumber
    END IF

    Temp_K = KW_TEMPERATURE_K
    perm = Water_Permittivity_Turner_2016(Frequency * 1.0d9, & ! Input
                                          Temp_K) ! Input

    perm_re = REAL(REAL(perm))  ! Double REAL is required to avoud issues with GNU Fortran
    perm_im = REAL(AIMAG(perm))
    kw = (perm - ONE )/(perm + TWO)
    Kw_2 = ABS(kw)**TWO      
    P1 = (M6_MM6 * Wavelength_m**4.0_fp) / (PI**5.0_fp * Kw_2)

    ! Calculate transmittance along the (slant) path, with the
    ! multiple-scattering treatment selected in Options
    CALL Path_Optical_Depth(AtmOptics, GeometryInfo, Options, Optical_Depth)
    DO k = 1, n_layers
       temp_sum = SUM(Optical_Depth(1:k))
       Transmittance(k) = EXP(-TWO * temp_sum)
    END DO
    
    ! Divide by dZ_m as water conntent calculations are based on
    ! vertical thickness of layers
    Reflectivity =  P1 * (AtmOptics%Backscat_Coefficient / dZ_m) ! mm^6 m^-3
    Reflectivity_Attenuated = Transmittance * Reflectivity ! mm^6 m^-3
    ! Wide-angle multiple scattering (TDTS) adds to the attenuated return
    IF (MS_Algorithm_Id(Options) == ACTIVE_MS_TDTS) THEN
       CALL TDTS_Reflectivity(AtmOptics, GeometryInfo, Options, Height, P1, Z_ms)
       Reflectivity_Attenuated = Reflectivity_Attenuated + Z_ms
    END IF
    ! Convert the unit to dBz
    ! Linear Z is defined everywhere (0 means no echo); dBZ is missing below threshold
    RTSolution%ReflectivityLinear = MM6_TO_CM6 * MAX(Reflectivity, ZERO)  ! cm^6m^-3
    WHERE (Reflectivity .GT.  REFLECTIVITY_THRESHOLD)
        RTSolution%Reflectivity       = TEN * LOG10(Reflectivity) ! [dBZ]
    ELSE WHERE
        RTSolution%Reflectivity       = MISSING_REFL
    END WHERE

    ! Convert the unit to dBz
    ! Note that Reflectivity can be greater than zero but Reflectivity_Attenuated
    ! can be still zero if Transmittance is zero
    RTSolution%Reflectivity_AttenuatedLinear = MM6_TO_CM6 * MAX(Reflectivity_Attenuated, ZERO) ! cm^6m^-3
    WHERE (Reflectivity_Attenuated .GT.  REFLECTIVITY_THRESHOLD)
        RTSolution%Reflectivity_Attenuated = TEN * LOG10(Reflectivity_Attenuated)      ! dBZ
    ELSE WHERE
        RTSolution%Reflectivity_Attenuated       = MISSING_REFL
    END WHERE
    
  END SUBROUTINE CRTM_Compute_Reflectivity

SUBROUTINE CRTM_Compute_Reflectivity_TL(Atm, &
                                       AtmOptics, &
                                       AtmOptics_TL, &
                                       GeometryInfo, &
                                       SensorIndex, &
                                       ChannelIndex, &
                                       RTSolution_TL, &
                                       Options)
  TYPE(CRTM_Atmosphere_type),     INTENT(IN)     :: Atm
  TYPE(CRTM_AtmOptics_type),      INTENT(IN)     :: AtmOptics, AtmOptics_TL
  TYPE(CRTM_GeometryInfo_type),   INTENT(IN)     :: GeometryInfo
  INTEGER,                        INTENT(IN)     :: SensorIndex, ChannelIndex
  TYPE(CRTM_RTSolution_type),     INTENT(IN OUT) :: RTSolution_TL
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)  :: Options

  ! Locals
  REAL(fp) :: Frequency, Wavenumber, Wavelength_m
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Reflectivity, Reflectivity_Attenuated
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Reflectivity_TL, Reflectivity_Attenuated_TL
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ReflectivityLinear, Reflectivity_AttenuatedLinear
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ReflectivityLinear_TL, Reflectivity_AttenuatedLinear_TL
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Transmittance, Transmittance_TL
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Optical_Depth, Optical_Depth_TL
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: P1, Kw_2, perm_re, perm_im
  REAL(fp), DIMENSION(0:AtmOptics%n_Layers) :: Height
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: dZ_m
  COMPLEX :: perm(AtmOptics%n_Layers)
  INTEGER :: k, n_Layers
  REAL(fp) :: temp_sum
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Z_ms, Z_ms_TL

  n_Layers = Atm%n_Layers
  CALL Set_TLAD_Limits()

  ! Compute layer interface heights if missing
  IF (ALL(Atm%Height .LT. EPSILON_FP)) THEN
     Height = Calculate_Height(Atm)
  ELSE
     Height = Atm%Height
  END IF

  dZ_m = Height(0:n_Layers-1) - Height(1:n_Layers)

  ! Compute wavelength
  IF (SpcCoeff_IsMicrowaveSensor(SC(SensorIndex))) THEN
     Frequency = SC(SensorIndex)%Frequency(ChannelIndex)
     Wavelength_m = POINT_01 / GHz_to_Inverse_cm(Frequency)
  ELSE IF (SpcCoeff_IsInfraredSensor(SC(SensorIndex))) THEN
     Wavenumber = SC(SensorIndex)%Wavenumber(ChannelIndex)
     Wavelength_m = POINT_01 / Wavenumber
  END IF

  ! Compute permittivity and related factor
  perm = Water_Permittivity_Turner_2016(Frequency * 1.0d9, &
                                        SPREAD(KW_TEMPERATURE_K, DIM=1, NCOPIES=n_Layers))
  perm_re = REAL(perm)
  perm_im = AIMAG(perm)
  Kw_2 = ABS((perm - ONE) / (perm + TWO)) ** TWO
  P1 = (M6_MM6 * Wavelength_m**4.0_fp) / (PI**5.0_fp * Kw_2)

  ! Forward (base state) calculations
  CALL Path_Optical_Depth(AtmOptics, GeometryInfo, Options, Optical_Depth)
  DO k = 1, n_Layers
     temp_sum = SUM(Optical_Depth(1:k))
     Transmittance(k) = EXP(-TWO * temp_sum)
  END DO

  Reflectivity = P1 * (AtmOptics%Backscat_Coefficient / dZ_m)
  Reflectivity_Attenuated = Transmittance * Reflectivity
  Z_ms_TL = ZERO
  IF (MS_Algorithm_Id(Options) == ACTIVE_MS_TDTS) THEN
     CALL TDTS_Reflectivity_TL(AtmOptics, AtmOptics_TL, GeometryInfo, Options, Height, P1, Z_ms, Z_ms_TL)
     Reflectivity_Attenuated = Reflectivity_Attenuated + Z_ms
  END IF

  ! Tangent-linear calculations
  CALL Path_Optical_Depth_TL(AtmOptics, AtmOptics_TL, GeometryInfo, Options, Optical_Depth_TL)
  DO k = 1, n_Layers
     temp_sum = SUM(Optical_Depth_TL(1:k))
     Transmittance_TL(k) = -TWO * Transmittance(k) * temp_sum
  END DO

  Reflectivity_TL = P1 * (AtmOptics_TL%Backscat_Coefficient / dZ_m)
  Reflectivity_Attenuated_TL = Transmittance * Reflectivity_TL + Reflectivity * Transmittance_TL + Z_ms_TL

  ! Linear reflectivity TL
  ReflectivityLinear = MM6_TO_CM6 * Reflectivity
  ReflectivityLinear_TL = MM6_TO_CM6 * Reflectivity_TL
  Reflectivity_AttenuatedLinear = MM6_TO_CM6 * Reflectivity_Attenuated
  Reflectivity_AttenuatedLinear_TL = MM6_TO_CM6 * Reflectivity_Attenuated_TL

  ! Convert to dBZ tangent linear if above threshold
  ! Linear outputs are smooth in Z (no log), so their TL is not masked
  RTSolution_TL%ReflectivityLinear = ReflectivityLinear_TL
  WHERE (Reflectivity > REFLECTIVITY_TLAD_MIN)
     Reflectivity_TL = TEN * Reflectivity_TL / (MAX(Reflectivity, REFLECTIVITY_TLAD_CAP) * LOG(TEN))
     RTSolution_TL%Reflectivity = Reflectivity_TL
  ELSEWHERE
     RTSolution_TL%Reflectivity = ZERO
  END WHERE

  RTSolution_TL%Reflectivity_AttenuatedLinear = Reflectivity_AttenuatedLinear_TL
  WHERE (Reflectivity_Attenuated > REFLECTIVITY_TLAD_MIN)
     Reflectivity_Attenuated_TL = TEN * Reflectivity_Attenuated_TL / &
                                  (MAX(Reflectivity_Attenuated, REFLECTIVITY_TLAD_CAP) * LOG(TEN))
     RTSolution_TL%Reflectivity_Attenuated = Reflectivity_Attenuated_TL
  ELSEWHERE
     RTSolution_TL%Reflectivity_Attenuated = ZERO
  END WHERE

END SUBROUTINE CRTM_Compute_Reflectivity_TL

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_Reflectivity_AD
!
! PURPOSE:
!       Subroutine to calculate the adjoint reflectivity for active instruments.
!       Includes both logarithmic (dBZ) and linear reflectivity adjoints:
!       ReflectivityLinear_AD and Reflectivity_AttenuatedLinear_AD.
!
! CALLING SEQUENCE:
!        CALL CRTM_Compute_Reflectivity_AD(Atm, &
!                                          AtmOptics     , &  ! Input
!                                          RTSolution    , &  ! Input
!                                          GeometryInfo  , &  ! Input
!                                          SensorIndex   , &  ! Input
!                                          ChannelIndex  , &  ! Input
!                                          AtmOptics_AD  , &  ! Input/Output
!                                          RTSolution_AD )    ! Input/Output
!
! INPUT ARGUMENTS:
!       Atm:            Structure containing the atmospheric state data.
!                       TYPE:       CRTM_Atmosphere_type
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:      Structure containing the combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       TYPE:       CRTM_AtmOptics_type
!                       ATTRIBUTES: INTENT(IN)
!
!       RTSolution:     Structure containing the solution to the RT equation
!                       for the given inputs.
!                       TYPE:       CRTM_RTSolution_type
!                       ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       TYPE:       CRTM_GeometryInfo_type
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:    Sensor index id used to access shared coefficient data.
!                       TYPE:       INTEGER
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:   Channel index id used to access shared coefficient data.
!                       TYPE:       INTEGER
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       RTSolution_AD:  Structure containing the RT solution adjoint inputs.
!                       TYPE:       CRTM_RTSolution_type
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       AtmOptics_AD:   Structure containing the adjoint combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols, including:
!                       - ReflectivityLinear_AD
!                       - Reflectivity_AttenuatedLinear_AD
!                       TYPE:       CRTM_AtmOptics_type
!                       ATTRIBUTES: INTENT(IN OUT)
!
!--------------------------------------------------------------------------------
SUBROUTINE CRTM_Compute_Reflectivity_AD(Atm, &
                                       AtmOptics, &
                                       RTSolution, &
                                       GeometryInfo, &
                                       SensorIndex, &
                                       ChannelIndex, &
                                       AtmOptics_AD, &
                                       RTSolution_AD, &
                                       Options)
  !-----------------------------------------------------------------------
  !  Compute adjoint sensitivities of reflectivity. The dBZ and linear
  !  adjoint inputs in RTSolution_AD are all accumulated, so a caller may
  !  seed any combination of them.
  !
  !  Supported adjoint targets:
  !    - RTSolution_AD%Reflectivity_Attenuated        (logarithmic dBZ)
  !    - RTSolution_AD%Reflectivity_AttenuatedLinear  (linear cm^6/m^3)
  !    - RTSolution_AD%Reflectivity                   (logarithmic dBZ)
  !    - RTSolution_AD%ReflectivityLinear             (linear cm^6/m^3)
  !
  !-----------------------------------------------------------------------

  ! Arguments
  TYPE(CRTM_Atmosphere_type),      INTENT(IN)     :: Atm
  TYPE(CRTM_GeometryInfo_type),    INTENT(IN)     :: GeometryInfo
  INTEGER,                         INTENT(IN)     :: SensorIndex
  INTEGER,                         INTENT(IN)     :: ChannelIndex
  TYPE(CRTM_AtmOptics_type),       INTENT(IN)     :: AtmOptics
  TYPE(CRTM_RTSolution_type),      TARGET, INTENT(IN) :: RTSolution
  TYPE(CRTM_AtmOptics_type),       INTENT(IN OUT) :: AtmOptics_AD
  TYPE(CRTM_RTSolution_type),      INTENT(IN OUT) :: RTSolution_AD
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)   :: Options

  ! Locals
  INTEGER :: k, j, n_Layers
  REAL(fp) :: Frequency, Wavenumber, Wavelength_m
  REAL(fp) :: temp_sum

  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Reflectivity, Reflectivity_Attenuated
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ReflectivityLinear, Reflectivity_AttenuatedLinear
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Transmittance, Optical_Depth
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: P1, Kw_2, perm_re, perm_im
  REAL(fp), DIMENSION(0:AtmOptics%n_Layers) :: Height
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: dZ_m
  COMPLEX :: perm(AtmOptics%n_Layers)

  ! Adjoint variables
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Transmittance_AD, Reflectivity_AD, Reflectivity_Attenuated_AD, Optical_Depth_AD
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ReflectivityLinear_AD, Reflectivity_AttenuatedLinear_AD
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Z_ms

  !-----------------------------------------------------------------------
  ! Initialization
  !-----------------------------------------------------------------------
  n_Layers = Atm%n_Layers
  CALL Set_TLAD_Limits()

  Transmittance_AD           = ZERO
  Reflectivity_AD            = ZERO
  Reflectivity_Attenuated_AD = ZERO
  ReflectivityLinear_AD            = ZERO
  Reflectivity_AttenuatedLinear_AD = ZERO
  Optical_Depth_AD           = ZERO

  !-----------------------------------------------------------------------
  ! Heights and geometry
  !-----------------------------------------------------------------------
  IF (ALL(Atm%Height .LT. EPSILON_FP)) THEN
      Height = Calculate_Height(Atm)
  ELSE
      Height = Atm%Height
  END IF

  dZ_m = Height(0:n_Layers-1) - Height(1:n_Layers)

  !-----------------------------------------------------------------------
  ! Frequency / wavelength
  !-----------------------------------------------------------------------
  IF (SpcCoeff_IsMicrowaveSensor(SC(SensorIndex))) THEN
      Frequency = SC(SensorIndex)%Frequency(ChannelIndex)
      Wavelength_m = POINT_01 / GHz_to_Inverse_cm(Frequency)
  ELSE IF (SpcCoeff_IsInfraredSensor(SC(SensorIndex))) THEN
      Wavenumber = SC(SensorIndex)%Wavenumber(ChannelIndex)
      Wavelength_m = POINT_01 / Wavenumber
  END IF

  !-----------------------------------------------------------------------
  ! Permittivity and forward reflectivity computation
  !-----------------------------------------------------------------------
  perm = Water_Permittivity_Turner_2016(Frequency * 1.0d9, &
                                        SPREAD(KW_TEMPERATURE_K, DIM=1, NCOPIES=n_Layers))
  perm_re = REAL(perm)
  perm_im = AIMAG(perm)
  Kw_2 = ABS((perm - ONE) / (perm + TWO))**TWO

  CALL Path_Optical_Depth(AtmOptics, GeometryInfo, Options, Optical_Depth)

  DO k = 1, n_Layers
      temp_sum = SUM(Optical_Depth(1:k))
      Transmittance(k) = EXP(-TWO * temp_sum)
  END DO

  P1 = (M6_MM6 * Wavelength_m**4.0_fp) / (PI**5.0_fp * Kw_2)
  Reflectivity = P1 * (AtmOptics%Backscat_Coefficient / dZ_m)
  Reflectivity_Attenuated = Transmittance * Reflectivity
  IF (MS_Algorithm_Id(Options) == ACTIVE_MS_TDTS) THEN
     CALL TDTS_Reflectivity(AtmOptics, GeometryInfo, Options, Height, P1, Z_ms)
     Reflectivity_Attenuated = Reflectivity_Attenuated + Z_ms
  END IF

  ! Linear reflectivities (cm^6/m^3)
  ReflectivityLinear = MM6_TO_CM6 * Reflectivity
  Reflectivity_AttenuatedLinear = MM6_TO_CM6 * Reflectivity_Attenuated

  !-----------------------------------------------------------------------
  ! ADJOINT PROPAGATION: determine which output adjoint is active
  !-----------------------------------------------------------------------

  ! === Case 1: Logarithmic attenuated reflectivity (dBZ)
  IF (ANY(RTSolution_AD%Reflectivity_Attenuated /= ZERO)) THEN
      WHERE (Reflectivity_Attenuated > REFLECTIVITY_TLAD_MIN)
        Reflectivity_Attenuated_AD = Reflectivity_Attenuated_AD + &
             TEN * RTSolution_AD%Reflectivity_Attenuated / &
             (MAX(Reflectivity_Attenuated, REFLECTIVITY_TLAD_CAP) * LOG(TEN))
      END WHERE

  END IF

  ! === Case 2: Linear attenuated reflectivity (cm^6/m^3)
  IF (ANY(RTSolution_AD%Reflectivity_AttenuatedLinear /= ZERO)) THEN
      Reflectivity_AttenuatedLinear_AD = Reflectivity_AttenuatedLinear_AD + RTSolution_AD%Reflectivity_AttenuatedLinear
      Reflectivity_Attenuated_AD = Reflectivity_Attenuated_AD + MM6_TO_CM6 * Reflectivity_AttenuatedLinear_AD
  END IF


  ! === Case 3: Logarithmic reflectivity (dBZ)
  IF (ANY(RTSolution_AD%Reflectivity /= ZERO)) THEN
      WHERE (Reflectivity > REFLECTIVITY_TLAD_MIN)
        Reflectivity_AD = Reflectivity_AD + &
             TEN * RTSolution_AD%Reflectivity / (MAX(Reflectivity, REFLECTIVITY_TLAD_CAP) * LOG(TEN))
      END WHERE

  END IF

  ! === Case 4: Linear reflectivity (cm^6/m^3)
  IF (ANY(RTSolution_AD%ReflectivityLinear /= ZERO)) THEN
      ReflectivityLinear_AD = ReflectivityLinear_AD + RTSolution_AD%ReflectivityLinear
      Reflectivity_AD = Reflectivity_AD + MM6_TO_CM6 * ReflectivityLinear_AD
  END IF

  !-----------------------------------------------------------------------
  ! Common adjoint propagation for reflectivity relationships
  !-----------------------------------------------------------------------
  Transmittance_AD = Transmittance_AD + Reflectivity * Reflectivity_Attenuated_AD
  Reflectivity_AD  = Reflectivity_AD  + Transmittance * Reflectivity_Attenuated_AD
  ! Wide-angle multiple scattering term
  IF (MS_Algorithm_Id(Options) == ACTIVE_MS_TDTS) &
     CALL TDTS_Reflectivity_AD(AtmOptics, GeometryInfo, Options, Height, P1, &
                               Reflectivity_Attenuated_AD, AtmOptics_AD)

  ! Backscatter coefficient adjoint
  AtmOptics_AD%Backscat_Coefficient = AtmOptics_AD%Backscat_Coefficient + Reflectivity_AD * P1 / dZ_m  

  !-----------------------------------------------------------------------
  ! Adjoint of transmittance → optical depth
  !-----------------------------------------------------------------------
  DO k = n_Layers, 1, -1
      DO j = 1, k
          Optical_Depth_AD(j) = Optical_Depth_AD(j) - TWO * Transmittance(k) * Transmittance_AD(k)
      END DO
  END DO

  !-----------------------------------------------------------------------
  ! Adjoint of the path optical depth (slant path + MS treatment)
  !-----------------------------------------------------------------------
  CALL Path_Optical_Depth_AD(AtmOptics, GeometryInfo, Options, Optical_Depth_AD, AtmOptics_AD)

END SUBROUTINE CRTM_Compute_Reflectivity_AD

!--------------------------------------------------------------------------------
!
! NAME:
!       Get_MS_Settings
!
! PURPOSE:
!       Return the active-sensor multiple-scattering algorithm id and the
!       forward-scattering fraction (eta) from the optional Options input.
!       Without Options the single-scattering treatment is used.
!
!--------------------------------------------------------------------------------
SUBROUTINE Get_MS_Settings(Options, ms_algorithm, eta)
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)  :: Options
  INTEGER,                           INTENT(OUT) :: ms_algorithm
  REAL(fp),                          INTENT(OUT) :: eta

  ms_algorithm = ACTIVE_MS_NONE
  eta          = ZERO
  IF (PRESENT(Options)) THEN
     ms_algorithm = INT(Options%Active_MS_Algorithm)
     eta          = MIN(MAX(REAL(Options%Active_MS_Eta, fp), ZERO), ONE)
  END IF

  SELECT CASE (ms_algorithm)
  CASE (ACTIVE_MS_NONE, ACTIVE_MS_EFFECTIVE_ATTENUATION, ACTIVE_MS_TDTS)
     ! supported
  CASE DEFAULT
     ms_algorithm = ACTIVE_MS_NONE
  END SELECT
END SUBROUTINE Get_MS_Settings

!--------------------------------------------------------------------------------
!
! NAME:
!       Path_Optical_Depth, Path_Optical_Depth_TL, Path_Optical_Depth_AD
!
! PURPOSE:
!       Layer optical depth along the radar path used for the two-way
!       attenuation, tau_path(k) = tau(k) / cos(theta).
!
!       ACTIVE_MS_NONE:
!         tau is AtmOptics%Optical_Depth (already delta-truncated in
!         CRTM_AtmOptics_Combine), i.e. the current CRTM behaviour.
!
!       ACTIVE_MS_EFFECTIVE_ATTENUATION:
!         A fraction eta of the scattered power is assumed to stay in the
!         radar beam (forward scattering), so only (1 - eta) of the
!         scattering optical depth attenuates the return:
!            tau_eff(k) = tau(k) * (1 - eta * w(k))
!         with w = AtmOptics%Single_Scatter_Albedo. eta = 0 recovers single
!         scattering; eta = 1 removes all scattering from the attenuation.
!
!--------------------------------------------------------------------------------
SUBROUTINE Path_Optical_Depth(AtmOptics, GeometryInfo, Options, tau)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN)  :: AtmOptics
  TYPE(CRTM_GeometryInfo_type),      INTENT(IN)  :: GeometryInfo
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)  :: Options
  REAL(fp),                          INTENT(OUT) :: tau(:)
  INTEGER  :: ms_algorithm, n
  REAL(fp) :: eta

  REAL(fp) :: w(AtmOptics%n_Layers), g(AtmOptics%n_Layers)

  CALL Get_MS_Settings(Options, ms_algorithm, eta)
  n = AtmOptics%n_Layers
  tau(1:n) = AtmOptics%Optical_Depth(1:n)
  IF (ms_algorithm == ACTIVE_MS_EFFECTIVE_ATTENUATION) THEN
     tau(1:n) = tau(1:n) * (ONE - eta * AtmOptics%Single_Scatter_Albedo(1:n))
  ELSE IF (ms_algorithm == ACTIVE_MS_TDTS) THEN
     CALL Unscaled_Optics(AtmOptics, tau(1:n), w, g)
  END IF
  tau(1:n) = tau(1:n) / GeometryInfo%Cosine_Sensor_Zenith
END SUBROUTINE Path_Optical_Depth

SUBROUTINE Path_Optical_Depth_TL(AtmOptics, AtmOptics_TL, GeometryInfo, Options, tau_TL)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN)  :: AtmOptics, AtmOptics_TL
  TYPE(CRTM_GeometryInfo_type),      INTENT(IN)  :: GeometryInfo
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)  :: Options
  REAL(fp),                          INTENT(OUT) :: tau_TL(:)
  INTEGER  :: ms_algorithm, n
  REAL(fp) :: eta

  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: w_TL, g_TL

  CALL Get_MS_Settings(Options, ms_algorithm, eta)
  n = AtmOptics%n_Layers
  tau_TL(1:n) = AtmOptics_TL%Optical_Depth(1:n)
  IF (ms_algorithm == ACTIVE_MS_EFFECTIVE_ATTENUATION) THEN
     tau_TL(1:n) = AtmOptics_TL%Optical_Depth(1:n) * (ONE - eta * AtmOptics%Single_Scatter_Albedo(1:n)) &
                 - eta * AtmOptics%Optical_Depth(1:n) * AtmOptics_TL%Single_Scatter_Albedo(1:n)
  ELSE IF (ms_algorithm == ACTIVE_MS_TDTS) THEN
     CALL Unscaled_Optics_TL(AtmOptics, AtmOptics_TL, tau_TL(1:n), w_TL, g_TL)
  END IF
  tau_TL(1:n) = tau_TL(1:n) / GeometryInfo%Cosine_Sensor_Zenith
END SUBROUTINE Path_Optical_Depth_TL

SUBROUTINE Path_Optical_Depth_AD(AtmOptics, GeometryInfo, Options, tau_AD, AtmOptics_AD)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN)     :: AtmOptics
  TYPE(CRTM_GeometryInfo_type),      INTENT(IN)     :: GeometryInfo
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options
  REAL(fp),                          INTENT(IN)     :: tau_AD(:)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN OUT) :: AtmOptics_AD
  INTEGER  :: ms_algorithm, n
  REAL(fp) :: eta
  REAL(fp) :: t_AD(SIZE(tau_AD))

  REAL(fp), DIMENSION(SIZE(tau_AD)) :: w_AD, g_AD

  CALL Get_MS_Settings(Options, ms_algorithm, eta)
  n = AtmOptics%n_Layers
  t_AD(1:n) = tau_AD(1:n) / GeometryInfo%Cosine_Sensor_Zenith
  IF (ms_algorithm == ACTIVE_MS_EFFECTIVE_ATTENUATION) THEN
     AtmOptics_AD%Single_Scatter_Albedo(1:n) = AtmOptics_AD%Single_Scatter_Albedo(1:n) &
          - eta * AtmOptics%Optical_Depth(1:n) * t_AD(1:n)
     t_AD(1:n) = t_AD(1:n) * (ONE - eta * AtmOptics%Single_Scatter_Albedo(1:n))
     AtmOptics_AD%Optical_Depth(1:n) = AtmOptics_AD%Optical_Depth(1:n) + t_AD(1:n)
  ELSE IF (ms_algorithm == ACTIVE_MS_TDTS) THEN
     w_AD = ZERO; g_AD = ZERO
     CALL Unscaled_Optics_AD(AtmOptics, t_AD(1:n), w_AD(1:n), g_AD(1:n), AtmOptics_AD)
  ELSE
     AtmOptics_AD%Optical_Depth(1:n) = AtmOptics_AD%Optical_Depth(1:n) + t_AD(1:n)
  END IF
END SUBROUTINE Path_Optical_Depth_AD

!--------------------------------------------------------------------------------
!
! NAME:
!       MS_Algorithm_Id
!
! PURPOSE:
!       Multiple-scattering algorithm selected in the optional Options input.
!
!--------------------------------------------------------------------------------
INTEGER FUNCTION MS_Algorithm_Id(Options)
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN) :: Options
  REAL(fp) :: eta
  CALL Get_MS_Settings(Options, MS_Algorithm_Id, eta)
END FUNCTION MS_Algorithm_Id

!--------------------------------------------------------------------------------
!
! NAME:
!       Unscaled_Optics, Unscaled_Optics_TL, Unscaled_Optics_AD
!
! PURPOSE:
!       Undo the delta truncation applied in CRTM_AtmOptics_Combine to get
!       the layer optical depth, single-scatter albedo and asymmetry factor
!       of the full phase function. With f = AtmOptics%Delta_Truncation,
!       tau' and w' the truncated values:
!          D   = 1 - f + f*w'
!          w   = w'/D
!          tau = tau'*D/(1 - f)
!          g   = c1/1.5, c1 = AtmOptics%Phase_Coefficient(1,1,k)
!       (the Legendre coefficients are normalised so that c0 = 1/2).
!       Layers without significant scattering keep tau' and get w = g = 0.
!
!--------------------------------------------------------------------------------
LOGICAL FUNCTION Has_Scattering(AtmOptics, k)
  TYPE(CRTM_AtmOptics_type), INTENT(IN) :: AtmOptics
  INTEGER,                   INTENT(IN) :: k
  Has_Scattering = AtmOptics%Include_Scattering .AND. AtmOptics%n_Legendre_Terms > 0 .AND. &
                   AtmOptics%Single_Scatter_Albedo(k) > BS_THRESHOLD_MS
END FUNCTION Has_Scattering

SUBROUTINE Unscaled_Optics(AtmOptics, tau, w, g)
  TYPE(CRTM_AtmOptics_type), INTENT(IN)  :: AtmOptics
  REAL(fp),                  INTENT(OUT) :: tau(:), w(:), g(:)
  INTEGER  :: k
  REAL(fp) :: f, wp, d
  DO k = 1, AtmOptics%n_Layers
     tau(k) = AtmOptics%Optical_Depth(k)
     w(k) = ZERO
     g(k) = ZERO
     IF (.NOT. Has_Scattering(AtmOptics, k)) CYCLE
     f  = AtmOptics%Delta_Truncation(k)
     wp = AtmOptics%Single_Scatter_Albedo(k)
     d  = ONE - f + f*wp
     w(k)   = wp/d
     tau(k) = AtmOptics%Optical_Depth(k)*d/(ONE - f)
     g(k)   = AtmOptics%Phase_Coefficient(1,1,k)/1.5_fp
  END DO
END SUBROUTINE Unscaled_Optics

SUBROUTINE Unscaled_Optics_TL(AtmOptics, AtmOptics_TL, tau_TL, w_TL, g_TL)
  TYPE(CRTM_AtmOptics_type), INTENT(IN)  :: AtmOptics, AtmOptics_TL
  REAL(fp),                  INTENT(OUT) :: tau_TL(:), w_TL(:), g_TL(:)
  INTEGER  :: k
  REAL(fp) :: f, wp, tp, d, f_TL, wp_TL, tp_TL, d_TL
  DO k = 1, AtmOptics%n_Layers
     tau_TL(k) = AtmOptics_TL%Optical_Depth(k)
     w_TL(k) = ZERO
     g_TL(k) = ZERO
     IF (.NOT. Has_Scattering(AtmOptics, k)) CYCLE
     f  = AtmOptics%Delta_Truncation(k);    f_TL  = AtmOptics_TL%Delta_Truncation(k)
     wp = AtmOptics%Single_Scatter_Albedo(k); wp_TL = AtmOptics_TL%Single_Scatter_Albedo(k)
     tp = AtmOptics%Optical_Depth(k);       tp_TL = AtmOptics_TL%Optical_Depth(k)
     d    = ONE - f + f*wp
     d_TL = f_TL*(wp - ONE) + f*wp_TL
     w_TL(k)   = wp_TL/d - wp*d_TL/(d*d)
     tau_TL(k) = (tp_TL*d + tp*d_TL)/(ONE - f) + tp*d*f_TL/((ONE - f)*(ONE - f))
     g_TL(k)   = AtmOptics_TL%Phase_Coefficient(1,1,k)/1.5_fp
  END DO
END SUBROUTINE Unscaled_Optics_TL

SUBROUTINE Unscaled_Optics_AD(AtmOptics, tau_AD, w_AD, g_AD, AtmOptics_AD)
  TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: AtmOptics
  REAL(fp),                  INTENT(IN)     :: tau_AD(:), w_AD(:), g_AD(:)
  TYPE(CRTM_AtmOptics_type), INTENT(IN OUT) :: AtmOptics_AD
  INTEGER  :: k
  REAL(fp) :: f, wp, tp, d, d_AD
  DO k = 1, AtmOptics%n_Layers
     IF (.NOT. Has_Scattering(AtmOptics, k)) THEN
        AtmOptics_AD%Optical_Depth(k) = AtmOptics_AD%Optical_Depth(k) + tau_AD(k)
        CYCLE
     END IF
     f  = AtmOptics%Delta_Truncation(k)
     wp = AtmOptics%Single_Scatter_Albedo(k)
     tp = AtmOptics%Optical_Depth(k)
     d  = ONE - f + f*wp
     AtmOptics_AD%Phase_Coefficient(1,1,k) = AtmOptics_AD%Phase_Coefficient(1,1,k) + g_AD(k)/1.5_fp
     ! tau = tp*d/(1-f)
     AtmOptics_AD%Optical_Depth(k) = AtmOptics_AD%Optical_Depth(k) + tau_AD(k)*d/(ONE - f)
     d_AD = tau_AD(k)*tp/(ONE - f)
     AtmOptics_AD%Delta_Truncation(k) = AtmOptics_AD%Delta_Truncation(k) &
                                      + tau_AD(k)*tp*d/((ONE - f)*(ONE - f))
     ! w = wp/d
     AtmOptics_AD%Single_Scatter_Albedo(k) = AtmOptics_AD%Single_Scatter_Albedo(k) + w_AD(k)/d
     d_AD = d_AD - w_AD(k)*wp/(d*d)
     ! d = 1 - f + f*wp
     AtmOptics_AD%Delta_Truncation(k) = AtmOptics_AD%Delta_Truncation(k) + d_AD*(wp - ONE)
     AtmOptics_AD%Single_Scatter_Albedo(k) = AtmOptics_AD%Single_Scatter_Albedo(k) + d_AD*f
  END DO
END SUBROUTINE Unscaled_Optics_AD

!--------------------------------------------------------------------------------
!
! NAME:
!       TDTS_Setup
!
! PURPOSE:
!       Radar geometry and unscaled optics for the TDTS calculation:
!         edge(0:n)  distance from the radar to the layer boundaries along
!                    the beam (m), flat-earth approximation
!         ext(n)     extinction coefficient (m-1) = tau/dz
!         rho        1/e half-width of the antenna pattern (rad), from the
!                    full 3-dB beamwidth as rho = theta_3dB/(2 sqrt(ln 2))
!
!--------------------------------------------------------------------------------
SUBROUTINE TDTS_Setup(AtmOptics, GeometryInfo, Options, Height, edge, ext, w, g, rho, cbe)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN)  :: AtmOptics
  TYPE(CRTM_GeometryInfo_type),      INTENT(IN)  :: GeometryInfo
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)  :: Options
  REAL(fp),                          INTENT(IN)  :: Height(0:)
  REAL(fp),                          INTENT(OUT) :: edge(0:), ext(:), w(:), g(:), rho, cbe
  REAL(fp) :: altitude, beamwidth, tau(AtmOptics%n_Layers)
  INTEGER  :: n
  TYPE(CRTM_Options_type) :: Default_Options

  IF (PRESENT(Options)) THEN
     altitude  = REAL(Options%Active_MS_Altitude, fp)
     beamwidth = REAL(Options%Active_MS_Beamwidth, fp)
     cbe       = REAL(Options%Active_MS_CBE, fp)
  ELSE
     altitude  = REAL(Default_Options%Active_MS_Altitude, fp)
     beamwidth = REAL(Default_Options%Active_MS_Beamwidth, fp)
     cbe       = REAL(Default_Options%Active_MS_CBE, fp)
  END IF
  rho = (beamwidth*PI/180.0_fp)/(TWO*SQRT(LOG(TWO)))
  n = AtmOptics%n_Layers
  edge(0:n) = (altitude - Height(0:n))/GeometryInfo%Cosine_Sensor_Zenith
  CALL Unscaled_Optics(AtmOptics, tau, w(1:n), g(1:n))
  ext(1:n) = tau/(Height(0:n-1) - Height(1:n))
END SUBROUTINE TDTS_Setup

!--------------------------------------------------------------------------------
!
! NAME:
!       TDTS_Reflectivity, TDTS_Reflectivity_TL, TDTS_Reflectivity_AD
!
! PURPOSE:
!       Apparent reflectivity (mm^6 m^-3) from wide-angle multiple scattering,
!       Z_ms = P1 * 4*pi * beta_ms, where beta_ms (m-1 sr-1) is the TDTS
!       apparent backscatter and P1*4*pi converts backscatter to Z in the same
!       way as P1 converts AtmOptics%Backscat_Coefficient/dz (Hogan and
!       Battaglia 2008, Eq. 1).
!
!--------------------------------------------------------------------------------
SUBROUTINE TDTS_Reflectivity(AtmOptics, GeometryInfo, Options, Height, P1, Z_ms)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN)  :: AtmOptics
  TYPE(CRTM_GeometryInfo_type),      INTENT(IN)  :: GeometryInfo
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)  :: Options
  REAL(fp),                          INTENT(IN)  :: Height(0:), P1(:)
  REAL(fp),                          INTENT(OUT) :: Z_ms(:)
  INTEGER  :: n
  REAL(fp) :: rho, cbe
  REAL(fp) :: edge(0:AtmOptics%n_Layers)
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ext, w, g, bscat

  n = AtmOptics%n_Layers
  CALL TDTS_Setup(AtmOptics, GeometryInfo, Options, Height, edge, ext, w, g, rho, cbe)
  CALL TDTS_Wide_Angle(n, edge, ext, w, g, rho, rho, cbe, bscat)
  Z_ms(1:n) = P1(1:n)*FOUR_PI*bscat
END SUBROUTINE TDTS_Reflectivity

SUBROUTINE TDTS_Reflectivity_TL(AtmOptics, AtmOptics_TL, GeometryInfo, Options, Height, P1, Z_ms, Z_ms_TL)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN)  :: AtmOptics, AtmOptics_TL
  TYPE(CRTM_GeometryInfo_type),      INTENT(IN)  :: GeometryInfo
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)  :: Options
  REAL(fp),                          INTENT(IN)  :: Height(0:), P1(:)
  REAL(fp),                          INTENT(OUT) :: Z_ms(:), Z_ms_TL(:)
  INTEGER  :: n
  REAL(fp) :: rho, cbe
  REAL(fp) :: edge(0:AtmOptics%n_Layers)
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ext, w, g, bscat, tau_TL, ext_TL, w_TL, g_TL, bscat_TL

  n = AtmOptics%n_Layers
  CALL TDTS_Setup(AtmOptics, GeometryInfo, Options, Height, edge, ext, w, g, rho, cbe)
  CALL Unscaled_Optics_TL(AtmOptics, AtmOptics_TL, tau_TL, w_TL, g_TL)
  ext_TL = tau_TL/(Height(0:n-1) - Height(1:n))
  CALL TDTS_Wide_Angle_TL(n, edge, ext, w, g, rho, rho, cbe, ext_TL, w_TL, g_TL, bscat, bscat_TL)
  Z_ms(1:n)    = P1(1:n)*FOUR_PI*bscat
  Z_ms_TL(1:n) = P1(1:n)*FOUR_PI*bscat_TL
END SUBROUTINE TDTS_Reflectivity_TL

SUBROUTINE TDTS_Reflectivity_AD(AtmOptics, GeometryInfo, Options, Height, P1, Z_ms_AD, AtmOptics_AD)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN)     :: AtmOptics
  TYPE(CRTM_GeometryInfo_type),      INTENT(IN)     :: GeometryInfo
  TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options
  REAL(fp),                          INTENT(IN)     :: Height(0:), P1(:), Z_ms_AD(:)
  TYPE(CRTM_AtmOptics_type),         INTENT(IN OUT) :: AtmOptics_AD
  INTEGER  :: n
  REAL(fp) :: rho, cbe
  REAL(fp) :: edge(0:AtmOptics%n_Layers)
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ext, w, g, bscat_AD, ext_AD, w_AD, g_AD

  n = AtmOptics%n_Layers
  CALL TDTS_Setup(AtmOptics, GeometryInfo, Options, Height, edge, ext, w, g, rho, cbe)
  bscat_AD = P1(1:n)*FOUR_PI*Z_ms_AD(1:n)
  ext_AD = ZERO; w_AD = ZERO; g_AD = ZERO
  CALL TDTS_Wide_Angle_AD(n, edge, ext, w, g, rho, rho, cbe, bscat_AD, ext_AD, w_AD, g_AD)
  CALL Unscaled_Optics_AD(AtmOptics, ext_AD/(Height(0:n-1) - Height(1:n)), w_AD, g_AD, AtmOptics_AD)
END SUBROUTINE TDTS_Reflectivity_AD

! Read the dBZ TL/AD limits from the environment (once)
SUBROUTINE Set_TLAD_Limits()
  CHARACTER(LEN=32) :: env_val
  INTEGER :: env_len, env_stat, io_stat
  REAL(fp) :: dbz

  IF (TLAD_LIMITS_SET) RETURN
!$OMP CRITICAL (crtm_active_tlad_limits)
  IF (.NOT. TLAD_LIMITS_SET) THEN
  CALL GET_ENVIRONMENT_VARIABLE('CRTM_ACTIVE_TLAD_MIN_DBZ', env_val, env_len, env_stat)
  IF (env_stat == 0 .AND. env_len > 0) THEN
     READ(env_val, *, IOSTAT=io_stat) dbz
     IF (io_stat == 0) REFLECTIVITY_TLAD_MIN = TEN**(dbz / TEN)
  END IF
  REFLECTIVITY_TLAD_CAP = REFLECTIVITY_TLAD_MIN
  CALL GET_ENVIRONMENT_VARIABLE('CRTM_ACTIVE_TLAD_CAP_DBZ', env_val, env_len, env_stat)
  IF (env_stat == 0 .AND. env_len > 0) THEN
     READ(env_val, *, IOSTAT=io_stat) dbz
     IF (io_stat == 0) REFLECTIVITY_TLAD_CAP = MAX(TEN**(dbz / TEN), REFLECTIVITY_TLAD_MIN)
  END IF
  TLAD_LIMITS_SET = .TRUE.
  END IF
!$OMP END CRITICAL (crtm_active_tlad_limits)
END SUBROUTINE Set_TLAD_Limits

END MODULE CRTM_Active_Sensor
