!
! radar_kinds
!
! Kinds and physical constants for the standalone radar operator.
!
MODULE radar_kinds
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER, PUBLIC :: fp = SELECTED_REAL_KIND(15, 300)

  REAL(fp), PARAMETER, PUBLIC :: ZERO = 0.0_fp, ONE = 1.0_fp, TWO = 2.0_fp, TEN = 10.0_fp
  REAL(fp), PARAMETER, PUBLIC :: HALF = 0.5_fp
  REAL(fp), PARAMETER, PUBLIC :: PI = 3.14159265358979323846_fp
  REAL(fp), PARAMETER, PUBLIC :: FOUR_PI = 4.0_fp*PI
  ! Speed of light (m/s)
  REAL(fp), PARAMETER, PUBLIC :: C_LIGHT = 2.99792458e8_fp
  ! Gravity (m/s2)
  REAL(fp), PARAMETER, PUBLIC :: GRAV = 9.80665_fp
  ! Gas constant of dry air (J/kg/K) and Rv/Rd
  REAL(fp), PARAMETER, PUBLIC :: R_DRY = 287.05_fp
  REAL(fp), PARAMETER, PUBLIC :: EPS_WV = 0.621970_fp
  ! dB <-> neper: 1 dB/km of power = DB_TO_NP neper/km of optical depth
  REAL(fp), PARAMETER, PUBLIC :: DB_TO_NP = 0.230258509299404568_fp   ! ln(10)/10
  ! m^6 -> mm^6
  REAL(fp), PARAMETER, PUBLIC :: M6_TO_MM6 = 1.0e18_fp
  ! Complex-step size used for derivatives of local (single layer) functions
  REAL(fp), PARAMETER, PUBLIC :: CSTEP = 1.0e-30_fp
END MODULE radar_kinds
