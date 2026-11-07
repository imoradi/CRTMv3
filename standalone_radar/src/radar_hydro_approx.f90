!
! radar_hydro_approx
!
! Approximate specific attenuation by hydrometeors, using formulas commonly
! used in the radar and propagation communities. They need only the water
! content (and temperature), not the scattering properties:
!
!   APPROX_CLOUD_LIQUID : Rayleigh absorption by liquid droplets,
!                         Recommendation ITU-R P.840-8 (double-Debye water
!                         permittivity, Eqs. 2-11).
!   APPROX_ICE          : Rayleigh absorption by ice (cloud ice, snow, graupel,
!                         hail), same form with the ice permittivity of
!                         Hufford (1991) as given by Maetzler (2006), mass
!                         converted to volume with a solid ice density of
!                         917 kg m-3. Scattering by large snow particles is
!                         not included, so at W band this underestimates the
!                         attenuation in heavy snow.
!   APPROX_RAIN         : Recommendation ITU-R P.838-3, gamma = k R^alpha
!                         (dB/km), with circular-polarisation coefficients
!                         k=(kH+kV)/2, alpha=(kH aH+kV aV)/(2k), and the rain
!                         rate from the rain water content with the
!                         Marshall-Palmer relation M = 0.072 R^0.88
!                         (M in g m-3, R in mm h-1).
!
! Written in complex arithmetic so derivatives are taken with a complex step.
! Inputs: f (GHz), water content rho (kg m-3), temperature t (K).
! Output: one-way absorption coefficient (Np m-1).
!
MODULE radar_hydro_approx
  USE radar_kinds, ONLY: fp, DB_TO_NP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: APPROX_NONE, APPROX_CLOUD_LIQUID, APPROX_ICE, APPROX_RAIN
  PUBLIC :: hydro_approx_attenuation_c
  PUBLIC :: water_eps_parts, ice_permittivity_c
  PUBLIC :: rain_coefficients_p838

  INTEGER, PARAMETER :: APPROX_NONE         = 0
  INTEGER, PARAMETER :: APPROX_CLOUD_LIQUID = 1
  INTEGER, PARAMETER :: APPROX_ICE          = 2
  INTEGER, PARAMETER :: APPROX_RAIN         = 3

  REAL(fp), PARAMETER :: RHO_ICE = 0.917_fp   ! g cm-3

CONTAINS

  ! ITU-R P.840-8 double-Debye permittivity of liquid water (real, imaginary parts)
  PURE SUBROUTINE water_eps_parts(f, t, ep, epp)
    REAL(fp),    INTENT(IN)  :: f
    COMPLEX(fp), INTENT(IN)  :: t
    COMPLEX(fp), INTENT(OUT) :: ep, epp
    COMPLEX(fp) :: th, e0, e1, fp1, fs
    REAL(fp), PARAMETER :: e2 = 3.52_fp
    th  = 300.0_fp/t
    e0  = 77.66_fp + 103.3_fp*(th - 1.0_fp)
    e1  = 0.0671_fp*e0
    fp1 = 20.20_fp - 146.0_fp*(th - 1.0_fp) + 316.0_fp*(th - 1.0_fp)**2
    fs  = 39.8_fp*fp1
    ep  = (e0 - e1)/(1.0_fp + (f/fp1)**2) + (e1 - e2)/(1.0_fp + (f/fs)**2) + e2
    epp = f*(e0 - e1)/(fp1*(1.0_fp + (f/fp1)**2)) + f*(e1 - e2)/(fs*(1.0_fp + (f/fs)**2))
  END SUBROUTINE water_eps_parts

  ! Ice permittivity, Hufford (1991) / Maetzler (2006). Returns real and
  ! imaginary parts as two complex-step carriers.
  PURE FUNCTION ice_permittivity_c(f, t) RESULT(eps)
    REAL(fp),    INTENT(IN) :: f
    COMPLEX(fp), INTENT(IN) :: t
    COMPLEX(fp) :: eps(2)
    COMPLEX(fp) :: th, a, b, x
    th = 300.0_fp/t - 1.0_fp
    a  = (0.00504_fp + 0.0062_fp*th)*EXP(-22.1_fp*th)
    x  = EXP(335.0_fp/t)
    b  = (0.0207_fp/t)*x/(x - 1.0_fp)**2 + 1.16e-11_fp*f**2 + EXP(-9.963_fp + 0.0372_fp*(t - 273.16_fp))
    eps(1) = 3.1884_fp + 9.1e-4_fp*(t - 273.16_fp)
    eps(2) = a/f + b*f
  END FUNCTION ice_permittivity_c

  ! ITU-R P.838-3 coefficients for circular polarisation
  PURE SUBROUTINE rain_coefficients_p838(f, k, alpha)
    REAL(fp), INTENT(IN)  :: f
    REAL(fp), INTENT(OUT) :: k, alpha
    REAL(fp), PARAMETER :: kh_a(4) = [-5.33980_fp, -0.35351_fp, -0.23789_fp, -0.94158_fp]
    REAL(fp), PARAMETER :: kh_b(4) = [-0.10008_fp,  1.26970_fp,  0.86036_fp,  0.64552_fp]
    REAL(fp), PARAMETER :: kh_c(4) = [ 1.13098_fp,  0.45400_fp,  0.15354_fp,  0.16817_fp]
    REAL(fp), PARAMETER :: kh_m = -0.18961_fp, kh_k = 0.71147_fp
    REAL(fp), PARAMETER :: kv_a(4) = [-3.80595_fp, -3.44965_fp, -0.39902_fp,  0.50167_fp]
    REAL(fp), PARAMETER :: kv_b(4) = [ 0.56934_fp, -0.22911_fp,  0.73042_fp,  1.07319_fp]
    REAL(fp), PARAMETER :: kv_c(4) = [ 0.81061_fp,  0.51059_fp,  0.11899_fp,  0.27195_fp]
    REAL(fp), PARAMETER :: kv_m = -0.16398_fp, kv_k = 0.63297_fp
    REAL(fp), PARAMETER :: ah_a(5) = [-0.14318_fp, 0.29591_fp, 0.32177_fp, -5.37610_fp, 16.1721_fp]
    REAL(fp), PARAMETER :: ah_b(5) = [ 1.82442_fp, 0.77564_fp, 0.63773_fp, -0.96230_fp, -3.29980_fp]
    REAL(fp), PARAMETER :: ah_c(5) = [-0.55187_fp, 0.19822_fp, 0.13164_fp,  1.47828_fp,  3.43990_fp]
    REAL(fp), PARAMETER :: ah_m = 0.67849_fp, ah_k = -1.95537_fp
    REAL(fp), PARAMETER :: av_a(5) = [-0.07771_fp, 0.56727_fp, -0.20238_fp, -48.2991_fp, 48.5833_fp]
    REAL(fp), PARAMETER :: av_b(5) = [ 2.33840_fp, 0.95545_fp,  1.14520_fp,  0.791669_fp, 0.791459_fp]
    REAL(fp), PARAMETER :: av_c(5) = [-0.76284_fp, 0.54039_fp,  0.26809_fp,  0.116226_fp, 0.116479_fp]
    REAL(fp), PARAMETER :: av_m = -0.053739_fp, av_k = 0.83433_fp
    REAL(fp) :: lf, kh, kv, ah, av
    lf = LOG10(f)
    kh = 10.0_fp**(SUM(kh_a*EXP(-((lf - kh_b)/kh_c)**2)) + kh_m*lf + kh_k)
    kv = 10.0_fp**(SUM(kv_a*EXP(-((lf - kv_b)/kv_c)**2)) + kv_m*lf + kv_k)
    ah = SUM(ah_a*EXP(-((lf - ah_b)/ah_c)**2)) + ah_m*lf + ah_k
    av = SUM(av_a*EXP(-((lf - av_b)/av_c)**2)) + av_m*lf + av_k
    k = 0.5_fp*(kh + kv)
    alpha = (kh*ah + kv*av)/(2.0_fp*k)
  END SUBROUTINE rain_coefficients_p838

  !
  ! Specific attenuation (Np m-1) of hydrometeor class 'kind'
  !
  PURE FUNCTION hydro_approx_attenuation_c(kind, f, rho, t) RESULT(alpha)
    INTEGER,     INTENT(IN) :: kind
    REAL(fp),    INTENT(IN) :: f
    COMPLEX(fp), INTENT(IN) :: rho, t
    COMPLEX(fp) :: alpha
    COMPLEX(fp) :: ep, epp, eps(2), m_gm3, rr
    REAL(fp) :: k, a
    m_gm3 = 1000.0_fp*rho
    alpha = (0.0_fp, 0.0_fp)
    IF (REAL(rho) <= 0.0_fp) RETURN
    SELECT CASE (kind)
    CASE (APPROX_CLOUD_LIQUID)
      CALL water_eps_parts(f, t, ep, epp)
      alpha = 0.819_fp*f*epp/((2.0_fp + ep)**2 + epp**2)*m_gm3
    CASE (APPROX_ICE)
      eps = ice_permittivity_c(f, t)
      alpha = 0.819_fp*f*eps(2)/((2.0_fp + eps(1))**2 + eps(2)**2)/RHO_ICE*m_gm3
    CASE (APPROX_RAIN)
      CALL rain_coefficients_p838(f, k, a)
      rr = (m_gm3/0.072_fp)**(1.0_fp/0.88_fp)
      alpha = k*rr**a
    END SELECT
    alpha = alpha*DB_TO_NP*1.0e-3_fp
  END FUNCTION hydro_approx_attenuation_c

END MODULE radar_hydro_approx
