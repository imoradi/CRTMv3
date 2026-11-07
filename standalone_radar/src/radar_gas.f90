!
! radar_gas
!
! Specific attenuation of microwaves by atmospheric gases (O2, N2 and H2O).
!
!   GAS_MODEL_ITU_LBL    : Recommendation ITU-R P.676-12 (08/2019), Annex 1.
!                          Line-by-line sum of 44 oxygen and 35 water vapour
!                          lines plus the dry-air (Debye + N2 pressure-induced)
!                          continuum. Valid 1-1000 GHz.
!   GAS_MODEL_ITU_APPROX : Recommendation ITU-R P.676-10 (09/2013), Annex 2.
!                          Closed-form approximations of the dry-air and water
!                          vapour attenuation commonly used in radar
!                          applications. Valid 1-350 GHz.
!
! The routines are written in complex arithmetic so that the derivatives
! with respect to temperature and humidity are obtained with a complex step
! (exact to machine precision). Inputs:
!   f    frequency (GHz)
!   p    total pressure (hPa)
!   t    temperature (K)
!   q    specific humidity (kg/kg)
! Output: one-way power absorption coefficient (Np/m).
!
MODULE radar_gas
  USE radar_kinds, ONLY: fp, DB_TO_NP, EPS_WV, CSTEP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GAS_MODEL_NONE, GAS_MODEL_ITU_LBL, GAS_MODEL_ITU_APPROX, GAS_MODEL_USER
  PUBLIC :: gas_absorption
  PUBLIC :: gas_absorption_c
  PUBLIC :: gas_absorption_deriv

  INTEGER, PARAMETER :: GAS_MODEL_NONE       = 0
  INTEGER, PARAMETER :: GAS_MODEL_ITU_LBL    = 1
  INTEGER, PARAMETER :: GAS_MODEL_ITU_APPROX = 2
  INTEGER, PARAMETER :: GAS_MODEL_USER       = 3   ! layer optical depth supplied by the caller

  ! ITU-R P.676-12 Tables 1 and 2: f0, a1..a6 (oxygen) and f0, b1..b6 (water vapour)
  INTEGER, PARAMETER :: N_O2 = 44
  REAL(fp), PARAMETER :: O2_C0(N_O2) = [ &
       50.474214_fp, 50.987745_fp, 51.50336_fp, 52.021429_fp, 52.542418_fp, &
       53.066934_fp, 53.595775_fp, 54.130025_fp, 54.67118_fp, 55.221384_fp, &
       55.783815_fp, 56.264774_fp, 56.363399_fp, 56.968211_fp, 57.612486_fp, &
       58.323877_fp, 58.446588_fp, 59.164204_fp, 59.590983_fp, 60.306056_fp, &
       60.434778_fp, 61.150562_fp, 61.800158_fp, 62.41122_fp, 62.486253_fp, &
       62.997984_fp, 63.568526_fp, 64.127775_fp, 64.67891_fp, 65.224078_fp, &
       65.764779_fp, 66.302096_fp, 66.836834_fp, 67.369601_fp, 67.900868_fp, &
       68.431006_fp, 68.960312_fp, 118.750334_fp, 368.498246_fp, 424.76302_fp, &
       487.249273_fp, 715.392902_fp, 773.83949_fp, 834.145546_fp ]
  REAL(fp), PARAMETER :: O2_C1(N_O2) = [ &
       0.975_fp, 2.529_fp, 6.193_fp, 14.32_fp, 31.24_fp, 64.29_fp, 124.6_fp, &
       227.3_fp, 389.7_fp, 627.1_fp, 945.3_fp, 543.4_fp, 1331.8_fp, 1746.6_fp, &
       2120.1_fp, 2363.7_fp, 1442.1_fp, 2379.9_fp, 2090.7_fp, 2103.4_fp, &
       2438.0_fp, 2479.5_fp, 2275.9_fp, 1915.4_fp, 1503.0_fp, 1490.2_fp, &
       1078.0_fp, 728.7_fp, 461.3_fp, 274.0_fp, 153.0_fp, 80.4_fp, 39.8_fp, &
       18.56_fp, 8.172_fp, 3.397_fp, 1.334_fp, 940.3_fp, 67.4_fp, 637.7_fp, &
       237.4_fp, 98.1_fp, 572.3_fp, 183.1_fp ]
  REAL(fp), PARAMETER :: O2_C2(N_O2) = [ &
       9.651_fp, 8.653_fp, 7.709_fp, 6.819_fp, 5.983_fp, 5.201_fp, 4.474_fp, &
       3.8_fp, 3.182_fp, 2.618_fp, 2.109_fp, 0.014_fp, 1.654_fp, 1.255_fp, &
       0.91_fp, 0.621_fp, 0.083_fp, 0.387_fp, 0.207_fp, 0.207_fp, 0.386_fp, &
       0.621_fp, 0.91_fp, 1.255_fp, 0.083_fp, 1.654_fp, 2.108_fp, 2.617_fp, &
       3.181_fp, 3.8_fp, 4.473_fp, 5.2_fp, 5.982_fp, 6.818_fp, 7.708_fp, &
       8.652_fp, 9.65_fp, 0.01_fp, 0.048_fp, 0.044_fp, 0.049_fp, 0.145_fp, &
       0.141_fp, 0.145_fp ]
  REAL(fp), PARAMETER :: O2_C3(N_O2) = [ &
       6.69_fp, 7.17_fp, 7.64_fp, 8.11_fp, 8.58_fp, 9.06_fp, 9.55_fp, 9.96_fp, &
       10.37_fp, 10.89_fp, 11.34_fp, 17.03_fp, 11.89_fp, 12.23_fp, 12.62_fp, &
       12.95_fp, 14.91_fp, 13.53_fp, 14.08_fp, 14.15_fp, 13.39_fp, 12.92_fp, &
       12.63_fp, 12.17_fp, 15.13_fp, 11.74_fp, 11.34_fp, 10.88_fp, 10.38_fp, &
       9.96_fp, 9.55_fp, 9.06_fp, 8.58_fp, 8.11_fp, 7.64_fp, 7.17_fp, 6.69_fp, &
       16.64_fp, 16.4_fp, 16.4_fp, 16.0_fp, 16.0_fp, 16.2_fp, 14.7_fp ]
  REAL(fp), PARAMETER :: O2_C4(N_O2) = [ &
       0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
       0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
       0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
       0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
       0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp ]
  REAL(fp), PARAMETER :: O2_C5(N_O2) = [ &
       2.566_fp, 2.246_fp, 1.947_fp, 1.667_fp, 1.388_fp, 1.349_fp, 2.227_fp, &
       3.17_fp, 3.558_fp, 2.56_fp, -1.172_fp, 3.525_fp, -2.378_fp, -3.545_fp, &
       -5.416_fp, -1.932_fp, 6.768_fp, -6.561_fp, 6.957_fp, -6.395_fp, &
       6.342_fp, 1.014_fp, 5.014_fp, 3.029_fp, -4.499_fp, 1.856_fp, 0.658_fp, &
       -3.036_fp, -3.968_fp, -3.528_fp, -2.548_fp, -1.66_fp, -1.68_fp, &
       -1.956_fp, -2.216_fp, -2.492_fp, -2.773_fp, -0.439_fp, 0.0_fp, 0.0_fp, &
       0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp ]
  REAL(fp), PARAMETER :: O2_C6(N_O2) = [ &
       6.85_fp, 6.8_fp, 6.729_fp, 6.64_fp, 6.526_fp, 6.206_fp, 5.085_fp, &
       3.75_fp, 2.654_fp, 2.952_fp, 6.135_fp, -0.978_fp, 6.547_fp, 6.451_fp, &
       6.056_fp, 0.436_fp, -1.273_fp, 2.309_fp, -0.776_fp, 0.699_fp, &
       -2.825_fp, -0.584_fp, -6.619_fp, -6.759_fp, 0.844_fp, -6.675_fp, &
       -6.139_fp, -2.895_fp, -2.59_fp, -3.68_fp, -5.002_fp, -6.091_fp, &
       -6.393_fp, -6.475_fp, -6.545_fp, -6.6_fp, -6.65_fp, 0.079_fp, 0.0_fp, &
       0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp ]
  INTEGER, PARAMETER :: N_H2O = 35
  REAL(fp), PARAMETER :: H2O_C0(N_H2O) = [ &
       22.23508_fp, 67.80396_fp, 119.99594_fp, 183.310087_fp, 321.22563_fp, &
       325.152888_fp, 336.227764_fp, 380.197353_fp, 390.134508_fp, &
       437.346667_fp, 439.150807_fp, 443.018343_fp, 448.001085_fp, &
       470.888999_fp, 474.689092_fp, 488.490108_fp, 503.568532_fp, &
       504.482692_fp, 547.67644_fp, 552.02096_fp, 556.935985_fp, &
       620.700807_fp, 645.766085_fp, 658.00528_fp, 752.033113_fp, &
       841.051732_fp, 859.965698_fp, 899.303175_fp, 902.611085_fp, &
       906.205957_fp, 916.171582_fp, 923.112692_fp, 970.315022_fp, &
       987.926764_fp, 1780.0_fp ]
  REAL(fp), PARAMETER :: H2O_C1(N_H2O) = [ &
       0.1079_fp, 0.0011_fp, 0.0007_fp, 2.273_fp, 0.047_fp, 1.514_fp, &
       0.001_fp, 11.67_fp, 0.0045_fp, 0.0632_fp, 0.9098_fp, 0.192_fp, &
       10.41_fp, 0.3254_fp, 1.26_fp, 0.2529_fp, 0.0372_fp, 0.0124_fp, &
       0.9785_fp, 0.184_fp, 497.0_fp, 5.015_fp, 0.0067_fp, 0.2732_fp, &
       243.4_fp, 0.0134_fp, 0.1325_fp, 0.0547_fp, 0.0386_fp, 0.1836_fp, &
       8.4_fp, 0.0079_fp, 9.009_fp, 134.6_fp, 17506.0_fp ]
  REAL(fp), PARAMETER :: H2O_C2(N_H2O) = [ &
       2.144_fp, 8.732_fp, 8.353_fp, 0.668_fp, 6.179_fp, 1.541_fp, 9.825_fp, &
       1.048_fp, 7.347_fp, 5.048_fp, 3.595_fp, 5.048_fp, 1.405_fp, 3.597_fp, &
       2.379_fp, 2.852_fp, 6.731_fp, 6.731_fp, 0.158_fp, 0.158_fp, 0.159_fp, &
       2.391_fp, 8.633_fp, 7.816_fp, 0.396_fp, 8.177_fp, 8.055_fp, 7.914_fp, &
       8.429_fp, 5.11_fp, 1.441_fp, 10.293_fp, 1.919_fp, 0.257_fp, 0.952_fp ]
  REAL(fp), PARAMETER :: H2O_C3(N_H2O) = [ &
       26.38_fp, 28.58_fp, 29.48_fp, 29.06_fp, 24.04_fp, 28.23_fp, 26.93_fp, &
       28.11_fp, 21.52_fp, 18.45_fp, 20.07_fp, 15.55_fp, 25.64_fp, 21.34_fp, &
       23.2_fp, 25.86_fp, 16.12_fp, 16.12_fp, 26.0_fp, 26.0_fp, 30.86_fp, &
       24.38_fp, 18.0_fp, 32.1_fp, 30.86_fp, 15.9_fp, 30.6_fp, 29.85_fp, &
       28.65_fp, 24.08_fp, 26.73_fp, 29.0_fp, 25.5_fp, 29.85_fp, 196.3_fp ]
  REAL(fp), PARAMETER :: H2O_C4(N_H2O) = [ &
       0.76_fp, 0.69_fp, 0.7_fp, 0.77_fp, 0.67_fp, 0.64_fp, 0.69_fp, 0.54_fp, &
       0.63_fp, 0.6_fp, 0.63_fp, 0.6_fp, 0.66_fp, 0.66_fp, 0.65_fp, 0.69_fp, &
       0.61_fp, 0.61_fp, 0.7_fp, 0.7_fp, 0.69_fp, 0.71_fp, 0.6_fp, 0.69_fp, &
       0.68_fp, 0.33_fp, 0.68_fp, 0.68_fp, 0.7_fp, 0.7_fp, 0.7_fp, 0.7_fp, &
       0.64_fp, 0.68_fp, 2.0_fp ]
  REAL(fp), PARAMETER :: H2O_C5(N_H2O) = [ &
       5.087_fp, 4.93_fp, 4.78_fp, 5.022_fp, 4.398_fp, 4.893_fp, 4.74_fp, &
       5.063_fp, 4.81_fp, 4.23_fp, 4.483_fp, 5.083_fp, 5.028_fp, 4.506_fp, &
       4.804_fp, 5.201_fp, 3.98_fp, 4.01_fp, 4.5_fp, 4.5_fp, 4.552_fp, &
       4.856_fp, 4.0_fp, 4.14_fp, 4.352_fp, 5.76_fp, 4.09_fp, 4.53_fp, 5.1_fp, &
       4.7_fp, 5.15_fp, 5.0_fp, 4.94_fp, 4.55_fp, 24.15_fp ]
  REAL(fp), PARAMETER :: H2O_C6(N_H2O) = [ &
       1.0_fp, 0.82_fp, 0.79_fp, 0.85_fp, 0.54_fp, 0.74_fp, 0.61_fp, 0.89_fp, &
       0.55_fp, 0.48_fp, 0.52_fp, 0.5_fp, 0.67_fp, 0.65_fp, 0.64_fp, 0.72_fp, &
       0.43_fp, 0.45_fp, 1.0_fp, 1.0_fp, 1.0_fp, 0.68_fp, 0.5_fp, 1.0_fp, &
       0.84_fp, 0.45_fp, 0.84_fp, 0.9_fp, 0.95_fp, 0.53_fp, 0.78_fp, 0.8_fp, &
       0.67_fp, 0.9_fp, 5.0_fp ]

CONTAINS

  ! Vapour partial pressure (hPa) and dry pressure from total pressure and q
  ELEMENTAL SUBROUTINE vapour_pressure(p, q, e, pd)
    COMPLEX(fp), INTENT(IN)  :: p, q
    COMPLEX(fp), INTENT(OUT) :: e, pd
    e  = q*p/(EPS_WV + (1.0_fp - EPS_WV)*q)
    pd = p - e
  END SUBROUTINE vapour_pressure

  ! ITU-R P.676-12 Annex 1, dry air (oxygen lines + dry continuum), dB/km
  PURE FUNCTION gamma_o_lbl(f, pd, e, t) RESULT(gam)
    REAL(fp),    INTENT(IN) :: f
    COMPLEX(fp), INTENT(IN) :: pd, e, t
    COMPLEX(fp) :: gam
    COMPLEX(fp) :: th, s, df, dl, ff, d, nd, nsum
    INTEGER :: i
    th = 300.0_fp/t
    nsum = (0.0_fp, 0.0_fp)
    DO i = 1, N_O2
      s  = O2_C1(i)*1.0e-7_fp*pd*th**3*EXP(O2_C2(i)*(1.0_fp - th))
      df = O2_C3(i)*1.0e-4_fp*(pd*th**(0.8_fp - O2_C4(i)) + 1.1_fp*e*th)
      df = SQRT(df*df + 2.25e-6_fp)
      dl = (O2_C5(i) + O2_C6(i)*th)*1.0e-4_fp*(pd + e)*th**0.8_fp
      ff = f/O2_C0(i)*((df - dl*(O2_C0(i) - f))/((O2_C0(i) - f)**2 + df*df) + &
                       (df - dl*(O2_C0(i) + f))/((O2_C0(i) + f)**2 + df*df))
      nsum = nsum + s*ff
    END DO
    d  = 5.6e-4_fp*(pd + e)*th**0.8_fp
    nd = f*pd*th**2*(6.14e-5_fp/(d*(1.0_fp + (f/d)**2)) + &
                     1.4e-12_fp*pd*th**1.5_fp/(1.0_fp + 1.9e-5_fp*f**1.5_fp))
    gam = 0.1820_fp*f*(nsum + nd)
  END FUNCTION gamma_o_lbl

  ! ITU-R P.676-12 Annex 1, water vapour lines, dB/km
  PURE FUNCTION gamma_w_lbl(f, pd, e, t) RESULT(gam)
    REAL(fp),    INTENT(IN) :: f
    COMPLEX(fp), INTENT(IN) :: pd, e, t
    COMPLEX(fp) :: gam
    COMPLEX(fp) :: th, s, df, ff, nsum
    INTEGER :: i
    th = 300.0_fp/t
    nsum = (0.0_fp, 0.0_fp)
    DO i = 1, N_H2O
      s  = H2O_C1(i)*1.0e-1_fp*e*th**3.5_fp*EXP(H2O_C2(i)*(1.0_fp - th))
      df = H2O_C3(i)*1.0e-4_fp*(pd*th**H2O_C4(i) + H2O_C5(i)*e*th**H2O_C6(i))
      df = 0.535_fp*df + SQRT(0.217_fp*df*df + 2.1316e-12_fp*H2O_C0(i)**2/th)
      ff = f/H2O_C0(i)*(df/((H2O_C0(i) - f)**2 + df*df) + df/((H2O_C0(i) + f)**2 + df*df))
      nsum = nsum + s*ff
    END DO
    gam = 0.1820_fp*f*nsum
  END FUNCTION gamma_w_lbl

  ! ITU-R P.676-10 Annex 2, dry air, dB/km (pd: dry pressure hPa)
  PURE FUNCTION gamma_o_approx(f, pd, t) RESULT(gam)
    REAL(fp),    INTENT(IN) :: f
    COMPLEX(fp), INTENT(IN) :: pd, t
    COMPLEX(fp) :: gam
    COMPLEX(fp) :: rp, rt, delta, xi1, xi2, xi3, xi4, xi5, xi6, xi7
    COMPLEX(fp) :: g54, g58, g60, g62, g64, g66
    rp = pd/1013.0_fp
    rt = 288.0_fp/t
    delta = -0.00306_fp*phi(3.211_fp, -14.94_fp, 1.583_fp, -16.37_fp)
    xi1 = phi( 0.0717_fp, -1.8132_fp,  0.0156_fp, -1.6515_fp)
    xi2 = phi( 0.5146_fp, -4.6368_fp, -0.1921_fp, -5.7416_fp)
    xi3 = phi( 0.3414_fp, -6.5851_fp,  0.2130_fp, -8.5854_fp)
    xi4 = phi(-0.0112_fp,  0.0092_fp, -0.1033_fp, -0.0009_fp)
    xi5 = phi( 0.2705_fp, -2.7192_fp, -0.3016_fp, -4.1033_fp)
    xi6 = phi( 0.2445_fp, -5.9191_fp,  0.0422_fp, -8.0719_fp)
    xi7 = phi(-0.1833_fp,  6.5589_fp, -0.2402_fp,  6.131_fp)
    g54 = 2.192_fp*phi(1.8286_fp, -1.9487_fp, 0.4051_fp, -2.8509_fp)
    g58 = 12.59_fp*phi(1.0045_fp,  3.5610_fp, 0.1588_fp,  1.2834_fp)
    g60 = 15.00_fp*phi(0.9003_fp,  4.1335_fp, 0.0427_fp,  1.6088_fp)
    g62 = 14.28_fp*phi(0.9886_fp,  3.4176_fp, 0.1827_fp,  1.3429_fp)
    g64 = 6.819_fp*phi(1.4320_fp,  0.6258_fp, 0.3177_fp, -0.5914_fp)
    g66 = 1.908_fp*phi(2.0717_fp, -4.1404_fp, 0.4910_fp, -4.8718_fp)
    IF (f <= 54.0_fp) THEN
      gam = ((7.2_fp*rt**2.8_fp)/(f**2 + 0.34_fp*rp**2*rt**1.6_fp) + &
             (0.62_fp*xi3)/((54.0_fp - f)**(1.16_fp*xi1) + 0.83_fp*xi2))*f**2*rp**2*1.0e-3_fp
    ELSE IF (f <= 60.0_fp) THEN
      gam = EXP(LOG(g54)/24.0_fp*(f - 58.0_fp)*(f - 60.0_fp) - &
                LOG(g58)/8.0_fp*(f - 54.0_fp)*(f - 60.0_fp) + &
                LOG(g60)/12.0_fp*(f - 54.0_fp)*(f - 58.0_fp))
    ELSE IF (f <= 62.0_fp) THEN
      gam = g60 + (g62 - g60)*(f - 60.0_fp)/2.0_fp
    ELSE IF (f <= 66.0_fp) THEN
      gam = EXP(LOG(g62)/8.0_fp*(f - 64.0_fp)*(f - 66.0_fp) - &
                LOG(g64)/4.0_fp*(f - 62.0_fp)*(f - 66.0_fp) + &
                LOG(g66)/8.0_fp*(f - 62.0_fp)*(f - 64.0_fp))
    ELSE IF (f <= 120.0_fp) THEN
      gam = (3.02e-4_fp*rt**3.5_fp + (0.283_fp*rt**3.8_fp)/((f - 118.75_fp)**2 + 2.91_fp*rp**2*rt**1.6_fp) + &
             (0.502_fp*xi6*(1.0_fp - 0.0163_fp*xi7*(f - 66.0_fp)))/((f - 66.0_fp)**(1.4346_fp*xi4) + 1.15_fp*xi5)) &
            *f**2*rp**2*1.0e-3_fp
    ELSE
      gam = ((3.02e-4_fp)/(1.0_fp + 1.9e-5_fp*f**1.5_fp) + &
             (0.283_fp*rt**0.3_fp)/((f - 118.75_fp)**2 + 2.91_fp*rp**2*rt**1.6_fp))*f**2*rp**2*rt**3.5_fp*1.0e-3_fp &
            + delta
    END IF
  CONTAINS
    PURE COMPLEX(fp) FUNCTION phi(a, b, c, d)
      REAL(fp), INTENT(IN) :: a, b, c, d
      phi = rp**a*rt**b*EXP(c*(1.0_fp - rp) + d*(1.0_fp - rt))
    END FUNCTION phi
  END FUNCTION gamma_o_approx

  ! ITU-R P.676-10 Annex 2, water vapour, dB/km (pd: dry pressure hPa, rho: g/m3)
  PURE FUNCTION gamma_w_approx(f, pd, rho, t) RESULT(gam)
    REAL(fp),    INTENT(IN) :: f
    COMPLEX(fp), INTENT(IN) :: pd, rho, t
    COMPLEX(fp) :: gam
    COMPLEX(fp) :: rp, rt, e1, e2
    rp = pd/1013.0_fp
    rt = 288.0_fp/t
    e1 = 0.955_fp*rp*rt**0.68_fp + 0.006_fp*rho
    e2 = 0.735_fp*rp*rt**0.50_fp + 0.0353_fp*rt**4*rho
    gam = ( (3.98_fp*e1*EXP(2.23_fp*(1.0_fp - rt)))/((f - 22.235_fp)**2 + 9.42_fp*e1**2)*g(22.0_fp) + &
            (11.96_fp*e1*EXP(0.70_fp*(1.0_fp - rt)))/((f - 183.310_fp)**2 + 11.14_fp*e1**2) + &
            (0.081_fp*e1*EXP(6.44_fp*(1.0_fp - rt)))/((f - 321.226_fp)**2 + 6.29_fp*e1**2) + &
            (3.660_fp*e1*EXP(1.60_fp*(1.0_fp - rt)))/((f - 325.153_fp)**2 + 9.22_fp*e1**2) + &
            (25.37_fp*e1*EXP(1.09_fp*(1.0_fp - rt)))/((f - 380.000_fp)**2) + &
            (17.40_fp*e1*EXP(1.46_fp*(1.0_fp - rt)))/((f - 448.000_fp)**2) + &
            (844.6_fp*e1*EXP(0.17_fp*(1.0_fp - rt)))/((f - 557.000_fp)**2)*g(557.0_fp) + &
            (290.0_fp*e1*EXP(0.41_fp*(1.0_fp - rt)))/((f - 752.000_fp)**2)*g(752.0_fp) + &
            (8.3328e4_fp*e2*EXP(0.99_fp*(1.0_fp - rt)))/((f - 1780.00_fp)**2)*g(1780.0_fp) ) &
          *f**2*rt**2.5_fp*rho*1.0e-4_fp
  CONTAINS
    PURE REAL(fp) FUNCTION g(fi)
      REAL(fp), INTENT(IN) :: fi
      g = 1.0_fp + ((f - fi)/(f + fi))**2
    END FUNCTION g
  END FUNCTION gamma_w_approx

  !
  ! gas_absorption_c: absorption coefficient (Np/m) in complex arithmetic
  !
  PURE FUNCTION gas_absorption_c(model, f, p, t, q) RESULT(alpha)
    INTEGER,     INTENT(IN) :: model
    REAL(fp),    INTENT(IN) :: f
    COMPLEX(fp), INTENT(IN) :: p, t, q
    COMPLEX(fp) :: alpha
    COMPLEX(fp) :: e, pd, rho
    CALL vapour_pressure(p, q, e, pd)
    SELECT CASE (model)
    CASE (GAS_MODEL_ITU_LBL)
      alpha = gamma_o_lbl(f, pd, e, t) + gamma_w_lbl(f, pd, e, t)
    CASE (GAS_MODEL_ITU_APPROX)
      rho = e*216.7_fp/t                 ! g/m3
      alpha = gamma_o_approx(f, pd, t) + gamma_w_approx(f, pd, rho, t)
    CASE DEFAULT
      alpha = (0.0_fp, 0.0_fp)
    END SELECT
    alpha = alpha*DB_TO_NP*1.0e-3_fp     ! dB/km -> Np/m
  END FUNCTION gas_absorption_c

  ! Real-valued absorption coefficient (Np/m)
  ELEMENTAL FUNCTION gas_absorption(model, f, p, t, q) RESULT(alpha)
    INTEGER,  INTENT(IN) :: model
    REAL(fp), INTENT(IN) :: f, p, t, q
    REAL(fp) :: alpha
    alpha = REAL(gas_absorption_c(model, f, CMPLX(p, 0.0_fp, fp), CMPLX(t, 0.0_fp, fp), &
                                  CMPLX(q, 0.0_fp, fp)))
  END FUNCTION gas_absorption

  ! Absorption coefficient and its derivatives with respect to t and q
  ELEMENTAL SUBROUTINE gas_absorption_deriv(model, f, p, t, q, alpha, dalpha_dt, dalpha_dq)
    INTEGER,  INTENT(IN)  :: model
    REAL(fp), INTENT(IN)  :: f, p, t, q
    REAL(fp), INTENT(OUT) :: alpha, dalpha_dt, dalpha_dq
    COMPLEX(fp) :: a
    REAL(fp) :: ht, hq
    ht = CSTEP*t
    hq = CSTEP
    a = gas_absorption_c(model, f, CMPLX(p, 0.0_fp, fp), CMPLX(t, ht, fp), CMPLX(q, 0.0_fp, fp))
    alpha = REAL(a)
    dalpha_dt = AIMAG(a)/ht
    a = gas_absorption_c(model, f, CMPLX(p, 0.0_fp, fp), CMPLX(t, 0.0_fp, fp), CMPLX(q, hq, fp))
    dalpha_dq = AIMAG(a)/hq
  END SUBROUTINE gas_absorption_deriv

END MODULE radar_gas
