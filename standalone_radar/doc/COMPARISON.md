# Standalone radar operator vs CRTM

Data: `collocated_ifs-geer_atms_npp_20170907_18_00.nc`, all 60048 IFS
profiles, 136 layers, looking straight down (zenith 0°; a second run at
15° on every 4th profile gives the same statistics). Radars: GPM DPR Ku
(13.6 GHz) and Ka (35.55 GHz) with `dpr_gpm`, CloudSat CPR (94.05 GHz) with
`cpr_cloudsat`. Hydrometeors: cloud liquid, cloud ice, rain and snow from
the CloudCoeff DDA Thompson 2008 table. CRTM run with 4 streams, overcast
layers, no extra layers. Driver: `compare/compare_crtm.f90`; analysis:
`compare/analyze.py`.

## Variants

| | Gas | Hydrometeor attenuation | CRTM rules |
| --- | --- | --- | --- |
| A | CRTM ODPS optical depth passed in | full, delta-truncated | yes |
| B | ITU-R P.676 line by line | full, delta-truncated | yes |
| C | ITU-R P.676 approximation | approximate formulas | yes |
| D | ITU-R P.676 line by line | full extinction, no truncation | no (cloud backscatter on) |
| E | ITU-R P.676 line by line | NuRTM DDA single particles, temperature-dependent default PSDs (`species_default_ssdb`), delta-truncated | n/a |
| F | ITU-R P.676 line by line | as E, full extinction | n/a |

## Reflectivity (1.34 million gates with echo per channel)

Attenuated Z, standalone minus CRTM (dB):

| Channel | A bias / max | B bias / sd / max | C bias / sd / max | D bias / sd |
| --- | --- | --- | --- | --- |
| Ku 13.6 | 2e-7 / 4e-6 | +0.002 / 0.003 / 0.015 | −0.03 / 0.11 / 2.5 | +3.7 / 7.6 |
| Ka 35.5 | 2e-7 / 4e-6 | +0.008 / 0.012 / 0.07 | −0.16 / 0.32 / 6.2 | +3.7 / 7.6 |
| W 94.05 | 4e-7 / 8e-6 | −0.008 / 0.048 / 0.31 | −0.08 / 1.39 / 32 | +3.6 / 7.5 |

Unattenuated Z agrees to 4e-6 dB in A, B and C (the rest is CRTM's single
precision |K|²). Every gate with an echo in CRTM has one in A, B and C (except 64 of 1.34 million W-band gates in C, where the approximate attenuation pushes Z below the cutoff).

![attenuated Z scatter](fig/zatt_scatter.png)
![mean difference with height](fig/zatt_diff_profile.png)

## Gas absorption

Median ratio of column gas optical depth to CRTM ODPS: line by line 0.98
(Ku), 0.99 (Ka), 1.01 (W); approximation 1.02, 1.06, 1.10. The line by line
model and ODPS agree to about 2%, which moves attenuated Z by at most
0.3 dB at 94 GHz.

![gas optical depth](fig/gas_od.png)

## What CRTM leaves out (variant D)

* **Cloud liquid and cloud ice backscatter is zero in CRTM.** These tables
  have a single-scatter albedo of 0, and CRTM's quality rule sets
  backscatter to zero when w ≤ 0. Turning it on (D) adds echo in 980 000
  more gates. For gates with only cloud liquid, the table gives a median of
  −6 dBZ against −54 dBZ from the single-particle calculation (E), so the
  cloud-liquid entries are about 50 dB too high and should not be used for
  backscatter. The cloud-ice entries give −17 dBZ, which is plausible.
* **Delta truncation of the extinction** (CRTM default) removes the forward
  peak from the attenuation. Using the full extinction raises the two-way
  PIA at the surface by about 0.5% at Ku and Ka and 1.5% at W, and by up to
  20 dB in the heaviest W-band columns (183 → 204 dB).
* CRTM's LUT phase coefficients seem to follow `(2l+1)/2·χ_l`
  (`pc1 = 1.5 g` checked). If that is so, the truncation fraction
  `f = pc(n)` is overstated by `(2n+1)/2`, and CRTM attenuation is slightly
  too low. Worth checking in CRTM.

![CRTM omissions](fig/crtm_omissions.png)

## Approximate attenuation (variant C)

The ITU-R formulas track the full calculation at Ku (sd 0.1 dB) and are
usable at Ka (sd 0.3 dB). At 94 GHz the Rayleigh formulas for snow do not
hold, and attenuated Z differs by up to 32 dB in heavy snow and rain, so the
full calculation should be used at W band.

## Jacobians

CRTM K-matrix (seeding attenuated Z one gate at a time) against the
standalone K for the 5 profiles with the most rain and snow, 136 × 136 per
species and channel. They agree to 1e-14 of the largest value everywhere
except layers whose water content is below the lowest density in the table
(1e-6 kg m⁻³). There the table is clamped, so the forward model is flat in
water content. The standalone operator gives the derivative of that flat
function (consistent with finite differences), while CRTM keeps the slope of
the interpolating polynomial at the table edge. This affects a few layers
with very small water content (up to 16% of the rain Jacobian at W band in
those layers).

The standalone TL, AD and K also pass their own tests (TL against finite
differences to 1e-8, AD against TL to 1e-15, K against TL to 1e-16) in 24
configurations: `test/test_radar_operator.f90` (30 cases including the
particle paths).

## Single particles from NuRTM with temperature-dependent PSDs (variants E, F)

The species are built from the NuRTM single-scattering databases
(`ssdb_nurtm/mw`, converted from the ARTS DDA database) with
`species_from_ssdb`:

| Species | Habit | PSD |
| --- | --- | --- |
| cloud liquid | LiquidSphere | gamma, N0 = 3e35, μ = 4, on d_veq (NuRTM validation default) |
| cloud ice | IceSphere | Thompson 2008 ice: exponential, m = π/6·890 D³, Ni from Cooper (1986) capped at 250 L⁻¹, Thompson slope limits |
| rain | LiquidSphere | Thompson 2008 rain (N0 = 8e6), on d_veq |
| snow | 6BulletRosette | Field et al. (2007) mid-latitude at the layer temperature |
| graupel, hail (not in the IFS file) | GemGraupel, GemHail | Field et al. (2007) mid-latitude at the layer temperature |

All cross sections are interpolated to the layer temperature, so all six
species depend on temperature; the frozen PSDs do as well.

In this run CRTM and the table variants (A to D) also use the 6-bullet
rosette for snow (`compare_crtm ... 6BulletRosette`), so the habit is the same
everywhere. As in NuRTM, the cross sections are interpolated linearly in
frequency and temperature, the PSD is evaluated on the database sizes and
renormalised to the database particle masses, and the asymmetry parameter
comes from F11. TL, AD and K pass the same tests as the other paths.

Median Z (dBZ) on gates with only one species, all 60048 profiles:

| Species | Ku: CRTM / E | Ka: CRTM / E | W: CRTM / E |
| --- | --- | --- | --- |
| rain | 2.2 / 2.3 | 2.9 / 2.9 | 2.7 / 2.6 |
| snow (6-bullet rosette in both) | −8.9 / −1.0 | −8.7 / −1.2 | −8.5 / −3.1 |
| cloud liquid | (0) / −54 | (0) / −54 | (0) / −54 |
| cloud ice | (0) / −26 | (0) / −26 | (0) / −25 |

* **Rain agrees with the CRTM table to 0.1 dB** at all three frequencies.
* **Snow is 5 to 8 dB higher than CRTM even with the same habit** (mean
  +7.7, +7.2, +5.2 dB), so the difference is in the PSD. (With the Thompson
  Field05 snow PSD used earlier it was 4 to 7 dB.) The scripts that built the CRTM DDA table
  (pyWeather `dda/gen_dda_coef`) use Thompson 2008 only for liquid (cloud
  water and rain). All frozen habits use the Field et al. (2007)
  mid-latitude PSD at a fixed 235 K, with the habit mass-size relation.
  With that PSD (`PSD_F07_MIDLAT`, `t_fix = 235`, `RADAR_SNOW_PSD=f07m_235` in
  the driver) the standalone snow agrees with CRTM to within 0.35 to 0.47 dB
  (sd 0.1 dB) at all three frequencies. The remaining offset is consistent
  with the mass-size relation (database masses here, the ARTS α, β fit in
  the table) and the temperature of the cross sections.
* **Fixed 235 K vs layer temperature for the Field (2007) snow PSD**
  (default vs `RADAR_SNOW_PSD=f07m_235`, all 60048 profiles). Using the layer
  temperature raises snow Z by 7.3 dB at Ku, 6.8 dB at Ka and 4.9 dB at W on
  snow-only gates, and the change grows with temperature:

  | Layer T (K) | Ku | Ka | W |
  | --- | --- | --- | --- |
  | 235–245 | +1.3 | +1.3 | +1.1 |
  | 245–255 | +3.4 | +3.3 | +2.7 |
  | 255–265 | +5.5 | +5.3 | +4.3 |
  | 265–273 | +7.6 | +7.1 | +5.0 |

  Warmer snow has fewer, larger particles in Field (2007), and most snow in
  this file is between 265 and 273 K. Over all gates with snow, attenuated Z
  changes by +1.0 (Ku), +0.9 (Ka) and +0.3 dB (W). The two-way PIA at the
  surface changes little at Ku and Ka (at most 2.4 dB), but at W it goes up
  by 0.3 dB on average and by up to 34 dB in heavy snow. So the fixed 235 K
  in the CRTM tables makes snow several dB too weak at the temperatures
  where most snow is. The layer-temperature PSD is therefore the default.
* **Cloud liquid** gives about −54 dBZ, realistic for these small water
  contents.
* **Cloud ice** with the Thompson PSD gives a median of −26 dBZ (90th
  percentile −8, 99th +7). The high tail is warm cloud ice: Cooper gives few
  crystals near 0 °C, so the IFS ice mass goes into large particles even
  with Thompson's 700 µm slope limit. A model with prognostic Ni can pass it
  (`ni` argument) instead.
* Surface two-way PIA, E vs B: means are nearly the same (W 5.49 vs 5.47 dB);
  the 1st to 99th percentile of the difference is −0.1 to +0.4 dB at Ka and
  −1.7 to +6.1 dB at W, and the W-band maximum rises from 183 to 200 dB.
* Full extinction (F) instead of truncation changes the surface PIA by at
  most 0.2 dB at Ka and 14.7 dB at W (mean 0.06 dB).

See [TECHNICAL_NOTE.md](TECHNICAL_NOTE.md) for the full write-up.
