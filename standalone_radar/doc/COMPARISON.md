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
  more gates and raises Z by several dB, up to about 18 dBZ near cloud top.
  That is too high for small cloud particles, so I would not use the cloud
  entries of the table for backscatter. Cloud backscatter needs proper
  single-particle data (for example the DDA database) through
  `species_from_particles`.
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
configurations: `test/test_radar_operator.f90`.
