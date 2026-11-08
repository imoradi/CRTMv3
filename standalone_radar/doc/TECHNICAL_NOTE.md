# Standalone radar operator: technical note

This note records what the standalone radar operator does, how it was
built and tested, and what the comparison with CRTM showed. It covers the
work on branch `feature/standalone-radar` up to November 2026. The detailed
comparison numbers are in [COMPARISON.md](COMPARISON.md), and the API is in
the [README](../README.md).

## 1. Purpose and scope

The goal was to take the CRTM radar operator out of CRTM so that it can:

* reproduce the CRTM active-sensor result for the same inputs;
* run for spaceborne and ground-based radars;
* give reflectivity with and without attenuation, in mm⁶ m⁻³ and dBZ;
* give zenith-independent (vertical) layer and cumulative optical depths,
  so the attenuation can be applied later for any angle, plus the two-way
  transmittance and path-integrated attenuation (PIA) for the given angle;
* use either a bulk lookup table or single particles with a particle size
  distribution (PSD), or bulk optics given directly;
* offer full and approximate attenuation;
* provide tangent linear (TL), adjoint (AD) and K-matrix (Jacobian) code;
* offer a short interface that a model such as GEOS can call.

The UFO/JEDI interface was dropped at the user's request. The operator has
no CRTM dependency; CRTM is used only by the comparison driver.

## 2. Code layout

| File | Contents |
| --- | --- |
| `src/radar_kinds.f90` | kinds and constants |
| `src/radar_interp.f90` | CRTM lookup-table interpolation (4-point blended Lagrange) |
| `src/radar_gas.f90` | gas absorption: ITU-R P.676-12 Annex 1 (line by line), Annex 2 (approximation), or user optical depth |
| `src/radar_hydro_approx.f90` | approximate hydrometeor attenuation |
| `src/radar_species.f90` | hydrometeor species: CloudCoeff LUT, particles, NuRTM SSDB, PSDs |
| `src/radar_tdts.F90` | TDTS wide-angle multiple scattering (Hogan and Battaglia 2008, as in CRTM) |
| `src/radar_operator.f90` | operator: forward, TL, AD, K |
| `src/radar_simple.f90` | short interface for models (GEOS) |
| `test/test_radar_operator.f90` | derivative tests |
| `compare/compare_crtm.f90`, `compare/analyze.py` | comparison with CRTM on the collocated IFS file |

Build with CMake (`-DRADAR_BUILD_TESTS=ON`); only netCDF-Fortran is needed.

## 3. Physics

### 3.1 Reflectivity

For each layer (range gate), the backscatter coefficient β (m⁻¹, summed over
species and any bulk input) gives

  Z = C·β,  C = 10¹⁸ λ⁴ / (π⁵ |K_w|²)  (mm⁶ m⁻³, λ in m),

with |K_w|² from Turner et al. (2016) at 273.15 K, as CRTM (or a fixed value
with `KW_FIXED`). dBZ = 10 log₁₀ Z, and −9999 marks Z = 0, as in CRTM.

### 3.2 Attenuation

Each layer has a vertical one-way optical depth τ = (α_gas + α_hydro)·Δz.
The cumulative optical depth to gate k counts the layers between the radar
and the gate plus a fraction `gate_fraction` of the gate's own layer
(1 = whole layer, as CRTM; 0.5 = to the layer centre). For a spaceborne radar
the sum runs from the top, for a ground radar from the bottom. Then

  trans2 = exp(−2 τ_cum / cos θ),  Z_att = Z·trans2,  PIA = 10 log₁₀(1/trans2).

Because `tau` and `tau_cum` are vertical, the attenuation for any other
zenith angle is `z_lin·exp(−2·tau_cum/cos θ)` without rerunning.

**Gas.** `GAS_MODEL_ITU_LBL` is ITU-R P.676-12 Annex 1 (O₂ and H₂O lines
plus continuum). `GAS_MODEL_ITU_APPROX` is the Annex 2 approximation.
`GAS_MODEL_USER` takes an optical depth from the caller; the comparison
driver passes the CRTM ODPS value this way.

**Hydrometeors, full.** The extinction from the scattering properties. Three
treatments of the forward-scattered part:

* `SCALING_DELTA` (CRTM default): τ′ = τ − f·σ_s with
  f = the phase coefficient at order `n_streams`. For particles, f comes
  from a Henyey-Greenstein expansion, f·σ_s = (2n+1)/2·gⁿ·σ_s, matching the
  CRTM LUT convention (Section 7.3).
* `SCALING_NONE`: full extinction, the usual single-scattering radar
  convention.
* `SCALING_EFFECTIVE`: the CRTM effective-attenuation option with η.

**Hydrometeors, approximate** (`ATTEN_APPROX`), formulas used in the radar
community:

| Species | Formula |
| --- | --- |
| cloud liquid | ITU-R P.840 (Rayleigh, double-Debye water permittivity) |
| cloud ice, snow, graupel, hail | Rayleigh absorption with the Mätzler ice permittivity |
| rain | ITU-R P.838-3, k·R^α, with rain rate from M = 0.072 R^0.88 |

### 3.3 Multiple scattering

The Hogan and Battaglia (2008) time-dependent two-stream (TDTS) wide-angle
term, ported from CRTM, is available for spaceborne radars
(`multiple_scattering = .true.`, with `altitude`, `beamwidth`, `cbe`). It
adds `z_ms` to the reflectivity.

### 3.4 CRTM-compatible rules

With `crtm_compat = .true.` (default for the LUT path) the operator follows
CRTM's quality rules: no backscatter when the single-scatter albedo is zero,
layers with water path below `wc_threshold` (1e-6 kg m⁻²) skipped, and the
CRTM dBZ TL/AD limits (`tlad_min_dbz`, `tlad_cap_dbz`).

## 4. Hydrometeor optics

Each species gives, per layer, extinction, scattering and backscatter
coefficients, the asymmetry parameter and the truncated fraction. There are
four ways to define them.

### 4.1 CRTM CloudCoeff table (`species_from_lut`)

Reads a CRTM CloudCoeff netCDF file and uses CRTM's interpolation in
frequency, density and (liquid) temperature. Names follow CRTM
(`cloud_liquid`, `rain`, `cloud_ice`, `snow`, `graupel`, `hail`), or any
habit in the file through `habit` (for example `SixBulletRosette`).

### 4.2 Single particles with a gamma PSD (`species_from_particles`)

Single-particle cross sections on a size grid plus a gamma PSD
N(D) = N₀ D^μ e^(−λD). λ is solved from the water content in each layer and
its derivative is carried through TL and AD.

### 4.3 NuRTM single-scattering databases (`species_from_ssdb`)

Reads a NuRTM SSDB netCDF file (one habit per file, converted from the ARTS
DDA database, Zenodo 4646605): frequency, temperature, d_veq, d_max, mass,
and ext, abs, bsc and the phase matrix on (frequency, temperature, size).
Following NuRTM:

* cross sections are interpolated linearly to the radar frequency and kept
  on the database temperature grid; in each layer they are interpolated
  linearly in temperature (clamped to the database range);
* backscatter is `bsc`, or 4π·F₁₁(180°) if `bsc` is missing;
* the asymmetry parameter is integrated from F₁₁;
* the PSD is evaluated on the database sizes with trapezoid bin widths and
  renormalised so that Σ N·m equals the water content with the database
  particle masses.

Databases used (`ssdb_nurtm/mw` on the Mac):

| File | Temperatures (K) |
| --- | --- |
| LiquidSphere | 230 to 310 |
| IceSphere, 6BulletRosette, GemGraupel, GemHail | 190, 230, 270 |

### 4.4 Bulk optics per layer

`prof%bulk_ext`, `bulk_sca`, `bulk_bsc`, `bulk_g` with `has_bulk = .true.`
are added to the species.

### 4.5 PSDs

| Constant | Definition |
| --- | --- |
| `PSD_GAMMA_FIXED_N0` | gamma, fixed N₀ and μ, λ from the water content |
| `PSD_GAMMA` | gamma with a log-log mass-size fit to the database |
| `PSD_T08_RAIN` | Thompson et al. (2008) rain, N₀ = 8e6 m⁻⁴ |
| `PSD_T08_SNOW` | exponential, Thompson (2004) intercept N₀ = min(2e8, 2e6·e^(−0.12·min(−0.001, T−273.15))) |
| `PSD_T08_SNOW_FIELD05` | Thompson 2008 snow: two-mode Field et al. (2005), moments from temperature |
| `PSD_T08_GRAUPEL` | N₀ = max(1e4, min(200/q, 5e6)) |
| `PSD_T08_ICE` | Thompson cloud ice: exponential, m = π/6·890·D³, N_i given or Cooper (1986) 5·e^(0.304(273.15−T)) m⁻³ capped at 250 L⁻¹, slope limited to [1/700 µm, 1e6 m⁻¹] |
| `PSD_F07_MIDLAT`, `PSD_F07_TROPICAL` | Field et al. (2007) moment-scaled PSD, Tc = T−273 clamped to [−70, −0.1] °C; layer temperature, or `t_fix` |

All PSDs depend on the layer temperature except `PSD_GAMMA*` and
`PSD_T08_RAIN`, and the temperature derivative of the PSD is carried
through TL and AD.

### 4.6 Default species (`species_default_ssdb`, `radar_simple_init_ssdb`)

The default set for the particle path uses temperature-dependent cross
sections for all six species, and temperature-dependent PSDs for all frozen
species:

| Species | Database | PSD | T dependence |
| --- | --- | --- | --- |
| cloud liquid | LiquidSphere (d_veq) | gamma, N₀ = 3e35, μ = 4 (NuRTM default) | cross sections |
| rain | LiquidSphere (d_veq) | Thompson 2008 rain | cross sections |
| cloud ice | IceSphere | Thompson 2008 ice, Cooper N_i | cross sections and PSD |
| snow | 6BulletRosette | Field 2007 mid-latitude, layer T | cross sections and PSD |
| graupel | GemGraupel | Field 2007 mid-latitude, layer T | cross sections and PSD |
| hail | GemHail | Field 2007 mid-latitude, layer T | cross sections and PSD |

The liquid PSDs (cloud liquid and rain) have no temperature dependence in
their definitions, so for them the temperature enters through the
permittivity in the cross sections. Snow, graupel and hail use the same
Field (2007) PSD as the CRTM DDA tables, but at the layer temperature
instead of the fixed 235 K (Section 6.4).

```fortran
call radar_simple_init_ssdb(h, 35.5_fp, '/path/ssdb_nurtm/mw', &
     ['cloud_liquid', 'cloud_ice   ', 'rain        ', 'snow        ', 'graupel     ', 'hail        '], &
     platform=PLATFORM_SPACEBORNE)
```

## 5. Derivatives and tests

TL, AD and K are written by hand for the operator. Inside the species and
gas code, derivatives with respect to water content and temperature use the
complex-step method (step 1e-30), so all PSDs, interpolation and the gas
model share one code path for value and derivative. `min`/`max` on complex
numbers compare the real parts.

`test/test_radar_operator.f90` checks, for every configuration at 13.6, 35.5
and 94 GHz:

* TL against finite differences (relative error about 1e-9);
* AD against TL (dot-product test, about 1e-15);
* K against TL (about 1e-15).

There are 30 cases: the LUT path in 8 configurations (spaceborne and ground,
full/approximate attenuation, delta/none/effective scaling, TDTS, user gas,
bulk input), single particles, and two SSDB cases (the default
temperature-dependent set, and a ground radar with Thompson Field05 snow,
Thompson graupel and a gamma ice PSD). All pass.

## 6. Comparison with CRTM

Data: `collocated_ifs-geer_atms_npp_20170907_18_00.nc`, all 60048 IFS
profiles with 136 layers, nadir. Radars: GPM DPR Ku (13.6 GHz) and Ka
(35.55 GHz) with `dpr_gpm`, CloudSat CPR (94.05 GHz) with `cpr_cloudsat`
(coefficients from fix_REL-3.1.2.0). CRTM uses the CloudCoeff DDA
"Thompson2008PSD" table with the 6-bullet rosette for snow, 4 streams,
overcast layers.

| Variant | Gas | Hydrometeors |
| --- | --- | --- |
| A | CRTM ODPS optical depth | CRTM table, delta truncation, CRTM rules |
| B | ITU-R line by line | as A |
| C | ITU-R approximation | approximate attenuation formulas |
| D | ITU-R line by line | table, full extinction, cloud backscatter on |
| E | ITU-R line by line | NuRTM particles, default T-dependent PSDs (4.6), delta truncation |
| F | ITU-R line by line | as E, full extinction |

### 6.1 Reproducing CRTM (A, B, C)

* **A matches CRTM to 4e-6 dB** in attenuated Z on 1.34 million gates per
  channel. The remainder is CRTM's single-precision |K_w|².
* **B** (own gas model): bias below 0.01 dB, at most 0.3 dB at W band.
  Column gas optical depth is within about 2% of ODPS.
* **C** (approximate attenuation): sd 0.1 dB (Ku), 0.3 dB (Ka), 1.4 dB (W);
  up to 32 dB at W band in heavy snow and rain, where the Rayleigh formulas
  fail. Use full attenuation at W band.

### 6.2 Jacobians

The CRTM K-matrix and the standalone K agree to 1e-14 of the largest value
for the 5 profiles with the most rain and snow (136 × 136 per species and
channel), except in layers whose water content is below the lowest table
density (1e-6 kg m⁻³). There the table is clamped; the standalone K is the
derivative of the clamped function (consistent with finite differences),
while CRTM keeps the slope of the interpolating polynomial.

### 6.3 Single particles: rain, cloud liquid, cloud ice

* Rain from LiquidSphere with the Thompson PSD agrees with the CRTM table
  to 0.1 dB at all three frequencies.
* Cloud liquid from particles gives about −54 dBZ, where the CRTM table
  gives about −6 dBZ (Section 7.2).
* Cloud ice: the NuRTM gamma PSD gave about 0 dBZ, too high. With the
  Thompson ice PSD the median is −26 dBZ (90th percentile −8, 99th +7). The
  high tail is warm ice, where Cooper gives few crystals. A model with
  prognostic N_i can pass it with `ni`.

### 6.4 Snow: habit, PSD and temperature

1. With the same habit (6-bullet rosette) in CRTM and the standalone
   operator, Thompson Field05 snow was 4 to 7 dB above the CRTM table, so
   the difference is in the PSD, not the habit.
2. The scripts that built the CRTM DDA table (pyWeather
   `dda/gen_dda_coef`: `arts_ssd.py`, `crtm_calc.py`, `psd_field.py`,
   `psd_thompson.py`) show that the "Thompson2008PSD" table uses Thompson
   only for cloud water and rain. All frozen habits (ice, snow, graupel,
   hail) use the Field et al. (2007) mid-latitude PSD at a fixed 235 K, with
   the ARTS α, β mass-size relation.
3. With that PSD at 235 K (`PSD_F07_MIDLAT`, `t_fix = 235`) the standalone
   snow agrees with CRTM to 0.35 to 0.47 dB (sd 0.1 dB). The remaining offset
   is from the mass-size relation (database masses vs the α, β fit) and the
   cross-section temperature.
4. Using the layer temperature instead of 235 K raises snow Z on snow-only
   gates by 7.3 dB (Ku), 6.8 dB (Ka) and 4.9 dB (W). The change grows with
   temperature, because warmer Field (2007) snow has fewer, larger
   particles:

   | Layer T (K) | Ku | Ka | W |
   | --- | --- | --- | --- |
   | 235–245 | +1.3 | +1.3 | +1.1 |
   | 245–255 | +3.4 | +3.3 | +2.7 |
   | 255–265 | +5.5 | +5.3 | +4.3 |
   | 265–273 | +7.6 | +7.1 | +5.0 |

   Over all gates with snow, attenuated Z rises by 1.0 (Ku), 0.9 (Ka) and
   0.3 dB (W). The W-band surface PIA rises by 0.3 dB on average and by up
   to 34 dB in heavy snow. Most snow in this file is between 265 and 273 K,
   so the fixed 235 K in the CRTM tables makes snow several dB too weak.

The layer-temperature Field (2007) PSD is now the default for snow, graupel
and hail (4.6).

### 6.5 Results with the temperature-dependent defaults (E, F)

Full run on all 60048 profiles with the defaults of 4.6 (the IFS file has
no graupel or hail, so those two are covered only by the derivative tests).

Median Z (dBZ) on gates with only one species:

| Species | gates | Ku: CRTM / E | Ka: CRTM / E | W: CRTM / E |
| --- | --- | --- | --- | --- |
| rain | 34284 | 2.2 / 2.3 | 2.9 / 2.9 | 2.7 / 2.6 |
| snow | 6078 | −8.9 / −1.0 | −8.7 / −1.2 | −8.5 / −3.1 |
| cloud liquid | 4551 | (0) / −54 | (0) / −54 | (0) / −54 |
| cloud ice | 254588 | (0) / −26 | (0) / −26 | (0) / −25 |

* Rain is within 0.1 dB of the CRTM table (mean difference +0.08, +0.01,
  −0.04 dB).
* Snow is higher than the CRTM table by 7.7 (Ku), 7.2 (Ka) and 5.2 dB (W) on
  average, almost all from the layer temperature in the Field (2007) PSD
  (6.4); at 235 K the two agree to within 0.5 dB.
* Cloud liquid and cloud ice are unchanged from 6.3.

Over all 1.34 million gates with echo per channel, attenuated Z from E minus
CRTM has a mean of +1.3 (Ku), +1.1 (Ka) and +0.4 dB (W), sd 6 to 7.5 dB. Most
of the spread is the cloud echo that CRTM leaves out (980 000 more gates
with echo in E).

Two-way PIA at the surface, E against B (same gas, CRTM table):

| Channel | B mean / max (dB) | E mean / max (dB) | E − B, 1st / 99th percentile |
| --- | --- | --- | --- |
| Ku | 0.25 / 14.5 | 0.25 / 14.7 | 0.0 / 0.0 |
| Ka | 1.36 / 82.4 | 1.37 / 83.5 | −0.1 / +0.4 |
| W | 5.47 / 183 | 5.49 / 200 | −1.7 / +6.1 |

Full extinction (F) instead of delta truncation (E) changes the surface PIA
by at most 0.2 dB at Ka and 14.7 dB at W (mean 0.06 dB at W); at Ku the
difference is negligible.

## 7. Findings about CRTM

1. **Cloud fraction.** `CRTM_Forward` zeroes cloud water content where
   `Cloud_Fraction < 1e-6`. Drivers must set the cloud fraction (we use 1),
   or every Z is missing.
2. **Cloud backscatter.** The CRTM water-cloud and ice-cloud tables have
   single-scatter albedo 0, so CRTM's quality rule sets their backscatter to
   zero. Turning it on (variant D) adds echo in 980 000 more gates and +3.7 dB
   on average. But the cloud-liquid table backscatter is about 50 dB too
   high (−6 dBZ vs −54 dBZ from particles) and should not be used. The
   cloud-ice entries (−17 dBZ) are plausible.
3. **Phase coefficients.** The CRTM table phase coefficients appear to
   follow (2l+1)/2·χ_l (checked: pc₁ = 1.5 g). If so, the delta-truncation
   fraction f = pc(n) is overstated by (2n+1)/2 and CRTM attenuation is a
   little too low. Worth checking in CRTM.
4. **Delta truncation.** Removing the forward peak lowers the attenuation:
   full extinction raises the surface two-way PIA by about 0.5% at Ku and Ka
   and 1.5% at W, up to 20 dB in the heaviest W-band columns.
5. **Fixed 235 K frozen PSD** in the DDA tables (6.4).
6. **Small water contents.** Below the table minimum the CRTM Jacobian keeps
   the polynomial slope while the forward model is flat (6.2).

## 8. Recommendations

* For GEOS: `radar_simple_init_ssdb` with the default species, full
  attenuation, `SCALING_NONE`, and the ITU line-by-line gas model. Use
  `tau_cum` to apply the attenuation for any zenith angle.
* To reproduce CRTM exactly: `radar_simple_init` with the CloudCoeff file
  and the default configuration.
* Approximate attenuation is fine at Ku, acceptable at Ka, not at W.
* Cloud ice is the most uncertain species; pass N_i when the model has it.

## 9. Open items

* Report the CRTM phase-coefficient convention and the cloud-liquid
  backscatter values to the CRTM team.
* Rebuild the CRTM DDA tables with a temperature-dependent frozen PSD.
* More habits from the DDA database (Zenodo 4646605) can be converted to the
  NuRTM SSDB format and used directly.
