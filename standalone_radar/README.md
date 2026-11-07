# Standalone radar operator

A radar reflectivity operator (forward, tangent linear, adjoint and K-matrix)
that runs without CRTM. With its default settings it reproduces the CRTM
active-sensor result (`CRTM_Active_Sensor`) for the same inputs, and it adds:

* spaceborne (looking down) and ground-based (looking up) radars;
* reflectivity with and without attenuation, in linear units (mm^6 m^-3) and dBZ;
* vertical layer and cumulative optical depths that do not depend on the
  zenith angle, so the attenuation can be applied later for any angle,
  plus the two-way transmittance and path-integrated attenuation (PIA) for
  the angle given;
* full or approximate gas and hydrometeor attenuation;
* three ways to give the hydrometeor optics: the CRTM bulk lookup table
  (CloudCoeff netCDF), single-particle properties with a gamma particle
  size distribution, or bulk optical properties per layer;
* optional wide-angle multiple scattering (Hogan and Battaglia 2008 TDTS, as
  in CRTM) for spaceborne radars.

## Building

```
mkdir build && cd build
cmake .. -DRADAR_BUILD_TESTS=ON      # needs netCDF-Fortran
make
./test_radar_operator /path/to/CloudCoeff_DDA_*.nc4
```

The test checks the TL against finite differences, the AD against the TL
(dot product test) and the K-matrix against the TL for 8 configurations at
13.6, 35.5 and 94 GHz.

## Source files

| File | Contents |
| --- | --- |
| `src/radar_kinds.f90` | kinds and constants |
| `src/radar_interp.f90` | the CRTM lookup-table interpolation (4-point blended Lagrange) |
| `src/radar_gas.f90` | gas absorption: ITU-R P.676-12 Annex 1 line by line, P.676 Annex 2 approximation, or a user optical depth |
| `src/radar_hydro_approx.f90` | approximate hydrometeor attenuation (ITU-R P.840 cloud, Rayleigh ice/snow, ITU-R P.838-3 rain) |
| `src/radar_species.f90` | hydrometeor species: CloudCoeff LUT, or particles plus a PSD |
| `src/radar_tdts.F90` | TDTS wide-angle multiple scattering |
| `src/radar_operator.f90` | the operator: forward, TL, AD, K |
| `src/radar_simple.f90` | a short interface for models such as GEOS |
| `test/test_radar_operator.f90` | derivative tests |
| `compare/compare_crtm.f90` | comparison with CRTM on a collocated IFS file |

## Simple interface (GEOS)

```fortran
use radar_kinds,  only: fp
use radar_simple
use radar_operator, only: PLATFORM_SPACEBORNE, PLATFORM_GROUND

integer :: h
call radar_simple_init(h, 13.6_fp, 'CloudCoeff_DDA_Thompson.nc4', &
     ['cloud_liquid', 'cloud_ice   ', 'rain        ', 'snow        ', 'graupel     '], &
     platform=PLATFORM_SPACEBORNE)

! per column; layer 1 is the top
wc(:,s) = radar_mixing_ratio_to_content(qx(:,s), p, t, q)     ! kg m-3
call radar_simple_run(h, zenith, p, t, q, z_int, wc, z_dbz, zatt_dbz, &
     z_lin=z_lin, zatt_lin=zatt_lin, tau=tau, tau_cum=tau_cum, pia=pia)

call radar_simple_destroy(h)
```

Inputs per layer: pressure `p` (hPa), temperature `t` (K), specific humidity
`q` (kg/kg), water content `wc(n_layers, n_species)` (kg m-3), and interface
heights `z_int(0:n_layers)` (m, top first).

Outputs per layer (range gate):

| Output | Meaning |
| --- | --- |
| `z_dbz`, `z_lin` | reflectivity without attenuation |
| `zatt_dbz`, `zatt_lin` | reflectivity attenuated along the beam at `zenith` |
| `tau` | vertical one-way optical depth of the layer (gas plus hydrometeors) |
| `tau_cum` | vertical one-way optical depth from the radar to the gate |
| `pia` | two-way path-integrated attenuation at `zenith` (dB) |

To attenuate for another zenith angle θ without rerunning:
`zatt_lin = z_lin * exp(-2 * tau_cum / cos θ)`. For a ground radar give
`zenith = 90 - elevation`. Missing dBZ (Z = 0) is -9999, as in CRTM.

`radar_simple_operator(h)` returns a pointer to the full operator, for
the TL, AD and K-matrix routines.

## Full interface

```fortran
type(radar_config_type)   :: cfg
type(radar_operator_type) :: op
type(radar_species_type)  :: sp
type(radar_profile_type)  :: prof, prof_tl, prof_ad
type(radar_output_type)   :: out, out_tl, out_ad
type(radar_k_type)        :: k

cfg%frequency = 94.05_fp
cfg%platform  = PLATFORM_SPACEBORNE
call radar_init(op, cfg, n_species=2)
call species_from_lut(sp, 'rain', 'CloudCoeff.nc4', cfg%frequency, cfg%n_streams)
call radar_add_species(op, 1, sp)
call species_from_particles(sp, 'snow', PHASE_FROZEN, cfg%frequency, d, mass, &
                            cext, csca, cbsc, asym, n0, mu)
call radar_add_species(op, 2, sp)

call radar_profile_create(prof, n, 2)     ! fill prof%p, t, q, z_int, rho, zenith
call radar_output_create(out, n)
call radar_forward(op, prof, out)
call radar_tl(op, prof, prof_tl, out_tl)
call radar_ad(op, prof, out_ad, prof_ad)   ! seeds in out_ad, accumulates in prof_ad
call radar_k(op, prof, OUT_ZATT_DBZ, k)    ! k%t(gate,layer), k%q, k%rho(gate,layer,species)
```

### Species

`species_from_lut(sp, name, file, frequency, n_streams [, habit])` reads a
CRTM CloudCoeff netCDF file. Names follow CRTM: `cloud_liquid`, `rain`
(liquid table, interpolated in frequency, density and temperature),
`cloud_ice`, `snow`, `graupel`, `hail` (the habit CRTM uses for each), or
any habit name in the file via `habit`.

`species_from_particles` takes single-particle properties on a size grid:
size `d` (m), mass (kg), extinction, scattering and backscattering cross
sections (m^2, backscatter = 4π dσ/dΩ at 180°) and asymmetry parameter, and
a gamma PSD `N(D) = n0 D^mu exp(-λ D)`. The slope λ is solved from the
water content for each layer, and its derivative is carried through the TL
and AD. Single-particle data from the DDA database (Zenodo record 4646605)
can be fed in here.

Bulk optical properties can be given directly per layer in
`prof%bulk_ext`, `bulk_sca`, `bulk_bsc` (m^-1) and `bulk_g` with
`prof%has_bulk = .true.`. They are added to the species.

### Options (`radar_config_type`)

| Field | Default | Choices |
| --- | --- | --- |
| `frequency` | 13.6 | GHz |
| `platform` | `PLATFORM_SPACEBORNE` | `PLATFORM_GROUND` |
| `gas_model` | `GAS_MODEL_ITU_LBL` | `GAS_MODEL_NONE`, `GAS_MODEL_ITU_APPROX`, `GAS_MODEL_USER` (`prof%tau_gas`) |
| `hydro_attenuation` | `ATTEN_FULL` | `ATTEN_NONE`, `ATTEN_APPROX` |
| `scaling` | `SCALING_DELTA` | `SCALING_NONE` (full extinction), `SCALING_EFFECTIVE` (CRTM effective attenuation, with `eta`) |
| `n_streams` | 4 | delta truncation order, as CRTM `n_Streams` |
| `multiple_scattering` | false | TDTS, spaceborne only (`altitude`, `beamwidth`, `cbe`) |
| `kw_model`, `kw2` | Turner 273 K | `KW_FIXED` with `kw2` |
| `gate_fraction` | 1.0 | part of the gate's own layer in its attenuation (CRTM uses the whole layer; 0.5 = layer centre) |
| `crtm_compat` | true | CRTM quality rules (no backscatter when the single-scatter albedo is zero) |
| `wc_threshold` | 1e-6 | water path (kg m^-2) below which a layer is skipped, as CRTM |
| `tlad_min_dbz`, `tlad_cap_dbz` | 0, 0 | dBZ TL/AD limits, as CRTM |

### Attenuation choices

Full attenuation uses the extinction from the scattering properties. With
`SCALING_DELTA` (CRTM default) the scattered part in the forward peak is
removed (`τ' = τ − f·σ_s`), which is what CRTM does; `SCALING_NONE` uses the
full extinction, which is the usual single-scattering radar convention.

Approximate attenuation follows the formulas used by the radar community:
ITU-R P.676 for gases, ITU-R P.840 (Rayleigh, double-Debye water) for cloud
liquid, Rayleigh absorption with the Mätzler ice permittivity for ice,
snow and graupel, and ITU-R P.838-3 `k R^α` for rain (rain rate from water
content with `M = 0.072 R^0.88`).

See [doc/COMPARISON.md](doc/COMPARISON.md) for the comparison with CRTM.
