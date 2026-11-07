"""Statistics and figures for the CRTM vs standalone radar comparison."""
import sys, json, numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fn, outdir = sys.argv[1], sys.argv[2]
d = nc.Dataset(fn)
cz, cza = d['crtm_z'][:].filled(-9999), d['crtm_zatt'][:].filled(-9999)
sz, sza = d['sa_z'][:].filled(-9999), d['sa_zatt'][:].filled(-9999)
ctg = d['crtm_tau_gas'][:]
stc = d['sa_tau_cum'][:]; sth = d['sa_tau_hydro'][:]
h = d['height'][:]
chans = ['DPR Ku 13.6 GHz', 'DPR Ka 35.5 GHz', 'CPR 94.05 GHz']
vars_ = ['A: CRTM-compatible, CRTM gas', 'B: CRTM-compatible, ITU LBL gas',
         'C: approximate gas + hydro', 'D: physical (full extinction)']
print('dims', cz.shape, sz.shape, ctg.shape, stc.shape)
stats = {}
for ic, cn in enumerate(chans):
    for v, vn in enumerate(vars_):
        r = {}
        for key, a, b in (('Z', sz[v, ic], cz[ic]), ('Zatt', sza[v, ic], cza[ic])):
            m = (a > -100) & (b > -100)
            df = (a - b)[m]
            only = int(((a > -100) != (b > -100)).sum())
            r[key] = dict(n=int(m.sum()), mean=float(df.mean()), std=float(df.std()),
                          maxabs=float(np.abs(df).max()), p99=float(np.percentile(np.abs(df), 99)),
                          mismatch_echo=only)
        stats[f'{ic}:{v}'] = r
        print(cn, vn, json.dumps(r))
json.dump(stats, open(f'{outdir}/stats.json', 'w'), indent=1)

# Figure 1: Zatt standalone vs CRTM per variant and channel (2D hist)
fig, ax = plt.subplots(3, 4, figsize=(16, 11.5), constrained_layout=True)
for ic in range(3):
    for v in range(4):
        a, b = sza[v, ic].ravel(), cza[ic].ravel()
        m = (a > -100) & (b > -100) & (b > -30)
        H = ax[ic, v].hexbin(b[m], a[m], gridsize=80, bins='log', cmap='viridis', mincnt=1,
                             extent=(-30, 60, -30, 60))
        ax[ic, v].plot([-30, 60], [-30, 60], 'k--', lw=.7)
        s = stats[f'{ic}:{v}']['Zatt']
        ax[ic, v].set_title(f'{chans[ic]}\n{vars_[v]}\nbias {s["mean"]:.3f} dB, sd {s["std"]:.3f} dB', fontsize=9)
        ax[ic, v].set_xlabel('CRTM attenuated Z (dBZ)'); ax[ic, v].set_ylabel('standalone (dBZ)')
fig.savefig(f'{outdir}/zatt_scatter.png', dpi=110)

# Figure 2: mean difference profile vs height (attenuated)
fig, ax = plt.subplots(1, 3, figsize=(15, 5.5), constrained_layout=True, sharey=True)
hb = np.arange(0, 16001, 500)
for ic in range(3):
    for v in range(3):
        a, b = sza[v, ic], cza[ic]
        m = (a > -100) & (b > -100)
        df = (a - b)[m]; hh = h[m] if h.shape == a.shape else np.broadcast_to(h, a.shape)[m]
        idx = np.digitize(hh, hb)
        mu = [df[idx == k].mean() if (idx == k).sum() > 20 else np.nan for k in range(1, len(hb))]
        ax[ic].plot(mu, (hb[:-1] + 250) / 1000, label=vars_[v])
    ax[ic].axvline(0, color='k', lw=.6); ax[ic].set_title(chans[ic])
    ax[ic].set_xlabel('mean standalone − CRTM attenuated Z (dB)')
ax[0].set_ylabel('height (km)'); ax[0].legend(fontsize=8)
fig.savefig(f'{outdir}/zatt_diff_profile.png', dpi=110)

# Figure 3: total gas optical depth (column) ODPS vs ITU LBL vs approx
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
for ic in range(3):
    cg = ctg[ic].sum(axis=-1)
    for v, lab in ((1, 'ITU-R P.676 line by line'), (2, 'ITU-R P.676 approx')):
        sg = stc[v, ic][:, -1] - sth[v, ic].sum(axis=-1)
        ax[ic].scatter(cg, sg, s=2, label=lab)
        print(chans[ic], lab, 'column gas OD ratio sa/crtm median %.3f' % np.median(sg / np.maximum(cg, 1e-12)))
    lim = [0, max(cg.max(), 1e-3) * 1.1]
    ax[ic].plot(lim, lim, 'k--', lw=.7); ax[ic].set_title(chans[ic])
    ax[ic].set_xlabel('CRTM ODPS column gas OD'); ax[ic].set_ylabel('standalone')
    ax[ic].legend(fontsize=8)
fig.savefig(f'{outdir}/gas_od.png', dpi=110)

# Figure 4: what CRTM leaves out (variant D): cloud liquid/ice backscatter, and
# full extinction vs delta-truncated extinction in the two-way PIA at the surface
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
for ic in range(3):
    a, b = sz[3, ic], cz[ic]
    m = (a > -100) & (b > -100)
    df = (a - b)[m]; hh = h[m]; idx = np.digitize(hh, hb)
    mu = [df[idx == k].mean() if (idx == k).sum() > 20 else np.nan for k in range(1, len(hb))]
    ax[0].plot(mu, (hb[:-1] + 250) / 1000, label=chans[ic])
    pb = (sz[1, ic] - sza[1, ic])[:, -1]; pdd = (sz[3, ic] - sza[3, ic])[:, -1]
    ok = (sz[1, ic][:, -1] > -100) & (sza[3, ic][:, -1] > -100)
    ax[1].scatter(pb[ok], pdd[ok], s=2, label=chans[ic])
ax[0].set_xlabel('mean Z with cloud backscatter − CRTM Z (dB)'); ax[0].set_ylabel('height (km)')
ax[0].set_title('Cloud liquid and ice backscatter (zero in CRTM)'); ax[0].legend(fontsize=8)
ax[1].plot([0, 200], [0, 200], 'k--', lw=.7); ax[1].set_xscale('symlog', linthresh=1); ax[1].set_yscale('symlog', linthresh=1)
ax[1].set_xlabel('surface PIA, delta-truncated extinction (dB)'); ax[1].set_ylabel('surface PIA, full extinction (dB)')
ax[1].set_title('Two-way PIA at the surface'); ax[1].legend(fontsize=8)
fig.savefig(f'{outdir}/crtm_omissions.png', dpi=110)
print('done')
