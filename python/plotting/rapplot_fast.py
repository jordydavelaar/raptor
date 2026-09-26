"""Fast drop-in replacement for rapplot.

The original rapplot calls pcolormesh (and h5py row reads) once per AMR block,
which is very slow for 1e4-1e5 blocks.  Here the blocks are instead
rasterised once onto a single uniform grid at the finest AMR pitch, and every
panel is drawn with one imshow (or one contour) call.  Coarser blocks are
replicated onto the fine grid, finest blocks are written last so they win in
any overlap.  Cells not covered by any block are NaN (transparent).

Same function names and signatures as rapplot, so
    import rapplot_fast as rapplot
is all that is needed.  Differences: contours are computed on the merged grid
(continuous across block boundaries, instead of per block), and overlay
contours are drawn with a single call.

For deeply refined images (e.g. SMR on the critical curve), a grid at the
finest pitch over the whole image does not fit in memory. All plot functions
therefore also accept
    window=(x0, x1, y0, y1)   rasterise only this region (rg, plotted
                              coordinates x = alpha, y = -beta)
    pitch=p                   grid pitch (rg), rounded up to dmin * 2^k;
                              blocks finer than that are area-averaged
Without them the grid is exactly as before: all blocks overlapping
+-halfrange, at the finest pitch present.

Assumes the layout produced by RAPTOR: each block is n x n cell-centred pixels
flattened to n*n, alpha varying along the first block index and beta along
the second, with power-of-two pitch ratios between blocks.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py

from rapplot import (G, MSUN, SPEED_OF_LIGHT, KPC, SEC_IN_DAY, MAS_IN_DEG,
                     read_data_id)

_cache = {}   # id(image) -> {'fields': {key: (N,n*n) array}, 'grids': {halfrange: _Grid}}


def _state(image):
    st = _cache.get(id(image))
    if st is None or st['image'] is not image:
        st = {'image': image, 'fields': {}, 'grids': {}}
        _cache[id(image)] = st
    return st


def _field(image, key):
    """Whole dataset in one h5py read (cached)."""
    f = _state(image)['fields']
    if key not in f:
        f[key] = np.asarray(image[key][:])
    return f[key]


class _Grid:
    """Geometry of the merged uniform grid (built once per image/window/pitch).

    Cell pitch is dmin * 2^k (dmin: finest block pitch among the blocks used),
    on the lattice of the camera corner, so every block at least as coarse as
    the grid covers whole cells (each pixel replicated r x r). Blocks finer
    than the grid are area-averaged into cells (NaN ignored); by the quadtree
    structure each cell is covered either by one coarser block or entirely by
    finer ones. Arrays are indexed internally as [alpha, beta].
    """

    max_cells = 2e9   # refuse grids larger than this (~8 GB per float32 field)

    def __init__(self, image, halfrange, window=None, pitch=None):
        alpha = _field(image, 'alpha')
        beta = _field(image, 'beta')
        N, npx2 = alpha.shape
        n = int(round(np.sqrt(npx2)))
        A = alpha.reshape(N, n, n)
        B = beta.reshape(N, n, n)
        bp = np.abs(A[:, 1, 0] - A[:, 0, 0])   # block pixel pitch
        x0 = A[:, 0, 0] - 0.5 * bp              # block lower edges
        y0 = B[:, 0, 0] - 0.5 * bp              # in (alpha, beta)
        x1 = x0 + n * bp
        y1 = y0 + n * bp
        cx, cy = x0.min(), y0.min()             # camera corner: lattice origin

        if window is not None:                  # plotted (x, y=-beta) -> (alpha, beta)
            wx0, wx1 = window[0], window[1]
            wy0, wy1 = -window[3], -window[2]
        elif halfrange is not None:
            wx0 = wy0 = -halfrange
            wx1 = wy1 = halfrange
        else:
            wx0, wx1, wy0, wy1 = -np.inf, np.inf, -np.inf, np.inf
        keep = (x1 > wx0) & (x0 < wx1) & (y1 > wy0) & (y0 < wy1)
        idx = np.nonzero(keep)[0]
        if len(idx) == 0:
            raise ValueError("no AMR blocks inside the requested window")

        self.dmin_data = bp[idx].min()
        k = 0
        if pitch is not None and pitch > self.dmin_data:
            k = int(np.ceil(np.log2(pitch / self.dmin_data) - 1e-9))
        gp = self.dmin_data * 2 ** k
        self.dmin = gp                           # grid cell pitch (rg)

        if window is None:                       # union of the blocks used, as before
            gx0, gx1 = x0[idx].min(), x1[idx].max()
            gy0, gy1 = y0[idx].min(), y1[idx].max()
        else:
            gx0, gx1, gy0, gy1 = wx0, wx1, wy0, wy1
        # snap outward to the lattice
        self.xmin = cx + np.floor((gx0 - cx) / gp + 1e-9) * gp
        self.ymin = cy + np.floor((gy0 - cy) / gp + 1e-9) * gp
        self.nx = int(np.ceil((gx1 - self.xmin) / gp - 1e-9))
        self.ny = int(np.ceil((gy1 - self.ymin) / gp - 1e-9))
        if self.nx * self.ny > self.max_cells:
            raise MemoryError(
                "grid of %d x %d cells at pitch %.3g rg; pass window= and/or a "
                "coarser pitch=" % (self.nx, self.ny, gp))
        self.n = n
        self.idx = idx

        # per-level (coarse -> fine) block lists with their grid offsets:
        #   ('rep', sel, r, ix0, iy0): pixel -> r x r cells
        #   ('avg', sel, p, x0, y0):   pixels averaged into cells
        self.levels = []
        lvl = np.round(np.log2(bp[idx] / gp)).astype(int)
        for L in sorted(np.unique(lvl), reverse=True):
            sel = idx[lvl == L]
            if L >= 0:
                r = 2 ** int(L)
                ix0 = np.round((x0[sel] - self.xmin) / gp).astype(np.int64)
                iy0 = np.round((y0[sel] - self.ymin) / gp).astype(np.int64)
                self.levels.append(('rep', sel, r, ix0, iy0))
            else:
                self.levels.append(('avg', sel, bp[sel[0]], x0[sel], y0[sel]))
        # imshow extent in plotted coordinates (x = alpha, y = -beta)
        self.extent_rg = (self.xmin, self.xmin + self.nx * gp,
                          -(self.ymin + self.ny * gp), -self.ymin)

    def rasterize(self, flat, dtype=np.float32):
        """flat: (N, n*n) values -> (ny, nx) image with y = -beta increasing upward."""
        n, nx, ny, gp = self.n, self.nx, self.ny, self.dmin
        G = np.full((nx, ny), np.nan, dtype=dtype)
        S = W = None
        for kind, sel, a, b, c in self.levels:
            blk = flat[sel].reshape(-1, n, n)
            if kind == 'rep':
                r, ix0, iy0 = a, b, c
                m = n * r
                if m <= 512:
                    ar = np.arange(m)
                    step = max(1, int(2e7 // (m * m)))
                    for i in range(0, len(sel), step):
                        ix = ix0[i:i + step, None] + ar           # (c, m)
                        iy = iy0[i:i + step, None] + ar
                        vb = blk[i:i + step].astype(dtype, copy=False)
                        if r > 1:
                            vb = vb.repeat(r, axis=1).repeat(r, axis=2)
                        if ix.min() >= 0 and ix.max() < nx and iy.min() >= 0 and iy.max() < ny:
                            G[ix[:, :, None], iy[:, None, :]] = vb
                        else:
                            M = (((ix >= 0) & (ix < nx))[:, :, None]
                                 & ((iy >= 0) & (iy < ny))[:, None, :])
                            IX = np.broadcast_to(ix[:, :, None], M.shape)[M]
                            IY = np.broadcast_to(iy[:, None, :], M.shape)[M]
                            G[IX, IY] = vb[M]
                else:   # large blocks (coarse block in a fine zoom): slice per block
                    for j in range(len(sel)):
                        ax0, ax1 = max(ix0[j], 0), min(ix0[j] + m, nx)
                        ay0, ay1 = max(iy0[j], 0), min(iy0[j] + m, ny)
                        if ax0 >= ax1 or ay0 >= ay1:
                            continue
                        px = (np.arange(ax0, ax1) - ix0[j]) // r
                        py = (np.arange(ay0, ay1) - iy0[j]) // r
                        G[ax0:ax1, ay0:ay1] = blk[j][np.ix_(px, py)]
            else:
                p, bx0, by0 = a, b, c
                if S is None:
                    S = np.zeros(nx * ny)
                    W = np.zeros(nx * ny)
                off = (np.arange(n) + 0.5) * p
                step = 200000
                for i in range(0, len(sel), step):
                    ci = np.floor((bx0[i:i + step, None] + off - self.xmin) / gp).astype(np.int64)
                    cj = np.floor((by0[i:i + step, None] + off - self.ymin) / gp).astype(np.int64)
                    ci = np.broadcast_to(ci[:, :, None], (len(ci), n, n))
                    cj = np.broadcast_to(cj[:, None, :], (len(cj), n, n))
                    v = blk[i:i + step].astype(np.float64)
                    ok = (ci >= 0) & (ci < nx) & (cj >= 0) & (cj < ny) & np.isfinite(v)
                    k = ci[ok] * ny + cj[ok]            # = index into G.ravel()
                    # area-weighted: levels finer than the grid may share a cell
                    S += np.bincount(k, weights=v[ok] * (p * p), minlength=nx * ny)
                    W += np.bincount(k, minlength=nx * ny) * (p * p)
        if S is not None:
            M = W > 0
            G.ravel()[M] = (S[M] / W[M]).astype(dtype)
        return G.T[::-1]      # rows = y=-beta ascending upward (origin='lower')

    def extent(self, mas):
        return tuple(e * mas for e in self.extent_rg)

    def centers(self, mas):
        """1D cell-centre coordinates (x, y=-beta), ascending, for contour."""
        x = (self.xmin + (np.arange(self.nx) + 0.5) * self.dmin) * mas
        y = (-(self.ymin + self.ny * self.dmin) + (np.arange(self.ny) + 0.5) * self.dmin) * mas
        return x, y


def _grid(image, halfrange, window=None, pitch=None):
    st = _state(image)['grids']
    key = (halfrange, None if window is None else tuple(window), pitch)
    if key not in st:
        st[key] = _Grid(image, halfrange, window, pitch)
    return st[key]


def _show(image, values, fig, ax, halfrange, mas, label, cmap, vmin, vmax, colorbar=True,
          window=None, pitch=None):
    g = _grid(image, halfrange, window, pitch)
    im = ax.imshow(g.rasterize(values), origin='lower', extent=g.extent(mas),
                   cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest',
                   aspect='equal')
    if colorbar:
        fig.colorbar(im, label=label, ax=ax)
    if window is None:
        ax.set_xlim(-halfrange * mas, halfrange * mas)
        ax.set_ylim(-halfrange * mas, halfrange * mas)
    else:
        ax.set_xlim(window[0] * mas, window[1] * mas)
        ax.set_ylim(window[2] * mas, window[3] * mas)
    return im


def read_data(folder, ind, data_id):
    """Same return as rapplot.read_data, but with whole-array reads.

    Note: like the original, min/max are initialised at -/+100, so max[j] is
    never below 100 (and min[j] never above -100).
    """
    file_name = folder + '/img_data_%d.h5' % ind
    print("Reading in: ", file_name)
    images = h5py.File(file_name, 'r')
    print(len(data_id))
    print(data_id)
    min_ = [-100., -100., -100., -100.]
    max_ = [100., 100., 100., 100.]
    for j in range(4):
        a = _field(images, data_id[j])
        max_[j] = max(max_[j], a.max())
        min_[j] = min(min_[j], a.min())
    return min_, max_, images


def plot_data_tau(image, data_id, ind, fig, ax, halfrange=40, mas=1, label="Stokes",
                  cmap="CMRmap", vmin=-8, vmax=2, **gridkw):
    v = np.log10(_field(image, data_id[ind]) + 1e-10)
    _show(image, v, fig, ax, halfrange, mas, label, cmap, vmin, vmax, **gridkw)


def plot_data_stokes(image, min, max, stokes_ind, data_id, fig, ax, halfrange=40, mas=1,
                     label="Stokes", cmap="afmhot", **gridkw):
    a = _field(image, data_id[stokes_ind])
    with np.errstate(invalid='ignore'):
        if stokes_ind == 0:
            _show(image, (a / max[stokes_ind]) ** 0.5, fig, ax, halfrange, mas, label, cmap, 0, 1, **gridkw)
        else:
            _show(image, a / max[stokes_ind], fig, ax, halfrange, mas, label, cmap, -1, 1, **gridkw)


def plot_data_norder(image, fig, ax, halfrange=40, mas=1, label="n (winding number)",
                     cmap="RdBu_r", vmin=-3, vmax=3, **gridkw):
    _show(image, _field(image, 'norder'), fig, ax, halfrange, mas, label, cmap, vmin, vmax, **gridkw)


def plot_data_norder_mino(image, fig, ax, halfrange=40, mas=1,
                          label="n (Mino-time half-orbit count)", cmap="RdBu_r", vmin=-3, vmax=3, **gridkw):
    # Sentinel -999 marks vortical geodesics (eta<=0): masked, not plotted.
    a = np.where(_field(image, 'norder_mino') < -100, np.nan, _field(image, 'norder_mino'))
    _show(image, a, fig, ax, halfrange, mas, label, cmap, vmin, vmax, **gridkw)


def plot_data_ncross(image, fig, ax, halfrange=40, mas=1, label="equatorial crossings",
                     cmap="tab10", vmin=0, vmax=9, **gridkw):
    _show(image, _field(image, 'ncross'), fig, ax, halfrange, mas, label, cmap, vmin, vmax, **gridkw)


def plot_data_dchi_grav(image, fig, ax, halfrange=40, mas=1, label=r"$\Delta\chi_{\rm grav}$ [deg]",
                        cmap="twilight_shifted", vmin=-90, vmax=90, **gridkw):
    # EVPA rotation from the sky-projected spin axis, source at infinity ->
    # camera (radians in the file, plotted in degrees; cyclic colormap since
    # it is defined mod 180 deg). NaN (captured / truncated rays) is
    # transparent.
    _show(image, np.degrees(_field(image, 'dchi_grav')), fig, ax, halfrange, mas, label,
          cmap, vmin, vmax, **gridkw)


def _contour(image, values, ax, levels, mas, halfrange, window=None, pitch=None, **kw):
    g = _grid(image, halfrange, window, pitch)
    x, y = g.centers(mas)
    Z = np.ma.masked_invalid(g.rasterize(values, dtype=np.float64))
    return ax.contour(x, y, Z, levels=levels, linewidths=0.5, **kw)


def overlay_norder_contours(image, ax, levels=(0.75, 1.25, 2.0, 3.0), mas=1, halfrange=40,
                            **gridkw):
    _contour(image, np.abs(_field(image, 'norder')), ax, levels, mas, halfrange, colors='cyan',
             **gridkw)


def overlay_norder_mino_contours(image, ax, levels=(1, 2, 3, 4, 5), mas=1, colors='cyan',
                                 label=True, halfrange=40, **gridkw):
    # norder_mino is a half-orbit count, so contour directly at integer
    # levels; sentinel -999 (vortical geodesics) is masked.
    a = np.where(_field(image, 'norder_mino') < -100, np.nan, _field(image, 'norder_mino'))
    cs = _contour(image, a, ax, levels, mas, halfrange, colors=colors, **gridkw)
    if label:
        ax.clabel(cs, levels, inline=True, fontsize=7, fmt='n=%d')


def critical_curve(a, inc_deg, n=4000):
    """Kerr critical curve (photon-shell image) in RAPTOR camera coordinates.

    Returns closed arrays (alpha, beta) in rg for a distant observer at
    inclination inc_deg (Bardeen 1973; Gralla & Lupsasca 2020, Eqs. 38-40):
      lambda = a + r/a (r - 2 Delta/(r-1)),  eta = r^3/a^2 (4 Delta/(r-1)^2 - r)
      alpha = -lambda/sin(i),  beta = +-sqrt(eta + a^2 cos^2 i - lambda^2 cot^2 i)
    over the photon-shell radii r where beta^2 >= 0. Same convention as
    initialize_photon() in metric.c (lambda = -alpha sin i, p_theta = beta),
    so plot it as (alpha, -beta) like the image, although it is up-down
    symmetric anyway. a < 0 is handled by mirroring alpha.
    """
    th = np.radians(inc_deg)
    s = np.sign(a) if a != 0 else 1.
    a = abs(a)
    if a < 1e-6:                       # Schwarzschild: circle of radius sqrt(27)
        phi = np.linspace(0., 2. * np.pi, n)
        return np.sqrt(27.) * np.cos(phi), np.sqrt(27.) * np.sin(phi)
    def lam_b2(r):
        D = r * r - 2. * r + a * a
        lam = a + r / a * (r - 2. * D / (r - 1.))
        eta = r ** 3 / a ** 2 * (4. * D / (r - 1.) ** 2 - r)
        return lam, eta + a * a * np.cos(th) ** 2 - lam ** 2 / np.tan(th) ** 2

    # beta^2 >= 0 on [r1, r2], inside the photon shell [r_pro, r_retro];
    # beta^2 is single-humped there, so bisect for the two roots.
    rpro = 2. * (1. + np.cos(2. / 3. * np.arccos(-a)))
    rret = 2. * (1. + np.cos(2. / 3. * np.arccos(a)))
    rs = np.linspace(rpro, rret, 2001)
    rpk = rs[np.argmax(lam_b2(rs)[1])]
    ends = []
    for lo, hi in ((rpro, rpk), (rret, rpk)):      # lo: b2 < 0 side, hi: b2 > 0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if lam_b2(mid)[1] < 0. else (lo, mid)
        ends.append(hi)
    # cosine spacing clusters samples at the roots, where beta ~ sqrt(r - r_root)
    r = ends[0] + (ends[1] - ends[0]) * 0.5 * (1. - np.cos(np.linspace(0., np.pi, n)))
    lam, b2 = lam_b2(r)
    alpha = -s * lam / np.sin(th)
    beta = np.sqrt(np.clip(b2, 0., None))
    # upper half r_min -> r_max, lower half back, closed
    return (np.concatenate([alpha, alpha[::-1], alpha[:1]]),
            np.concatenate([beta, -beta[::-1], beta[:1]]))


def overlay_critical_curve(ax, a, inc_deg, mas=1, color='w', lw=0.6, ls='--', alpha=0.9,
                           **kw):
    """Draw the critical curve on ax in plotted coordinates (x=alpha, y=-beta)."""
    xa, yb = critical_curve(a, inc_deg)
    return ax.plot(xa * mas, -yb * mas, color=color, lw=lw, ls=ls, alpha=alpha, **kw)


def plot_data_polfrac(image, max, data_id, fig, ax, halfrange=10, mas=1, label="|m|", cmap="afmhot", **gridkw):
    I = _field(image, data_id[0])
    Q = _field(image, data_id[1])
    U = _field(image, data_id[2])
    with np.errstate(invalid='ignore', divide='ignore'):
        array = np.sqrt(Q ** 2. + U ** 2) / I
    array[I / max[0] < 1e-7] = 0
    _show(image, array, fig, ax, halfrange, mas, label, cmap, 0, 1, **gridkw)


def plot_data_RM(image, min, max, ind_1, ind_2, data_id, lam1, lam2, fig, ax, halfrange=10,
                 mas=1, label="RM", cmap="RdBu", **gridkw):
    I1 = _field(image, data_id[ind_1])
    Q1 = _field(image, data_id[ind_1 + 20])
    U1 = _field(image, data_id[ind_1 + 40])
    I2 = _field(image, data_id[ind_2])
    Q2 = _field(image, data_id[ind_2 + 20])
    U2 = _field(image, data_id[ind_2 + 40])

    EVPA_1 = 0.5 * np.angle(Q1 + 1j * U1)
    EVPA_2 = 0.5 * np.angle(Q2 + 1j * U2)
    EVPA_2[EVPA_1 == 0.0] = 0
    EVPA_1[EVPA_2 == 0.0] = 0

    d = EVPA_2 - EVPA_1
    wrap = np.abs(d) > 0.75 * np.pi
    d[wrap] -= np.sign(d[wrap]) * np.pi     # vectorised version of the original double loop
    d[I1 / max[0] < 1e-6] = 0
    d[I2 / max[0] < 1e-6] = 0

    _show(image, np.sign(d) * np.abs(d) ** 0.25, fig, ax, halfrange, mas, label, cmap, -1, 1, **gridkw)
