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
    """Geometry of the merged uniform grid (built once per image/window)."""

    def __init__(self, image, halfrange):
        alpha = _field(image, 'alpha')
        beta = _field(image, 'beta')
        N, npx2 = alpha.shape
        n = int(round(np.sqrt(npx2)))
        A = alpha.reshape(N, n, n)
        B = beta.reshape(N, n, n)
        pitch = np.abs(A[:, 1, 0] - A[:, 0, 0])
        x0 = A[:, 0, 0] - 0.5 * pitch          # block lower edges
        y0 = B[:, 0, 0] - 0.5 * pitch          # in (alpha, beta)
        x1 = x0 + n * pitch
        y1 = y0 + n * pitch

        if halfrange is not None:
            keep = (x1 > -halfrange) & (x0 < halfrange) & (y1 > -halfrange) & (y0 < halfrange)
        else:
            keep = np.ones(N, bool)
        idx = np.nonzero(keep)[0]
        if len(idx) == 0:
            raise ValueError("no AMR blocks inside the requested window")

        dmin = pitch[idx].min()
        self.xmin, self.ymin = x0[idx].min(), y0[idx].min()
        self.dmin = dmin
        self.nx = int(round((x1[idx].max() - self.xmin) / dmin))
        self.ny = int(round((y1[idx].max() - self.ymin) / dmin))
        self.n = n
        self.idx = idx
        # per-level (coarse -> fine) block lists with their fine-grid offsets
        self.levels = []
        lvl = np.round(np.log2(pitch[idx] / dmin)).astype(int)
        for L in sorted(np.unique(lvl), reverse=True):
            sel = idx[lvl == L]
            r = 2 ** int(L)
            ix0 = np.round((x0[sel] - self.xmin) / dmin).astype(int)
            iy0 = np.round((y0[sel] - self.ymin) / dmin).astype(int)
            m = n * r
            ar = np.arange(m)
            ix = (ix0[:, None] + ar)[:, :, None]
            iy = (iy0[:, None] + ar)[:, None, :]
            self.levels.append((sel, r, ix, iy))
        # imshow extent in plotted coordinates (x = alpha, y = -beta)
        self.extent_rg = (self.xmin, self.xmin + self.nx * dmin,
                          -(self.ymin + self.ny * dmin), -self.ymin)

    def rasterize(self, flat, dtype=np.float32):
        """flat: (N, n*n) values -> (ny, nx) image with y = -beta increasing upward."""
        n = self.n
        G = np.full((self.nx, self.ny), np.nan, dtype=dtype)
        for sel, r, ix, iy in self.levels:
            blk = flat[sel].reshape(-1, n, n).astype(dtype, copy=False)
            if r > 1:
                blk = blk.repeat(r, axis=1).repeat(r, axis=2)
            G[ix, iy] = blk
        return G.T[::-1]      # rows = y=-beta ascending upward (origin='lower')

    def extent(self, mas):
        return tuple(e * mas for e in self.extent_rg)

    def centers(self, mas):
        """1D cell-centre coordinates (x, y=-beta), ascending, for contour."""
        x = (self.xmin + (np.arange(self.nx) + 0.5) * self.dmin) * mas
        y = (-(self.ymin + self.ny * self.dmin) + (np.arange(self.ny) + 0.5) * self.dmin) * mas
        return x, y


def _grid(image, halfrange):
    st = _state(image)['grids']
    if halfrange not in st:
        st[halfrange] = _Grid(image, halfrange)
    return st[halfrange]


def _show(image, values, fig, ax, halfrange, mas, label, cmap, vmin, vmax, colorbar=True):
    g = _grid(image, halfrange)
    im = ax.imshow(g.rasterize(values), origin='lower', extent=g.extent(mas),
                   cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest',
                   aspect='equal')
    if colorbar:
        fig.colorbar(im, label=label, ax=ax)
    ax.set_xlim(-halfrange * mas, halfrange * mas)
    ax.set_ylim(-halfrange * mas, halfrange * mas)
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
                  cmap="CMRmap", vmin=-8, vmax=2):
    v = np.log10(_field(image, data_id[ind]) + 1e-10)
    _show(image, v, fig, ax, halfrange, mas, label, cmap, vmin, vmax)


def plot_data_stokes(image, min, max, stokes_ind, data_id, fig, ax, halfrange=40, mas=1,
                     label="Stokes", cmap="afmhot"):
    a = _field(image, data_id[stokes_ind])
    with np.errstate(invalid='ignore'):
        if stokes_ind == 0:
            _show(image, (a / max[stokes_ind]) ** 0.5, fig, ax, halfrange, mas, label, cmap, 0, 1)
        else:
            _show(image, a / max[stokes_ind], fig, ax, halfrange, mas, label, cmap, -1, 1)


def plot_data_norder(image, fig, ax, halfrange=40, mas=1, label="n (winding number)",
                     cmap="RdBu_r", vmin=-3, vmax=3):
    _show(image, _field(image, 'norder'), fig, ax, halfrange, mas, label, cmap, vmin, vmax)


def plot_data_norder_mino(image, fig, ax, halfrange=40, mas=1,
                          label="n (Mino-time half-orbit count)", cmap="RdBu_r", vmin=-3, vmax=3):
    # Sentinel -999 marks vortical geodesics (eta<=0): masked, not plotted.
    a = np.where(_field(image, 'norder_mino') < -100, np.nan, _field(image, 'norder_mino'))
    _show(image, a, fig, ax, halfrange, mas, label, cmap, vmin, vmax)


def plot_data_ncross(image, fig, ax, halfrange=40, mas=1, label="equatorial crossings",
                     cmap="tab10", vmin=0, vmax=9):
    _show(image, _field(image, 'ncross'), fig, ax, halfrange, mas, label, cmap, vmin, vmax)


def plot_data_dchi_grav(image, fig, ax, halfrange=40, mas=1, label=r"$\Delta\chi_{\rm grav}$ [deg]",
                        cmap="twilight_shifted", vmin=-90, vmax=90):
    # EVPA rotation from the sky-projected spin axis, source at infinity ->
    # camera (radians in the file, plotted in degrees; cyclic colormap since
    # it is defined mod 180 deg). NaN (captured / truncated rays) is
    # transparent.
    _show(image, np.degrees(_field(image, 'dchi_grav')), fig, ax, halfrange, mas, label,
          cmap, vmin, vmax)


def _contour(image, values, ax, levels, mas, halfrange, **kw):
    g = _grid(image, halfrange)
    x, y = g.centers(mas)
    Z = np.ma.masked_invalid(g.rasterize(values, dtype=np.float64))
    return ax.contour(x, y, Z, levels=levels, linewidths=0.5, **kw)


def overlay_norder_contours(image, ax, levels=(0.75, 1.25, 2.0, 3.0), mas=1, halfrange=40):
    _contour(image, np.abs(_field(image, 'norder')), ax, levels, mas, halfrange, colors='cyan')


def overlay_norder_mino_contours(image, ax, levels=(1, 2, 3, 4, 5), mas=1, colors='cyan',
                                 label=True, halfrange=40):
    # norder_mino is a half-orbit count, so contour directly at integer
    # levels; sentinel -999 (vortical geodesics) is masked.
    a = np.where(_field(image, 'norder_mino') < -100, np.nan, _field(image, 'norder_mino'))
    cs = _contour(image, a, ax, levels, mas, halfrange, colors=colors)
    if label:
        ax.clabel(cs, levels, inline=True, fontsize=7, fmt='n=%d')


def plot_data_polfrac(image, max, data_id, fig, ax, halfrange=10, mas=1, label="|m|", cmap="afmhot"):
    I = _field(image, data_id[0])
    Q = _field(image, data_id[1])
    U = _field(image, data_id[2])
    with np.errstate(invalid='ignore', divide='ignore'):
        array = np.sqrt(Q ** 2. + U ** 2) / I
    array[I / max[0] < 1e-7] = 0
    _show(image, array, fig, ax, halfrange, mas, label, cmap, 0, 1)


def plot_data_RM(image, min, max, ind_1, ind_2, data_id, lam1, lam2, fig, ax, halfrange=10,
                 mas=1, label="RM", cmap="RdBu"):
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

    _show(image, np.sign(d) * np.abs(d) ** 0.25, fig, ax, halfrange, mas, label, cmap, -1, 1)
