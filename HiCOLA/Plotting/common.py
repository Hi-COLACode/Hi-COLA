"""Arguments, style, table reading and ratio bands shared by the hicola_plot_* tools."""
import argparse
from contextlib import contextmanager
import logging
from pathlib import Path
import re

import numpy as np

# Categorical colours, assigned in this order; the first two also mark run and reference.
COLORS = ('tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan')
RUN_COLOR, REFERENCE_COLOR = COLORS[:2]
INK = '#52514e'
# Inner (narrower) band is darker, as in the Hi-COLA papers' grey accuracy bands.
BAND_COLORS = ('#c9c8c3', '#e6e5e1')
MAX_BANDS = 2
DEFAULT_BANDS = (0.01, 0.025)


def parser(description, input_help):
    """Return the argument parser shared by every tool."""
    tool = argparse.ArgumentParser(description=description)
    tool.add_argument('input', type=Path, help=input_help)
    tool.add_argument('--reference', type=Path,
                      help='The same kind of folder from a reference run (for example GR); '
                           'adds panels showing this run divided by the reference')
    tool.add_argument('--bands', type=float, nargs='*', default=list(DEFAULT_BANDS), metavar='FRACTION',
                      help=f'Up to {MAX_BANDS} fractional bands drawn around a ratio of 1 '
                           '(default: 0.01 0.025, i.e. ±1%% and ±2.5%%); give no values for none')
    tool.add_argument('--label', help='Run name used in titles (default: taken from the folder name)')
    tool.add_argument('--reference-label', help='Reference name used in labels (default: taken from its folder name)')
    tool.add_argument('--out', type=Path, help='Folder for the figure (default: plots/ beside an "output" folder, otherwise inside the input folder)')
    tool.add_argument('--format', default='png', choices=('png', 'pdf', 'svg'))
    return tool


def run_name(folder):
    """Name a run after its folder, skipping a trailing 'output' folder."""
    folder = Path(folder).resolve()
    return folder.parent.name if folder.name == 'output' else folder.name


def parse(tool, argv=None):
    args = tool.parse_args(argv)
    if len(args.bands) > MAX_BANDS:
        tool.error(f'--bands accepts at most {MAX_BANDS} values')
    if any(band <= 0 for band in args.bands):
        tool.error('--bands values must be positive fractions, for example 0.01')
    args.bands = sorted(args.bands)
    args.label = args.label or run_name(args.input)
    if args.reference is not None:
        args.reference_label = args.reference_label or run_name(args.reference)
    folder = Path(args.input).resolve()
    args.out = args.out or (folder.parent if folder.name == 'output' else folder)/'plots'
    logging.debug('Plot arguments: %s', args)
    return args


def pyplot():
    # The tools only write files, so they never need an interactive backend.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


@contextmanager
def style():
    """Apply the shared Hi-COLA figure style inside the block only."""
    plt = pyplot()
    from cycler import cycler
    with plt.rc_context({
            'axes.prop_cycle': cycler(color=COLORS), 'axes.grid': True, 'grid.alpha': .25,
            'axes.titlesize': 'medium', 'axes.labelsize': 'medium', 'font.size': 10,
            'legend.fontsize': 'small', 'legend.frameon': False, 'lines.linewidth': 1.8,
            'figure.dpi': 110, 'figure.constrained_layout.use': True, 'savefig.dpi': 110}):
        yield plt


def sequential(count):
    """Colours on one hue, light to dark, for an ordered quantity such as redshift or k."""
    return pyplot().get_cmap('Blues')(np.linspace(.4, 1, count)) if count else []


def read_table(path, names=None):
    """Read a whitespace table, naming columns from a '#' header or from ``names``.

    Returns a dict mapping column name to values, in file order.
    """
    data = np.loadtxt(path, ndmin=2)
    if data.size == 0:
        raise ValueError(f'Empty table: {path}')
    labels = list(names or [])
    if not labels:
        with open(path) as stream:
            header = stream.readline()
        if header.startswith('#'):
            # Power spectrum headers end in a 'ShotnoiseSubtracted = true' flag;
            # '(unit)' tokens belong to the preceding name.
            for token in re.sub(r'\S+\s*=.*$', '', header[1:]).split():
                if token.startswith('(') and labels:
                    labels[-1] += f' {token}'
                else:
                    labels.append(token)
    if len(labels) < data.shape[1]:
        labels += [f'column {column}' for column in range(len(labels), data.shape[1])]
    return dict(zip(labels, data.T))


def log_if_wide(axis, values, which='y'):
    """Use a log axis for positive quantities spanning two decades or more."""
    values = np.asarray(values, dtype=float)
    positive = values[values > 0]
    if positive.size == values.size and positive.size and positive.max()/positive.min() >= 100:
        (axis.set_yscale if which == 'y' else axis.set_xscale)('log')


def ratio_axis(axis, bands, label):
    """Prepare a ratio panel: a line at 1 and up to two fractional bands around it."""
    # Widest band first, so the narrower, darker band is drawn on top of it.
    for band, color in zip(reversed(bands), reversed(BAND_COLORS[:len(bands)])):
        axis.axhspan(1-band, 1+band, color=color, lw=0, zorder=0)
    axis.axhline(1, color=INK, lw=1, zorder=1)
    axis.set_ylabel(label)


def overlay(axis, method, x, y, reference_x, reference_y):
    """Draw this run (solid) and the reference run (dashed) on one axis."""
    getattr(axis, method)(x, y, color=RUN_COLOR)
    getattr(axis, method)(reference_x, reference_y, color=REFERENCE_COLOR, ls='--')


def comparison_legend(figure, args):
    """One figure-level legend naming the two runs and the ratio bands."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Line2D([], [], color=RUN_COLOR, label=args.label),
               Line2D([], [], color=REFERENCE_COLOR, ls='--', label=args.reference_label)]
    handles += [Patch(color=color, label=f'±{100*band:g}%') for band, color in zip(args.bands, BAND_COLORS)]
    figure.legend(handles=handles, loc='outside lower center', ncol=len(handles))


def grid(plt, count, rows_per_cell, height_ratios, cell_size, max_columns=3):
    """Small multiples: ``count`` cells of ``rows_per_cell`` stacked panels, sharing x.

    Returns the figure and a list of per-cell lists of axes; unused cells are hidden.
    """
    columns = min(max_columns, count)
    rows = -(-count//columns)
    figure, axes = plt.subplots(rows*rows_per_cell, columns, squeeze=False, sharex=True,
                                figsize=(cell_size[0]*columns, cell_size[1]*rows),
                                height_ratios=list(height_ratios)*rows)
    cells = []
    for index in range(rows*columns):
        row, column = divmod(index, columns)
        cell = [axes[row*rows_per_cell+offset, column] for offset in range(rows_per_cell)]
        if index < count:
            cells.append(cell)
        else:
            for axis in cell:
                axis.set_visible(False)
    return figure, cells


def interpolate_on(x, reference_x, reference_y):
    """Evaluate a reference column on another run's grid, interpolating in log x."""
    order = np.argsort(reference_x)
    return np.interp(np.log(x), np.log(reference_x[order]), reference_y[order])


def up_to_today(table):
    """Keep rows with a <= 1; backend tables continue past today only for interpolation."""
    keep = table['a'] <= 1+1e-12
    return {name: values[keep] for name, values in table.items()}


def save(figure, args, name):
    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out/f'{name}.{args.format}'
    figure.savefig(path)
    print(path)
    return path
