"""hicola_plot_growth: linear growth D1(a, k) and G_eff/G of one run, for each k.

Reads the backend's gravitymodel_*_k*.txt files. With --reference, each k the two
runs share gets its own panels, overlaying this run (solid) and the reference run
(dashed): D1, then D1/D1_ref with up to two fractional bands around 1, then G_eff/G.
"""
import logging
import re
import sys

from HiCOLA.Plotting import common

WAVENUMBER = re.compile(r'_k(\d+(?:\.\d+)?)\.txt$')
D1, GEFF = 'D1(a,k)', 'GeffG(a,k)'


def tables(folder):
    """Map k (as written in the file name) to its gravitymodel file, in increasing k."""
    found = {}
    for path in folder.glob('gravitymodel_*_k*.txt'):
        match = WAVENUMBER.search(path.name)
        if match:
            found[match[1]] = path
    return dict(sorted(found.items(), key=lambda item: float(item[0])))


def read(path):
    return common.up_to_today(common.read_table(path))


def single(plt, runs, args):
    """D1 and G_eff/G of one run, one line per k, darker for larger k."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for (k, path), color in zip(runs.items(), common.sequential(len(runs))):
        run = read(path)
        label = f'k = {float(k):g} h/Mpc'
        axes[0].loglog(run['a'], run[D1], color=color, label=label)
        axes[1].semilogx(run['a'], run[GEFF], color=color, label=label)
    axes[0].set(xlabel='a', ylabel='D$_1$(a, k)')
    axes[1].set(xlabel='a', ylabel='G$_{\\rm eff}$ / G')
    axes[1].legend()
    figure.suptitle(f'{args.label}: linear growth')
    return figure


def compared(plt, runs, references, shared, args):
    """One column of D1, D1 ratio and G_eff/G panels per shared k."""
    figure, cells = common.grid(plt, len(shared), 3, (2, 1, 1.4), (4.2, 7.5), max_columns=4)
    for index, (k, (top, middle, bottom)) in enumerate(zip(shared, cells)):
        run, reference = read(runs[k]), read(references[k])
        common.overlay(top, 'loglog', run['a'], run[D1], reference['a'], reference[D1])
        middle.semilogx(run['a'], run[D1]/common.interpolate_on(run['a'], reference['a'], reference[D1]),
                        color=common.RUN_COLOR)
        first = index % 4 == 0
        common.ratio_axis(middle, args.bands, f'D$_1$ / D$_{{1,\\rm {args.reference_label}}}$' if first else '')
        common.overlay(bottom, 'semilogx', run['a'], run[GEFF], reference['a'], reference[GEFF])
        top.set_title(f'k = {float(k):g} h/Mpc')
        if first:
            top.set_ylabel('D$_1$(a, k)')
            bottom.set_ylabel('G$_{\\rm eff}$ / G')
        bottom.set_xlabel('a')
    common.comparison_legend(figure, args)
    figure.suptitle(f'{args.label} and {args.reference_label}: linear growth')
    return figure


def main(argv=None):
    tool = common.parser(__doc__.splitlines()[0], 'Simulation output folder containing gravitymodel_*_k*.txt files')
    tool.add_argument('--wavenumbers', type=float, nargs='+', metavar='K',
                      help='Only plot these k values in h/Mpc (default: all)')
    args = common.parse(tool, argv)

    def chosen(found):
        if args.wavenumbers is None:
            return found
        return {k: path for k, path in found.items()
                if any(abs(float(k)-value) <= 1e-6*max(value, 1e-12) for value in args.wavenumbers)}

    runs = chosen(tables(args.input))
    if not runs:
        tool.error(f'No gravitymodel_*_k*.txt files in {args.input}')
    with common.style() as plt:
        if args.reference:
            references = chosen(tables(args.reference))
            shared = [k for k in runs if k in references]
            if not shared:
                tool.error(f'{args.reference} has no growth tables at the k values of {args.input}')
            for k in sorted(set(runs)-set(shared), key=float):
                logging.warning('No reference growth table at k = %s; it is not plotted', k)
            figure = compared(plt, runs, references, shared, args)
        else:
            figure = single(plt, runs, args)
        common.save(figure, args, 'growth')
        plt.close(figure)
    return 0


if __name__ == '__main__':
    sys.exit(main())
