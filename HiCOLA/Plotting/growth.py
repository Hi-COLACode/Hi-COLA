"""hicola_plot_growth: linear growth D1(a, k) and G_eff/G of one run, for each k.

Reads the backend's gravitymodel_*_k*.txt files. With --reference, a panel under D1
shows D1/D1_ref for each k the two runs share, with up to two fractional bands around 1.
"""
import logging
import re
import sys

from HiCOLA.Plotting import common

WAVENUMBER = re.compile(r'_k(\d+(?:\.\d+)?)\.txt$')


def tables(folder):
    """Map k (as written in the file name) to its gravitymodel file, in increasing k."""
    found = {}
    for path in folder.glob('gravitymodel_*_k*.txt'):
        match = WAVENUMBER.search(path.name)
        if match:
            found[match[1]] = path
    return dict(sorted(found.items(), key=lambda item: float(item[0])))


def main(argv=None):
    tool = common.parser(__doc__.splitlines()[0], 'Simulation output folder containing gravitymodel_*_k*.txt files')
    args = common.parse(tool, argv)
    runs = tables(args.input)
    if not runs:
        tool.error(f'No gravitymodel_*_k*.txt files in {args.input}')
    references = tables(args.reference) if args.reference else {}
    if args.reference and not set(runs) & set(references):
        tool.error(f'{args.reference} has no growth tables at the k values of {args.input}')

    with common.style() as plt:
        layout = [['D1', 'Geff']]
        if args.reference:
            layout.append(['ratio', 'Geff'])
        figure, axes = plt.subplot_mosaic(layout, figsize=(12, 4.5 if not args.reference else 6),
                                          height_ratios=(2, 1) if args.reference else None)
        for (k, path), color in zip(runs.items(), common.sequential(len(runs))):
            run = common.up_to_today(common.read_table(path))
            label = f'k = {float(k):g} h/Mpc'
            axes['D1'].loglog(run['a'], run['D1(a,k)'], color=color, label=label)
            axes['Geff'].semilogx(run['a'], run['GeffG(a,k)'], color=color, label=label)
            if k in references:
                reference = common.read_table(references[k])
                axes['ratio'].semilogx(run['a'], run['D1(a,k)']/common.interpolate_on(
                    run['a'], reference['a'], reference['D1(a,k)']), color=color)
            elif args.reference:
                logging.warning('No reference growth table at k = %s', k)
        axes['D1'].set_ylabel('D$_1$(a, k)')
        axes['Geff'].set_ylabel('G$_{\\rm eff}$ / G')
        axes['Geff'].legend()
        if args.reference:
            common.ratio_axis(axes['ratio'], args.bands, f'D$_1$ / D$_{{1,\\rm {args.reference_label}}}$')
            if args.bands:
                axes['ratio'].legend()
            axes['D1'].sharex(axes['ratio'])
            axes['D1'].tick_params(labelbottom=False)
        for name in ('ratio' if args.reference else 'D1', 'Geff'):
            axes[name].set_xlabel('a')
        figure.suptitle(f'{args.label}: linear growth')
        common.save(figure, args, 'growth')
        plt.close(figure)
    return 0


if __name__ == '__main__':
    sys.exit(main())
