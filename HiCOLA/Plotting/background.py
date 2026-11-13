"""hicola_plot_background: the expansion history and density fractions of one run.

Reads the backend's cosmology_*.txt. With --reference, every panel overlays this run
(solid) and the reference run (dashed), and a panel under H/H0 shows H/H_ref with up
to two fractional bands around 1. Density fractions keep one colour per component,
so there the two runs differ by line style only.
"""
import sys

import numpy as np

from HiCOLA.Plotting import common

SUBSCRIPTS = {'M': 'm', 'R': 'r', 'Nu': '\\nu', 'MNu': 'm\\nu', 'Lambda': '\\Lambda'}


def table(folder, tool):
    matches = sorted(folder.glob('cosmology_*.txt'))
    if len(matches) != 1:
        tool.error(f'Expected one cosmology_*.txt in {folder}; found {len(matches)}')
    return common.up_to_today(common.read_table(matches[0]))


def fractions(table):
    # Constant-zero fractions (for example massless neutrinos) only add clutter.
    return [name for name, values in table.items()
            if name.startswith('Omega') and name != 'OmegaNu_exact' and np.any(values)]


def main(argv=None):
    tool = common.parser(__doc__.splitlines()[0], 'Simulation output folder containing cosmology_*.txt')
    args = common.parse(tool, argv)
    run = table(args.input, tool)
    reference = table(args.reference, tool) if args.reference else None
    a = run['a']
    with common.style() as plt:
        layout = [['H', 'dlogH', 'omega']]
        if reference:
            layout.append(['ratio', 'dlogH', 'omega'])
        figure, axes = plt.subplot_mosaic(layout, figsize=(14, 6 if reference else 4.5),
                                          height_ratios=(2, 1) if reference else None)
        if reference:
            common.overlay(axes['H'], 'loglog', a, run['H/H0'], reference['a'], reference['H/H0'])
            common.overlay(axes['dlogH'], 'semilogx', a, run['dlogH/dloga'], reference['a'], reference['dlogH/dloga'])
        else:
            axes['H'].loglog(a, run['H/H0'])
            axes['dlogH'].semilogx(a, run['dlogH/dloga'])
        axes['H'].set_ylabel('H / H$_0$')
        axes['dlogH'].set_ylabel('d ln H / d ln a')
        names = fractions(run)
        for name, color in zip(names, common.COLORS[2:] if reference else common.COLORS):
            label = f'$\\Omega_{{\\rm {SUBSCRIPTS.get(name[5:], name[5:])}}}$'
            axes['omega'].semilogx(a, run[name], color=color, label=label)
            if reference and name in reference:
                axes['omega'].semilogx(reference['a'], reference[name], color=color, ls='--')
        axes['omega'].set_ylabel('density fraction')
        axes['omega'].legend(title=f'solid: {args.label}\ndashed: {args.reference_label}' if reference else None)
        if reference:
            axes['ratio'].semilogx(a, run['H/H0']/common.interpolate_on(a, reference['a'], reference['H/H0']),
                                   color=common.RUN_COLOR)
            common.ratio_axis(axes['ratio'], args.bands, f'H / H$_{{\\rm {args.reference_label}}}$')
            axes['H'].sharex(axes['ratio'])
            axes['H'].tick_params(labelbottom=False)
            common.comparison_legend(figure, args)
            figure.suptitle(f'{args.label} and {args.reference_label}: background')
        else:
            figure.suptitle(f'{args.label}: background')
        for name in ('ratio' if reference else 'H', 'dlogH', 'omega'):
            axes[name].set_xlabel('a')
        common.save(figure, args, 'background')
        plt.close(figure)
    return 0


if __name__ == '__main__':
    sys.exit(main())
