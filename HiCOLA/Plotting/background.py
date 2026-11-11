"""hicola_plot_background: the expansion history and density fractions of one run.

Reads the backend's cosmology_*.txt. With --reference, a panel under H/H0 shows
H/H_ref against the scale factor, with up to two fractional bands around 1.
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


def main(argv=None):
    tool = common.parser(__doc__.splitlines()[0], 'Simulation output folder containing cosmology_*.txt')
    args = common.parse(tool, argv)
    run = table(args.input, tool)
    a = run['a']
    with common.style() as plt:
        layout = [['H', 'dlogH', 'omega']]
        if args.reference:
            layout.append(['ratio', 'dlogH', 'omega'])
        figure, axes = plt.subplot_mosaic(layout, figsize=(14, 4.5 if not args.reference else 6),
                                          height_ratios=(2, 1) if args.reference else None)
        axes['H'].loglog(a, run['H/H0'])
        axes['H'].set_ylabel('H / H$_0$')
        axes['dlogH'].semilogx(a, run['dlogH/dloga'])
        axes['dlogH'].set_ylabel('d ln H / d ln a')
        for name, values in run.items():
            # Constant-zero fractions (for example massless neutrinos) only add clutter.
            if name.startswith('Omega') and name != 'OmegaNu_exact' and np.any(values):
                axes['omega'].semilogx(a, values, label=f'$\\Omega_{{\\rm {SUBSCRIPTS.get(name[5:], name[5:])}}}$')
        axes['omega'].set_ylabel('density fraction')
        axes['omega'].legend()
        if args.reference:
            reference = table(args.reference, tool)
            axes['ratio'].semilogx(a, run['H/H0']/common.interpolate_on(a, reference['a'], reference['H/H0']))
            common.ratio_axis(axes['ratio'], args.bands, f'H / H$_{{\\rm {args.reference_label}}}$')
            if args.bands:
                axes['ratio'].legend()
            axes['H'].sharex(axes['ratio'])
            axes['H'].tick_params(labelbottom=False)
        for name in ('ratio' if args.reference else 'H', 'dlogH', 'omega'):
            axes[name].set_xlabel('a')
        figure.suptitle(f'{args.label}: background')
        common.save(figure, args, 'background')
        plt.close(figure)
    return 0


if __name__ == '__main__':
    sys.exit(main())
