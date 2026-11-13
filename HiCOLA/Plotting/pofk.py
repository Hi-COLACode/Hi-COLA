"""hicola_plot_pofk: the matter power spectrum of one run at every output redshift.

With --reference, each redshift the two runs share gets its own panel, overlaying
this run (solid) and the reference run (dashed), with the boost P(k)/P_ref(k) below
and up to two fractional bands around 1.
"""
import logging
import re
import sys

from HiCOLA.Plotting import common

REDSHIFT = re.compile(r'_z(\d+\.\d+)\.txt$')
K_LABEL = 'k  [h/Mpc]'
P_LABEL = 'P(k)  [(Mpc/h)$^3$]'


def spectra(folder, redshifts=None):
    """Map redshift to power spectrum file for the top-level pofk*_z*.txt files in ``folder``."""
    found = {}
    for path in sorted(folder.glob('pofk*_z*.txt')):
        match = REDSHIFT.search(path.name)
        if match and (redshifts is None or any(abs(float(match[1])-z) < 5e-4 for z in redshifts)):
            found[match[1]] = path
    return dict(sorted(found.items(), key=lambda item: -float(item[0])))


def columns(path):
    table = list(common.read_table(path).values())
    return table[0], table[1]


def single(plt, runs, args):
    """All redshifts of one run on one axis, darker for later outputs."""
    figure, axis = plt.subplots(figsize=(7, 5))
    for (z, path), color in zip(runs.items(), common.sequential(len(runs))):
        axis.loglog(*columns(path), color=color, label=f'z = {z}')
    axis.set(xlabel=K_LABEL, ylabel=P_LABEL, title=f'{args.label}: matter power spectrum')
    axis.legend(title='Output', loc='upper left', bbox_to_anchor=(1.01, 1))
    return figure


def compared(plt, runs, references, shared, args):
    """One panel pair per shared redshift: both spectra above, their ratio below."""
    figure, cells = common.grid(plt, len(shared), 2, (2, 1), (4.6, 5))
    for index, (z, (top, bottom)) in enumerate(zip(shared, cells)):
        k, power = columns(runs[z])
        reference_k, reference_power = columns(references[z])
        common.overlay(top, 'loglog', k, power, reference_k, reference_power)
        bottom.semilogx(k, power/common.interpolate_on(k, reference_k, reference_power), color=common.RUN_COLOR)
        common.ratio_axis(bottom, args.bands, f'P / P$_{{\\rm {args.reference_label}}}$' if index % 3 == 0 else '')
        top.set_title(f'z = {z}')
        if index % 3 == 0:
            top.set_ylabel(P_LABEL)
        bottom.set_xlabel(K_LABEL)
    common.comparison_legend(figure, args)
    figure.suptitle(f'{args.label} and {args.reference_label}: matter power spectrum')
    return figure


def main(argv=None):
    tool = common.parser(__doc__.splitlines()[0], 'Simulation output folder containing pofk*_z*.txt files')
    tool.add_argument('--redshifts', type=float, nargs='+', metavar='Z',
                      help='Only plot these output redshifts (default: all)')
    args = common.parse(tool, argv)
    runs = spectra(args.input, args.redshifts)
    if not runs:
        tool.error(f'No pofk*_z*.txt files in {args.input}')
    with common.style() as plt:
        if args.reference:
            references = spectra(args.reference, args.redshifts)
            shared = [z for z in runs if z in references]
            if not shared:
                tool.error(f'{args.reference} has no power spectra at the redshifts of {args.input}')
            for z in sorted(set(runs)-set(shared)):
                logging.warning('No reference power spectrum at z = %s; it is not plotted', z)
            figure = compared(plt, runs, references, shared, args)
        else:
            figure = single(plt, runs, args)
        common.save(figure, args, 'pofk')
        plt.close(figure)
    return 0


if __name__ == '__main__':
    sys.exit(main())
