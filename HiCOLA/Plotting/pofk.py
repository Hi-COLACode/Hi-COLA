"""hicola_plot_pofk: the matter power spectrum of one run at every output redshift.

With --reference, a lower panel shows the boost P(k)/P_ref(k) at each redshift the
two runs share, with up to two fractional bands around 1.
"""
import logging
import re
import sys

from HiCOLA.Plotting import common

REDSHIFT = re.compile(r'_z(\d+\.\d+)\.txt$')


def spectra(folder, redshifts=None):
    """Map redshift to power spectrum file for the top-level pofk*_z*.txt files in ``folder``."""
    found = {}
    for path in sorted(folder.glob('pofk*_z*.txt')):
        match = REDSHIFT.search(path.name)
        if match and (redshifts is None or any(abs(float(match[1])-z) < 5e-4 for z in redshifts)):
            found[match[1]] = path
    return dict(sorted(found.items(), key=lambda item: -float(item[0])))


def main(argv=None):
    tool = common.parser(__doc__.splitlines()[0], 'Simulation output folder containing pofk*_z*.txt files')
    tool.add_argument('--redshifts', type=float, nargs='+', metavar='Z',
                      help='Only plot these output redshifts (default: all)')
    args = common.parse(tool, argv)
    runs = spectra(args.input, args.redshifts)
    if not runs:
        tool.error(f'No pofk*_z*.txt files in {args.input}')
    references = spectra(args.reference, args.redshifts) if args.reference else {}
    shared = [z for z in runs if z in references]
    if args.reference and not shared:
        tool.error(f'{args.reference} has no power spectra at the redshifts of {args.input}')
    for z in sorted(set(runs) - set(shared)) if args.reference else ():
        logging.warning('No reference power spectrum at z = %s; plotting it without a ratio', z)

    with common.style() as plt:
        if args.reference:
            figure, (top, bottom) = plt.subplots(2, 1, figsize=(7, 7), sharex=True, height_ratios=(2, 1))
        else:
            figure, top = plt.subplots(figsize=(7, 5))
        for (z, path), color in zip(runs.items(), common.sequential(len(runs))):
            table = list(common.read_table(path).values())
            k, power = table[0], table[1]
            top.loglog(k, power, color=color, label=f'z = {z}')
            if z in shared:
                reference = list(common.read_table(references[z]).values())
                bottom.semilogx(k, power/common.interpolate_on(k, reference[0], reference[1]), color=color)
        top.set_ylabel('P(k)  [(Mpc/h)$^3$]')
        top.set_title(f'{args.label}: matter power spectrum')
        top.legend(title='Output', loc='upper left', bbox_to_anchor=(1.01, 1))
        if args.reference:
            common.ratio_axis(bottom, args.bands, f'P / P$_{{\\rm {args.reference_label}}}$')
            if args.bands:
                bottom.legend(loc='upper left', bbox_to_anchor=(1.01, 1))
            bottom.set_xlabel('k  [h/Mpc]')
        else:
            top.set_xlabel('k  [h/Mpc]')
        common.save(figure, args, 'pofk')
        plt.close(figure)
    return 0


if __name__ == '__main__':
    sys.exit(main())
