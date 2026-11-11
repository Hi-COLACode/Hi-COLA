"""hicola_plot_frontend: the expansion and force tables written by the Hi-COLA frontend.

Reads <model>_*_expansion.txt (a, E, E'/E) and <model>_*_force.txt (a, screening
coefficient, coupling and, for Horndeski models, G_G4/G_N). With --reference, a panel
under E shows E/E_ref, with up to two fractional bands around 1.
"""
import sys

from HiCOLA.Plotting import common

# The frontend writes these tables without headers; the column order is fixed by its writers.
EXPANSION = ('a', 'E', "E'/E")
FORCE = ('a', 'screening coefficient', 'coupling', 'G$_{G4}$ / G$_N$')


def table(folder, kind, names, tool):
    matches = sorted(folder.glob(f'*_{kind}.txt'))
    if len(matches) != 1:
        tool.error(f'Expected one *_{kind}.txt in {folder}; found {len(matches)}')
    return common.read_table(matches[0], names)


def main(argv=None):
    tool = common.parser(__doc__.splitlines()[0], 'Frontend output folder containing *_expansion.txt and *_force.txt')
    args = common.parse(tool, argv)
    expansion = table(args.input, 'expansion', EXPANSION, tool)
    force = table(args.input, 'force', FORCE, tool)
    force_names = [name for name in FORCE[1:] if name in force]
    with common.style() as plt:
        layout = [['E', 'dE', *force_names]]
        if args.reference:
            layout.append(['ratio', 'dE', *force_names])
        figure, axes = plt.subplot_mosaic(layout, figsize=(4.2*len(layout[0]), 4.5 if not args.reference else 6),
                                          height_ratios=(2, 1) if args.reference else None)
        a = expansion['a']
        axes['E'].loglog(a, expansion['E'])
        axes['E'].set_ylabel('E = H / H$_0$')
        axes['dE'].semilogx(a, expansion["E'/E"])
        axes['dE'].set_ylabel("E' / E")
        for name in force_names:
            axes[name].semilogx(force['a'], force[name])
            common.log_if_wide(axes[name], force[name])
            axes[name].set_ylabel(name)
        if args.reference:
            reference = table(args.reference, 'expansion', EXPANSION, tool)
            axes['ratio'].semilogx(a, expansion['E']/common.interpolate_on(a, reference['a'], reference['E']))
            common.ratio_axis(axes['ratio'], args.bands, f'E / E$_{{\\rm {args.reference_label}}}$')
            if args.bands:
                axes['ratio'].legend()
            axes['E'].sharex(axes['ratio'])
            axes['E'].tick_params(labelbottom=False)
        for name in ('ratio' if args.reference else 'E', 'dE', *force_names):
            axes[name].set_xlabel('a')
        figure.suptitle(f'{args.label}: frontend expansion and force tables')
        common.save(figure, args, 'frontend')
        plt.close(figure)
    return 0


if __name__ == '__main__':
    sys.exit(main())
