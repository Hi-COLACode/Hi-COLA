"""Measure how close each group of summary outputs came to its tolerance.

The per-run figures come from the hicola_plot_* tools during the simulations;
this adds the two regression-only figures: velocity statistics and the overview.
"""
import logging

from compare import numeric, particles, velocity_statistics
from plots import overview_plot, velocity_plot, worst_ratio


def _group(name):
    for prefix, group in (('pofk', 'power spectra'), ('gravitymodel', 'growth tables'),
                          ('cosmology', 'expansion table')):
        if name.startswith(prefix):
            return group
    return 'other tables'


def summarise(work, cases, ranks, out):
    """Write the overview and velocity figures into ``out``.

    Returns the worst |difference| / tolerance for each group of checks and any
    output that could not be measured; problems are recorded, never raised.
    """
    checks, errors = {}, {}

    def attempt(key, action):
        try:
            action()
        except Exception as error:
            logging.debug('Could not summarise %s', key, exc_info=True)
            errors[key] = f'{type(error).__name__}: {error}'

    def record(group, candidate, baseline, kind):
        ratio = worst_ratio(numeric(candidate), numeric(baseline), kind)
        checks[group] = max(ratio, checks.get(group, 0.))

    frontend = work/'baseline'/'frontend'
    for path in sorted(frontend.glob('*/*.txt')):
        name = path.relative_to(frontend)
        attempt(f'frontend/{name}', lambda: record(
            f'frontend / {name.parts[0]}', work/'candidate'/'frontend'/name, path, 'tables'))

    for key in (f'{case}-np{n}' for n in ranks for case in cases):
        output = work/'baseline'/key/'output'
        for path in sorted(output.rglob('*.txt')):
            name = path.relative_to(output)
            attempt(f'{key}/{name}', lambda: record(
                f'{key} / {_group(path.name)}', work/'candidate'/key/'output'/name, path,
                'power' if 'pofk' in path.name else 'tables'))

        def velocities(key=key):
            actual, reference = (velocity_statistics(particles(sorted(
                (work/label/key/'output').rglob('gadget_z0.000.*')))[2]) for label in ('candidate', 'baseline'))
            checks[f'{key} / velocity statistics'] = worst_ratio(actual, reference, 'velocity_statistics')
            velocity_plot(actual, reference, f'{key}: velocity statistics at z = 0',
                          out/key/'velocity_statistics.png')
        attempt(f'{key}/velocity_statistics', velocities)

    if checks:
        attempt('overview', lambda: overview_plot(checks, out/'overview.png'))
    return dict(checks=checks, errors=errors)
