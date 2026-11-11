"""Regression-only figures: the checks against their tolerances, and velocity statistics.

Differences are measured in units of the comparison allowance atol + rtol*|baseline|,
so a value of 1 is the pass limit of the summary comparison.
"""
from pathlib import Path

import numpy as np
from compare import TOLERANCES
from HiCOLA.Plotting import common

EXCEEDED_COLOR = '#e34948'
# Ratios below this, including exact reproduction (zero), are shown as upper limits on the log axis.
OVERVIEW_FLOOR = 1e-3


def tolerance_ratio(actual, reference, kind):
    """Return the signed difference in units of the comparison's allowance."""
    rtol, atol = TOLERANCES[kind]
    return (actual-reference)/(atol+rtol*np.abs(reference))


def worst_ratio(actual, reference, kind):
    """Return the largest |difference| / allowance; above 1 the comparison fails."""
    actual, reference = np.asarray(actual), np.asarray(reference)
    if actual.shape != reference.shape or actual.size == 0:
        raise ValueError(f'Cannot compare shapes {actual.shape} and {reference.shape}')
    return float(np.max(np.abs(tolerance_ratio(actual, reference, kind))))


def velocity_plot(actual, reference, title, path):
    """Compare mean velocity components and RMS speed from the z = 0 snapshot."""
    names = ('mean $v_x$', 'mean $v_y$', 'mean $v_z$', 'RMS speed')
    actual, reference = np.asarray(actual), np.asarray(reference)
    ratio = tolerance_ratio(actual, reference, 'velocity_statistics')
    worst = float(np.max(np.abs(ratio)))
    x = np.arange(len(names))
    with common.style() as plt:
        figure, (top, bottom) = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True, height_ratios=(2, 1))
        top.plot(x, reference, ls='none', marker='o', ms=9, color=common.COLORS[0], label='baseline')
        top.plot(x, actual, ls='none', marker='x', ms=9, mew=2, color=common.COLORS[1], label='candidate')
        top.set_ylabel('km/s')
        top.legend()
        top.set_title(f'{title}   (max |Δ|/tol {worst:.3g}{"  EXCEEDS TOLERANCE" if worst > 1 else ""})')
        bottom.plot(x, ratio, marker='o', color=EXCEEDED_COLOR if worst > 1 else common.INK)
        for limit in (-1, 1):
            bottom.axhline(limit, color='black', ls=':', lw=1)
        bottom.set_ylim(*(np.array([-1, 1])*max(1.25, 1.1*worst)))
        bottom.set_ylabel('Δ / tolerance')
        bottom.set_xticks(x, names)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path)
        plt.close(figure)
    return Path(path)


def overview_plot(checks, path):
    """Show how close each group of checks came to its tolerance, worst per group."""
    names = list(checks)[::-1]
    values = np.array([checks[name] for name in names])
    exceeded = values > 1
    below = values < OVERVIEW_FLOOR
    within = ~exceeded & ~below
    y = np.arange(len(names))
    with common.style() as plt:
        figure, axis = plt.subplots(figsize=(8, 1.5+.26*len(names)))
        axis.scatter(values[within], y[within], s=50, color=common.COLORS[0], label='within tolerance')
        if below.any():
            # Upper limits: a left-pointing arrow from the floor, as the true value lies below it.
            axis.errorbar(np.full(below.sum(), OVERVIEW_FLOOR), y[below], xerr=OVERVIEW_FLOOR*.6,
                          xuplims=True, fmt='o', ms=6, color=common.COLORS[0], capsize=0)
        if exceeded.any():
            axis.scatter(values[exceeded], y[exceeded], s=60, marker='X', color=EXCEEDED_COLOR,
                         label='exceeds tolerance')
        axis.axvline(1, color='black', ls=':', lw=1)
        axis.set_xscale('log')
        axis.set_xlim(OVERVIEW_FLOOR/4, max(5, 2*float(values.max())))
        axis.set_yticks(y, names)
        axis.set_ylim(-1, len(names))
        axis.grid(axis='y', visible=False)
        axis.set_xlabel('worst |candidate − baseline| / tolerance')
        axis.set_title('Summary checks against their tolerances')
        axis.legend(loc='upper left', bbox_to_anchor=(1.01, 1))
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path)
        plt.close(figure)
    return Path(path)
