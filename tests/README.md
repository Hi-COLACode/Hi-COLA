# Hi-COLA tests

From the repository root, with Python 3 and the dependencies installed:

```sh
python -m pip install -e '.[tests]'
MPLBACKEND=Agg python -m unittest discover -s tests -v
```

To measure Python statement and branch coverage:

```sh
MPLBACKEND=Agg python -m coverage run -m unittest discover -s tests -v
python -m coverage report
python -m coverage html
```

Open `htmlcov/index.html` for uncovered lines and branches. CI publishes a summary
and uploads HTML/XML reports as `python-coverage`, including when tests fail.
Coverage includes untested Python modules under `HiCOLA`, excludes FML, and does
not measure C++ execution. No minimum percentage is enforced yet.

After a successful default-branch push, CI publishes the generated badge value
to `coverage.json` on the dedicated `coverage-badges` branch. This branch contains
only badge data; it is created automatically on the first publication. Shields.io
reads that file for the README badge, with a short caching delay. Until the first
publication, the badge will be unavailable. Pull requests upload reports without
publishing. The badge shows the latest successful default-branch run; older runs
cannot overwrite a newer default-branch revision's result.

Only the separate publishing job has `contents: write`, using the built-in
`GITHUB_TOKEN`; no additional account or secret is needed. Repository rules must
allow that token to update `coverage-badges`. The README endpoint targets the
canonical `Hi-COLACode/Hi-COLA` repository; forks can change that owner/name in the
badge URL to display their own published value. The initial local public-branch
baseline was 41% on 21 October 2026.

Discovery is restricted to `tests/`; neither `FML/` nor `HiCOLA/FML/` is
collected or modified. The base branch is `Hi-COLA_public`. Its C++ code is
entirely under FML, so this suite tests Python and the actual SciPy integration
and root-solving boundaries. No C++ simulator coverage is claimed.

The suite checks scalar and NumPy numeric types, list/array inputs where
supported, malformed types and configurations, callback argument order,
normal and threshold solver paths, and exported table layout and precision.
Independent references include exact pure-component expansion histories,
finite-difference density derivatives, hand-derived polynomial derivatives,
an exact linear closure root, and an integrated Einstein-de Sitter solution.
Tolerances account for finite differences, ODE error, and five-digit ASCII
output. Inputs outside physical domains do not yet have consistent validation;
these tests do not imply that all such inputs raise sensible errors.

Two tests are explicitly marked `expectedFailure` to record existing defects:

- `test_hubble_closure_argument_substitution`: closure index 0 references an
  undefined `f_phi`, rather than passing the documented callback arguments.
- `test_g3_phi_derivative`: the `phi` selector constructs a symbol rather than
  taking the derivative.

Unexpected success fails the suite, so remove the marker when fixing either
implementation. Production code is unchanged in this branch.

Coverage is a starting baseline, not every function. Remaining work includes
independent references for the full Horndeski expressions and lambdified model
builder, complete scanner/legacy-reader paths, ESS Cartesian export and
tracking/root helpers, broader full-solver models, and consistent domain/error
contracts. Executable generation and plotting scripts need subprocess tests
because importing them executes their workflows.

If Hi-COLA-specific backend classes under `HiCOLA/FML/COLASolver` are later
included in scope, use tiny deterministic simulations through a subprocess
wrapper: check exit status and error diagnostics, parse the generated expansion
and force tables, compare against analytic GR limits, then compare power
spectra statistically with explicit tolerances. Do not compare stochastic
simulation outputs byte-for-byte or replace the simulator with a mock.
