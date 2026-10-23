# Hi-COLA tests

From the repository root, with Python 3 and the dependencies installed:

```sh
python -m pip install -r requirements-test.txt
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

CI uploads `coverage.xml` to Codecov after generating the local reports. Before
using this integration, enable `Hi-COLACode/Hi-COLA` in Codecov and add its upload
token as the repository's GitHub Actions secret `CODECOV_TOKEN`. If the token is
absent, including on fork pull requests, the Codecov upload is skipped. Unit tests
and local coverage reporting still run.

With a token configured, upload failures fail the job; the HTML/XML artifacts
remain available. Codecov
analyses coverage produced by GitHub Actions; it does not execute the unit tests.
The README badge is deferred until collaborator approval and account setup.
The initial local public-branch coverage baseline was 41% on 21 October 2026.

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

### Python lint checks

The Ruff workflow checks pull requests and pushes to the default branch. To run
the same check locally:

```sh
python -m pip install ruff==0.17.0
ruff check .
```

This checks syntax and core correctness rules, without requiring formatting or
local Git hooks. FML and scratch work are excluded. The two frontend modules
that create SymPy symbols through `exec(declare_symbols())` have an exception
for undefined-name checks (`F821`); other selected rules still apply.
