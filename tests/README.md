# Hi-COLA tests

From the repository root, with Python 3 and the dependencies installed:

```sh
python -m pip install -r requirements-test.txt
MPLBACKEND=Agg python -m unittest discover -s tests -v
```

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
