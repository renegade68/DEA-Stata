# Reference-weight and result-saving regression tests

Run with Python 3, Git, and a licensed Stata 16+ batch executable. No additional
Python packages or network access are needed. Choose a new output directory:

```sh
python3 tests/test_saveref.py \
  --stata /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp \
  --output /tmp/dea-saveref-validation \
  --baseline-ref 41ea1dccda72b4a49c0964e358a0f76196709994
```

The optional baseline argument reproduces the old defects and compares 32
solver/specification combinations against that committed version. Omit it to
test the current working files alone. Each case runs in a fresh Stata batch
session with isolated ado paths. The runner returns a nonzero status if any
case fails; it checks a completion marker because Stata can exit with status
zero even when a do-file fails.

Coverage includes:

- Simplex and IPM; CRS, VRS, NIRS and NDRS; input/output orientation; stages 1/2.
- Reference variables equal to the returned matrix, primal feasibility,
  convexity restrictions, and unchanged baseline scores and weights.
- Sample alignment, missing values, names, metadata, collisions and automatic IDs.
- New and existing files, write failures, internal-variable exclusion and
  preservation of unrelated variables in memory.
- Automatic retries, multiplier weights, radial super-efficiency, partially
  infeasible VRS super-efficiency, selected other saving paths, and 120-DMU samples.

Retry tests use a separate generated copy that invalidates the initial simplex
attempt before the display/return step. The real IPM solver remains unchanged.
These tests establish control-flow behavior, not the frequency of numerical
failures in empirical data. VRS super-efficiency additionally exercises natural
infeasibility; an infeasible DMU remaining missing is an expected model outcome.

Generated source copies, `manifest.json`, `results.json`, and each case's
`test.do`, `test.log`, and datasets remain in the specified output directory.
Instrumented copies must not be distributed as the package. Successful old-defect
tests mean the defect was reproduced, not that the old version was correct.

These checks do not establish correctness of every DEA model or performance on
very large samples. Cross-engine weights need not be identical when optima are
nonunique; baseline comparisons use the same engine and the same data.

## Help examples and public syntax

```sh
python3 tests/test_help_syntax.py \
  --stata /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp \
  --output /tmp/dea-help-syntax-validation
```

This runs the marked example blocks directly from all 15 tracked help files,
including their original replication counts and `preserve`/`restore` sequences.
It decodes SMCL brace/indentation escapes and records the help source line for
each command in `help_line_map.json`; it does not silently insert `dmu()` or
modify example commands to make them pass. The fuzzy help's additional alpha=1
check in Remarks is also executed.

The additional cases compare both syntax forms, `if`/`in` with numeric and
string equality, solver choices, reference/saving combinations, documented
abbreviations and aliases, all 16 core models, bootstrap, wildcard/range
varlists, legacy `using`, and expected rejection of conflicting options.
Use `--only REGEX` to rerun selected cases in a new output directory.
