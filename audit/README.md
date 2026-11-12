# Scientific audit archive

The [independent audit of DOSY Lab 2.0.0](SCIENTIFIC_AUDIT_V2.0.md) records the findings for tag `v2.0.0`. The original [benchmark JSON](results/benchmark.json), [fit comparison](results/fit_and_ilt.png), [map comparison](results/map_fidelity.png) and executable audit scripts are preserved here. The unchanged historical source is in audit commit `42895f1`.

Use the [2.1 remediation plan](../docs/SCIENTIFIC_REMEDIATION_V2.1.md) to trace the corrections. The [QC study](../docs/QC_AND_MAPS_2.1.md), [map restoration comparison](../docs/ILT_MAP_RESTORATION_2.1.md) and [regularization study](../docs/ILT_REGULARIZATION_2.1.md) document later validation and remaining limitations. Reproducible 2.1 scripts and their results are in [benchmarks](../benchmarks/).

## Reproduce the historical audit

Run the audit from a detached checkout of the historical audit commit. The script checks the application source against the 2.0.0 tag and writes results inside that checkout. From the current repository, after installing the application dependencies:

~~~powershell
$py = (Resolve-Path '.\.venv\Scripts\python.exe').Path
git worktree add --detach ../DOSY-audit-reproduction 42895f1
Push-Location ../DOSY-audit-reproduction
$env:OPENBLAS_NUM_THREADS = '1'
& $py audit/run_audit.py
Pop-Location
~~~

The optional plot script also requires Matplotlib. Version 2.1 tests run separately with `tests/run.py`. The strychnine comparisons which use a saved user project require that local `.proj` and its source JDF. The tracked reference JDFs and synthetic fixtures are sufficient for the automated tests.
