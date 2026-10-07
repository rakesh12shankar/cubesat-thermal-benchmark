# Reproduction

## Lightweight review

Install `requirements.txt`, run `check_inputs.py`, then run `render_results.py` and `check_release.py` from the repository root. Plotting uses the included compact tables and vector-extracted references. The resulting JSON reports per-part RMSE, maximum absolute difference and mean bias on the original 10 s phase grid.

The compact `last_orbit.csv` files retain both absolute solver time and phase time. They contain volume averages, stored energy and diagnostic input/output power; they are not complete nodal solutions. Full solver datasets remain local because they are large.

## Fresh MAPDL run

Use `python run_benchmark.py --ansys "PATH_TO_MAPDL_EXECUTABLE" --epsilon 0.5 --h 15 --duration 44640 --warm-orbits 6`. A typical Windows executable is `C:/Program Files/ANSYS Inc/v241/ansys/bin/winx64/ANSYS241.exe`. Other versions require their own compatibility checks.

The runner generates mesh and input, solves, checks solver/export error counts, then computes thermal averages, energy diagnostics and periodicity. It uses two shared-memory cores. It refuses to overwrite a nonempty case folder. Set `ANSYS_EXECUTABLE` as an alternative to `--ansys`.

The generated run has its own `applied_flux.csv`; postprocessing reads that snapshot before falling back to the declared baseline history. Inspect `solution_audit.json`, the solver log, and at least two matching final cycles. If the maximum nodal cycle change exceeds 1 K, more cycles are needed. For stronger settling, use a smaller tolerance and report it explicitly.

For epsilon=0 use at least 10 total cycles with the last two at 10 s; for epsilon=1 the recorded initial run used 5 cycles. These are starting points, not guaranteed periodicity for changed geometry or properties. Run each configuration in an isolated checkout or new case directory.

## Recorded histories

The original emissivity cases used three 60 s warm-up cycles, then 10 s cycles and additional continuation to reach 10,7,5 total cycles for epsilon=0,0.5,1 respectively. The coarse qualification run added an eighth epsilon=0.5 orbit. The fine mesh used five coarse warm-up cycles followed by two 10 s cycles. The step/radiation sensitivity cases started from the same original seventh-orbit nodal state and solved one additional cycle with the altered numerical setting.

The included baseline generator provides a fresh complete simulation rather than recreating the historical restart sequence. Small differences from recorded compact outputs can arise from warm-up histories. Use matching phases and actual periodicity diagnostics when comparing results.

The separate equation-literal load case has a captured complete APDL input with the seeded initial nodal values embedded. If included, rerun it with `python run_saved_case.py data/ansys/captured_equation_literal --ansys "PATH_TO_MAPDL_EXECUTABLE"`. This case uses its own applied-load snapshot and changes only the heat-input interpretation. Inspect its result audit before comparing.

The initial literal case was extended by two additional refined cycles because its cycle change was 1.024 K. The captured input represents those final two cycles, seeded from the preceding exact nodal state. The included provenance records this staged history. The runner prepends the captured initial nodal values to exported data for phase-zero periodicity checks.

The recommended generator additionally uses internal view-factor balancing and tighter nonlinear heat convergence. The controlled `captured_balanced_strict` input, when present, starts from the original coarse eighth-orbit state and solves two refined cycles with these settings. Run it with `python run_saved_case.py data/ansys/captured_balanced_strict --ansys "PATH_TO_MAPDL_EXECUTABLE"`. Historical output tables remain available for comparison.

## Analytical verification

Run `python verify_radiative_cooling.py --ansys "PATH_TO_MAPDL_EXECUTABLE"`. It solves an independently specified uniform radiating cube at 10 s and 1 s steps. The exact solution is T(t)=[T0^-3+3 epsilon sigma A t/(rho c V)]^-1/3 for a 0 K sink. The check requires a maximum error below 0.02 K at 1 s and improvement from the coarser step. It protects existing verification folders.

Run `python verify_closed_enclosure.py --ansys "PATH_TO_MAPDL_EXECUTABLE"` for the independent isolated gray enclosure. This checks invariant total stored energy/mean temperature, temperature bounds and decreasing variance under internal radiation exchange. It uses tight heat convergence and closed-enclosure view-factor adjustment. Earlier failed diagnostic variants are retained locally and their numerical summaries are included so the reason for tightening the solver settings is visible.

## CAD

The STEP and Parasolid files can be opened independently. The native part was built with SolidWorks 2019; the optional PowerShell COM scripts expect its interop assembly and template at the installed default paths. Adjust these paths for another installation. Use Windows PowerShell 5.1 for that workflow. CAD rebuilding is not required for the analytic mesh generator.

## Source curves

Journal Figures 9 and 11 were extracted from vector paths, calibrated to their own axis grid coordinates and interpolated at 10 s. The PDF's SHA-256 is retained in `config/source_hashes.json` and `results/journal_validation.json`; the paper is not redistributed. The temperature curves agree with the supplemental thesis curves to approximately 0.008 K after interpolation, which reflects their very similar plotted datasets and calibration, not independent physical validation.

Reference PDFs are available through the DOI/publisher and the thesis repository links in the README. PyMuPDF is an optional extraction dependency listed separately in `requirements-reference.txt`. Published reference samples resolve plot data, not the authors' unrounded solution histories.

To repeat the extraction, install `requirements-reference.txt` and run `python extract_journal_curves.py --pdf "PATH_TO_DOWNLOADED_PAPER" --output extracted_reference`. This protects existing extraction folders and checks the supplied PDF hash before using fixed vector paths. A differently packaged publisher PDF requires visual recalibration.

The packaged Python checks used Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.11.2. `requirements-lock.txt` records the main tested versions; `requirements.txt` allows a broader compatible range. GitHub Actions runs the lightweight checks; a licensed commercial solver is not run in CI. Workflow actions follow the official [checkout](https://github.com/actions/checkout) and [setup-python](https://github.com/actions/setup-python) examples.
