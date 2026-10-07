# Conservation and nonlinear solver settings

An independent isolated radiation enclosure exposed drift with default nonlinear heat convergence. Enforcing view-factor closure/reciprocity alone reduced the drift modestly; tightening heat convergence brought it below one microkelvin in mean temperature over 500 s. The initial tolerance criterion (0.005 K) was retained. Diagnostic results were preserved rather than replacing the acceptance threshold.

The CubeSat comparison keeps all physical inputs fixed. Both updated cases start from the same original coarse eighth-orbit state and solve two additional 10 s cycles. The balanced case adds internal view-factor closure and reciprocity; balanced_strict additionally sets relative nonlinear heat convergence to 1e-8. The open exterior enclosure is not forced closed.

| Settings | Final-cycle nodal change (K) | Energy RMS residual (W) | Residual / throughput (%) |
|---|---:|---:|---:|
| original | 0.0493 | 0.055399 | 0.32848 |
| balanced | 0.3014 | 0.009713 | 0.05756 |
| balanced_strict | 0.3085 | 0.000001 | 0.00001 |

The largest main-part average change from the original trajectory is 0.3284 K for balanced_strict. Its pooled equal-part-weight journal RMSE is 4.28 K. The conservation improvement is therefore distinguished from resolving the paper-comparison residual.

| Part | Original journal RMSE (K) | Balanced/strict journal RMSE (K) |
|---|---:|---:|
| panel1 | 6.93 | 6.87 |
| panel2 | 3.78 | 3.73 |
| panel3 | 6.63 | 6.62 |
| panel4 | 5.82 | 5.74 |
| panel5 | 3.76 | 3.84 |
| panel6 | 2.97 | 3.05 |
| pcb4 | 2.45 | 2.44 |
| pcb3 | 3.64 | 3.64 |
| pcb2 | 3.19 | 3.11 |
| pcb1 | 2.02 | 2.02 |
| battery | 2.67 | 2.62 |

## Release settings

The current generator uses closed internal-enclosure balancing and strict nonlinear heat convergence. The controlled improved epsilon=0.5 input and its initial nodal values are captured in the public package. Earlier sensitivity/emissivity studies retain their original settings and are labeled historical. They do not establish formal convergence of every updated case.

The energy diagnostic uses volume-integrated nodal temperature averages and the radiation surface-temperature approximation described in the analysis script; it is not an exact integration of every solver reaction. The analytical closed-enclosure check provides an independent conservation test.

The solver reports generic warnings for zero temperature offset and the combination of balancing with a space temperature. The model uses Kelvin with TOFFST=0. Balancing is scoped to closed internal enclosure 1; the 2.7 K space temperature is scoped to open exterior enclosure 2. No exterior row sum is forced to unity. Solver/export logs report zero errors.

The remaining journal discrepancy is explicitly retained. The repository presents a verified model and numerical comparison study, with reproducible evidence and limitations, rather than a claim of exact author-model recovery.

[ANSYS VFSM documentation](https://mapdl.docs.pyansys.com/version/stable/mapdl_commands/aux12/_autosummary/ansys.mapdl.core._commands.aux12.radiosity_solver.RadiositySolver.vfsm.html) describes closure/reciprocity enforcement; [RADOPT documentation](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/ans_cmd/Hlp_C_RADOPT.html) describes radiosity and multipass flux convergence. The commands were exercised with MAPDL 2024 R1; documentation versions differ.
