# Qualification overview

Full details, plots and assumptions are in [the study report](Numerical_Uncertainty_Report.md).

| Numerical family | Final adjacent average-curve change (K) | Fixed screening targets |
|---|---:|---|
| spatial | 0.0693 | PASS |
| temporal | 0.2116 | PASS |
| angular | 0.0016 | PASS |

Finest tested bonded case: `convergence_h7.5_features3`. Temporal and angular sequences use the original solid mesh. These separate tests do not establish a bound for every combination of mesh, time step and angular resolution.

Input bands are conditional on the six assumed engineering intervals. Their sensitivities were computed on the original mesh; their numerical-resolution limitations are reported separately. Interface variants retain their original meshes. Source geometry details and measured input-variability data remain necessary for stronger validation claims.
