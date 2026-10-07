# Incident heat-input interpretation study

Journal equations 20–23 list solar, albedo and Earth infrared inputs without absorption factors; Table 2 nevertheless lists solar absorptivity 0.77 and external emissivity 0.72. The original reconstruction applies 0.77 to solar/albedo and 0.72 to Earth infrared. A separately labeled equation-literal case applies the complete incident sum as heat, leaving emissivities and all other physical/numerical inputs unchanged. This tests a source ambiguity and is not a fitted correction.

| Part | Spectral-absorption RMSE (K) | Equation-literal RMSE (K) | Literal mean bias (K) |
|---|---:|---:|---:|
| panel1 | 6.93 | 15.72 | 15.59 |
| panel2 | 3.78 | 16.62 | 16.33 |
| panel3 | 6.63 | 16.44 | 16.13 |
| panel4 | 5.82 | 16.79 | 16.13 |
| panel5 | 3.76 | 21.41 | 20.94 |
| panel6 | 2.97 | 20.59 | 20.14 |
| pcb4 | 2.45 | 17.79 | 17.27 |
| pcb3 | 3.64 | 17.29 | 16.74 |
| pcb2 | 3.19 | 17.52 | 17.46 |
| pcb1 | 2.02 | 18.19 | 17.78 |
| battery | 2.67 | 17.71 | 17.56 |

Equal-part-weight pooled RMSE changes from 4.31 K to 17.90 K. The equation-literal reading worsens aggregate agreement. This does not recover the author's actual heat-input implementation.

## Solver checks

- spectral_absorption: final-orbit maximum nodal cycle change 0.049 K; RMS energy residual 0.0554 W (0.328% of heat throughput).
- equation_literal: final-orbit maximum nodal cycle change 0.003 K; RMS energy residual 0.0793 W (0.360% of heat throughput).

The literal case starts from the exported original seventh-orbit baseline state at equivalent phase zero, then solves five 60 s warm-up cycles and two 10 s cycles. Its 1.024 K periodicity difference prompted two additional refined cycles seeded from its exact final nodal state. The settled results are used when available. Solver and export logs report zero errors. Its input and initial nodal values are retained for reproduction. The independent analytical cooling test verifies the external radiation law to less than 0.001 K at a 1 s step.

## Interpretation

The ambiguity is now quantified rather than left as an untested explanation. The spectral baseline and literal case are both retained; neither is silently substituted or described as a confirmed author model. Geometry/contact reconstruction and the source's approximation of partially obstructed radiation views remain unresolved. Their effects have not been independently isolated. The project is suitable to present as a verified orbital thermal model and documented numerical benchmark/discrepancy study, with the limits stated alongside the results.

See `load_interpretation.json` for the numerical metrics. The public figure `load_interpretation.png` can be regenerated from compact output tables.
