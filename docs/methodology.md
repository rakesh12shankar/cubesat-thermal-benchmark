# Methodology and limitations

## Governing problem

In each solid, rho c dT/dt = div(k grad T). All interfaces use perfect thermal continuity through shared mesh nodes. Vacuum regions are omitted; internal exposed surfaces exchange diffuse-gray radiation through MAPDL's radiosity enclosure formulation, with the PCBs and battery obstructing views. Exterior surfaces radiate to a 2.7 K sink with emissivity 0.72.

The baseline absorbs solar/albedo irradiance with alpha=0.77 and Earth infrared with epsilon=0.72. Its applied heat flux is q_abs = 0.77(q_solar + q_albedo) + 0.72 q_EarthIR. This is a declared spectral-absorption interpretation. Journal equations 20–23 show incoming terms without these factors, despite listing absorptivity in Table 2. A separate `equation_literal` case applies the incident sum directly and preserves all other inputs. Neither interpretation is silently presented as confirmed author code.

## Inputs from the paper

| Part | Density kg/m³ | Specific heat J/(kg K) | Conductivity W/(m K) |
|---|---:|---:|---:|
| Solar panels | 2325 | 1103 | 1.03 |
| Structure/supports | 2810 | 948 | 140 |
| PCBs | 2120 | 975 | 0.64 |
| Battery | 2247 | 1110 | 23 |

Internal emissivity is 0, 0.5 or 1 according to the case. There is no baseline internal heat generation. Thermal properties are constant. Coordinates used to construct the geometry are mm; MAPDL coordinates and thermal quantities use SI and Kelvin.

## Geometry reconstruction

The external envelope is 100 × 100 × 100 mm. Panels are 2 mm thick. Published full panel prisms would overlap at the edges, so this model partitions edge ownership: Z panels, then X panels, then Y panels. This retains six external faces of 0.01 m² each but slightly changes individual panel volumes from six full prisms.

A 5 mm square frame is positioned just inside the panels. Four 5 mm square supports run between the Z panels, at inferred X/Y offsets 10–15 and 85–90 mm. Four 90 × 90 × 2 mm PCBs have bottom surfaces at Z=20,37,61,78 mm. The 60 × 60 × 9 mm battery is centered above PCB2, at Z=39–48 mm. Corner/frame intrusions and support intersections are partitioned to avoid overlapping solids.

The PCB end dimensions are interpreted from the outer envelope, leaving 18 mm clearances inside the end panels. Frame/support positions, panel-edge ownership and the battery centering are reconstruction assumptions. The author CAD and exact contact definition were not available. The included interface table measures the modeled contact areas; it does not verify their exact correspondence to the author model.

Saved CAD verification checked all 351 body pairs for overlap. Total solid volume is 0.000247936 m³; modeled mass is 0.57848424 kg; heat capacity is 605.3914104 J/K. Perfect contact at sharp shared edges also depends on the chosen conformal representation.

## Orbit and heat inputs

The declared approximation is a circular, 431 km, zero-beta, nadir orbit with period 5,580 s. The paper lists 51.6° inclination and 0° ascending node; this project does not reproduce a complete date-dependent orbit propagator. Face 2 points toward Earth. Solar irradiance is 1,367 W/m², Earth infrared is 237 W/m², and the albedo coefficient is 0.3. The albedo calculation uses the stated cosine phase approximation. Earth view factors are integrated over the finite apparent disk using 80-point polar quadrature and 240 azimuth directions.

The analytic histories are compared with journal Figure 9 before absorption factors are applied. Eclipse transitions remain on the original time axis; both full-orbit errors and errors outside ±30 s transition bands are reported. No phase fitting was used.

## Solver and numerical evidence

ANSYS MAPDL 2024 R1 uses 8-node SOLID278 elements, transient conduction, and RDSF radiosity surface conditions. The initial in-plane mesh has 2,192 elements / 4,093 nodes; the finer mesh has 3,124 elements / 5,901 nodes. Thin regions are partitioned by geometry; h=15 or 10 mm describes the maximum in-plane subdivision, not a uniform cube edge length.

The current generator adds `VFSM,DEFINE,1,2,1000,1E-8` for closure and reciprocity of the closed internal enclosure, and `CNVTOL,HEAT,,1E-8,2,1E-6` for nonlinear heat convergence. Closure is not imposed on the open exterior enclosure. The earlier recorded emissivity/mesh/time-step studies use unadjusted view factors and default heat convergence; a separate balanced/strict case measures the update without modifying geometry or properties. [ANSYS VFSM documentation](https://mapdl.docs.pyansys.com/version/stable/mapdl_commands/aux12/_autosummary/ansys.mapdl.core._commands.aux12.radiosity_solver.RadiositySolver.vfsm.html) describes the closure/reciprocity options.

Thermal runs begin at 300 K or use an explicitly captured nodal state at an equivalent orbital phase. Final orbital cycles use a maximum 10 s time step. The emissivity-0.5 sensitivity cases check a 5 s step and hemisphere divisions 40 instead of 20.

The largest average-temperature differences are 0.815 K for the mesh, 0.398 K for the step, and 0.020 K for radiation divisions. These are changes between selected cases, not formal error bounds. The fine-mesh case's remaining nodal periodicity change is approximately 0.995 K; it is compared with coarse orbit 8, whose change is approximately 0.049 K. The mesh comparison therefore includes warm-up drift.

## What is and is not validated

The analytical cooling case checks the external radiation law and integration against a closed-form solution. The CAD/mesh checks establish internal consistency. Energy and periodicity checks establish conservation and settling within the reported diagnostics. Direct journal comparison establishes the actual remaining discrepancy.

An independent isolated gray enclosure has six separate equal-capacity walls, one initially at 400 K and the other five at 300 K. Its exact capacity-weighted mean is 316.6666667 K, with constant total stored energy. Default heat convergence produces approximately 0.006 K mean drift over 500 s; view-factor balancing alone leaves approximately 0.0053 K. With balancing and strict heat convergence, mean drift is below 0.000001 K and energy drift below 0.00012 J, while wall-temperature variance decreases. These checks distinguish conservation from agreement with the paper.

The spectral-absorption baseline has approximately 2–7 K volume-average RMSE and up to approximately 15 K maximum error. Numerical sensitivities below 1 K do not explain those residuals by themselves. The source boundary-condition ambiguity, geometry/contact choices and different treatment of partial radiation obstructions remain material limitations. The contribution of each has not been isolated. No arbitrary contact coefficient, conductivity adjustment or time shift was fitted to improve agreement.

The comparison is between two numerical models, not with experiment. This is a portfolio benchmark and discrepancy study, not a qualified spacecraft design tool.
