# Spatial-Epidemiological-Model-Simulations

This repository contains MATLAB simulation codes for spatial and epidemiological modeling projects. The code is organized so that different modeling approaches remain separate. In particular, the DDE-BIFTOOL continuation scripts are not combined with the PDE simulation codes.

## Repository organization

| Folder | Purpose |
|---|---|
| `dde-bifurcation/` | MATLAB scripts for a delayed differential equation model analyzed with DDE-BIFTOOL. This part focuses on steady states, Hopf bifurcation detection, periodic-orbit continuation, and Floquet multiplier stability. |
| PDE simulation files | PDE-based spatial epidemiological simulations should remain separate from the DDE-BIFTOOL files. These scripts should be documented as their own model when they are added or reorganized. |

## DDE-BIFTOOL continuation codes

The folder `dde-bifurcation/` contains the delay differential equation continuation experiment. This is a separate numerical study from the PDE simulations.

| File | Purpose |
|---|---|
| `dde-bifurcation/Function_def6.m` | Defines the symbolic delayed model and generates the DDE-BIFTOOL-compatible function file `sym_congestionControlModel.m`. |
| `dde-bifurcation/sym_congestionControlModel.m` | Automatically generated MATLAB function containing the right-hand side and derivatives required by DDE-BIFTOOL. This file should be regenerated from `Function_def6.m` when needed. |
| `dde-bifurcation/Main_code6.m` | Runs the numerical continuation experiment, computes steady-state stability, locates Hopf bifurcations, branches off periodic solutions, and plots Floquet multiplier information. |

### DDE model used in the continuation code

The symbolic generator defines the delayed system

```text
dw/dt = delay*15*w*(1 - w/k) - delay*w*q/(5 + w),

dq/dt = delay*exp(-0.1*delay)*w(t - tau)*q(t - tau)/(5 + w(t - tau)) - delay*0.1*q.
```

Here `k` is the carrying-capacity-type parameter, `tau` represents the delay index used by DDE-BIFTOOL, and `delay` is the continuation parameter used in the bifurcation experiment. This model is treated as a DDE model only and should not be described as part of the PDE simulations.

### Running the DDE code

1. Install MATLAB and DDE-BIFTOOL.
2. Open MATLAB in the `dde-bifurcation/` folder.
3. Make sure the `ddebiftoolpath` variable in both MATLAB scripts points to your local DDE-BIFTOOL installation.
4. Run:

```matlab
Function_def6
```

This generates `sym_congestionControlModel.m`.

5. Then run:

```matlab
Main_code6
```

This performs the continuation and produces the stability and bifurcation plots.

### What the DDE main script does

`dde-bifurcation/Main_code6.m` performs the following steps:

1. Loads the required DDE-BIFTOOL folders.
2. Defines model parameters and a positive steady-state initial guess.
3. Corrects the steady state using DDE-BIFTOOL.
4. Computes the eigenvalue stability of the corrected steady state.
5. Continues the steady-state branch with respect to the delay parameter.
6. Locates Hopf bifurcation points along the branch.
7. Computes Hopf normal-form information.
8. Branches off periodic solutions from the Hopf point.
9. Computes Floquet multipliers along the periodic-orbit branch.
10. Plots the maximum Floquet multiplier magnitude against the delay parameter.

## PDE simulation codes

The PDE simulation codes should be documented separately from the DDE continuation codes. When PDE scripts are added, this section should describe the spatial domain, state variables, diffusion or movement terms, infection terms, boundary conditions, time-stepping method, and output plots.

A suggested organization is:

```text
pde-simulations/
  main simulation scripts
  helper functions
  plotting scripts
```

This keeps the PDE implementation independent from the DDE-BIFTOOL continuation workflow.

## Dependencies

For the DDE-BIFTOOL folder:

- MATLAB
- Symbolic Math Toolbox
- DDE-BIFTOOL
- DDE-BIFTOOL symbolic extension
- DDE-BIFTOOL periodic-orbit and normal-form utilities

For the PDE simulation files, dependencies should be listed separately after the PDE scripts are finalized.
