# Spatial-Epidemiological-Model-Simulations

This repository contains MATLAB codes for studying delay-driven dynamics in a spatial or epidemiological modeling setting. The current scripts use DDE-BIFTOOL to define a delayed two-variable model, compute steady states, continue equilibrium branches, detect Hopf bifurcations, branch off periodic orbits, and examine the stability of the resulting periodic solutions through Floquet multipliers.

## Project overview

The model has two state variables, denoted in the code by `w` and `q`. The script treats the delay as a continuation parameter and studies how changes in the delay affect stability and oscillatory behavior. This makes the repository useful for numerical experiments involving delayed population, epidemic, or spatial interaction models where the main goal is to detect transitions from stable equilibria to periodic behavior.

The workflow is organized around three MATLAB files:

| File | Purpose |
|---|---|
| `Function_def6.m` | Defines the symbolic delayed model and generates the DDE-BIFTOOL-compatible function file `sym_congestionControlModel.m`. |
| `sym_congestionControlModel.m` | Automatically generated MATLAB function containing the right-hand side and derivatives required by DDE-BIFTOOL. This file should be regenerated from `Function_def6.m` when needed. |
| `Main_code6.m` | Runs the numerical continuation experiment, computes steady-state stability, locates Hopf bifurcations, branches off periodic solutions, and plots Floquet multiplier information. |

## Mathematical model used in the MATLAB code

The symbolic generator defines the delayed system

```text
dw/dt = delay*15*w*(1 - w/k) - delay*w*q/(5 + w),

dq/dt = delay*exp(-0.1*delay)*w(t - tau)*q(t - tau)/(5 + w(t - tau)) - delay*0.1*q.
```

Here `k` is the carrying-capacity-type parameter, `tau` represents the delay index used by DDE-BIFTOOL, and `delay` is the continuation parameter used in the simulations. The first equation combines logistic growth with nonlinear interaction, while the second equation includes delayed feedback through `w(t - tau)` and `q(t - tau)`.

## How to run the code

1. Install MATLAB and DDE-BIFTOOL.
2. Make sure the `ddebiftoolpath` variable in both MATLAB scripts points to your local DDE-BIFTOOL installation.
3. Run:

```matlab
Function_def6
```

This generates `sym_congestionControlModel.m`.

4. Then run:

```matlab
Main_code6
```

This performs the continuation and produces the stability and bifurcation plots.

## What the main script does

`Main_code6.m` performs the following steps:

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

## Dependencies

- MATLAB
- Symbolic Math Toolbox
- DDE-BIFTOOL
- DDE-BIFTOOL symbolic extension
- DDE-BIFTOOL periodic-orbit and normal-form utilities

## Notes

The file `sym_congestionControlModel.m` is generated automatically. It may be long because it contains derivative information up to the order requested by DDE-BIFTOOL. For readability, the cleaner source of the model is `Function_def6.m`, while `Main_code6.m` contains the continuation and bifurcation experiment.
