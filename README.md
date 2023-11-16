This repository contains the R scripts for the following paper:

Loh, W. W. (In press). Estimating curvilinear time-varying treatment effects: combining g-estimation of structural nested mean models with time-varying effect models for longitudinal causal inference. *Psychological Methods.*

There are two R scripts with the main functions for carrying out the estimation procedures described in the paper:
- **time_varying-helper_funs.R**: standard g-estimation focusing on an end-of-study outcome
- **time_varying-lagged_effect-wrapper.R**: g-estimation using TVEM for multiple waves of outcome simultaneously

There are two folders implementing these procedures:
- **simulation-studies**: Simulation studies 1 and 2 as described in the paper
- **data-illustration-covid19_sleep**: Illustration using a real-world dataset (https://osf.io/gpxwa/)
