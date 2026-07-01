# Dual-Pulse Gruneisen-Relaxation Photoacoustic Simulation Code

This repository contains the custom MATLAB code used to support the simulation and application-validation results in the manuscript:

**Mechanism-to-Application Co-Optimization of Dual-Pulse Gruneisen-Relaxation Photoacoustics via a Hybrid Computational Modeling Framework**

The code includes source-term parameter analysis, PAWS validation, GR-PAM resolution validation, and weak-nonlinearity compensation simulations.

## Repository structure

```text
dual-pulse-GR-PA-code/
│
├─ README.md
├─ LICENSE
├─ Code_Availability_Statement.txt
├─ .gitignore
│
├─ Fig2_Fig4_parameter_analysis/
│   ├─ main_generate_manuscript_Fig2_Fig4.m
│   ├─ run_all_figures.m
│   ├─ Best_Impulse_Response_Avg.mat
│   └─ outputs/
│
├─ Fig5_PAWS/
│   ├─ run_PAWS_dV_noisy.m
│   └─ outputs/
│
├─ Fig5_GRPAM_resolution/
│   ├─ run_GRPAM_resolution_Bessel_validation.m
│   └─ outputs/
│
├─ Fig6_compensation/
│   ├─ run_Fig6_CW_preheating.m
│   ├─ run_Fig6_multi_pulse_heat_accumulation.m
│   └─ outputs/
│
└─ outputs/
```

## Code description

### Fig2_Fig4_parameter_analysis

`main_generate_manuscript_Fig2_Fig4.m` generates the parameter-analysis results for the dual-pulse GR nonlinear photoacoustic model, including fluence-dependent scaling, absorption-coefficient-dependent scaling, pulse-interval and thermal-memory effects, coupled parameter maps, sensitivity analysis, and detector-response-related analysis.

`run_all_figures.m` is a lightweight wrapper for the Fig.2--Fig.4 script.

The empirical impulse-response file `Best_Impulse_Response_Avg.mat` is included in this folder and should contain `h_t_final`.

### Fig5_PAWS

`run_PAWS_dV_noisy.m` performs photoacoustic wavefront shaping simulation using a two-stage optimization workflow:

1. linear pre-focusing;
2. nonlinear GR-feedback-based optimization.

The simulation compares representative operating regimes, including global best, thermal mismatch, and weak nonlinearity.

### Fig5_GRPAM_resolution

`run_GRPAM_resolution_Bessel_validation.m` supports the GR-PAM resolution-validation analysis, including lateral PSF narrowing, axial optical sectioning, noise-limited effective resolution, and Bessel sidelobe-suppression validation.

### Fig6_compensation

`run_Fig6_CW_preheating.m` simulates continuous-wave preheating compensation under weak nonlinear conditions.

`run_Fig6_multi_pulse_heat_accumulation.m` simulates multi-pulse thermal accumulation compensation under weak nonlinear conditions.

## Requirements

The code was developed in MATLAB.

Recommended environment:

```text
MATLAB R2022b or later
Statistics and Machine Learning Toolbox
Image Processing Toolbox
k-Wave Toolbox
```

Some scripts use k-Wave for acoustic reference simulations. Please install k-Wave and add it to the MATLAB path before running these scripts.

## How to run

Run each script independently from its own folder.

For Fig.2--Fig.4 parameter analysis:

```matlab
cd Fig2_Fig4_parameter_analysis
run_all_figures
```

For PAWS validation:

```matlab
cd Fig5_PAWS
run_PAWS_dV_noisy
```

For GR-PAM resolution validation:

```matlab
cd Fig5_GRPAM_resolution
run_GRPAM_resolution_Bessel_validation
```

For weak-nonlinearity compensation:

```matlab
cd Fig6_compensation
run_Fig6_CW_preheating
run_Fig6_multi_pulse_heat_accumulation
```

Generated figures and source data are saved to the corresponding `outputs` folders.

## Notes

The scripts are intended to reproduce the main computational trends and figure-generation results reported in the manuscript. Randomized simulations, such as PAWS optimization, may show small run-to-run variations depending on the MATLAB version, random seed, and system configuration.

## License

This code is released under the MIT License.
