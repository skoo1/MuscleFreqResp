# MuscleFreqResp (MFR)

This repository contains MATLAB implementations of two numerical Hill-type muscle models and scripts for frequency response analysis as described in our submitted manuscript.

> M. Kim., S. Yoon., M. Roh., S. Koo. (2025). *Frequency Response Characterization of Numerical Hill-Type Muscle Models.* (Submitted)

The repository includes:
- MATLAB code, tested on **Windows 11** using **MATLAB R2025b**
- Implementations of the Thelen (2003) and Millard (2012) muscle models  
- Scripts for generating Bode plots for both active and passive state analyses  
- Supplementary data for model verification

---

# Muscle Models

## **1. Thelen Muscle Model (TMM)**  
MATLAB implementation based on:

> Thelen, D. G. (2003). *Adjustment of muscle mechanics model parameters to simulate dynamic contractions in older adults.* Journal of Biomechanical Engineering, 125(1), 70–77.

**Implementation Notes**

- The governing equations described in the paper were implemented directly in MATLAB.

- To ensure correctness, core curve functions (F–L, F–V, passive elements, activation dynamics) were reproduced from the original formulation.

- For transparency, we include comparison plots for F–L, F–V, and combined F–L–V relationships in the `supplementary/` directory.
  
- These checks confirm that the implemented equations match the expected behavior of the published model.

## **2. Millard Muscle Model (MMM)**  
The Millard model is implemented using the official MATLAB port distributed by the original author:

https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort

Reference:  
> Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). Flexing computational muscle: modeling and simulation of musculotendon dynamics. *Journal of Biomechanical Engineering*, 135(2), 021005.

**Implementation Notes**

- This project provides wrapper scripts that configure and test three variants:
  - `MMM` (classic model)
  - `MMM-DEq` (damped-equilibrium version)
  - `MMM-Rigid` (rigid-tendon condition)

---

# Main Scripts

The main MATLAB scripts for frequency response analysis and Bode plot generation were developed by the present authors for this study.

## **1. Structures**

<img width="1163" height="627" alt="image" src="https://github.com/user-attachments/assets/191a3087-b5ef-4843-be54-fc5e27ef3e9c" />

## **2. How to run**

Please follow the steps below to run the simulation:

(...)

---

# Other Informations

## Author  
- Minseung Kim, Ph.D. Student (KAIST)

## Revision Notes
- August 2025: Implementation
- August 2025: Edit Readme, main.m, etc.
- November 2025: Debugging a code run error
