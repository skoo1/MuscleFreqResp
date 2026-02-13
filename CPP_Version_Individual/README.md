
# MuscleFreqResp (MFR)

This repository contains MATLAB implementations of two numerical Hill-type muscle models and scripts for frequency response analysis as described in our submitted manuscript.

The repository includes:
- C++ code, tested on **Windows 11** using **Visual Studio 2022**
- Implementations of the Thelen (2003) and Millard (2012) muscle models
- Scripts for generating Bode plots for both active and passive state analyses

---

# Installation & Usage

### 1. Prerequisite (Crucial Step)

---

### 2. Running the C++ code

---

# Muscle Models

## 1. Thelen Muscle Model (TMM)
MATLAB implementation based on:

> Thelen, D. G. (2003). *Adjustment of muscle mechanics model parameters to simulate dynamic contractions in older adults.* Journal of Biomechanical Engineering, 125(1), 70–77.

**Implementation Notes**
- The governing equations described in the paper were implemented directly in MATLAB.
- To ensure correctness, core curve functions (F–L, F–V, passive elements, activation dynamics) were reproduced from the original formulation.

## 2. Millard Muscle Model (MMM)
The Millard model is implemented using the official [MATLAB port](https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort) distributed by the original author.

Reference:
> Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). Flexing computational muscle: modeling and simulation of musculotendon dynamics. *Journal of Biomechanical Engineering*, 135(2), 021005.

**Implementation Notes**
- This project provides wrapper scripts that configure and test three variants:
  - `MMM` (classic model)
  - `MMM-DEq` (damped-equilibrium version)
  - `MMM-Rigid` (rigid-tendon condition)

---

# Result Images

## Thelen Muscle Model Active Test
#### >> TMM_active
<img src="./images/TMM_active_result.png" alt="Thelen Active Result" width="400" />

## Thelen Muscle Model Passive Test
#### >> TMM_passive
<img src="./images/TMM_passive_result.png" alt="Thelen Passive Result" width="400" />

## Millard Muscle Model Active Test
#### >> MMM_active
<img src="./images/MMM_active_result.png" alt="Millard Active Result" width="400" />

## Millard Muscle Model Passive Test
#### >> MMM_passive
<img src="./images/MMM_passive_result.png" alt="Millard Passive Result" width="400" />

---

# Other Information

## Author
- **Minseung Kim**, Ph.D. Student (KAIST)
- **Seungbum Koo**, Ph.D. (KAIST)

## Revision Notes
- **August 2025**: Initial implementation
- **August 2025**: README and main script update
- **November 2025**: Debugged runtime errors
- **January 2026**: Minor bug fix and refactoring
- **February 2026**: C++ implementation