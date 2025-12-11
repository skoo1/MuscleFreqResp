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

## **1. File structure**

<img width="1116" height="575" alt="image" src="https://github.com/user-attachments/assets/cf4929f5-9179-462b-a1a8-6566013458b5" />


## **2. How to run**

### Please follow the steps below to run the simulation.
> **Note:** Make sure that the repository paths are correctly set before running the code.
>
> In each model execution file, the variables `saveFolder1` and `saveFolder2` must be set to your data storage directory.
>
> For example:
> 'saveFolder1 = 'C:\Users\user\Desktop\MuscleFreqResp-main\MMM test\MMM\src\MMM_result\beforeFFT';'

#### 1. Run `MFR_main.m`
#### 2. Choose the simulation mode and adjust the simulation configurations
- Default muscle parameter values correspond to the human soleus muscle (based on the OpenSim Gait2392 model)

#### 3. After the simulation finishes, all Bode plot data points are saved in:  
- `(TMM or MMM)_result/afterFFT`

#### 4. To visualize the Bode plot under various conditions, press the 'Open Bode viewer' button on the GUI:  
- Multiple result files can be selected
- The selected files are passed to `MFR_Bode_viewer.m` to generate the corresponding Bode plots

### GUI Sample

<p align="center">
  <img src="https://github.com/user-attachments/assets/54e6c0e2-122b-405d-8b79-0c89e7f92fd5" width="700">
</p>

---

# Other Informations

## Author  
- Minseung Kim, Ph.D. Student (KAIST)

## Revision Notes
- August 2025: Initial implementation
- August 2025: README and main script update
- November 2025: Debugged runtime errors
