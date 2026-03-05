# MATLAB Implementation: Muscle Frequency Response (Individual)

This directory contains standalone MATLAB scripts for analyzing the frequency response of the **Thelen (2003)** and **Millard (2012)** muscle models. This version is ideal for quick testing, parameter tuning, and individual model verification without the need for C++ compilation.

---

## 🛠 Prerequisites

* **MATLAB R2025b** or later is recommended.

---

## 📂 Source Files Description

### 1. Muscle Test Scripts
These scripts perform the core simulation for each model and state:
* `TMM_active.m` / `TMM_passive.m`: Frequency response tests for the **Thelen Muscle Model (TMM)**.
* `MMM_active.m` / `MMM_passive.m`: Frequency response tests for the **Millard Muscle Model (MMM)**.

### 2. Visualization Tool
* `MFR_Bode_viewer.m`: A dedicated utility script to load simulation results and generate **Bode plots** (Magnitude and Phase) for comparison.

---

## 🚀 How to Use

1. **Open MATLAB** and navigate to this directory (`Matlab_Version_Individual`).
2. **Run a Test Script**:
   Select a model you wish to test (e.g., `TMM_active.m`) and run it.
   ```matlab
   run('TMM_active.m')
   ```
  
This will execute the frequency sweep simulation and store the response data in the workspace or as a temporary file.

3. **Visualize Results**:  
   Use the viewer script to generate professional plots:
   ```
   run('MFR_Bode_viewer.m')
   ```

---

## 📊 Key Features
* Individual Control: Each script is self-contained, allowing for easy modification of frequency ranges, amplitudes, or muscle parameters.
* Direct Comparison: Use the viewer to overlay results from different models (TMM vs. MMM) to analyze dynamic differences.
* Rapid Prototyping: Ideal for researchers who want to quickly test a specific hypothesis before running large-scale C++ simulations.

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

