# MuscleFreqResp (MFR)

This repository contains the code and data for analyzing the frequency response of Thelen (2003) and Millard (2012) muscle models during active and passive tests, as presented in our manuscript (link will be provided upon publication).

The repository includes:
- C++ OpenSim Implementations: Active and passive tests of Thelen (2003) and Millard (2012) models. (Tested on Windows 11 with Visual Studio 2022)
- MATLAB Implementations: The same tests conducted using MATLAB-based muscle models. (Tested on Windows 11 with MATLAB R2025b)
- Thelen (2003) Model Source: MATLAB implementation of the Thelen muscle model.
- Analysis Scripts: Scripts for generating Bode plots for both active and passive state analyses.

---

# Installation & Usage

### 1. Prerequisite (Crucial Step)
To use the C++ code, you must include the OpenSim-Core library.
1. Download the source code from the official [OpenSim-Core] (https://github.com/opensim-org/opensim-core) repository.
2. Compile and install the OpenSim-Core.

To use the Matlab code, you must include the original Matlab implementation code of the Millard Muscle Model.
1. Create a folder named **`MMM`** in the root directory of this repository.
2. Download the code from the official [Millard2012EquilibriumMuscleMatlabPort](https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort) repository.
3. Place all the downloaded files into the `MMM` folder.
* It is already included in this repository for your convenience.

### 2. Running the code
C++ Code - Compile it along with OpenSim-Core and run the muscle_app.exe.
- Modify CMakeLists.txt if necessary.
- The code should be compiled and work on all OS if OpenSim-Core has been compiled. 

Matlab Code - Simply download this repository and run the .m files in MATLAB. No complex installation is required.
- Open MATLAB and navigate to the repository folder.
- Run the main.m which runs active and passive tests on the muscles. 

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

---

# Result Images

## Active Tests on Thelen (2003) and Millard (2012)
<img src="./images/active_results.png" alt="Active Results" width="400" />

## Passive Tests on Thelen (2003) and Millard (2012)
<img src="./images/passive_results.png" alt="Passive Results" width="400" />

---

# Other Information

## Author
- **Minseung Kim**, Ph.D. Student (KAIST)
- **Seungwoo Yoon**, Ph.D. (KAIST)
- **Seungbum Koo**, Ph.D. (KAIST)

## Revision Notes
- **August 2025**: Initial implementation
- **August 2025**: README and main script update
- **November 2025**: Debugged runtime errors
- **January 2026**: Minor bug fix and refactoring
- **February 2026**: C++ implementation of tests using OpenSim functions