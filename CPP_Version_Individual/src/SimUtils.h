// By Minseung Kim, Seungwoo Yoon and Seungbum Koo
// KAIST, Daejeon, South Korea
// February 23, 2026

// SimUtils.h
#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// Common utility function (declared as inline to prevent multiple definitions)
inline std::vector<double> logspace(double startExp, double endExp, int num) {
    std::vector<double> result;
    if (num <= 1) { result.push_back(std::pow(10, endExp)); return result; }
    double step = (endExp - startExp) / (num - 1);
    for (int i = 0; i < num; ++i) result.push_back(std::pow(10, startExp + i * step));
    return result;
}

inline double sinwave(double freq, double t, double a0, double a, double b) {
    double val = (a0 - a) / b;
    if (val > 1.0) val = 1.0; if (val < -1.0) val = -1.0;
    double phase = std::asin(val);
    return a + std::sin(2.0 * M_PI * freq * t + phase) * b;
}

// Declaration of functions for each simulation run
int run_MMM_active(std::string MMM_type, double L_mn_input);
int run_MMM_passive(std::string MMM_type);
int run_TMM_active(double L_mn_input);
int run_TMM_passive();