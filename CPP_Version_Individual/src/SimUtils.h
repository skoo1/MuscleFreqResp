// SimUtils.h
#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// 공통 유틸리티 함수 (inline으로 선언하여 중복 정의 방지)
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

// 각 시뮬레이션 실행 함수 선언
int run_MMM_active();
int run_MMM_passive();
int run_TMM_active();
int run_TMM_passive();