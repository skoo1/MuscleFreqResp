#include <iostream>
#include "SimUtils.h"

int main() {
    try {
        std::cout << "========================================" << std::endl;
        std::cout << "   OpenSim Muscle Simulations Runner    " << std::endl;
        std::cout << "========================================" << std::endl;

        // 1. MMM Active
        run_MMM_active();
        std::cout << "\n[Done] MMM Active Simulation.\n" << std::endl;

        // 2. MMM Passive
        run_MMM_passive();
        std::cout << "\n[Done] MMM Passive Simulation.\n" << std::endl;

        // 3. TMM Active
        run_TMM_active();
        std::cout << "\n[Done] TMM Active Simulation.\n" << std::endl;

        // 4. TMM Passive
        run_TMM_passive();
        std::cout << "\n[Done] TMM Passive Simulation.\n" << std::endl;

        std::cout << "All simulations completed successfully." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "CRITICAL ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}