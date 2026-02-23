#include <iostream>
#include "SimUtils.h"

int main() {
    try {
        std::cout << "========================================" << std::endl;
        std::cout << "   OpenSim Muscle Simulations Runner    " << std::endl;
        std::cout << "========================================" << std::endl;

        // 1. MMM Active
        run_MMM_active("Classic", 0.6);
        std::cout << "\n[Done] MMM-Classic Active Simulation. (L_mn_input=0.6)\n" << std::endl;

        run_MMM_active("Classic", 0.8);
        std::cout << "\n[Done] MMM-Classic Active Simulation. (L_mn_input=0.8)\n" << std::endl;

        run_MMM_active("Classic", 1.0);
        std::cout << "\n[Done] MMM-Classic Active Simulation. (L_mn_input=1.0)\n" << std::endl;

        run_MMM_active("Classic", 1.2);
        std::cout << "\n[Done] MMM-Classic Active Simulation. (L_mn_input=1.2)\n" << std::endl;

        run_MMM_active("Classic", 1.4);
        std::cout << "\n[Done] MMM-Classic Active Simulation. (L_mn_input=1.4)\n" << std::endl;

        run_MMM_active("DEq", 0.6);
        std::cout << "\n[Done] MMM-DEq Active Simulation. (L_mn_input=0.6)\n" << std::endl;

        run_MMM_active("DEq", 0.8);
        std::cout << "\n[Done] MMM-DEq Active Simulation. (L_mn_input=0.8)\n" << std::endl;

        run_MMM_active("DEq", 1.0);
        std::cout << "\n[Done] MMM-DEq Active Simulation. (L_mn_input=1.0)\n" << std::endl;

        run_MMM_active("DEq", 1.2);
        std::cout << "\n[Done] MMM-DEq Active Simulation. (L_mn_input=1.2)\n" << std::endl;

        run_MMM_active("DEq", 1.4);
        std::cout << "\n[Done] MMM-DEq Active Simulation. (L_mn_input=1.4)\n" << std::endl;

        run_MMM_active("Rigid", 0.6);
        std::cout << "\n[Done] MMM-Rigid Active Simulation. (L_mn_input=0.6)\n" << std::endl;

        run_MMM_active("Rigid", 0.8);
        std::cout << "\n[Done] MMM-Rigid Active Simulation. (L_mn_input=0.8)\n" << std::endl;

        run_MMM_active("Rigid", 1.0);
        std::cout << "\n[Done] MMM-Rigid Active Simulation. (L_mn_input=1.0)\n" << std::endl;

        run_MMM_active("Rigid", 1.2);
        std::cout << "\n[Done] MMM-Rigid Active Simulation. (L_mn_input=1.2)\n" << std::endl;

        run_MMM_active("Rigid", 1.4);
        std::cout << "\n[Done] MMM-Rigid Active Simulation. (L_mn_input=1.4)\n" << std::endl;

        // 2. MMM Passive
        run_MMM_passive("Classic");
        std::cout << "\n[Done] MMM-Classic Passive Simulation.\n" << std::endl;

        run_MMM_passive("DEq");
        std::cout << "\n[Done] MMM-DEq Passive Simulation.\n" << std::endl;

        run_MMM_passive("Rigid");
        std::cout << "\n[Done] MMM-Rigid Passive Simulation.\n" << std::endl;

        // 3. TMM Active
        run_TMM_active(0.6);
        std::cout << "\n[Done] TMM Active Simulation. (L_mn_input=0.6)\n" << std::endl;

        run_TMM_active(0.8);
        std::cout << "\n[Done] TMM Active Simulation. (L_mn_input=0.8)\n" << std::endl;

        run_TMM_active(1.0);
        std::cout << "\n[Done] TMM Active Simulation. (L_mn_input=1.0)\n" << std::endl;

        run_TMM_active(1.2);
        std::cout << "\n[Done] TMM Active Simulation. (L_mn_input=1.2)\n" << std::endl;

        run_TMM_active(1.4);
        std::cout << "\n[Done] TMM Active Simulation. (L_mn_input=1.4)\n" << std::endl;

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