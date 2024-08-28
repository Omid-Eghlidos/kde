/**
 * @file main.cpp
 * @brief Entry point for the application that calculates distribution functions using
 * Kernel Density Estimation (KDE).
 *
 * This program reads input files, processes data using KDE, computes various
 * structural distribution functions, and outputs the results. It supports both
 * serial and parallel processing modes.
 *
 * @author Omid Eghlidos, Ph.D.
 * @date August 2024
 */


#include <fmt/format.h>
#include <Eigen/Dense>
#include <iostream>
#include <filesystem>
#include <vector>
#include <omp.h>
#include "string_tools.h"
#include "inputs.h"
#include "lmp_data.h"
#include "lmp_dump.h"
#include "cg_data.h"
#include "cg_dump.h"
#include "kde.h"
#include "process.h"
#include "outputs.h"


/**
 * @brief Main function to execute the KDE-based distribution calculations.
 *
 * This function initializes the input parameters, reads the LAMMPS data
 * and dump files, generates coarse-grained data, computes the distribution
 * functions, and writes the results to output files.
 *
 * @param argc The number of command-line arguments.
 * @param argv The command-line arguments.
 * @return int Returns 0 on successful execution.
 */
int main(int argc, char **argv) {
    //! Takes the path to the input file
    std::string input_file;
    if (argc == 1) {
        //! Checks for a default file if no arguments are provided
        if (std::filesystem::exists("kde.ini"))
            input_file = "kde.ini";
        else
            std::cerr << "Input file is not entered.\nTry ./kde <input_file>\n";
    }
    else if (argc == 2) {
        input_file = argv[1];
    }
    //! Stores the input settings
    Inputs inputs(input_file);
    //! Sets up parallel processing settings
    Eigen::initParallel();
    omp_set_num_threads(inputs.num_threads);
    omp_set_nested(1);
    //! Reads and stores the LAMMPS-formatted atomic data file
    LMP_data lmp_data(inputs);
    //! Reads and stores the LAMMPS-formatted atomic trajectory/dump file
    LMP_dump lmp_dump(inputs.lmp_dump);
    //! Generates the LAMMPS-formatted coarse-grained (CG) data file
    CG_data cg_data(inputs, lmp_data);
    //! Generates the LAMMPS-formatted coarse-grained (CG) trajectory/dump file
    CG_dump cg_dump(lmp_data, lmp_dump, cg_data);
    //! Computes the distributions using Kernel Density Estimation (KDE)
    KDE kde(inputs, lmp_data, lmp_dump, cg_data, cg_dump);
    //! Computes the distributions (parallel or serial)
    Process distributions(inputs, kde, cg_dump.cg_ts.size());
    //! Writes the distributions to output files
    Outputs outputs(inputs.tag, lmp_data, lmp_dump, distributions);

    //! Returns success
    return 0;
}

