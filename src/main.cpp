#include <fmt/format.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "string_tools.h"
#include "inputs.h"
#include "lammps_data.h"
#include "lammps_dump.h"
#include "cg_data.h"
#include "cg_dump.h"
#include "distributions.h"
#include "process.h"
#include "outputs.h"


int main(int argc, char **argv) {
    Eigen::initParallel();
    // Take the path to the input file
    std::string input_file;
    if (argc != 2) {
        std::cerr << "Input file is not entered.\nTry ./kde <input_file>\n";
    }
    input_file = argv[1];
    // Store the input settings
    Inputs inputs(input_file);
    // Read and store the lammps format data file
    LMP_data lmp_data(inputs);
    // Read and store the lammps format dump file
    LMP_dump lmp_dump(inputs.lmp_dump);
    // Generate the lammps format coarse-grained data file
    CG_data cg_data(inputs, lmp_data);
    // Generate the lammps format coarse-grained dump file
    CG_dump cg_dump(lmp_data, lmp_dump, cg_data);
    // Computing the distributions using Kernel Density Estimation (KDE) if true
    KDE kde(inputs, lmp_data, lmp_dump, cg_data, cg_dump);
    // Compute the distributions (parallel or serial)
    Process distributions(kde, cg_dump.cg_ts.size());
    if (inputs.process == "parallel") distributions.parallel();
    else distributions.serial();
    // Writing the results into output files
    Outputs outputs(inputs.tag, lmp_data, lmp_dump, distributions);
}

