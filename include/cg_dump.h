#pragma once
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "lammps_data.h"
#include "lammps_dump.h"
#include "cg_data.h"


using Eigen::MatrixXd;
using Eigen::Vector3d;

typedef std::vector<Vector3d> BeadsCoords;


// Stores specifications of the dumped timesteps (both AA and CG)
struct Timestep {
    int timestep;
    MatrixXd box;
    BeadsCoords coords;
};


// Class to store coarse grained timesteps
class CG_dump {
    public:
        CG_dump(const LMP_data &data, const LMP_dump &dump, const CG_data &cgd);
        // Store the coarse-grained timesteps
        std::vector<Timestep> cg_ts;
    private:
        // Compute the CG timesteps using the mapping matrix
        void _coarsen_timesteps();

        const LMP_data &_data;
        const LMP_dump &_dump;
        const CG_data &_cg_data;
};

