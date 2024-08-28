#pragma once
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "lmp_data.h"
#include "lmp_dump.h"
#include "cg_data.h"

using Eigen::MatrixXd;
using Eigen::Vector3d;

typedef std::vector<Vector3d> BeadsCoords;


/**
 * @struct Timestep
 * @brief Stores specifications of the dumped timesteps (both AA and CG).
 *
 * The Timestep structure holds information about each dumped timestep,
 * including the timestep index, the simulation box, and the coordinates of the beads.
 */
struct Timestep {
    //! Timestep index
    int timestep;
    //! Simulation box matrix
    MatrixXd box;
    //! Coordinates of the beads at this timestep
    BeadsCoords coords;
};


/**
 * @class CG_dump
 * @brief Manages coarse-grained timesteps for the simulation.
 *
 * The CG_dump class is responsible for storing and managing the coarse-grained
 * timesteps generated from the simulation.
 */
class CG_dump {
    public:
        /**
         * @brief Constructs a CG_dump object from LAMMPS data, dump, and coarse-grained data.
         *
         * @param data The LAMMPS data used for the simulation.
         * @param dump The dump file data containing the atomic coordinates.
         * @param cgd The coarse-grained data used to generate the timesteps.
         */
        CG_dump(const LMP_data &data, const LMP_dump &dump, const CG_data &cgd);

        //! Store the coarse-grained timesteps
        std::vector<Timestep> cg_ts;
    private:
        /**
         * @brief Computes the coarse-grained timesteps using the mapping matrix.
         *
         * This function processes the atomic data to generate the corresponding
         * coarse-grained timesteps.
         */
        void _coarsen_timesteps();

        //! Reference to the LAMMPS data
        const LMP_data &_data;
        //! Reference to the dump file data
        const LMP_dump &_dump;
        //! Reference to the coarse-grained data
        const CG_data &_cg_data;
};

