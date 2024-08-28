#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>

using Eigen::MatrixXd;
using Eigen::Vector3d;


/**
 * @class LMP_dump
 * @brief Reads and stores data from a LAMMPS dump file.
 *
 * The LMP_dump class is responsible for parsing and storing the data
 * from a LAMMPS dump file, including information on timesteps,
 * simulation box dimensions, and atom coordinates.
 */
class LMP_dump {
    public:
        /**
         * @brief Constructs an LMP_dump object and reads data from the specified dump file.
         *
         * @param path The path to the LAMMPS dump file.
         */
        LMP_dump(std::string path);

        //! The Number of atoms in the system
        int num_atoms;
        //! Stores the timestep numbers
        std::vector<int> timesteps;
        //! Stores the simulation box dimension in each timestep
        std::vector<MatrixXd> box;
        //! Stores the simulation box bounds of each timestep
        std::vector<MatrixXd> bounds;
        //! Indicates the type of the coordinates (e.g., 'x y z', 'xs ys zs') in the dump file
        std::string coords_type;
        //! Stores the atom coordinates of each timestep
        std::vector<MatrixXd> atom_coords;
};

