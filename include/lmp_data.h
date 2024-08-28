#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <map>
#include "inputs.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;


/**
 * @class LMP_data
 * @brief Reads and stores data from a LAMMPS data file.
 *
 * The LMP_data class handles the reading, processing, and storage of
 * data from a LAMMPS data file. It includes functionality for managing
 * simulation box information, atom coordinates, molecular chains, and bonds.
 */
class LMP_data {
    public:
        /**
         * @brief Constructs an LMP_data object and reads data from the provided inputs.
         *
         * @param inputs The input parameters for the LAMMPS data file.
         */
        LMP_data(const Inputs &inputs);

        /**
         * @brief Retrieves the mass of a specific atom type.
         *
         * @param i The index of the atom.
         * @return double The mass of the atom.
         */
        double get_atom_mass(int i) const {
            return atom_masses.at(atom_types[i]);
        }

        //! Stores the initial simulation box dimensions
        MatrixXd box = MatrixXd::Zero(3, 3);
        //! Stores the initial coordinates of the atoms
        MatrixXd atom_coords;
        //! Stores the molecule number associated with each atom
        std::vector<int> atom_mols;
        //! Stores the type associated with each atom
        std::vector<int> atom_types;
        //! Stores the charge associated with each atom
        std::vector<double> atom_charges;
        //! Stores the mass corresponding to each atom type
        std::map<int, double> atom_masses;
        //! Stores the bonds of atoms by using a bond table
        std::vector<std::vector<int>> bond_table;
        //! Stores the atoms within molecular chains
        std::vector<std::vector<int>> chains;
        //! Indicates if the simulation box is triclinic (true) or orthogonal (false)
        bool triclinic = false;
    private:
        /**
         * @brief Identifies molecular chains within the system and associates atoms with them.
         *
         * This function processes the atom data to identify molecules
         * and assigns atoms to their respective molecular chains.
         */
        void _find_molecules();

        //! Reference to the input parameters
        const Inputs &_inputs;
};

