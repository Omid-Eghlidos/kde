#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <map>
#include "inputs.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;


// Reads and stores LAMMPS data file.
class LMP_data {
    public:
        LMP_data(const Inputs &inputs);
        // Deletes bonds that cross the simulation cell (for periodic crystal).
        void remove_bonds_across_pbc();
        // Find the mass for each atom type
        double get_atom_mass(int i) const {
            return atom_masses.at(atom_types[i]);
        }

        // Store the initial simulation box
        MatrixXd box = MatrixXd::Zero(3, 3);
        // Store the initial coordinates of the atoms
        MatrixXd atom_coords;
        // Store the molecule number of each atom
        std::vector<int> atom_mols;
        // Store the types of each atom
        std::vector<int> atom_types;
        // Store the charge of each atom
        std::vector<double> atom_charges;
        // Store the mass of each atom type
        std::map<int, double> atom_masses;
        // Store the bonds of atoms by using a bond table
        std::vector<std::vector<int>> bond_table;
        // Store atoms of each chains
        std::vector<std::vector<int>> chains;
        // Simulation box type
        bool triclinic = false;
    private:
        // Find each molecule chain and each atoms
        void _find_molecules();
           
        // Input settings
        const Inputs &_inputs;
};

