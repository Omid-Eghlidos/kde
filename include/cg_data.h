#pragma once
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "inputs.h"
#include "lmp_data.h"
#include "lmp_dump.h"

using Eigen::MatrixXd;
using Eigen::Vector3d;

typedef std::vector<std::vector<MatrixXd>> BlockAtoms;


/**
 * @struct Bead
 * @brief Stores parameters of the beads in the coarse-grained system.
 *
 * The Bead structure holds information about each bead, including its ID,
 * type, chain number, charge, and the list of atoms it contains.
 */
struct Bead {
    //! Index of the bead.
    int id;
    //! Chain number that the bead belongs to.
    int chain;
    //! Type of the bead.
    int type;
    //! Bead's alphabetic name.
    char name;
    //! Charge of the bead as sum of all the atom charges inside of it.
    double charge = 0.0;
    //! List of atoms inside the bead.
    std::vector<int> atoms;
    //! Bead mapping weights.
    std::vector<double> atom_weights;
    //! Coordinates of the bead (can be unwrapped, scaled, or unscaled).
    Vector3d coords = Vector3d::Zero(3);
};


/**
 * @class CG_data
 * @brief Manages the coarse-grained (CG) data for the simulation.
 *
 * The CG_data class is responsible for handling the specifications of the
 * coarse-grained system, including bead types, bonds, angles, torsions, and impropers.
 */
class CG_data {
public:
    /**
     * @brief Constructs a CG_data object from inputs and LAMMPS data.
     *
     * @param inputs The input parameters for the CG system.
     * @param data The LAMMPS data used to define the CG system.
     */
    CG_data(const Inputs &inputs, const LMP_data &data);

    /**
     * @brief Determines the bead type considering any existing symmetry in the mapping.
     *
     * @param type The original bead type.
     * @return The mapped bead type.
     */
    int get_bead_type(int type) const;

    /**
     * @brief Finds the bond type for two beads among the existing bond types.
     *
     * @param i Index of the first bead.
     * @param j Index of the second bead.
     * @return The bond type.
     */
    int get_bond_type(int i, int j) const;

    /**
     * @brief Finds the angle type for three beads among the existing angle types.
     *
     * @param i Index of the first bead.
     * @param j Index of the second bead.
     * @param k Index of the third bead.
     * @return The angle type.
     */
    int get_angle_type(int i, int j, int k) const;

    /**
     * @brief Finds the torsion type for four beads among the existing torsion types.
     *
     * @param i Index of the first bead.
     * @param j Index of the second bead.
     * @param k Index of the third bead.
     * @param l Index of the fourth bead.
     * @return The torsion type.
     */
    int get_torsion_type(int i, int j, int k, int l) const;

    /**
     * @brief Finds the improper type for four beads among the existing improper types.
     *
     * @param i Index of the first bead.
     * @param j Index of the second bead.
     * @param k Index of the third bead.
     * @param l Index of the fourth bead.
     * @return The improper type.
     */
    int get_improper_type(int i, int j, int k, int l) const;

    //! Beads defined according to the mapping matrix.
    std::vector<Bead> beads;
    //! Store the unique bead types.
    std::vector<int> bead_types;
    //! Store the bead bonds.
    std::vector<std::vector<int>> bead_bonds;
    //! Store the unique bead bond types.
    std::vector<std::tuple<int, int>> bead_bond_types;
    //! Store the bead angles.
    std::vector<std::tuple<int, int, int>> bead_angles;
    //! Store the unique bead angle types.
    std::vector<std::tuple<int, int, int>> bead_angle_types;
    //! Store the bead torsions.
    std::vector<std::tuple<int, int, int, int>> bead_torsions;
    //! Store the unique bead torsion types.
    std::vector<std::tuple<int, int, int, int>> bead_torsion_types;
    //! Store the bead impropers.
    std::vector<std::tuple<int, int, int, int>> bead_impropers;
    //! Store the unique bead improper types.
    std::vector<std::tuple<int, int, int, int>> bead_improper_types;
private:
    /**
     * @brief Makes beads based on the input definition.
     */
    void _make_beads();

    /**
     * @brief Finds beads for each defined type.
     *
     * @param bead_name The name of the bead.
     * @param bead The bead type structure.
     */
    void _find_beads_of_type(char bead_name, BeadType bead);

    /**
     * @brief Finds the bead coordinates for each bead.
     *
     * @param atoms List of atoms in the bead.
     * @param weights Weights of the atoms in the bead.
     * @return Vector3d The calculated coordinates of the bead.
     */
    Vector3d _find_beads_coords(std::vector<int> atoms, std::vector<double> weights);

    /**
     * @brief Identifies bonds between the beads.
     */
    void _identify_beads_bonds();

    /**
     * @brief Determines unique bond types for bead bonds.
     *
     * @param i Index of the first bead.
     * @param j Index of the second bead.
     */
    void _bond_types(int i, int j);

    /**
     * @brief Determines if two beads are bonded together.
     *
     * @param bead1_atoms List of atoms in the first bead.
     * @param bead2_atoms List of atoms in the second bead.
     * @return true If the beads are bonded.
     * @return false If the beads are not bonded.
     */
    bool _are_bonded(std::vector<int> bead1_atoms, std::vector<int> bead2_atoms);

    /**
     * @brief Identifies angles between the beads.
     */
    void _identify_beads_angles();

    /**
     * @brief Determines unique angle types for bead angles.
     *
     * @param i Index of the first bead.
     * @param j Index of the second bead.
     * @param k Index of the third bead.
     */
    void _angle_types(int i, int j, int k);

    /**
     * @brief Identifies torsions between the beads.
     */
    void _identify_beads_torsions();

    /**
     * @brief Determines unique torsion types for bead torsions.
     *
     * @param i Index of the first bead.
     * @param j Index of the second bead.
     * @param k Index of the third bead.
     * @param l Index of the fourth bead.
     */
    void _torsion_types(int i, int j, int k, int l);

    /**
     * @brief Identifies impropers between the beads.
     */
    void _identify_beads_impropers();

    /**
     * @brief Determines unique improper types for bead torsions.
     *
     * @param i Index of the first bead.
     * @param j Index of the second bead.
     * @param k Index of the third bead.
     * @param l Index of the fourth bead.
     */
    void _improper_types(int i, int j, int k, int l);

    //! Reference to the input parameters.
    const Inputs &_inputs;
    //! Reference to the LAMMPS data.
    const LMP_data &_data;
};

