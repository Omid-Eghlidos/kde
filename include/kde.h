#pragma once
#include <fmt/format.h>
#include <Eigen/Dense>
#include <iostream>
#include <stdio.h>
#include <map>
#include "inputs.h"
#include "lmp_data.h"
#include "lmp_dump.h"
#include "cg_data.h"
#include "cg_dump.h"

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::ArrayXd;

//! Radial Distribution Function (RDF) mapped by bead type pairs
using RDF = std::map<std::tuple<char, char>, ArrayXd>;
//! Bond Distribution Function (BDF) mapped by bead type pairs
using BDF = std::map<std::tuple<char, char>, ArrayXd>;
//! Angle Distribution Function (ADF) mapped by bead type triples
using ADF = std::map<std::tuple<char, char, char>, ArrayXd>;
//! Torsion Distribution Function (TDF) mapped by bead type quadruples
using TDF = std::map<std::tuple<char, char, char, char>, ArrayXd>;
//! Improper Distribution Function (IDF) mapped by bead type quadruples
using IDF = std::map<std::tuple<char, char, char, char>, ArrayXd>;


/**
 * @class KDE
 * @brief Calculates various distribution functions of the system using
 * Kernel Density Estimation (KDE).
 *
 * The KDE class computes distributions such as RDF, BDF, ADF, TDF, and IDF
 * using Kernel Density Estimation based on the method described by McCabe in 2014.
 */
class KDE {
    public:
        /**
         * @brief Constructs a KDE object to compute distributions.
         *
         * @param inputs The input parameters for the system.
         * @param sys The LAMMPS data representing the system.
         * @param dump The LAMMPS dump containing atomic coordinates.
         * @param cg_data The coarse-grained data of the system.
         * @param cg_dump The coarse-grained dump data.
         */
        KDE(const Inputs &inputs, const LMP_data &sys, const LMP_dump &dump,
            const CG_data &cg_data, const CG_dump &cg_dump);

        /**
         * @brief Computes the radial distribution function (RDF) for all
         * bond types at a specific timestep using KDE.
         *
         * @param timestep The timestep at which to compute the RDF.
         */
        void calculate_rdf(int timestep);

        /**
         * @brief Computes the bond distribution function (BDF) for all
         * bond types at a specific timestep using KDE.
         *
         * @param timestep The timestep at which to compute the BDF.
         */
        void calculate_bdf(int timestep);

        /**
         * @brief Computes the angle distribution function (ADF) for all
         * bond types at a specific timestep using KDE.
         *
         * @param timestep The timestep at which to compute the ADF.
         */
        void calculate_adf(int timestep);

        /**
         * @brief Computes the torsion distribution function (TDF) for all
         * bond types at a specific timestep using KDE.
         *
         * @param timestep The timestep at which to compute the TDF.
         */
        void calculate_tdf(int timestep);

        /**
         * @brief Computes the improper distribution function (IDF) for all
         * bond types at a specific timestep using KDE.
         *
         * @param timestep The timestep at which to compute the IDF.
         */
        void calculate_idf(int timestep);

        //! Determines what distributions to compute
        std::map<std::string, bool> compute;
        //! Stores radial distances at which to compute RDF
        ArrayXd radius;
        //! Stores distances at which to compute BDF
        ArrayXd lengths;
        //! Stores angles at which to compute ADF
        ArrayXd angles;
        //! Stores torsion angles at which to compute TDF
        ArrayXd torsions;
        //! Stores improper angles at which to compute IDF
        ArrayXd impropers;
        //! Stores RDF for each timestep
        std::vector<RDF> rdfs;
        //! Stores BDF for each timestep
        std::vector<BDF> bdfs;
        //! Stores ADF for each timestep
        std::vector<ADF> adfs;
        //! Stores TDF for each timestep
        std::vector<TDF> tdfs;
        //! Stores IDF for each timestep
        std::vector<IDF> idfs;
    private:
        /**
         * @brief Determines if atoms i and j are bonded together.
         *
         * @param i Index of the first atom.
         * @param j Index of the second atom.
         * @return true If atoms i and j are first neighbors (bonded).
         * @return false If atoms i and j are not first neighbors.
         */
        bool _are_first_neighbors(int i , int j);

        /**
         * @brief Determines if atoms i and j are bonded to a shared atom (second neighbors).
         *
         * @param i Index of the first atom.
         * @param j Index of the second atom.
         * @return true If atoms i and j are second neighbors.
         * @return false If atoms i and j are not second neighbors.
         */
        bool _are_second_neighbors(int i, int j);

        /**
         * @brief Determines if atoms i and j are three bonds away from each
         * other (third neighbors).
         *
         * @param i Index of the first atom.
         * @param j Index of the second atom.
         * @return true If atoms i and j are third neighbors.
         * @return false If atoms i and j are not third neighbors.
         */
        bool _are_third_neighbors(int i , int j);

        /**
         * @brief Divides a range into equidistant segments.
         *
         * @param start The starting value of the range.
         * @param end The ending value of the range.
         * @param num_step The number of segments to divide the range into.
         * @return ArrayXd An array of equidistant segments within the range.
         */
        ArrayXd _linspace(double start, double end, int num_step);

        /**
         * @brief Finds the bond vector unwrapped from the periodic simulation box.
         *
         * @param timestep The timestep at which to compute the bond vector.
         * @param i Index of the first atom.
         * @param j Index of the second atom.
         * @return Vector3d The unwrapped bond vector between atoms i and j.
         */
        Vector3d _unwrapped_bond_vector(int timestep, int i, int j);

        /**
         * @brief Calculates the angle between two unit vectors.
         *
         * @param v1 The first unit vector.
         * @param v2 The second unit vector.
         * @return double The angle between the two vectors in radians.
         */
        double _theta_angle(const VectorXd &v1, const VectorXd &v2);

        /**
         * @brief Determines the bond type between atoms i and j.
         *
         * @param i Index of the first atom.
         * @param j Index of the second atom.
         * @return std::tuple<char, char> The bond type as a tuple of two characters.
         */
        std::tuple<char, char> _determine_bond_type(int i, int j);

        /**
         * @brief Determines the angle type between atoms i, j, and k.
         *
         * @param i Index of the first atom.
         * @param j Index of the second atom.
         * @param k Index of the third atom.
         * @return std::tuple<char, char, char> The angle type as a tuple of three characters.
         */
        std::tuple<char, char, char> _determine_angle_type(int i, int j, int k);

        /**
         * @brief Determines the torsion angle type between atoms i, j, k, and l.
         *
         * @param i Index of the first atom.
         * @param j Index of the second atom.
         * @param k Index of the third atom.
         * @param l Index of the fourth atom.
         * @return std::tuple<char, char, char, char> The torsion type as
         * a tuple of four characters.
         */
        std::tuple<char, char, char, char> _determine_torsion_type(int i, int j, int k, int l);

        /**
         * @brief Determines the improper angle type between atoms i, j, k, and l.
         *
         * @param i Index of the first atom.
         * @param j Index of the second atom.
         * @param k Index of the third atom.
         * @param l Index of the fourth atom.
         * @return std::tuple<char, char, char, char> The improper type as a
         * tuple of four characters.
         */
        std::tuple<char, char, char, char> _determine_improper_type(int i, int j, int k,int l);

        /**
         * @brief Merges local results from parallel processing into global results.
         *
         * This function consolidates local computations into a global to
         * avoid race conditions in parallel processing.
         * It updates the global distribution map and the global count map
         * with data from local threads.
         *
         * @tparam MapType The type of the distribution map (e.g., RDF, BDF, ADF, TDF, IDF).
         * @tparam CountMapType The type of the count map used for keeping track of types.
         * @param global_map The global map where results from all threads are accumulated.
         * @param local_map The local map containing results computed by a single thread.
         * @param global_count_map The global count map where counts from all
         * threads are accumulated.
         * @param local_count_map The local count map containing counts computed
         * by a single thread.
         */
        template <typename DF, typename TypeCount>
        void _merge_results(DF &df, const DF &df_, TypeCount &type_count,
                                             const TypeCount &type_count_);

        //! Reference to the input parameters
        const Inputs &_inputs;
        //! Reference to the LAMMPS data
        const LMP_data &_data;
        //! Reference to the LAMMPS dump data
        const LMP_dump &_dump;
        //! Reference to the coarse-grained data
        const CG_data &_cg_data;
        //! Reference to the coarse-grained trajectory/dump data
        const CG_dump &_cg_dump;
};

