#pragma once
#include <fmt/format.h>
#include <Eigen/Dense>
#include <iostream>
#include <stdio.h>
#include <map>
#include "inputs.h"
#include "lammps_data.h"
#include "lammps_dump.h"
#include "cg_data.h"
#include "cg_dump.h"


using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using RDF = std::map<std::tuple<char, char>, ArrayXd>;
using BDF = std::map<std::tuple<char, char>, ArrayXd>;
using ADF = std::map<std::tuple<char, char, char>, ArrayXd>;
using TDF = std::map<std::tuple<char, char, char, char>, ArrayXd>;
using IDF = std::map<std::tuple<char, char, char, char>, ArrayXd>;


// A class to calculate the distributions of the system using Kernell Density
// Estimation (KDE) - 2014 McCabe
class KDE {
    public:
        KDE(const Inputs &inputs, const LMP_data &sys, const LMP_dump &dump,
            const CG_data &cg_data, const CG_dump &cg_dump);
        // Calculates the distribution functions
        void calculate_rdf(int timestep);
        void calculate_bdf(int timestep);
        void calculate_adf(int timestep);
        void calculate_tdf(int timestep);
        void calculate_idf(int timestep);

        // Variables
        std::map<std::string, bool> compute;
        // Average distributions over all the timesteps
        ArrayXd radius;
        ArrayXd lengths;
        ArrayXd angles;
        ArrayXd torsions;
        ArrayXd impropers;
        // Store the distributions for each timestep
        std::vector<RDF> rdfs;
        std::vector<BDF> bdfs;
        std::vector<ADF> adfs;
        std::vector<TDF> tdfs;
        std::vector<IDF> idfs;
    private:
        // Helper functions
        // Check if the atoms are 1st, 2nd, or 3rd neighbors of each other
        bool _are_first_neighbors(int i , int j);
        bool _are_second_neighbors(int i, int j);
        bool _are_third_neighbors(int i , int j);
        // Divide a range into a equidistance segments
        ArrayXd _linspace(double start, double end, int num_step);
        // Find the bond vector unwrapped from the periodic simulation box
        Vector3d _unwrapped_bond_vector(int timestep, int i, int j);
        // Calculates the angle between two unit vectors
        double _theta_angle(const VectorXd &v1, const VectorXd &v2);

        // Variables
        const Inputs &_inputs;
        const LMP_data &_data;
        const LMP_dump &_dump;
        const CG_data &_cg_data;
        const CG_dump &_cg_dump;
};

