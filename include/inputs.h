#pragma once
#include <string>
#include <vector>
#include <map>
#include <thread>

/**
 * @struct BeadType
 * @brief Represents a bead type in the coarse-grained system.
 *
 * The BeadType structure defines a bead type by specifying the atom types
 * that make up the bead and the corresponding weights for coarse-graining.
 */
struct BeadType {
    //! Bead type
    int type;
    //! Atom types that make up the bead
    std::vector<int> atom_types;
    //! CG weights for each atom inside the bead
    std::vector<double> weights;
};


/**
 * @class Inputs
 * @brief Manages the input parameters for the coarse-grained simulation.
 *
 * The Inputs class stores all the input parameters required to run the
 * coarse-grained simulation, including file paths, type definitions, and
 * process settings.
 */
class Inputs {
    public:
        /**
         * @brief Constructs an Inputs object by reading parameters from a file.
         *
         * @param path The path to the file containing input parameters.
         */
        Inputs(std::string path);

        //! LAMMPS format data file
        std::string lmp_data;
        //! LAMMPS format dump file
        std::string lmp_dump;
        //! Output files tag (default: "CG")
        std::string tag = "CG";
        //! System phase (true if amorphous, false if crystalline)
        bool amorphous = true;
        //! Input file atom types definition
        std::map<std::string, int> atom_types;
        //! Bead types definition
        std::map<char, BeadType> bead_types;
        //! Enable or disable parallel processing (default: true)
        bool parallel = true;
        //! Number of threads to use for parallel processing (default: hardware concurrency)
        int num_threads = std::thread::hardware_concurrency();
        //! RDF settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> rdf;
        //! Special bonds settings : [1st neighbor] [2nd neighbor] [3rd neighbor]
        std::vector<int> special_bonds;
        //! BDF settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> bdf;
        //! ADF settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> adf;
        //! TDF settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> tdf;
        //! IDF settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> idf;
};

