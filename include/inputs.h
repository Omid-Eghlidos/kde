#pragma once
#include <string>
#include <vector>
#include <map>


struct BeadType {
    // Bead type
    int type;
    // Atom types that make up the bead
    std::vector<int> atom_types;
    // CG weights for each atom inside the bead
    std::vector<double> weights;
};


// Class to store the input parameters
class Inputs {
    public:
        Inputs(std::string path);
        // Inputs/outputs settings
        // LAMMPS format data file
        std::string lmp_data;
        // LAMMPS format dump file
        std::string lmp_dump;
        // Output files tag
        std::string tag = "CG";
        // Type definition settings
        // System phase
        bool amorphous = true;
        // Input file atom types definition
        std::map<std::string, int> atom_types;
        // Bead types definition
        std::map<char, BeadType> bead_types;
        // Distribution settings
        // Process method (using <thread> library)
        std::string process = "parallel";
        // RDF 4D vector settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> rdf;
        // Special bonds 3D vector : [1st neighbor] [second neighbor] [third neighbor]
        std::vector<int> special_bonds;
        // BDF 4D vector settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> bdf;
        // ADF 4D vector settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> adf;
        // TDF 4D vector settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> tdf;
        // IDF 4D vector settings : [min] [max] [steps] [kde_bandwidth]
        std::vector<double> idf;
};

