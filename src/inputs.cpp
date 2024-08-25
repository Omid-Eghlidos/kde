#include "inputs.h"
#include "string_tools.h"
#include <iostream>


Inputs::Inputs(std::string path) {
    std::fstream fid(path);
    while (fid) {
        auto row = split(read_line(fid));
        if (row.empty()) continue;
        // Inputs/ouputs setting
        if (row[0] == "lammps_data") {
            lmp_data = row[1];
        }
        else if (row[0] == "lammps_dump") {
            lmp_dump = row[1];
        }
        else if (row[0] == "tag") {
            tag = row[1];
        }
        // Type definitions
        else if (row[0] == "phase" && row[1] == "crystal") {
            amorphous = false;
        }
        else if (row[0] == "type") {
            atom_types[row[2]] = str2u32(row[1]);
        }
        else if (row[0] == "bead") {
            BeadType bead;
            for (size_t i=2; i<row.size(); i++) {
                bead.atom_types.push_back(atom_types.find(row[i])->second);
            }
            bead.type = bead_types.size();
            bead_types[row[1][0]] = bead;
        }
        else if (row[0] == "weight") {
            for (size_t i=2; i<row.size(); i++) {
                bead_types[row[1][0]].weights.push_back(str2dbl(row[i]));
            }
        }
        // Distribution settings
        else if (row[0] == "process") {
            process = row[1];
        }
        else if (row[0] == "rdf") {
            for (size_t i=1; i<row.size(); i++) {
                rdf.push_back(str2dbl(row[i]));
            }
        }
        else if (row[0] == "special_bonds") {
            for (size_t i=1; i<row.size(); i++) {
                special_bonds.push_back(str2u32(row[i]));
            }
        }
        else if (row[0] == "bdf") {
            for (size_t i=1; i<row.size(); i++) {
                bdf.push_back(str2dbl(row[i]));
            }
        }
        else if (row[0] == "adf") {
            for (size_t i=1; i<row.size(); i++) {
                adf.push_back(str2dbl(row[i]));
            }
        }
        else if (row[0] == "tdf") {
            for (size_t i=1; i<row.size(); i++) {
                tdf.push_back(str2dbl(row[i]));
            }
        }
        else if (row[0] == "idf") {
            for (size_t i=1; i<row.size(); i++) {
                idf.push_back(str2dbl(row[i]));
            }
        }
    }
}

