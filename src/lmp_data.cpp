#include "lmp_data.h"
#include "string_tools.h"
#include <fmt/format.h>
#include <iostream>


LMP_data::LMP_data(const Inputs &inputs) : _inputs(inputs) {
    fmt::print("\n#################### Initializing #######################\n");
    std::fstream fid;
    fid.open(_inputs.lmp_data, std::fstream::in);
    if (!fid) {
        std::cout << "Data file does not exist.\n";
        exit(0);
    }
    // Triclinic box (covers the cubical box as well)
    while (fid) {
        auto args = split(before(read_line(fid), "#"));
        if (args.empty()) {
            continue;
        }
        else if (args.size() == 2 && args.back() == "atoms") {
            int n = from_string<int>(args[0]);
            atom_coords = MatrixXd::Zero(n, 3);
            atom_mols.assign(n, 0);
            atom_types.assign(n, 0);
            atom_charges.assign(n, 0.0);
            bond_table.assign(n, {});
        }
        else if (args.size() == 4 && args[3] == "xhi") {
            box(0,0) = str2dbl(args[1]) - str2dbl(args[0]);
        }
        else if (args.size() == 4 && args[3] == "yhi") {
            box(1,1) = str2dbl(args[1]) - str2dbl(args[0]);
        }
        else if (args.size() == 4 && args[3] == "zhi") {
            box(2,2) = str2dbl(args[1]) - str2dbl(args[0]);
        }
        else if (args.size() == 6 && args[5] == "yz") {
            triclinic = true;
            box(0,1) = str2dbl(args[0]);
            box(0,2) = str2dbl(args[1]);
            box(1,2) = str2dbl(args[2]);
        }
        else if (args.size() == 1 && args[0] == "Masses") {
            while (fid) {
                auto pos = fid.tellp();
                args = split(before(read_line(fid), "#"));
                if (args.size() == 1) {
                    fid.seekp(pos);
                    break;
                }
                else if (args.size() == 2) {
                    auto t = from_string<int>(args[0]) - 1;
                    atom_masses[t] = from_string<double>(args[1]);
                }
            }
        }
        else if (args.size() == 1 && args[0] == "Atoms") {
            while (fid) {
                auto pos = fid.tellp();
                args = split(before(read_line(fid), "#"));
                if (args.size() == 1) {
                    fid.seekp(pos);
                    break;
                }
                if (args.size() > 6) {
                    auto i = from_string<int>(args[0]) - 1;
                    atom_mols[i] = from_string<int>(args[1]);
                    atom_types[i] = from_string<int>(args[2]);
                    atom_charges[i] = from_string<double>(args[3]);
                    atom_coords(i,0) = from_string<double>(args[4]);
                    atom_coords(i,1) = from_string<double>(args[5]);
                    atom_coords(i,2) = from_string<double>(args[6]);
                    Vector3d xs = (box.inverse() * atom_coords.row(i).transpose());
                    atom_coords.row(i) -= (box * xs.array().floor().matrix());
                }
            }
        }
        else if (args.size() == 1 && args[0] == "Bonds") {
            while (fid) {
                auto pos = fid.tellp();
                args = split(before(read_line(fid), "#"));
                if (args.size() == 1) {
                    fid.seekp(pos);
                    break;
                }
                if (args.size() == 4) {
                    int i = from_string<int>(args[2]) - 1;
                    int j = from_string<int>(args[3]) - 1;
                    bond_table[i].push_back(j);
                    bond_table[j].push_back(i);
                }
            }
        }
    }
    for (auto &b : bond_table) {
        std::sort(b.begin(), b.end());
    }
    fmt::print("Read {} atoms from {}.\n", atom_coords.rows(), _inputs.lmp_data);
    _find_molecules();
}


// Uses bond connectivity to group atoms into molecules.
void LMP_data::_find_molecules() {
    std::vector<bool> done(atom_types.size(), false);
    for (int i=0; i<(int)atom_types.size(); ++i) {
        if (done[i]) continue;
        chains.push_back({i});
        done[i] = true;
        std::vector<int> queue = bond_table[i];
        while (!queue.empty()) {
            auto j = queue.back();
            queue.pop_back();
            if (done[j]) continue;
            chains.back().push_back(j);
            done[j] = true;
            for (auto k: bond_table[j]) {
                queue.push_back(k);
            }
        }
    }
    for (auto &c: chains) {
        std::sort(c.begin(), c.end());
    }
    fmt::print("Found {} chains in this system.\n", chains.size());
}

