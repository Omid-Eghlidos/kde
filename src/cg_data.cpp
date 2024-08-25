#include "cg_data.h"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tuple>
#include <iostream>
#include <fstream>


// Class to store specifications of the coarse grained system
CG_data::CG_data(const Inputs &i, const LMP_data &d): _inputs(i), _data(d) {
    fmt::print("\n################## Coarse Graining ######################\n");
    _make_beads();
    _identify_beads_bonds();
    _identify_beads_angles();
    _identify_beads_torsions();
    _identify_beads_impropers();
}


void CG_data::_make_beads() {
    // Find the beads in each system for each defined bead type
    for (auto b : _inputs.bead_types) {
        _find_beads_of_type(b.first, b.second);
    }
}


void CG_data::_find_beads_of_type(char bead_name, BeadType bead) {
    fmt::print("Beads of type '{}': ", bead_name);
    // Go through all the atoms and find the ones that matches the bead definition
    int bead_type_count = 0;
    std::vector<bool> visited(_data.atom_types.size(), false);
    for (size_t i=0; i < _data.chains.size(); i++) {
        for (size_t j=0; j < _data.chains[i].size(); j++) {
            std::vector<int> potential_bead(bead.atom_types.size());
            // Continue if already visited this atom or the type of this atom
            // is not equal to the first one in the bead
            int a = _data.chains[i][j];
            if (visited[a] || (bead.atom_types[0] != _data.atom_types[a])) continue;
            potential_bead[0] = a;
            for (size_t k=1; k < bead.atom_types.size(); k++) {
                for (size_t l=0; l < _data.bond_table[a].size(); l++) {
                    int b = _data.bond_table[a][l];
                    if (!visited[b] && bead.atom_types[k] == _data.atom_types[b]) {
                        potential_bead[k] = b;
                        a = b;
                        visited[b] = true;
                    }
                }
            }
            if (potential_bead.size() == bead.atom_types.size()) {
                Bead bead_;
                bead_.id = beads.size();
                bead_.chain = i;
                bead_.type = bead.type;
                bead_.name = bead_name;
                bead_.atoms = potential_bead;
                for (auto b : potential_bead) {
                    bead_.charge += _data.atom_charges[b];
                }
                bead_.atom_weights = bead.weights;
                bead_.coords = _find_beads_coords(potential_bead, bead.weights);
                // Check it does not already exist
                if (!(std::find(bead_types.begin(), bead_types.end()
                                          , bead.type) != bead_types.end())) {
                    bead_types.push_back(bead.type);
                }
                beads.push_back(bead_);
                bead_type_count++;
            }
        }
    }
    fmt::print("Found {} and their coordinates.\n", bead_type_count);
}


Vector3d CG_data::_find_beads_coords(std::vector<int> atoms, std::vector<double> weights) {
    // Map atoms into the periodic cell w.r.t first atom of the bead
    Vector3d x0 = _data.atom_coords.row(atoms[0]);
    Vector3d coords = x0 * weights[0];
    for (size_t j = 1; j < atoms.size(); j++) {
        Vector3d x = _data.atom_coords.row(atoms[j]);
        Vector3d dx = x - x0;
        Vector3d ds = (_data.box.inverse() * dx).array().round().matrix();
        dx -= _data.box * ds;
        coords += (x0 + dx) * weights[j];
    }
    return coords;
}


// Two beads are bonded if any atoms in one is bonded to any atom in the other
void CG_data::_identify_beads_bonds() {
    bead_bonds.assign(beads.size(), {});
    for (auto b1 : beads) {
        for (auto b2 : beads) {
            if ((b1.chain == b2.chain) && (b1.id != b2.id)) {
                if (_are_bonded(b1.atoms, b2.atoms)) {
                    bead_bonds[b1.id].push_back(b2.id);
                    _bond_types(b1.type, b2.type);
                }
            }
        }
    }
    int num_bonds = 0;
    for (size_t i = 0; i < bead_bonds.size(); i++) {
        for (size_t j = 0; j < bead_bonds[i].size(); j++) {
            if ((int)i < bead_bonds[i][j]) {
                num_bonds++;
            }
        }
    }
    fmt::print("Found {} bonds and {} bond types.\n"
              , num_bonds, bead_bond_types.size());
}


bool CG_data::_are_bonded(std::vector<int> b1_atoms, std::vector<int> b2_atoms) {
    for (auto atom1 : b1_atoms) {
        for (auto bonded_to_atom1 : _data.bond_table[atom1]) {
            for (auto atom2 : b2_atoms) {
                if (bonded_to_atom1 == atom2) {
                    return true;
                }
            }
        }
    }
    return false;
}


void CG_data::_bond_types(int i, int j) {
    std::tuple<int, int> type;
    if (i < j) {
        type = std::make_tuple(i, j);
    }
    else {
        type = std::make_tuple(j, i);
    }
    // Check it does not already exist
    if (!(std::find(bead_bond_types.begin(), bead_bond_types.end(), type) !=
                                                      bead_bond_types.end())) {
        bead_bond_types.push_back(type);
    }
}


int CG_data::get_bond_type(int i, int j) const {
    size_t bbt;
    std::tuple<int, int> type;
    if (beads[i].type < beads[j].type) {
        type = std::make_tuple(beads[i].type, beads[j].type);
    }
    else {
        type = std::make_tuple(beads[j].type, beads[i].type);
    }
    for (bbt=0; bbt<bead_bond_types.size(); bbt++) {
        if (type == bead_bond_types[bbt]) {
            break;
        }
    }
    return bbt;
}


// Three linearly bonded beads (no loop) make an angle
void CG_data::_identify_beads_angles() {
    for (auto b : beads) {
        for (auto i : bead_bonds[b.id]) {
            for (auto j : bead_bonds[b.id]) {
                if (i < j) {
                    bead_angles.push_back(std::make_tuple(i, b.id, j));
                    _angle_types(i, b.id, j);
                }
            }
        }
    }
    fmt::print("Found {} bond angles and {} angle types.\n"
              , bead_angles.size(), bead_angle_types.size());
}


void CG_data::_angle_types(int i, int j, int k) {
    std::tuple<int, int, int> type;
    if (beads[i].type < beads[k].type) {
        type = std::make_tuple(beads[i].type, beads[j].type, beads[k].type);
    }
    else if (beads[i].type > beads[k].type) {
        type = std::make_tuple(beads[k].type, beads[j].type, beads[i].type);
    }
    // Check it does not already exist
    if (!(std::find(bead_angle_types.begin(), bead_angle_types.end(), type) !=
                                                      bead_angle_types.end())) {
        bead_angle_types.push_back(type);
    }
}


int CG_data::get_angle_type(int i, int j, int k) const {
    int at = 0;
    for (auto bat : bead_angle_types) {
        at++;
        if (std::make_tuple(beads[i].type, beads[j].type, beads[k].type) == bat) {
            break;
        }
    }
    return at;
}


// Four linearly bonded beads (no loop) make a torsion/dihedral angle
void CG_data::_identify_beads_torsions() {
    for (int j=0; j<(int)bead_bonds.size(); j++) {
        for (auto i : bead_bonds[j]) {
            for (auto k : bead_bonds[j]) {
                if ((k > j) && (k != i)) {
                    for (auto l : bead_bonds[k]) {
                        if ((l != j) && (l != i)){
                            bead_torsions.push_back(std::make_tuple(i, j, k, l));
                            _torsion_types(i, j, k, l);
                        }
                    }
                }
            }
        }
    }
    fmt::print("Found {} torsion angles and {} torsion angle types.\n"
              , bead_torsions.size(), bead_torsion_types.size());
}


void CG_data::_torsion_types(int i, int j, int k, int l) {
    std::tuple<int, int, int, int> type;
    if (beads[i].type > beads[l].type) {
        type = std::make_tuple(beads[l].type, beads[k].type,
                               beads[j].type, beads[i].type);
    }
    else if ((beads[i].type == beads[l].type) && (beads[j].type > beads[k].type)) {
        type = std::make_tuple(beads[l].type, beads[k].type,
                               beads[j].type, beads[i].type);
    }
    else {
        type = std::make_tuple(beads[i].type, beads[j].type,
                               beads[k].type, beads[l].type);
    }
    // Check it does not already exist
    if (!(std::find(bead_torsion_types.begin(), bead_torsion_types.end(), type)
                                                != bead_torsion_types.end())) {
        bead_torsion_types.push_back(type);
    }
}


int CG_data::get_torsion_type(int i, int j, int k, int l) const {
    int tt = 0;
    for (auto btt : bead_torsion_types) {
        tt++;
        if (std::make_tuple(beads[i].type, beads[j].type
                          , beads[k].type, beads[l].type) == btt) {
            break;
        }
    }
    return tt;
}


// Three beads connected to a bead with three of them in a plane and one in another
void CG_data::_identify_beads_impropers() {
    for (int j=0; j<(int)bead_bonds.size(); j++) {
        // Only beads with three bonds or more can form improper angle
        if (bead_bonds[j].size() > 2) {
            for (auto i : bead_bonds[j]) {
                for (auto k : bead_bonds[j]) {
                    for (auto l : bead_bonds[j]) {
                        if (i < k && k < l){
                            bead_impropers.push_back(std::make_tuple(i, j, k, l));
                            _improper_types(i, j, k, l);
                        }
                    }
                }
            }
        }
    }
    fmt::print("Found {} improper angles and {} improper angle types.\n"
              , bead_impropers.size(), bead_improper_types.size());
}


void CG_data::_improper_types(int i, int j, int k, int l) {
    std::tuple<int, int, int, int> type;
    // TODO: Needs to be double checked
    if (beads[i].type < beads[l].type) {
        type = std::make_tuple(beads[i].type, beads[j].type,
                               beads[k].type, beads[l].type);
    }
    if (beads[i].type == beads[j].type && beads[j].type == beads[k].type &&
        beads[j].type < beads[k].type) {
        type = std::make_tuple(beads[i].type, beads[j].type,
                               beads[k].type, beads[l].type);
    }
    // Check it does not already exist
    if (!(std::find(bead_improper_types.begin(), bead_improper_types.end(), type)
                                            != bead_improper_types.end())) {
        bead_improper_types.push_back(type);
    }
}


int CG_data::get_improper_type(int i, int j, int k, int l) const {
    int it = 0;
    for (auto bit : bead_improper_types) {
        it++;
        if (std::make_tuple(beads[i].type, beads[j].type
                          , beads[k].type, beads[l].type) == bit) {
            break;
        }
    }
    return it;
}

