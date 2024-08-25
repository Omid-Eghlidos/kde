#include "distributions.h"
#include <cmath>


KDE::KDE(const Inputs &inputs, const LMP_data &data, const LMP_dump &dump,
         const CG_data &cg_data, const CG_dump &cg_dump) : _inputs(inputs),
         _data(data), _dump(dump), _cg_data(cg_data), _cg_dump(cg_dump) {
    // Radius
    compute["rdf"] = false;
    if (!_inputs.rdf.empty()) {
        compute["rdf"] = true;
        radius = _linspace(_inputs.rdf[0], _inputs.rdf[1], _inputs.rdf[2]);
    }
    // Bond lengths
    compute["bdf"] = false;
    if (!_inputs.bdf.empty()) {
        compute["bdf"] = true;
        lengths = _linspace(_inputs.bdf[0], _inputs.bdf[1], _inputs.bdf[2]);
    }
    // Bond angles
    compute["adf"] = false;
    if (!_inputs.adf.empty()) {
        compute["adf"] = true;
        angles = _linspace(_inputs.adf[0], _inputs.adf[1], _inputs.adf[2]);
    }
    // Torsion angles
    compute["tdf"] = false;
    if (!_inputs.tdf.empty()) {
        compute["tdf"] = true;
        torsions = _linspace(_inputs.tdf[0], _inputs.tdf[1], _inputs.tdf[2]);
    }
    // Improper angles
    compute["idf"] = false;
    if (!_inputs.idf.empty()) {
        compute["idf"] = true;
        impropers = _linspace(_inputs.idf[0], _inputs.idf[1], _inputs.idf[2]);
    }
}


// Computes the radial distribution function across all bead types within the
// timestep using KDE.
void KDE::calculate_rdf(int timestep) {
    RDF rdf;
    // KDE Bandwidth
    auto w = _inputs.rdf[3];
    // Bead type counts
    std::map<char, int> bead_type_count;
    // Find different pair types and their distribution
    #pragma omp parallel for
    for (size_t i=0; i<_cg_data.beads.size(); i++) {
        char type_i = _cg_data.beads[i].name;
        if (!(bead_type_count.find(type_i) != bead_type_count.end())) {
            bead_type_count[type_i] = 0;
        }
        bead_type_count[type_i] += 1;
        for (size_t j=i+1; j<_cg_data.beads.size(); j++) {
            if ((_inputs.special_bonds[0] == 0) && (_are_first_neighbors(i, j))) {
                continue;
            }
            if ((_inputs.special_bonds[1] == 0) && (_are_second_neighbors(i, j))) {
                continue;
            }
            if ((_inputs.special_bonds[2] == 0) && (_are_third_neighbors(i, j))) {
                continue;
            }
            char type_j = _cg_data.beads[j].name;
            std::tuple<char, char> pair_type;
            // Prevent the reverse type to be counted separately
            if (type_i < type_j) {
                pair_type = std::make_tuple(type_i, type_j);
            }
            else {
                pair_type = std::make_tuple(type_j, type_i);
            }
            // Bond vector between two beads in system coordinates within PBC
            auto r = _unwrapped_bond_vector(timestep, i, j);
            // Initialize the count if the type already does not exist
            if (!(rdf.find(pair_type) != rdf.end())) {
                rdf[pair_type] = ArrayXd::Zero(radius.rows());
            }
            #pragma omp critical
            rdf[pair_type] += exp(-0.5 * ((radius - r.norm()) / w).square());
        }
    }
    auto z = w * sqrt(2.0 * M_PI);
    auto V = _dump.box[timestep].determinant();
    auto dr = radius[1] - radius[0];
    auto vs = 4.0 * M_PI / 3.0 * ((pow(radius + dr, 3) - pow(radius, 3)));
    for (auto pair_type : rdf) {
        int n1 = bead_type_count[std::get<0>(pair_type.first)];
        int n2 = bead_type_count[std::get<1>(pair_type.first)];
        ArrayXd s = V / n1 / n2 / vs / z * dr;
        // Handles symmetry in pair distance r_ij = r_ji when i = j
        if (std::get<0>(pair_type.first) == std::get<1>(pair_type.first)) {
            s *= 2.0;
        }
        rdf[pair_type.first] *= s;
    }
    rdfs.push_back(rdf);
}


// Computes the bond distribution function across all bond types within the
// timestep using KDE.
void KDE::calculate_bdf(int timestep) {
    // If there is no bond return
    if (_cg_data.bead_bonds.size() == 0) return;
    BDF bdf;
    // KDE Bandwidth
    auto w = _inputs.bdf[3];
    // Count and store number of each bond type
    std::map<std::tuple<char, char>, int> bond_type_count;
    // Find different bond types and their distribution
    #pragma omp parallel for
    for (size_t i=0; i<_cg_data.bead_bonds.size(); i++) {
        for (auto j : _cg_data.bead_bonds[i]) {
            char type_i = _cg_data.beads[i].name;
            char type_j = _cg_data.beads[j].name;
            std::tuple<char, char> bond_type;
            // Prevent the reverse type to be counted separately
            if (type_i < type_j) {
                bond_type = std::make_tuple(type_i, type_j);
            }
            else {
                bond_type = std::make_tuple(type_j, type_i);
            }
            // Initialize the count if the type already does not exist
            if (!(bond_type_count.find(bond_type) != bond_type_count.end())) {
                bond_type_count[bond_type] = 0;
            }
            bond_type_count[bond_type] += 1;
            // Bond vector between two beads in system coordinates within PBC
            Vector3d r = _unwrapped_bond_vector(timestep, i, j);
            if (!(bdf.find(bond_type) != bdf.end())) {
                bdf[bond_type] = ArrayXd::Zero(lengths.rows());
            }
            #pragma omp critical
            bdf[bond_type] += exp(-0.5 * ((lengths - r.norm()) / w).square());
        }
    }
    auto z = w * sqrt(2.0 * M_PI);
    for (auto bond_type : bdf) {
        bdf[bond_type.first] /= (z * bond_type_count[bond_type.first]);
    }
    bdfs.push_back(bdf);
}


// Computes the angle distribution function across all bond angle types
// within the timestep using KDE.
void KDE::calculate_adf(int timestep) {
    // If there is no angle return
    if (_cg_data.bead_angles.size() == 0) return;
    ADF adf;
    // KDE Bandwidth
    auto w = _inputs.adf[3];
    // Count and store number of each angle type
    std::map<std::tuple<char, char, char>, int> angle_type_count;
    // Find different angle types and their distribution
    #pragma omp parallel for
    for (auto a_type : _cg_data.bead_angles) {
        int i = std::get<0>(a_type);
        int j = std::get<1>(a_type);
        int k = std::get<2>(a_type);
        char type_i = _cg_data.beads[i].name;
        char type_j = _cg_data.beads[j].name;
        char type_k = _cg_data.beads[k].name;
        std::tuple<char, char, char> angle_type;
        // Prevent the reverse type to be counted separately
        if (type_i <= type_k) {
            angle_type = std::make_tuple(type_i, type_j, type_k);
        }
        else {
            angle_type = std::make_tuple(type_k, type_j, type_i);
        }
        // Initialize the count if the type already does not exist
        if (!(angle_type_count.find(angle_type) != angle_type_count.end())) {
            angle_type_count[angle_type] = 0;
        }
        angle_type_count[angle_type] += 1;
        // Bond vector between two beads in system coordinates within PBC
        auto r1 = _unwrapped_bond_vector(timestep, j, i);
        auto r2 = _unwrapped_bond_vector(timestep, j, k);
        auto theta = _theta_angle(r1, r2) * 180.0 / M_PI;
        if (!(adf.find(angle_type) != adf.end())) {
            adf[angle_type] = ArrayXd::Zero(angles.rows());
        }
        #pragma omp critical
        adf[angle_type] += exp(-0.5 * ((angles - theta) / w).square());
    }
    auto z = w * sqrt(2.0 * M_PI);
    for (auto angle_type : adf) {
        adf[angle_type.first] /= (z * angle_type_count[angle_type.first]);
    }
    adfs.push_back(adf);
}


// Computes the torsion distribution function across all torsion angle types
// within the timestep using KDE.
void KDE::calculate_tdf(int timestep) {
    // If there is no torsion return
    if (_cg_data.bead_torsions.size() == 0) return;
    TDF tdf;
    // KDE Bandwidth
    auto w = _inputs.tdf[3];
    // Count and store number of each torsion type
    std::map<std::tuple<char, char, char, char>, int> torsion_type_count;
    // Find different torsion types and their distribution
    #pragma omp parallel for
    for (auto t_type : _cg_data.bead_torsions) {
        int i = std::get<0>(t_type);
        int j = std::get<1>(t_type);
        int k = std::get<2>(t_type);
        int l = std::get<3>(t_type);
        char type_i = _cg_data.beads[i].name;
        char type_j = _cg_data.beads[j].name;
        char type_k = _cg_data.beads[k].name;
        char type_l = _cg_data.beads[l].name;
        std::tuple<char, char, char, char> torsion_type;
        // Prevent the reverse type to be counted separately
        if (type_i <= type_l) {
            torsion_type = std::make_tuple(type_i, type_j, type_k, type_l);
        }
        else {
            torsion_type = std::make_tuple(type_l, type_j, type_k, type_i);
        }
        // Initialize the count if the type already does not exist
        if (!(torsion_type_count.find(torsion_type) != torsion_type_count.end())) {
            torsion_type_count[torsion_type] = 0;
        }
        torsion_type_count[torsion_type] += 1;
        // Bond vector between two beads in system coordinates within PBC
        Vector3d v1 = _unwrapped_bond_vector(timestep, i, j);
        v1 /= v1.norm();
        Vector3d v2 = _unwrapped_bond_vector(timestep, j, k);
        v2 /= v2.norm();
        Vector3d v3 = _unwrapped_bond_vector(timestep, k, l);
        v3 /= v3.norm();
        auto v12 = v1.cross(v2);
        auto v23 = v2.cross(v3);
        auto y = v12.cross(v2).dot(v23);
        auto x = v12.dot(v23);
        auto phi = atan2(y, x) * 180.0 / M_PI;
        ArrayXd xx = torsions - phi;
        for (int i=0; i<xx.rows(); i++) {
            if (xx[i] < -180.0) xx[i] += 360;
            if (xx[i] > 180.0) xx[i] -= 360;
        }
        if (!(tdf.find(torsion_type) != tdf.end())) {
            tdf[torsion_type] = ArrayXd::Zero(torsions.rows());
        }
        #pragma omp critical
        tdf[torsion_type] += exp(-0.5 * (xx / w).square());
    }
    auto z = w * sqrt(2.0 * M_PI);
    for (auto torsion_type : tdf) {
        tdf[torsion_type.first] /= z * torsion_type_count[torsion_type.first];
    }
    tdfs.push_back(tdf);
}


// Computes the improper distribution function (IDF) across all improper types
// within the timestep using KDE.
void KDE::calculate_idf(int timestep) {
    // If there is no improper return
    if (_cg_data.bead_impropers.size() == 0) return;
    IDF idf;
    // KDE Bandwidth
    auto w = _inputs.idf[3];
    // Count and store number of each angle type
    std::map<std::tuple<char, char, char, char>, int> improper_type_count;
    // Find different angle types and their distribution
    #pragma omp parallel for
    for (auto i_type : _cg_data.bead_impropers) {
        int i = std::get<0>(i_type);
        int j = std::get<1>(i_type);
        int k = std::get<2>(i_type);
        int l = std::get<3>(i_type);
        char type_i = _cg_data.beads[i].name;
        char type_j = _cg_data.beads[j].name;
        char type_k = _cg_data.beads[k].name;
        char type_l = _cg_data.beads[l].name;
        std::tuple<char, char, char, char> improper_type;
        // Prevent the reverse type to be counted separately
        if (type_i <= type_l) {
            improper_type = std::make_tuple(type_i, type_j, type_k, type_l);
        }
        else {
            improper_type = std::make_tuple(type_l, type_j, type_k, type_i);
        }
        // Initialize the count if the type already does not exist
        if (!(improper_type_count.find(improper_type) != improper_type_count.end())) {
            improper_type_count[improper_type] = 0;
        }
        improper_type_count[improper_type] += 1;
        // Bond vector between two beads in system coordinates within PBC
        Vector3d v1 = _unwrapped_bond_vector(timestep, j, i);
        v1 /= v1.norm();
        Vector3d v2 = _unwrapped_bond_vector(timestep, j, k);
        v2 /= v2.norm();
        Vector3d v3 = _unwrapped_bond_vector(timestep, j, l);
        v3 /= v3.norm();
        auto v12 = v1.cross(v2);
        auto v23 = v2.cross(v3);
        auto psi = _theta_angle(v12, v23) * 180.0 / M_PI;
        if (!(idf.find(improper_type) != idf.end())) {
            idf[improper_type] = ArrayXd::Zero(impropers.rows());
        }
        #pragma omp critical
        idf[improper_type] += exp(-0.5 * ((impropers - psi) / w).square());
    }
    auto z = w * sqrt(2.0 * M_PI);
    for (auto improper_type : idf) {
        idf[improper_type.first] /= (z * improper_type_count[improper_type.first]);
    }
    idfs.push_back(idf);
}


// Helper functions
ArrayXd KDE::_linspace(double start, double end, int num_step) {
    ArrayXd linspaced = ArrayXd::Zero(num_step);
    if (num_step <= 0) {
        fmt::print("Number of segments cannot be zero or negative!");
        exit(0);
    }
    else if (num_step == 1) {
        linspaced(0) = start;
        linspaced(1) = end;
    }
    else {
        double delta = (end - start) / (num_step - 1);
        for (int i=0; i<num_step; i++) {
            linspaced(i) = start + i * delta;
        }
    }
    return linspaced;
}


bool KDE::_are_first_neighbors(int i, int j) {
    bool first_neighbors = false;
    for (auto b1 : _data.bond_table[i]) {
        if (b1 == j) {
            first_neighbors = true;
            break;
        }
    }
    return first_neighbors;
}


bool KDE::_are_second_neighbors(int i, int j) {
    bool second_neighbors = false;
    for (auto b1 : _data.bond_table[i]) {
        for (auto b2 : _data.bond_table[j]) {
            if (b1 == b2) {
                second_neighbors = true;
                break;
            }
        }
    }
    return second_neighbors;
}


bool KDE::_are_third_neighbors(int i, int j) {
    bool third_neighbors = false;
    for (auto b1 : _data.bond_table[i]) {
        for (auto b2 : _data.bond_table[b1]) {
            for (auto b3 : _data.bond_table[j]) {
                if (b2 == b3) {
                    third_neighbors = true;
                    break;
                }
            }
        }
    }
    return third_neighbors;
}


Vector3d KDE::_unwrapped_bond_vector(int ts, int i, int j) {
    Vector3d x_i = _cg_dump.cg_ts[ts].coords[i];
    Vector3d x_j = _cg_dump.cg_ts[ts].coords[j];
    Vector3d dx = x_j - x_i;
    // Solving linear equation of : box * dx = dr
    Vector3d ds = (_cg_dump.cg_ts[ts].box.inverse() * dx).array().round().matrix();
    return dx - _cg_dump.cg_ts[ts].box * ds;
}


double KDE::_theta_angle(const VectorXd &V1, const VectorXd &V2) {
    return std::acos(V1.dot(V2) / (V1.norm() * V2.norm()));
}

