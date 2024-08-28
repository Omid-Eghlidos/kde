#include "kde.h"
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


void KDE::calculate_rdf(int timestep) {
    RDF rdf;
    // Bead type counts
    std::map<char, int> bead_type_count;
    // KDE Bandwidth
    auto w = _inputs.rdf[3];
    #pragma omp parallel
    {
        // Local rdf to prevent race conditions
        RDF rdf_;
        // Local bead type count to prevent race conditions
        std::map<char, int> bead_type_count_;
        // Find different pair types and their distribution
        #pragma omp for
        for (size_t i = 0; i < _cg_data.beads.size(); i++) {
            char type_i = _cg_data.beads[i].name;
            if (!(bead_type_count_.find(type_i) != bead_type_count_.end())) {
                bead_type_count_[type_i] = 0;
            }
            bead_type_count_[type_i] += 1;
            for (size_t j = i+1; j < _cg_data.beads.size(); j++) {
                if ((_inputs.special_bonds[0] == 0) && (_are_first_neighbors(i, j)))
                    continue;
                if ((_inputs.special_bonds[1] == 0) && (_are_second_neighbors(i, j)))
                    continue;
                if ((_inputs.special_bonds[2] == 0) && (_are_third_neighbors(i, j)))
                    continue;
                // Bond vector between two beads in system coordinates within PBC
                auto r = _unwrapped_bond_vector(timestep, i, j);
                // Determine pair type
                char type_j = _cg_data.beads[j].name;
                std::tuple<char, char> pair_type = std::make_tuple(type_i, type_j);
                if (!(rdf_.find(pair_type) != rdf_.end())) {
                    rdf_[pair_type] = ArrayXd::Zero(radius.rows());
                }
                rdf_[pair_type] += exp(-0.5 * ((radius - r.norm()) / w).square());
            }
        }
        // Combine the local into the global results
        _merge_results(rdf, rdf_, bead_type_count, bead_type_count_);
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
    // Transfer ownership of the data rather than copy it for better memory use
    rdfs.push_back(std::move(rdf));
}


void KDE::calculate_bdf(int timestep) {
    // If there is no bond return
    if (_cg_data.bead_bonds.size() == 0)
        return;
    BDF bdf;
    // Count and store number of each bond type
    std::map<std::tuple<char, char>, int> bond_type_count;
    // KDE Bandwidth
    auto w = _inputs.bdf[3];
    // Find different bond types and their distribution
    #pragma omp parallel
    {
        // Local rdf to prevent race conditions
        BDF bdf_;
        // Local bead type count to prevent race conditions
        std::map<std::tuple<char, char>, int> bond_type_count_;
        #pragma omp for
        for (size_t i=0; i<_cg_data.bead_bonds.size(); i++) {
            for (auto j : _cg_data.bead_bonds[i]) {
                // Bond vector between two beads in system coordinates within PBC
                Vector3d r = _unwrapped_bond_vector(timestep, i, j);
                // Determine the bond type
                auto bond_type = _determine_bond_type(i, j);
                if (!(bdf_.find(bond_type) != bdf_.end())) {
                    bdf_[bond_type] = ArrayXd::Zero(lengths.rows());
                }
                bdf_[bond_type] += exp(-0.5 * ((lengths - r.norm()) / w).square());
            }
        }
        // Combine the local into the global results
        _merge_results(bdf, bdf_, bond_type_count, bond_type_count_);
    }
    auto z = w * sqrt(2.0 * M_PI);
    for (auto bond_type : bdf) {
        bdf[bond_type.first] /= (z * bond_type_count[bond_type.first]);
    }
    // Transfer ownership of the data rather than copy it for better memory use
    bdfs.push_back(std::move(bdf));
}


std::tuple<char, char> KDE::_determine_bond_type(int i, int j) {
    std::tuple<char, char> bond_type;
    char type_i = _cg_data.beads[i].name;
    char type_j = _cg_data.beads[j].name;
    // Prevent the reverse type to be counted separately
    if (type_i < type_j) {
        bond_type = std::make_tuple(type_i, type_j);
    }
    else {
        bond_type = std::make_tuple(type_j, type_i);
    }
    return bond_type;
}


void KDE::calculate_adf(int timestep) {
    // If there is no angle return
    if (_cg_data.bead_angles.size() == 0)
        return;
    // Final adf and angle type count after adding results from all threads
    ADF adf;
    std::map<std::tuple<char, char, char>, int> angle_type_count;
    // KDE Bandwidth
    auto w = _inputs.adf[3];
    // Find different angle types and their distribution in parallel
    #pragma omp parallel
    {
        // Local adf and angle type to prevent race conditions
        ADF adf_;
        std::map<std::tuple<char, char, char>, int> angle_type_count_;

        #pragma omp for
        for (auto a_type : _cg_data.bead_angles) {
            int i = std::get<0>(a_type);
            int j = std::get<1>(a_type);
            int k = std::get<2>(a_type);
            // Bond vector between two beads in system coordinates within PBC
            auto r1 = _unwrapped_bond_vector(timestep, j, i);
            auto r2 = _unwrapped_bond_vector(timestep, j, k);
            auto theta = _theta_angle(r1, r2) * 180.0 / M_PI;
            // Determine the angle type
            auto angle_type = _determine_angle_type(i, j, k);
            // Initialize the count if the type already does not exist
            if (!(angle_type_count_.find(angle_type) != angle_type_count_.end())) {
                angle_type_count_[angle_type] = 0;
            }
            angle_type_count_[angle_type] += 1;
            // Initialize local tdf and store the current type distribution
            if (!(adf_.find(angle_type) != adf_.end())) {
                adf_[angle_type] = ArrayXd::Zero(angles.rows());
            }
            adf_[angle_type] += exp(-0.5 * ((angles - theta) / w).square());
        }
        // Combine the local into the global results
        _merge_results(adf, adf_, angle_type_count, angle_type_count_);
    }
    auto z = w * sqrt(2.0 * M_PI);
    for (auto angle_type : adf) {
        adf[angle_type.first] /= (z * angle_type_count[angle_type.first]);
    }
    // Transfer ownership of the data rather than copy it for better memory use
    adfs.push_back(std::move(adf));
}


std::tuple<char, char, char> KDE::_determine_angle_type(int i, int j, int k) {
    std::tuple<char, char, char> angle_type;
    char type_i = _cg_data.beads[i].name;
    char type_j = _cg_data.beads[j].name;
    char type_k = _cg_data.beads[k].name;
    // Prevent the reverse type to be counted separately
    if (type_i <= type_k) {
        angle_type = std::make_tuple(type_i, type_j, type_k);
    }
    else {
        angle_type = std::make_tuple(type_k, type_j, type_i);
    }
    return angle_type;
}


void KDE::calculate_tdf(int timestep) {
    // If there is no torsion return
    if (_cg_data.bead_torsions.size() == 0)
        return;
    // Final tdf and torsion angle type count after adding results from all threads
    TDF tdf;
    std::map<std::tuple<char, char, char, char>, int> torsion_type_count;
    // KDE Bandwidth
    auto w = _inputs.tdf[3];
    // Find different torsion types and their distribution in parallel
    #pragma omp parallel
    {
        // Local tdf and torsion angle type to prevent race conditions
        TDF tdf_;
        std::map<std::tuple<char, char, char, char>, int> torsion_type_count_;

        #pragma omp for
        for (auto t_type : _cg_data.bead_torsions) {
            int i = std::get<0>(t_type);
            int j = std::get<1>(t_type);
            int k = std::get<2>(t_type);
            int l = std::get<3>(t_type);
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
            // Determine torsion type
            auto torsion_type = _determine_torsion_type(i, j, k, l);
            // Initialize the count if the type already does not exist
            if (!(torsion_type_count_.find(torsion_type) != torsion_type_count_.end())) {
                torsion_type_count_[torsion_type] = 0;
            }
            torsion_type_count_[torsion_type] += 1;
            // Initialize local tdf and store the current type distribution
            if (!(tdf_.find(torsion_type) != tdf_.end())) {
                tdf_[torsion_type] = ArrayXd::Zero(torsions.rows());
            }
            tdf_[torsion_type] += exp(-0.5 * (xx / w).square());
        }
        // Combine the local into the global results
        _merge_results(tdf, tdf_, torsion_type_count, torsion_type_count_);
    }
    auto z = w * sqrt(2.0 * M_PI);
    for (auto torsion_type : tdf) {
        tdf[torsion_type.first] /= z * torsion_type_count[torsion_type.first];
    }
    // Transfer ownership of the data rather than copy it for better memory use
    tdfs.push_back(std::move(tdf));
}


std::tuple<char, char, char, char> KDE::_determine_torsion_type(int i, int j, int k, int l) {
    std::tuple<char, char, char, char> torsion_type;
    char type_i = _cg_data.beads[i].name;
    char type_j = _cg_data.beads[j].name;
    char type_k = _cg_data.beads[k].name;
    char type_l = _cg_data.beads[l].name;
    // Prevent the reverse type to be counted separately
    if (type_i <= type_l) {
        torsion_type = std::make_tuple(type_i, type_j, type_k, type_l);
    }
    else {
        torsion_type = std::make_tuple(type_l, type_j, type_k, type_i);
    }
    return torsion_type;
}


void KDE::calculate_idf(int timestep) {
    // If there is no improper return
    if (_cg_data.bead_impropers.size() == 0)
        return;
    // Final idf and improper angle type count after adding results from all threads
    IDF idf;
    std::map<std::tuple<char, char, char, char>, int> improper_type_count;
    // KDE Bandwidth
    auto w = _inputs.idf[3];
    // Find different angle types and their distribution
    #pragma omp parallel
    {
        // Local idf and improper angle type to prevent race conditions
        IDF idf_;
        std::map<std::tuple<char, char, char, char>, int> improper_type_count_;

        for (auto i_type : _cg_data.bead_impropers) {
            int i = std::get<0>(i_type);
            int j = std::get<1>(i_type);
            int k = std::get<2>(i_type);
            int l = std::get<3>(i_type);
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
            // Determine improper angle type
            auto improper_type = _determine_improper_type(i, j, k, l);
            // Initialize the count if the type already does not exist
            if (!(improper_type_count_.find(improper_type) != improper_type_count_.end())) {
                improper_type_count_[improper_type] = 0;
            }
            improper_type_count_[improper_type] += 1;
            // Initialize local idf and store the current type distribution
            if (!(idf_.find(improper_type) != idf_.end())) {
                idf_[improper_type] = ArrayXd::Zero(impropers.rows());
            }
            idf_[improper_type] += exp(-0.5 * ((impropers - psi) / w).square());
        }
        // Combine the local into the global results
        _merge_results(idf, idf_, improper_type_count, improper_type_count_);
    }
    auto z = w * sqrt(2.0 * M_PI);
    for (auto improper_type : idf) {
        idf[improper_type.first] /= (z * improper_type_count[improper_type.first]);
    }
    // Transfer ownership of the data rather than copy it for better memory use
    idfs.push_back(std::move(idf));
}


std::tuple<char, char, char, char> KDE::_determine_improper_type(int i, int j, int k, int l) {
    std::tuple<char, char, char, char> improper_type;
    char type_i = _cg_data.beads[i].name;
    char type_j = _cg_data.beads[j].name;
    char type_k = _cg_data.beads[k].name;
    char type_l = _cg_data.beads[l].name;
    // Prevent the reverse type to be counted separately
    if (type_i <= type_l) {
        improper_type = std::make_tuple(type_i, type_j, type_k, type_l);
    }
    else {
        improper_type = std::make_tuple(type_l, type_j, type_k, type_i);
    }
    return improper_type;
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
    return std::find(_cg_data.bead_bonds[i].begin(), _cg_data.bead_bonds[i].end(), j) !=
                                                     _cg_data.bead_bonds[i].end();
}


bool KDE::_are_second_neighbors(int i, int j) {
    for (auto b1 : _cg_data.bead_bonds[i]) {
        if (_are_first_neighbors(b1, j)) {
                return true;
        }
    }
    return false;
}


bool KDE::_are_third_neighbors(int i, int j) {
    for (auto b1 : _cg_data.bead_bonds[i]) {
        for (auto b2 : _cg_data.bead_bonds[b1]) {
            if (_are_first_neighbors(b2, j)) {
                return true;
            }
        }
    }
    return false;
}


Vector3d KDE::_unwrapped_bond_vector(int ts, int i, int j) {
    Vector3d x_i = _cg_dump.cg_ts[ts].coords[i];
    Vector3d x_j = _cg_dump.cg_ts[ts].coords[j];
    // Displacement vector in Cartesian coordinates
    Vector3d dx = x_j - x_i;
    // Convert to fractional coordinates by solving: box^-1 * dx = dxs
    Vector3d dxs = (_cg_dump.cg_ts[ts].box.inverse() * dx).array().round().matrix();
    return dx - _cg_dump.cg_ts[ts].box * dxs;
}


double KDE::_theta_angle(const VectorXd &V1, const VectorXd &V2) {
    return std::acos(V1.dot(V2) / (V1.norm() * V2.norm()));
}


template <typename DF, typename TypeCount>
void KDE::_merge_results(DF &df, const DF &df_, TypeCount &type_count,
                                          const TypeCount &type_count_) {
    #pragma omp critical
    {
        for (const auto& item : df_) {
            if (df.find(item.first) == df.end()) {
                df[item.first] = item.second;
            } else {
                df[item.first] += item.second;
            }
        }

        for (const auto& item : type_count_) {
            type_count[item.first] += item.second;
        }
    }
}

