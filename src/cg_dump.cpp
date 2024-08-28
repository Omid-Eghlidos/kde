#include "cg_dump.h"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tuple>
#include <iostream>
#include <fstream>


CG_dump::CG_dump(const LMP_data &data, const LMP_dump &dump, const CG_data &cgd) :
    _data(data), _dump(dump), _cg_data(cgd) {
    _coarsen_timesteps();
}


void CG_dump::_coarsen_timesteps () {
    cg_ts = std::vector<Timestep> (_dump.timesteps.size());
    for (size_t ts = 0; ts < _dump.timesteps.size(); ts++) {
        cg_ts[ts].timestep = _dump.timesteps[ts];
        cg_ts[ts].box = _dump.box[ts];
        cg_ts[ts].coords = BeadsCoords (_cg_data.beads.size());
        for (auto bead : _cg_data.beads) {
            // Position of the first atom in the bead
            Vector3d x0 = _dump.atom_coords[ts].row(bead.atoms[0]);
            Vector3d bead_coords = x0 * bead.atom_weights[0];
            for (size_t i = 1; i < bead.atoms.size(); i++) {
                Vector3d x = _dump.atom_coords[ts].row(bead.atoms[i]);
                Vector3d dx = x - x0;
                Vector3d ds = (_dump.box[ts].inverse() * dx).array().round().matrix();
                dx -= _dump.box[ts] * ds;
                bead_coords += (x0 + dx) * bead.atom_weights[i];
            }
            cg_ts[ts].coords[bead.id] = bead_coords;
        }
    }
    fmt::print("Found coarse coordinates for {} timsteps.\n", cg_ts.size());
}

