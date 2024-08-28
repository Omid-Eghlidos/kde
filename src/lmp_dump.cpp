#include "lmp_dump.h"
#include "string_tools.h"
#include <fmt/format.h>
#include <iostream>


LMP_dump::LMP_dump(std::string path) {
    std::fstream fid;
    fid.open(path, std::fstream::in);
    if (!fid) {
        std::cout << "Dump file does not exist.\n";
        exit(0);
    }
    double xlo = 0, xhi, ylo = 0, yhi, zlo = 0, zhi;
    while(fid) {
        auto line = read_line(fid);
        auto args = split(line);
        if (args.empty())
            continue;
        else if (line == "ITEM: TIMESTEP") {
            int n = from_string<int>(read_line(fid));
            timesteps.push_back(n);
        }
        else if (line == "ITEM: NUMBER OF ATOMS") {
            num_atoms = from_string<int>(read_line(fid));
        }
        else if (startswith(line, "ITEM: BOX BOUNDS")) {
            MatrixXd timestep_box = MatrixXd::Zero(3, 3);
            MatrixXd timestep_bounds = MatrixXd::Zero(2, 3);
            bool triclinic = false;
            for (auto arg : args) {
                if (arg == "xy" || arg == "xz" || arg == "yz") {
                    triclinic = true;
                }
            }
            if (triclinic == true) {
                // Triclinic box : https://lammps.sandia.gov/doc/Howto_triclinic.html
                // Read xlo_bound and xhi_bound xy
                args = split(read_line(fid));
                auto xlo_bound = str2dbl(args[0]);
                auto xhi_bound = str2dbl(args[1]);
                auto xy = str2dbl(args[2]);
                // Read ylo_bound and yhi_bound yz
                args = split(read_line(fid));
                auto ylo_bound = str2dbl(args[0]);
                auto yhi_bound = str2dbl(args[1]);
                auto xz = str2dbl(args[2]);
                // Read zlo_bound and zhi_bound xy
                args = split(read_line(fid));
                auto zlo_bound = str2dbl(args[0]);
                auto zhi_bound = str2dbl(args[1]);
                auto yz = str2dbl(args[2]);

                // Box parameters xlo, xhi, ylo, yhi, zlo, zhi
                xlo = xlo_bound - std::min(0.0, std::min(xy, std::min(xz, xy+xz)));
                xhi = xhi_bound - std::max(0.0, std::min(xy, std::min(xz, xy+xz)));
                ylo = ylo_bound - std::min(0.0, yz);
                yhi = yhi_bound - std::max(0.0, yz);
                zlo = zlo_bound;
                zhi = zhi_bound;

                timestep_box(0, 0) = xhi - xlo;
                timestep_box(0, 1) = xy;
                timestep_box(0, 2) = xz;
                timestep_box(1, 0) = 0.0;
                timestep_box(1, 1) = yhi - ylo;
                timestep_box(1, 2) = yz;
                timestep_box(2, 0) = 0.0;
                timestep_box(2, 1) = 0.0;
                timestep_box(2, 2) = zhi - zlo;
            }
            else {
                args = split(read_line(fid));
                xlo = str2dbl(args[0]);
                xhi = str2dbl(args[1]);
                args = split(read_line(fid));
                ylo = str2dbl(args[0]);
                yhi = str2dbl(args[1]);
                args = split(read_line(fid));
                zlo = str2dbl(args[0]);
                zhi = str2dbl(args[1]);

                timestep_box(0, 0) = xhi - xlo;
                timestep_box(1, 1) = yhi - ylo;
                timestep_box(2, 2) = zhi - zlo;
            }
            box.push_back(timestep_box);
        }
        else if (startswith(line, "ITEM: ATOMS")) {
            if (std::find(args.begin(), args.end(), "xs") != args.end()) {
                coords_type = "scaled";
            }
            else if (std::find(args.begin(), args.end(), "x") != args.end()) {
                coords_type = "unscaled";
            }
            else if (std::find(args.begin(), args.end(), "xu") != args.end()) {
                coords_type = "unwrapped";
            }
            else {
                fmt::print("Unrecognized coordinate type.\n");
                fmt::print("Recognized coordinate types are xs, x, and xu.\n");
                exit(0);
            }
            MatrixXd timestep_coords(num_atoms, 3);
            int counter = 0;
            while(fid && counter < num_atoms) {
                line = read_line(fid);
                args = split(line);
                int num_args = args.size();
                int i = from_string<int>(args[0]) - 1;
                timestep_coords(i, 0) = str2dbl(args[num_args - 3]) - xlo;
                timestep_coords(i, 1) = str2dbl(args[num_args - 2]) - ylo;
                timestep_coords(i, 2) = str2dbl(args[num_args - 1]) - zlo;
                Vector3d xs = (box.back().inverse() * timestep_coords.row(i).transpose());
                timestep_coords.row(i) -= (box.back() * xs.array().floor().matrix());
                counter++;
            }
            atom_coords.push_back(timestep_coords);
        }
    }
    fmt::print("Read {} timesteps from {}.\n", timesteps.size(), path);
}

