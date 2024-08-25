#include "outputs.h"
#include <iostream>
#include <fstream>
#include <fmt/format.h>
#include "fmt/ostream.h"


Outputs::Outputs(std::string tag, const LMP_data &data, const LMP_dump &dump
        , const Process &results) : _tag(tag), _data(data), _dump(dump)
                                  , _results(results) {
    fmt::print("\n################ Writing the outputs ####################\n");
    // Write the distributions if true
    if (!_results.rdfs.empty()) {
        write_distributions("rdf");
    }
    if (!_results.bdfs.empty()) {
        write_distributions("bdf");
    }
    if (!_results.adfs.empty()) {
        write_distributions("adf");
    }
    if (!_results.tdfs.empty()) {
        write_distributions("tdf");
    }
    if (!_results.idfs.empty()) {
        write_distributions("idf");
    }
}


void Outputs::write_distributions(std::string distribution) const {
    if (distribution == "rdf") {
        // Write RDFs
        for (auto pair_type : _results.rdfs) {
            char type_i = std::get<0>(pair_type.first);
            char type_j = std::get<1>(pair_type.first);
            std::string path = "rdf-" + _tag + "_" + type_i + type_j + ".txt";
            std::fstream fid(path, std::fstream::out);
            fmt::print(fid, "Radius\t RDF \n");
            fmt::print(fid, "------\t-----\n");
            for (int i=0; i<pair_type.second.rows(); i++) {
                fmt::print(fid, "{:6.3f}\t{:6.3f}\n"
                              , _results._kde->radius(i), pair_type.second(i));
            }
        }
        fmt::print("RDFs are written into files.\n");
    }
    else if (distribution == "bdf") {
        // Write BDFs
        for (auto bond_type : _results.bdfs) {
            char type_i = std::get<0>(bond_type.first);
            char type_j = std::get<1>(bond_type.first);
            std::string path = distribution + "-" + _tag + "_"
                             + type_i + type_j + ".txt";
            std::fstream fid(path, std::fstream::out);
            fmt::print(fid, "Lengths  BDF \n");
            fmt::print(fid, "------\t-----\n");
            for (int i=0; i<bond_type.second.rows(); i++) {
                fmt::print(fid, "{:6.3f}\t{:6.3f}\n"
                              , _results._kde->lengths(i), bond_type.second(i));
            }
        }
        fmt::print("BDFs are written into files.\n");
    }
    else if (distribution == "adf") {
        // Write ADFs
        for (auto angle_type : _results.adfs) {
            char type_i = std::get<0>(angle_type.first);
            char type_j = std::get<1>(angle_type.first);
            char type_k = std::get<2>(angle_type.first);
            std::string path = distribution + "-" + _tag + "_"
                             + type_i + type_j + type_k + ".txt";
            std::fstream fid(path, std::fstream::out);
            fmt::print(fid, " Angles\t  ADF \n");
            fmt::print(fid, " ------\t -----\n");
            for (int i=0; i<angle_type.second.rows(); i++) {
                fmt::print(fid, "{:7.3f}\t{:7.3f}\n"
                              , _results._kde->angles(i), angle_type.second(i));
            }
        }
        fmt::print("ADFs are written into files.\n");
    }
    else if (distribution == "tdf") {
        // Write TDFs
        for (auto torsion_type : _results.tdfs) {
            char type_i = std::get<0>(torsion_type.first);
            char type_j = std::get<1>(torsion_type.first);
            char type_k = std::get<2>(torsion_type.first);
            char type_l = std::get<3>(torsion_type.first);
            std::string path = distribution + "-" + _tag + "_"
                             + type_i + type_j + type_k + type_l + ".txt";
            std::fstream fid(path, std::fstream::out);
            fmt::print(fid, "Angles\t\t TDF \n");
            fmt::print(fid, "------\t\t-----\n");
            for (int i=0; i<torsion_type.second.rows(); i++) {
                fmt::print(fid, "{:7.4f}\t{:7.4f}\n"
                              ,_results._kde->torsions(i), torsion_type.second(i));
            }
        }
        fmt::print("TDFs are written into files.\n");
    }
    else if (distribution == "idf") {
        // Write IDFs
        for (auto improper_type : _results.idfs) {
            char type_i = std::get<0>(improper_type.first);
            char type_j = std::get<1>(improper_type.first);
            char type_k = std::get<2>(improper_type.first);
            char type_l = std::get<3>(improper_type.first);
            std::string path = distribution + "-" + _tag + "_"
                             + type_i + type_j + type_k + type_l + ".txt";
            std::fstream fid(path, std::fstream::out);
            fmt::print(fid, "Angles\t\t IDF \n");
            fmt::print(fid, "------\t\t-----\n");
            for (int i=0; i<improper_type.second.rows(); i++) {
                fmt::print(fid, "{:7.4f}\t{:7.4f}\n"
                              ,_results._kde->impropers(i), improper_type.second(i));
            }
        }
        fmt::print("IDFs are written into files.\n");
    }
}

