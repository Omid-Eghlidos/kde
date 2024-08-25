#pragma once
#include "lammps_data.h"
#include "lammps_dump.h"
#include "distributions.h"
#include "process.h"


class Outputs {
    public:
        // Constructor with distribution results
        Outputs(std::string tag, const LMP_data &data, const LMP_dump &dump
                               , const Process &results);
        // Write the distributions
        void write_distributions(std::string distribution) const;
    private:
        std::string _tag;
        const LMP_data &_data;
        const LMP_dump &_dump;
        const Process &_results;
};

