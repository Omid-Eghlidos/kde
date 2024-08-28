#pragma once
#include "lmp_data.h"
#include "lmp_dump.h"
#include "kde.h"
#include "process.h"


/**
 * @class Outputs
 * @brief Handles the output of distribution results to files.
 *
 * The Outputs class is responsible for writing the results of the
 * coarse-grained simulation distributions to files, using the provided
 * tag, LAMMPS data, dump, and process results.
 */
class Outputs {
    public:
        /**
         * @brief Constructs an Outputs object with the given distribution results.
         *
         * @param tag The tag for the output files.
         * @param data Reference to the LMP_data object containing LAMMPS data.
         * @param dump Reference to the LMP_dump object containing dump data.
         * @param results Reference to the Process object containing the distribution results.
         */
        Outputs(std::string tag, const LMP_data &data, const LMP_dump &dump
                               , const Process &results);

        /**
         * @brief Writes the distribution results to files.
         *
         * This function writes the specified distribution results (e.g., RDF, BDF)
         * to output files using the provided tag and the results from the process.
         *
         * @param distribution The name of the distribution to write (e.g., "rdf", "bdf").
         */
        void write_distributions(std::string distribution) const;
    private:
        //! The tag used for naming output files
        std::string _tag;
        //! Reference to the LMP_data object
        const LMP_data &_data;
        //! Reference to the LMP_dump object
        const LMP_dump &_dump;
        //! Reference to the Process object containing the results
        const Process &_results;
};

