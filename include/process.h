#pragma once
#include <fmt/format.h>
#include <Eigen/Dense>
#include <condition_variable>
#include <queue>
#include <chrono>
#include "kde.h"

using RDF = std::map<std::tuple<char, char>, ArrayXd>;
using BDF = std::map<std::tuple<char, char>, ArrayXd>;
using ADF = std::map<std::tuple<char, char, char>, ArrayXd>;
using TDF = std::map<std::tuple<char, char, char, char>, ArrayXd>;
using IDF = std::map<std::tuple<char, char, char, char>, ArrayXd>;


/**
 * @class Process
 * @brief Handles the processing and averaging of distribution functions over
 * multiple timesteps.
 *
 * The Process class computes and averages the distribution functions (RDF, BDF, ADF, TDF, IDF)
 * over all timesteps in a simulation, using the inputs and Kernel Density Estimation (KDE)
 * methods.
 */
class Process {
    public:
        /**
         * @brief Constructs a Process object to manage distribution calculations.
         *
         * @param inputs The input parameters for the simulation.
         * @param kde Reference to the KDE object for calculating distributions.
         * @param num_timesteps The number of timesteps to process.
         */
        Process(const Inputs &inputs, KDE &kde, int num_timesteps);

        //! Stores the averaged Radial Distribution Function (RDF)
        RDF rdfs;
        //! Stores the averaged Bond Distribution Function (BDF)
        BDF bdfs;
        //! Stores the averaged Angle Distribution Function (ADF)
        ADF adfs;
        //! Stores the averaged Torsion Distribution Function (TDF)
        TDF tdfs;
        //! Stores the averaged Improper Distribution Function (IDF)
        IDF idfs;
    private:
        /**
         * @brief Runs the calculations for the radial distribution functions (RDFs)
         * over all timesteps.
         *
         * This function manages the calculation of the RDFs using KDE for each
         * timestep and stores the results.
         */
        void _compute_rdf();

        /**
         * @brief Runs the calculations for the bond distribution functions (BDFs)
         * over all timesteps.
         *
         * This function manages the calculation of the BDFs using KDE for each
         * timestep and stores the results.
         */
        void _compute_bdf();

        /**
         * @brief Runs the calculations for the bond angle distribution functions (ADFs)
         * over all timesteps.
         *
         * This function manages the calculation of the ADFs using KDE for each
         * timestep and stores the results.
         */
        void _compute_adf();

        /**
         * @brief Runs the calculations for the torsion angle distribution functions (TDFs)
         * over all timesteps.
         *
         * This function manages the calculation of the TDFs using KDE for each
         * timestep and stores the results.
         */
        void _compute_tdf();

        /**
         * @brief Runs the calculations for the improper angle distribution functions (IDFs)
         * over all timesteps.
         *
         * This function manages the calculation of the IDFs using KDE for each
         * timestep and stores the results.
         */
        void _compute_idf();

        /**
         * @brief Averages the distribution functions over all timesteps.
         *
         * @tparam T The type of the distribution function (RDF, BDF, ADF, TDF, IDF).
         * @param _pdfs A vector containing the distribution functions for each timestep.
         * @param pdfs The resulting averaged distribution function.
         */
        template <class T>
        void _timesteps_average(const std::vector<T> &_pdfs, T &pdfs);

        /**
         * @brief Calculates the wall time for processing of a function
         *
         * This function calculate the amount of wall time for running a function.
         *
         * @tparam Time The timepoint class
         * @param t_start Starting time of the computation right before a function call.
         */
        template <typename Time>
        void _compute_wall_time(const Time& t_start);

        //! Reference to the input parameters
        const Inputs &_inputs;
        //! Pointer to the KDE object used for distribution calculations
        KDE *_kde;
        //! The number of timesteps to process
        int _num_timesteps;

        //! Allows the Outputs class to access private members of this class
        friend class Outputs;
};

