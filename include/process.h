#pragma once
#include <fmt/format.h>
#include <Eigen/Dense>
#include <thread>
#include <condition_variable>
#include <queue>
#include <chrono>
#include "distributions.h"


using RDF = std::map<std::tuple<char, char>, ArrayXd>;
using BDF = std::map<std::tuple<char, char>, ArrayXd>;
using ADF = std::map<std::tuple<char, char, char>, ArrayXd>;
using TDF = std::map<std::tuple<char, char, char, char>, ArrayXd>;
using IDF = std::map<std::tuple<char, char, char, char>, ArrayXd>;


class Process {
    public:
        Process(KDE &kde, int num_timesteps);
        void serial();
        void parallel();
        // Average distributions over all the timesteps
        RDF rdfs;
        BDF bdfs;
        ADF adfs;
        TDF tdfs;
        IDF idfs;
    private:
        // Function to take the average over all timesteps for each distribution
        template <class T>
        void timesteps_average(const std::vector<T> &_pdfs, T &pdfs);
        // Variables
        KDE *_kde;
        int _num_threads;
        int _num_timesteps;

        friend class Outputs;
};

