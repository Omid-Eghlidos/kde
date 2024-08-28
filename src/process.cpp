#include "process.h"


Process::Process(const Inputs &inputs, KDE &kde, int num_ts) : _inputs(inputs),
                                         _kde(&kde), _num_timesteps(num_ts) {
    fmt::print("\n################### Distributions #######################\n");
    if (_inputs.parallel) {
        fmt::print("///////////////// Parallel Processing ///////////////////\n\n");
        fmt::print("Running on {} threads ...\n\n", _inputs.num_threads);
    }
    else {
        fmt::print("////////////////// Serial Processing ////////////////////\n\n");
    }
    // Determine runtime
    auto t_start = std::chrono::steady_clock::now();

    if (_kde->compute["rdf"])
        _compute_rdf();
    if (_kde->compute["bdf"])
        _compute_bdf();
    if (_kde->compute["adf"])
        _compute_adf();
    if (_kde->compute["tdf"])
        _compute_tdf();
    if (_kde->compute["idf"])
        _compute_idf();

    fmt::print("\nTotal Time: ");
    _compute_wall_time(t_start);
}


void Process::_compute_rdf() {
    auto t_start = std::chrono::steady_clock::now();
    for (int timestep = 0; timestep < _num_timesteps; timestep++) {
        _kde->calculate_rdf(timestep);
    }
    // Averages over all the timesteps
    _timesteps_average(_kde->rdfs, rdfs);
    fmt::print("Computed RDFs for {} pair types - ET: ", rdfs.size());
    _compute_wall_time(t_start);
}


void Process::_compute_bdf() {
    auto t_start = std::chrono::steady_clock::now();
    for (int timestep = 0; timestep < _num_timesteps; timestep++) {
        _kde->calculate_bdf(timestep);
    }
    if (!_kde->bdfs.empty()) {
        // Averages over all the timesteps
        _timesteps_average(_kde->bdfs, bdfs);
        fmt::print("Computed BDFs for {} bond types - ET: ", bdfs.size());
        _compute_wall_time(t_start);
    }
    else {
        fmt::print("No bonds exist - skip.\n");
    }
}


void Process::_compute_adf() {
    auto t_start = std::chrono::steady_clock::now();
    for (int timestep=0; timestep<_num_timesteps; timestep++) {
        _kde->calculate_adf(timestep);
    }
    if (!_kde->adfs.empty()) {
        // Averages over all the timesteps
        _timesteps_average(_kde->adfs, adfs);
        fmt::print("Computed ADFs for {} angle types - ET: ", adfs.size());
        _compute_wall_time(t_start);
    }
    else {
        fmt::print("No bond angles exist - skip.\n");
    }
}


void Process::_compute_tdf() {
    auto t_start = std::chrono::steady_clock::now();
    for (int timestep=0; timestep<_num_timesteps; timestep++) {
        _kde->calculate_tdf(timestep);
    }
    if (!_kde->tdfs.empty()) {
        // Averages over all the timesteps
        _timesteps_average(_kde->tdfs, tdfs);
        fmt::print("Computed TDFs for {} torsion types - ET: ", tdfs.size());
        _compute_wall_time(t_start);
    }
    else {
        fmt::print("No torsion angles exist - skip.\n");
    }
}


void Process::_compute_idf() {
    auto t_start = std::chrono::steady_clock::now();
    for (int timestep = 0; timestep < _num_timesteps; timestep++) {
        _kde->calculate_idf(timestep);
    }
    if (!_kde->idfs.empty()) {
        // Averages over all the timesteps
        _timesteps_average(_kde->idfs, idfs);
        fmt::print("Computed IDFs for {} improper types - ET: ", idfs.size());
        _compute_wall_time(t_start);
    }
    else {
        fmt::print("No improper angles exist - skip.\n");
    }
}


// Calculate the average distributions for timesteps
template <class T>
void Process::_timesteps_average(const std::vector<T> &_pdfs, T &pdfs) {
    // Timesteps distributions (tsd)
    for (auto ts_pdf : _pdfs) {
        for (auto type : ts_pdf) {
            if (!(pdfs.find(type.first) != pdfs.end())) {
                pdfs[type.first] = ArrayXd::Zero(type.second.rows());
            }
            pdfs[type.first] += type.second;
        }
    }
    for (auto type : pdfs) {
        pdfs[type.first] /= _num_timesteps;
    }
}


// Calculate wall time between two time points before and after a function call
template <typename Time>
void Process::_compute_wall_time(const Time& t_start) {
    // Time point after the call
    auto t_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> t_elapsed = t_end - t_start;
    int min = t_elapsed.count() / 60;
    int hr = min / 60;
    double sec = t_elapsed.count() - min * 60;
    fmt::print("{:02d} (hr) : {:02d} (min) : {:4.2f} (s)\n", hr, min, sec);
}

