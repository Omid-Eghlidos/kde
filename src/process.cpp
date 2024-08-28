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
    _run();
    // Wall time
    auto t_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> t_elapsed = t_end - t_start;
    int minutes = t_elapsed.count() / 60;
    int hours = minutes / 60;
    double seconds = t_elapsed.count() - minutes * 60;
    fmt::print("\nFinished in: {:02d} (hr) : {:02d} (min) : {:6.3f} (s).\n"
                            , hours       , minutes      , seconds);
}


void Process::_run() {
    // RDF
    if (_kde->compute["rdf"]) {
        for (int timestep = 0; timestep < _num_timesteps; timestep++) {
            _kde->calculate_rdf(timestep);
        }
        // Averages over all the timesteps
        _timesteps_average(_kde->rdfs, rdfs);
        fmt::print("Computed RDFs for {} pair types.\n", rdfs.size());
    }
    // BDF
    if (_kde->compute["bdf"]) {
        for (int timestep = 0; timestep < _num_timesteps; timestep++) {
            _kde->calculate_bdf(timestep);
        }
        if (!_kde->bdfs.empty()) {
            // Averages over all the timesteps
            _timesteps_average(_kde->bdfs, bdfs);
            fmt::print("Computed BDFs for {} bond types.\n", bdfs.size());
        }
        else {
            fmt::print("No bonds exist to compute their distribution.\n");
        }
    }
    // ADF
    if (_kde->compute["adf"]) {
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_adf(timestep);
        }
        if (!_kde->adfs.empty()) {
            // Averages over all the timesteps
            _timesteps_average(_kde->adfs, adfs);
            fmt::print("Computed ADFs for {} angle types.\n", adfs.size());
        }
        else {
            fmt::print("No bond angles exist to compute their distribution.\n");
        }
    }
    // TDF
    if (_kde->compute["tdf"]) {
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_tdf(timestep);
        }
        if (!_kde->tdfs.empty()) {
            // Averages over all the timesteps
            _timesteps_average(_kde->tdfs, tdfs);
            fmt::print("Computed TDFs for {} torsion types.\n", tdfs.size());
        }
        else {
            fmt::print("No torsion angles exist to compute their distribution.\n");
        }
    }
    // IDF
    if (_kde->compute["idf"]) {
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_idf(timestep);
        }
        if (!_kde->idfs.empty()) {
            // Averages over all the timesteps
            _timesteps_average(_kde->idfs, idfs);
            fmt::print("Computed IDFs for {} improper types.\n", idfs.size());
        }
        else {
            fmt::print("No improper angles exist to compute their distribution.\n");
        }
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

