#include "process.h"


Process::Process(KDE &kde, int num_ts) : _kde(&kde), _num_timesteps(num_ts) {
    _num_threads = std::thread::hardware_concurrency();
}


void Process::serial() {
    fmt::print("\n################### Distributions #######################\n");
    fmt::print("////////////////// Serial Processing ////////////////////\n\n");
    auto t_start = std::chrono::steady_clock::now();
    // RDF
    if (_kde->compute["rdf"]) {
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_rdf(timestep);
        }
        // Averages over all the timesteps
        timesteps_average(_kde->rdfs, rdfs);
        fmt::print("Computed RDFs for {} pair types.\n", rdfs.size());
    }
    // BDF
    if (_kde->compute["bdf"]) {
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_bdf(timestep);
        }
        // Averages over all the timesteps
        timesteps_average(_kde->bdfs, bdfs);
        fmt::print("Computed BDFs for {} bond types.\n", bdfs.size());
    }
    // ADF
    if (_kde->compute["adf"]) {
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_adf(timestep);
        }
        // Averages over all the timesteps
        timesteps_average(_kde->adfs, adfs);
        fmt::print("Computed ADFs for {} angle types.\n", adfs.size());
    }
    // TDF
    if (_kde->compute["tdf"]) {
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_tdf(timestep);
        }
        timesteps_average(_kde->tdfs, tdfs);
        fmt::print("Computed TDFs for {} torsion types.\n", tdfs.size());
    }
    // IDF
    if (_kde->compute["idf"]) {
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_idf(timestep);
        }
        timesteps_average(_kde->idfs, idfs);
        fmt::print("Computed IDFs for {} improper types.\n", idfs.size());
    }
    // Wall time
    auto t_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> t_elapsed = t_end - t_start;
    int minutes = t_elapsed.count() / 60;
    int hours = minutes / 60;
    double seconds = t_elapsed.count() - minutes * 60;
    fmt::print("\nFinished in {:02d} (hr) : {:02d} (min) : {:6.3f} (s).\n"
                            , hours       , minutes      , seconds);
}


void Process::parallel() {
    fmt::print("\n################### Distributions #######################\n");
    fmt::print("///////////////// Parallel Processing ///////////////////\n\n");
    auto t_start = std::chrono::steady_clock::now();
    fmt::print("Running on {} threads ...\n\n", _num_threads);
    // RDF
    if (_kde->compute["rdf"]) {
        #pragma omp parallel for
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_rdf(timestep);
        }
        // Averages over all the timesteps
        timesteps_average(_kde->rdfs, rdfs);
        fmt::print("Computed RDFs for {} pair types.\n", rdfs.size());
    }
    // BDF
    if (_kde->compute["bdf"]) {
        #pragma omp parallel for
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_bdf(timestep);
        }
        if (!_kde->bdfs.empty()) {
            // Averages over all the timesteps
            timesteps_average(_kde->bdfs, bdfs);
            fmt::print("Computed BDFs for {} bond types.\n", bdfs.size());
        }
        else {
            fmt::print("No bonds exist to compute their distribution.\n");
        }
    }
    // ADF
    if (_kde->compute["adf"]) {
        #pragma omp parallel for
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_adf(timestep);
        }
        if (!_kde->adfs.empty()) {
            // Averages over all the timesteps
            timesteps_average(_kde->adfs, adfs);
            fmt::print("Computed ADFs for {} angle types.\n", adfs.size());
        }
        else {
            fmt::print("No bond angles exist to compute their distribution.\n");
        }
    }
    // TDF
    if (_kde->compute["tdf"]) {
        #pragma omp parallel for
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_tdf(timestep);
        }
        if (!_kde->tdfs.empty()) {
            // Averages over all the timesteps
            timesteps_average(_kde->tdfs, tdfs);
            fmt::print("Computed TDFs for {} torsion types.\n", tdfs.size());
        }
        else {
            fmt::print("No torsion angles exist to compute their distribution.\n");
        }
    }
    // IDF
    if (_kde->compute["idf"]) {
        #pragma omp parallel for
        for (int timestep=0; timestep<_num_timesteps; timestep++) {
            _kde->calculate_idf(timestep);
        }
        if (!_kde->idfs.empty()) {
            // Averages over all the timesteps
            timesteps_average(_kde->idfs, idfs);
            fmt::print("Computed IDFs for {} improper types.\n", idfs.size());
        }
        else {
            fmt::print("No improper angles exist to compute their distribution.\n");
        }
    }
    // Wall time
    auto t_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> t_elapsed = t_end - t_start;
    int minutes = t_elapsed.count() / 60;
    int hours = minutes / 60;
    double seconds = t_elapsed.count() - minutes * 60;
    fmt::print("\nFinished in: {:02d} (hr) : {:02d} (min) : {:6.3f} (s).\n"
                            , hours       , minutes      , seconds);
}


// Calculate the average distributions for timesteps
template <class T>
void Process::timesteps_average(const std::vector<T> &_pdfs, T &pdfs) {
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

