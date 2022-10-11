
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
using namespace std;

#include <cblas.h>

#include <mpi.h>
#include <boost/mpi.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <chrono>
using namespace std::chrono;

#include "ReactionDiffusion.h"


int main(int argc, char* argv[]) {
    int rank;
    int size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // run on C cores = c^2 for integer 1 ≤ c ≤ 8.

    // accept command line arguments
    po::options_description opts("Available options.");
    opts.add_options()
        ("dt", po::value<double>()->default_value(0.001, "0.001"), "Time-step to use.")
        ("T", po::value<double>()->default_value(100, "100"), "Total integration time.")
        ("Nx", po::value<int>()->default_value(101, "101"), "Number of grid points in x.")
        ("Ny", po::value<int>()->default_value(101, "101"), "Number of grid points in y.")
        ("a", po::value<double>()->default_value(0.75, "0.75"), "Value of parameter a.")
        ("b", po::value<double>()->default_value(0.06, "0.06"), "Value of parameter b.")
        ("mu1", po::value<double>()->default_value(5, "5"), "Value of parameter mu1.")
        ("mu2", po::value<double>()->default_value(0, "0"), "Value of parameter mu2.")
        ("eps", po::value<double>()->default_value(50, "50"), "Value of parameter epsilon.")
        ("help", "Print help message.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (rank == 0) {
        if (vm.count("help")) {
            std::cout << opts << std::endl;
        }
    }


    // initialise parameters from command line or the default values
    double dt = vm["dt"].as<double>();
    double integration_time = vm["T"].as<double>();
    int Nx = vm["Nx"].as<int>();
    int Ny = vm["Ny"].as<int>();
    double a = vm["a"].as<double>();
    double b = vm["b"].as<double>();
    double mu_1 = vm["mu1"].as<double>();
    double mu_2 = vm["mu2"].as<double>();
    double eps = vm["eps"].as<double>();

    if (0 && rank == 0) {
        cout << "" << endl;
        cout << "dt: " << dt << endl;
        cout << "integration_time: " << integration_time << endl;
        cout << "Nx: " << Nx << endl;
        cout << "Ny: " << Ny << endl;
        cout << "a: " << a << endl;
        cout << "b: " << b << endl;
        cout << "mu_1: " << mu_1 << endl;
        cout << "mu_2: " << mu_2 << endl;
        cout << "eps: " << eps << endl;
        cout << "" << endl;
    }


    ReactionDiffusion barkley(dt, integration_time, Nx, 
            Ny, a, b, mu_1, 
            mu_2, eps, rank, size);
    
    // double myvar = 0;
    // if (rank == 0) cout << "\n myvar: " << myvar << endl;
    barkley.SetInitialConditions();

    barkley.TimeIntegrate();
    // if (rank == 0) cout << "\n myvar: " << myvar << endl;


    MPI_Finalize();
    return 0;
}

