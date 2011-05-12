#ifndef OPTIONS_H_INCLUDED
#define OPTIONS_H_INCLUDED

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

namespace bpo = boost::program_options;

void parse_options(int argc, char **argv, options* opt) {
    cout << "\nParse Program Options\n";
    bpo::options_description desc("Configuration ");
    desc.add_options()
        ("help", "show possible options")
        ("path", bpo::value<string>(), "path to write and load from")
        ("t_step", bpo::value<float>(), "simulation time step")
        ("total_steps", bpo::value<int>(), "total time steps to compute")
        ("buffer", bpo::value<int>(), "buffer size for correlation function")
        ("system_size", bpo::value<int>(), "lattice size s x s")
        ("krylov_dim", bpo::value<int>(), "dimension of krylov space")
        ("append", bpo::value<bool>(), "set 1 to append to previous simulation")
        ("job", bpo::value<char>(), "job to simulate d := dos, w := wavepacket, c := correlation function, v := play wave propagation")
        ("grid_w", bpo::value<int>(), "cluster grid width")
        ("grid_h", bpo::value<int>(), "cluster grid height")
        ("defects_A", bpo::value<float>(), "defects in sublattice A in %")
        ("defects_B", bpo::value<float>(), "defects in sublattice B in %")
        ("omp_threads", bpo::value<int>(), "number of omp threads")
        ("serial", bpo::value<bool>(), "serial or not")
        ("graphene", bpo::value<bool>(), "graphene or square lattice")
        ("debug", bpo::value<bool>(), "debug mode")
        ("seed_system", bpo::value<int>(), "system seed")
        ("seed_disorder", bpo::value<int>(), "disorder seed")
        ("dos_crop", bpo::value<float>(), "nimmt nicht alle zeitschritte mit")
    ;

    bpo::variables_map vm;
    ifstream in("sim.cfg");
    bpo::store(bpo::parse_command_line(argc, argv, desc), vm);
    bpo::store(bpo::parse_config_file(in, desc), vm);
    in.close();
    bpo::notify(vm);

    //help-----
    if (vm.count("help")) { cout << desc << "\n"; exit(1); }

    //job------
    if (vm.count("job")) opt->job = vm["job"].as<char>();
    else { cout << desc << "\n"; exit(1); }

    //path-----
    if (vm.count("path")) opt->path = vm["path"].as<string>();
    //else { cout << desc << "\n"; exit(1); }

    //load-----
    if (vm.count("append")) opt->append = vm["append"].as<bool>();
    if (vm.count("serial")) opt->serial = vm["serial"].as<bool>();
    if (vm.count("graphene")) opt->graphene = vm["graphene"].as<bool>();
    if (vm.count("debug")) opt->debug = vm["debug"].as<bool>();

    //other parameter
    if (vm.count("total_steps")) opt->T = vm["total_steps"].as<int>();
    if (vm.count("t_step")) opt->dt = vm["t_step"].as<float>();
    if (vm.count("buffer")) opt->N_buffer = vm["buffer"].as<int>();
    if (vm.count("system_size")) opt->k = vm["system_size"].as<int>();
    if (vm.count("krylov_dim")) opt->m = vm["krylov_dim"].as<int>();

    if (vm.count("grid_w")) opt->grid_w = vm["grid_w"].as<int>();
    if (vm.count("grid_h")) opt->grid_h = vm["grid_h"].as<int>();

    if (vm.count("seed_system")) opt->seed_system = vm["seed_system"].as<int>();
    if (vm.count("seed_disorder")) opt->seed_disorder = vm["seed_disorder"].as<int>();
    if (vm.count("dos_crop")) opt->dos_crop = vm["dos_crop"].as<float>();

    //specific defects
    if (vm.count("defects_A")) opt->def_A = vm["defects_A"].as<float>();
    if (vm.count("defects_B")) opt->def_B = vm["defects_B"].as<float>();

    //cout << "\nopt sys size : " << opt->k << "\n";
    //cout << "\nopt grid w : " << opt->grid_w << "\n";
    //cout << "\nopt grid h : " << opt->grid_h << "\n";
}


#endif // OPTIONS_H_INCLUDED
