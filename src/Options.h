#ifndef OPTIONS_H_INCLUDED
#define OPTIONS_H_INCLUDED

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

namespace bpo = boost::program_options;


struct head {
    //system, seeds, disorder
    int k;//system size kxk
    int seed_system;
    int seed_disorder;
    float dA;
    float dB;

    //simulation
    int m;//krylov dimension
    float dt;//timestep
    int T;//total timesteps to perform

    //final results
    int N;//korrelation funktion size

    //mpi variables
    int grid_id; //grid id, needed to resume or append
    int grid_w;
    int grid_h;

    string getPath() {//BS
        char path[200];
        sprintf(path, "%ix%ion%ix%i_ID%i_A%f\%_B%f\%_m%i_dt%f_T%i_%d.bin", k*grid_w, k*grid_h, grid_w, grid_h, grid_id, dA*100, dB*100, m, dt, T, (int)time(0));
        return string(path);
    }

    void writeHead(string path) {
        ofstream file(path.c_str(), fstream::in | fstream::out | fstream::binary);
        file.seekp(ios_base::beg);
        file.write((char*)this, sizeof(head));
        file.close();
    }

    void readHead(string path) {
        ifstream file(path.c_str(), fstream::in | fstream::binary);
        if (!file.is_open()) {
            cout << "\nError while loading state, no file " << path << "\n";
            return;
        }

        file.seekg(ios_base::beg);
        file.read((char*)this, sizeof(head));
        file.close();
    }

    void printHead() {
        cout << "\nPrint Head : ";
        cout << "\n k " << k;
        cout << "\n seed_disorder " << seed_disorder;
        cout << "\n seed_system " << seed_system;
        cout << "\n dA " << dA;
        cout << "\n dB " << dB;
        cout << "\n m " << m;
        cout << "\n dt " << dt;
        cout << "\n T " << T;
        cout << "\n N " << N;
        cout << "\n grid_id " << grid_id;
        cout << "\n grid_w " << grid_w;
        cout << "\n grid_h " << grid_h;
    }
};

class options : public head {
    public:
        char job;

        string path;//datapath
        int N_buffer;//buffer
        float dos_crop;

        bool append;
        bool serial;
        bool graphene;
        bool debug;

        //diffusion
        int R;
        char order;
        float t_skip;

        //recorder options
        int frame_skip;
        int frame_w;
        int frame_h;

        int argc;
        char** argv;

        static options* get() {
            static options* singleton_opt = 0;
            if (singleton_opt == 0) singleton_opt = new options();
            return singleton_opt;
        }

        void parse(int _argc, char** _argv) {
            argc = _argc;
            argv = _argv;

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
                ("serial", bpo::value<bool>(), "serial or not")
                ("graphene", bpo::value<bool>(), "graphene or square lattice")
                ("debug", bpo::value<bool>(), "debug mode")
                ("seed_system", bpo::value<int>(), "system seed")
                ("seed_disorder", bpo::value<int>(), "disorder seed")
                ("dos_crop", bpo::value<float>(), "crop the input data before fourier transforming")
                ("dif_r", bpo::value<int>(), "diffusion, max radius to store")
                ("dif_tskip", bpo::value<float>(), "diffusion, betrachtete zeiten")
                ("dif_order", bpo::value<char>(), "diffusion, order parameter to sort data [t]imewise or by [r]adius")
                ("frame_skip", bpo::value<int>(), "steps between two snapshots of the recorder")
                ("frame_w", bpo::value<int>(), "width of recorded frames")
                ("frame_h", bpo::value<int>(), "height of recorded frames")
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
            if (vm.count("job")) job = vm["job"].as<char>();
            else { cout << desc << "\n"; exit(1); }

            //path-----
            if (vm.count("path")) path = vm["path"].as<string>();
            //else { cout << desc << "\n"; exit(1); }

            //load-----
            if (vm.count("append")) append = vm["append"].as<bool>();
            if (vm.count("serial")) serial = vm["serial"].as<bool>();
            if (vm.count("graphene")) graphene = vm["graphene"].as<bool>();
            if (vm.count("debug")) debug = vm["debug"].as<bool>();

            //other parameter
            if (vm.count("total_steps")) T = vm["total_steps"].as<int>();
            if (vm.count("t_step")) dt = vm["t_step"].as<float>();
            if (vm.count("buffer")) N_buffer = vm["buffer"].as<int>();
            if (vm.count("system_size")) k = vm["system_size"].as<int>();
            if (vm.count("krylov_dim")) m = vm["krylov_dim"].as<int>();

            if (vm.count("grid_w")) grid_w = vm["grid_w"].as<int>();
            if (vm.count("grid_h")) grid_h = vm["grid_h"].as<int>();

            if (vm.count("seed_system")) seed_system = vm["seed_system"].as<int>();
            if (vm.count("seed_disorder")) seed_disorder = vm["seed_disorder"].as<int>();
            if (vm.count("dos_crop")) dos_crop = vm["dos_crop"].as<float>();

            //specific defects
            if (vm.count("defects_A")) dA = vm["defects_A"].as<float>();
            if (vm.count("defects_B")) dB = vm["defects_B"].as<float>();


            if (vm.count("dif_r")) R = vm["dif_r"].as<int>();
            if (vm.count("dif_order")) order = vm["dif_order"].as<char>();
            if (vm.count("dif_tskip")) t_skip = vm["dif_tskip"].as<float>();

            if (vm.count("frame_skip")) frame_skip = vm["frame_skip"].as<int>();
            if (vm.count("frame_w")) frame_w = vm["frame_w"].as<int>();
            if (vm.count("frame_h")) frame_h = vm["frame_h"].as<int>();
        }

};



#endif // OPTIONS_H_INCLUDED
