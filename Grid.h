#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "Options.h"

class mpi_grid {
    public:
        //MPI zeugs
        int process_N, process_i, name_len;
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        MPI_Datatype type;//datatype

        //grid breite und hoehe
        int W;
        int H;
        //nachbar nodes
        int west;
        int east;
        int north;
        int south;

        //eigene position im grid
        int gx, gy;

        int getNode(int i, int j) {
            if (i == W) i=0;
            if (i < 0 ) i=W-1;
            if (j == H) j=0;
            if (j < 0 ) j=H-1;
            return i+j*W;
        }

        mpi_grid() {
            options* opt = options::get();

            MPI_Init(&opt->argc, &opt->argv);
            MPI_Comm_size(MPI_COMM_WORLD, &process_N);
            MPI_Comm_rank(MPI_COMM_WORLD, &process_i);
            MPI_Get_processor_name(processor_name, &name_len);

            type = MPI_DOUBLE_COMPLEX;

            cout << "\nGrid " << process_i << " from " << process_N << " initiated\n";

            W = opt->grid_w;
            H = opt->grid_h;

            //row-wise
            gx = process_i%W;
            gy = process_i/W;

            west = getNode(gx-1,gy);
            east = getNode(gx+1,gy);
            north = getNode(gx,gy-1);
            south = getNode(gx,gy+1);

            //wait for all processes
            MPI_Barrier(MPI_COMM_WORLD);
        }

        ~mpi_grid() { ; }

    public:

        static mpi_grid* get() {
            static mpi_grid* singleton_grid = 0;
            options* opt = options::get();
            if (singleton_grid == 0 and opt->serial == 0) singleton_grid = new mpi_grid();
            return singleton_grid;
        }

        void endMPI() {
            MPI_Finalize();
        }

        bool iamroot() {
            return (process_i == 0);
        }

        int id() {
            return process_i;
        }

        void waitforuser() {
            if (iamroot()) cout << "\nPress a key to start...\n";
            if (iamroot()) getchar();
            MPI_Barrier(MPI_COMM_WORLD);
        }

        void getBounds(int k, cplx* northbound, cplx* southbound, cplx* westbound, cplx* eastbound) {//grid zero is in north west
            //cout << "\nProcess " << process_i << " send and receive bounds from " << west << " " << east << " " << north << " " << south << "\n";


            int tag = 0;
            MPI_Status status;
            MPI_Request request[4];
            MPI_Barrier(MPI_COMM_WORLD);

            //get and fill the bounds of the system
            //send and do not wait
            tag = 1;
            MPI_Isend (westbound+k, k, type, east, tag, MPI_COMM_WORLD, &request[0]);//send my west
            MPI_Recv (eastbound, k, type, west, tag, MPI_COMM_WORLD, &status);

            tag = 2;
            MPI_Isend (eastbound+k, k, type, west, tag, MPI_COMM_WORLD, &request[1]);//send my east
            MPI_Recv (westbound, k, type, east, tag, MPI_COMM_WORLD, &status);

            tag = 3;
            MPI_Isend (northbound+k, k, type, south, tag, MPI_COMM_WORLD, &request[2]);//send my north
            MPI_Recv (southbound, k, type, north, tag, MPI_COMM_WORLD, &status);

            tag = 4;
            MPI_Isend (southbound+k, k, type, north, tag, MPI_COMM_WORLD, &request[3]);//send my souths
            MPI_Recv (northbound, k, type, south, tag, MPI_COMM_WORLD, &status);

            //cout << "\nProcess " << process_i << " received all and waits\n";

            MPI_Barrier(MPI_COMM_WORLD);
            //cout << "\nProcess " << process_i << " done waiting\n";

            //check the request objects to free internal memory
            MPI_Waitall(4, request, MPI_STATUSES_IGNORE);
        }

        cplx gatherSum(cplx& c) {
            cplx rbuf[process_N];
            MPI_Gather(&c, 1, type, rbuf, 1,type,0,MPI_COMM_WORLD);

            cplx c_sum = 0;
            if (process_i == 0) for (int i=0;i<process_N;i++) c_sum += rbuf[i];

            MPI_Bcast(&c_sum, 1, type, 0, MPI_COMM_WORLD);
            return c_sum;
        }
};

#endif // GRID_H_INCLUDED
