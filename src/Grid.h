#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "Options.h"

typedef complex<double> cplx;

class mpi_grid {
    public:
        //MPI zeugs
        int process_N, process_i, name_len;
        string processor_name;

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

        int getNode(int i, int j);

        mpi_grid();
        ~mpi_grid();

    public:
        static mpi_grid* get();
        void endMPI();
        bool iamroot();
        int id();
        void waitforuser();
        void getBounds(int k, cplx* northbound, cplx* southbound, cplx* westbound, cplx* eastbound);
        cplx gatherSum(cplx& c);
        void barrier();
};

#endif // GRID_H_INCLUDED
