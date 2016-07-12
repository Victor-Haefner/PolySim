#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "Options.h"

typedef complex<double> cplx;

class mpi_grid {
    private:
        //MPI zeugs
        int process_N, process_i, name_len;
        string processor_name;

        //nachbar nodes
        int west;
        int east;
        int north;
        int south;


        int getNode(int i, int j);

        mpi_grid();
        ~mpi_grid();

    public:
        //grid breite und hoehe
        int W, H;

        //eigene position im grid
        int gx, gy;

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
