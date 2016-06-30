#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include <boost/timer.hpp>

#include <fftw3.h>

using namespace std;
using boost::timer;

//TODO
//graphene visualisierung
//negative zeiten fuer dos
//intervall der anzeige in der config angeben

typedef complex<double> cplx;

float simQ = 0.00001;//guete kriterium
float pi = 4*std::atan(1);

#define GRAPHENE
//#define THREADNUMBER 8//über 30% geschwindigkeit
#define THREADNUMBER 1//keine threads

#include "Simulator.h"
#include "Grid.h"
#include "Options.h"

//MEMORY : 256x256 array is 1 mb!

int main(int argc, char** argv) {
    options* opt = options::get();
    opt->parse(argc, argv);

    Simulator sim;
    sim.start();

    cout << "\nSimulation End\n";
    return 0;
}
