#ifndef SIMULATOR_H_INCLUDED
#define SIMULATOR_H_INCLUDED

#include "storage.h"
#include "Grid.h"
#include "krylov.h"
#include "Zustandsdichte.h"
#include "cimgVisual.h"
#include "Wavepacket.h"
#include "Timeevolution.h"
#include "Waveplayer.h"

class Simulator {
        options* opt;
        grid* gr;

        krylovRaum K;
        timeEvolution U;

        int id;

        void initRandomSystem(storage* s) {//needs to be checked!
            cout << "\nInit random system";

            K.set(s,gr);                cout << "\nK set\n";
            U.set(s);                   cout << "\nU set\n";

            //disorder
            s->distributeRandomDefects(id + opt->seed_disorder, opt->dA, opt->dB);

            //if the append option is active, the file in path will be loaded and new data appended
            if (opt->append) s->load(opt->path);
            else {//generate new random state
                s->krylov_basis[0].setRandom(id + opt->seed_system);

                K.normalize(s->krylov_basis[0]);
                s->krylov_basis[0].apply_mask(s->defects_mask);
                s->krylov_basis[0].copy(s->initial_state);
            }
        }

        void compute_correlation_function() {
            cout << "\nStart correlation function computation\n";

            storage* s = new storage(opt);
            initRandomSystem(s);

            vector<double>* stats = new vector<double>(5,0);
            timeline* time = 0;
            if (s->opt->debug) {
                time = new timeline(800,300,6);
                time->setVector(stats);
            }

            timer t0;                   cout << "\ntimer set\n";

            for (int i=0,j=0; i<s->T; i++) {

                //projection on root state
                cplx c = s->initial_state.mult(s->krylov_basis[0]);
                gr->gatherSum(c);
                s->dos.buffer()[j] = c;

                cout << "\nSimulation step " << i << " " << c << " on " << id << "\n";

                //compute timestep
                //construct krylov basis and construct hamiltonian in krylov space
                K.process();//MPI -> hamilton ist fuer alle processe gleich, deshalb zeitentwicklung identisch, kein problem
                cplx* Vk = U.evolve(K.getHamilton());//evolve with hamiltonian from krylov space and given timestep
                K.convert(Vk);//write new state back into world space and into the krylov basis


                //flush when needed
                j++;
                if (i>0 and j == opt->N_buffer) {
                    if (id == 0) s->append();//MPI?
                    j = 0;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            cout << "\nEnd correlation function computation\n";
        }

        void compute_dos() {
            cout << "\nStart DOS Simulation\n";

            storage* s = new storage(opt);

            Zustandsdichte dos;
            vector<string> paths = dos.getPaths(opt->path);//get all file names beginning with opt->path
            if (paths.size()==0) cout << "\nWarning! size of paths vector 0!\n";


            s->load(paths[0]);
            dos.set(s, paths.size());

            for (unsigned int i=0;i<paths.size();i++) {
                s->load(paths[i]);
                dos.process(opt->dos_crop);
            }

            dos.saveData("DOS");
            //cout << "\nStart DOS Visualization\n";

            //dos.draw(-3, 3, -0.1, 0.8);
            //dos.draw(-0.2, 0.4, -0.1, 0.2);
            //dos.savePlot();
        }

    public:
        Simulator() {
            opt = 0;
            gr = 0;

            id = 0;
            if (gr) id = gr->id();
        }

        void start(options* o, grid* _gr = 0) {
            opt = o;
            gr = _gr;

            switch (opt->job) {
                case 'c':
                    compute_correlation_function();//system s, number of steps
                    break;
                case 'd':
                    if (id == 0) compute_dos();
                    break;
            }

            gr->endMPI();
        }

};

#endif // SIMULATOR_H_INCLUDED
