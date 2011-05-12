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

        void initRandomSystem(systm* s) {
            cout << "\nInit random system";

            s->set(opt);
            if (gr) s->distributeRandomDefects(gr->id(), opt->def_A, opt->def_B);
            else s->distributeRandomDefects(0, opt->def_A, opt->def_B);

            if (opt->append) s->load(opt->path);

            K.set(s,gr);                cout << "\nK set\n";
            U.set(s);                   cout << "\nU set\n";

            if (!s->started) {//seed 1299107795
                if (gr) K.setRandom(gr->id());
                else K.setRandom(0);          cout << "\nK random state set\n";
                K.saveState();          cout << "\nK state saved\n";
            }
        }

        void compute_wavepacket() {
            cout << "\nStart wavepacket\n";
            systm* s = new systm;
            sequence* seq = new sequence;
            s->set(opt);

            s->distributeRandomDefects(gr->id(), opt->def_A, opt->def_B);

            seq->set(opt,s);

            krylovRaum K;
            timeEvolution U;
            wavepacket W;

            K.set(s,gr);
            U.set(s);

            //global wavepacket
            int wx = s->k/2 - gr->gx*s->k;
            int wy = s->k/2 - gr->gy*s->k;
            W.set(s->vsys, s->k, wx, wy, pi/2, pi/4, 4);//sys vector, sys size, x0, y0, impuls, phi, packet width

            vector<double>* stats = new vector<double>(5,0);
            timeline* time = 0;
            player* pl = 0;
            if (s->debug) {
                //timeline time(800,300,6);
                time = new timeline(800,300,6);
                time->setVector(stats);
                cout << "\ntimeline set\n";
                pl = new player;
                pl->set(seq);
                cout << "\nplayer set\n";
            }

            seq->save();
            seq->copyData();

            //Simulation
            for (int i=0;i<seq->lenght;i++) {
                cout << "\nSimulation step " << i << " on " << gr->id() << "\n";

                (*stats)[3] = K.process();//construct krylov basis and construct hamiltonian in krylov space
                cplx* Vk = U.evolve(K.getHamilton());//evolve with hamiltonian from krylov space and given timestep
                K.convert(Vk);//write new state back into world space and into the krylov basis

                seq->copyData();//write in buffer
                seq->append();//write to file

                if (s->debug) {
                    (*stats)[1] = U.getDiagQ();
                    (*stats)[0] = U.getQ();

                    time->draw();
                    time->update();

                    pl->draw(0);

                    if (i == 0) gr->waitforuser();
                    cplx tmp = s->footprint();
                    cout << "\nFootprint : " << gr->gatherSum(tmp)/4. << "   , press key to continue\n";
                }
                getchar();
            }

            cout << "\nSimulation end\n";

        }

        void propagate_wavepacket() {
            systm* s = new systm;
            sequence* seq = new sequence;
            s->set(opt);
            seq->set(opt,s);

            seq->load(opt->path);

            cout << "\nStart Wavepacket Propagation\n";
            player pl;
            pl.set(seq);

            pl.play();

            getchar();
            exit(0);
        }

        void compute_correlation_function() {
            cout << "\nStart correlation function computation\n";

            systm* s = new systm;
            initRandomSystem(s);

            vector<double>* stats = new vector<double>(5,0);
            timeline* time = 0;
            if (s->debug) {
                time = new timeline(800,300,6);
                time->setVector(stats);
            }

            timer t0;                   cout << "\ntimer set\n";

            for (int i=0,j=0; i<s->T; i++) {

                //projection on root state
                s->v0vt_buffer[j] = K.mult0();//MPI?
                if (gr) cout << "\nSimulation step " << i << " " << s->v0vt_buffer[j] << " on " << gr->id() << "\n";
                else cout << "\nSimulation step " << i << " " << s->v0vt_buffer[j] << " on " << 0 << "\n";

                //compute timestep
                //construct krylov basis and construct hamiltonian in krylov space
                (*stats)[3] = K.process();//MPI -> hamilton ist fuer alle processe gleich, deshalb zeitentwicklung identisch, kein problem
                cplx* Vk = U.evolve(K.getHamilton());//evolve with hamiltonian from krylov space and given timestep
                K.convert(Vk);//write new state back into world space and into the krylov basis


                //flush when needed
                j++;
                if (i>0 and j == s->N_buffer) {
                    if (gr->iamroot()) s->append();//MPI?
                    j = 0;
                }

                if (s->debug) {
                    (*stats)[1] = U.getDiagQ();
                    (*stats)[0] = U.getQ();

                    time->draw();
                    time->update();

                    cplx tmp = s->footprint();
                    cout << "\nFootprint : " << gr->gatherSum(tmp) << "\n";
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            cout << "\nEnd correlation function computation\n";
        }

        void compute_dos() {
            cout << "\nStart DOS Simulation\n";

            systm* s = new systm;
            s->set(opt);

            Zustandsdichte dos;
            vector<string> paths = dos.getPaths(opt->path);
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

        void compute_EV() {
            cout << "\nStart DC computation\n";

            systm* s = new systm;
            initRandomSystem(s);

            //512x512 states are 4 mb big, this means with 512 memory I can do easily 100 states
            // per node! and flush them at the end, or with a buffer inbetween

            for (int i=0,j=0; i<s->T; i++) {

                //stats
                if (gr) cout << "\nSimulation step " << i << " " << s->v0vt_buffer[j] << " on " << gr->id() << "\n";
                else cout << "\nSimulation step " << i << " " << s->v0vt_buffer[j] << " on " << 0 << "\n";

                //compute timestep
                //construct krylov basis and construct hamiltonian in krylov space
                K.process();//MPI -> hamilton ist fuer alle processe gleich, deshalb zeitentwicklung identisch, kein problem
                cplx* Vk = U.evolve(K.getHamilton());//evolve with hamiltonian from krylov space and given timestep
                K.convert(Vk);//write new state back into world space and into the krylov basis
            }
            MPI_Barrier(MPI_COMM_WORLD);
            cout << "\nEnd eigenvector computation\n";
        }

        void compute_DC() {
            cout << "\nStart DC computation\n";

            systm* s = new systm;
            initRandomSystem(s);

            //J_x on system -> vector from options?
            K.do_J_p0(1,0,0);

            for (int i=0,j=0; i<s->T; i++) {

                //projection on root state
                s->v0vt_buffer[j] = K.multE();//MPI?

                //stats
                if (gr) cout << "\nSimulation step " << i << " " << s->v0vt_buffer[j] << " on " << gr->id() << "\n";
                else cout << "\nSimulation step " << i << " " << s->v0vt_buffer[j] << " on " << 0 << "\n";

                //compute timestep
                //construct krylov basis and construct hamiltonian in krylov space
                (*stats)[3] = K.process();//MPI -> hamilton ist fuer alle processe gleich, deshalb zeitentwicklung identisch, kein problem
                cplx* Vk = U.evolve(K.getHamilton());//evolve with hamiltonian from krylov space and given timestep
                K.convert(Vk);//write new state back into world space and into the krylov basis


                //flush when needed
                j++;
                if (i>0 and j == s->N_buffer) {
                    if (gr->iamroot()) s->append();//MPI?
                    j = 0;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            cout << "\nEnd DC computation\n";
        }

    public:
        Simulator() {
            opt = 0;
            gr = 0;
        }

        void start(options* o, grid* _gr = 0) {
            opt = o;
            gr = _gr;

            switch (opt->job) {
                case 'w':
                    compute_wavepacket();
                    break;
                case 'v':
                    if (gr->iamroot()) propagate_wavepacket();
                    break;
                case 'c':
                    compute_correlation_function();//system s, number of steps
                    break;
                case 'd':
                    if (gr) {
                        if (gr->iamroot()) compute_dos();
                    } else compute_dos();
                    break;
                case 'e':
                    compute_EV();
                    break;
                case 't':
                    compute_DC();
                    break;
            }

            gr->endMPI();
        }

};

#endif // SIMULATOR_H_INCLUDED
