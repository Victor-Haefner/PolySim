#ifndef SIMULATOR_H_INCLUDED
#define SIMULATOR_H_INCLUDED

#include "storage.h"
#include "Grid.h"
#include "krylov.h"
#include "Zustandsdichte.h"
#include "Wavepacket.h"
#include "Timeevolution.h"
#include "recorder.h"
#include "Diffusion.h"

class Simulator {
        options* opt;
        mpi_grid* gr;

        krylovRaum K;
        timeEvolution U;

        int id;

        void initRandomSystem(storage* s) {//needs to be checked!
            cout << "\nInit random system\n";

            //if the append option is active, the file in path will be loaded and new data appended
            if (opt->append) s->load(opt->path);

            K.set(s,gr);                cout << "\nK set\n";
            U.set(s);                   cout << "\nU set\n";

            if (!opt->append) {//generate new random state
                s->allocate();
                s->krylov_basis[0].setRandom(id + opt->seed_system);//system
                s->distributeRandomDefects(id + opt->seed_disorder, opt->dA, opt->dB);//disorder

                s->krylov_basis[0].apply_mask(s->defects_mask);
                s->krylov_basis[0].normalize();
                s->initial_state.copy(s->krylov_basis[0]);
            }
        }

        void compute_correlation_function() {
            cout << "\nStart correlation function computation\n";

            storage* s = new storage(opt);
            initRandomSystem(s);

            timer t0;                   cout << "\ntimer set\n";

            for (int i=0,j=0; i<s->T; i++) {

                //projection on root state
                cplx c = s->initial_state.mult(s->krylov_basis[0]);
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

            if (gr) MPI_Barrier(MPI_COMM_WORLD);
            cout << "\nEnd correlation function computation after " << t0.elapsed() << " s\n";
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

            dos.saveDosE("DOS");
        }

        void compute_diffusion_constant() {

            storage* s = new storage(opt);
            K.set(s,gr);                cout << "\nK set\n";
            U.set(s);                   cout << "\nU set\n";

            s->allocate();

            diffusion diff;

            int x0 = s->k/2;

            //delta peak
            //s->krylov_basis[0].setDelta(r0);

            //gaus wave packet
            wavepacket w;
            w.set(s->krylov_basis[0].data(), s->k, x0, x0, 0, 0, 2);


            s->distributeRandomDefects(id + opt->seed_disorder, opt->dA, opt->dB);//disorder
            s->krylov_basis[0].apply_mask(s->defects_mask);
            s->krylov_basis[0].normalize();
            s->initial_state.copy(s->krylov_basis[0]);

            /*------------------end preparation-------------------*/

            diff.set(s, opt, x0, x0);

            //distance from wave origin
            for (int t=0;t<s->T;t++) {//propagate in time
                cout << "\nSim t " << t;

                //cout << " norm " << s->krylov_basis[0].norm();
                //cout << " sqsum " << s->krylov_basis[0].sqsum();

                //compute timestep
                //construct krylov basis and construct hamiltonian in krylov space
                K.process();//MPI -> hamilton ist fuer alle processe gleich, deshalb zeitentwicklung identisch, kein problem
                cplx* Vk = U.evolve(K.getHamilton());//evolve with hamiltonian from krylov space and given timestep
                K.convert(Vk);//write new state back into world space and into the krylov basis


                diff.process(t);
                //diff.getShape(t);
            }

            diff.save();
            //diff.saveShape();
        }

        void propagate_wavepacket() {
            cout << "\nPropagate Wavepacket\n";

            storage* s = new storage(opt);
            K.set(s,gr);                cout << "\nK set\n";
            U.set(s);                   cout << "\nU set\n";

            s->allocate();
            s->distributeRandomDefects(id + opt->seed_disorder, opt->dA, opt->dB);//disorder

            int x0 = s->k/2;
            int y0 = x0;

            //gaus wave packet
            wavepacket w;

            //-----------SCENARIOS-------------------------------------------------------------
            //zerfliessendes gaus packet in der mitte
            w.set(s->krylov_basis[0].data(), s->k, x0, y0, 0, 0, 2);




            //defekt in der mitte------------------------------------------------------------------
            //propagierendes packet
            //w.set(s->krylov_basis[0].data(), s->k, x0-15, y0, 1, 0, 2);
            //s->defects_mask[s->k*y0+x0] = cplx(0,0);
            //---------------------------------------------------------------------------------



            //double split---------------------------------------------------------------------
            /*w.set(s->krylov_basis[0].data(), s->k, x0-15, y0, 1, 0, 2);
            for (int i=0;i<50;i++)
                if (i!=3) { s->defects_mask[s->k*(y0+i)+x0] = cplx(0,0); s->defects_mask[s->k*(y0-i)+x0] = cplx(0,0); }
            for (int i=0;i<15;i++)
                { s->defects_mask[s->k*(y0+i)+x0-30] = cplx(0,0); s->defects_mask[s->k*(y0-i)+x0-30] = cplx(0,0); }
            for (int i=0;i<30;i++)
                { s->defects_mask[s->k*(y0+14)+x0-30+i] = cplx(0,0); s->defects_mask[s->k*(y0-14)+x0-30+i] = cplx(0,0); }*/
            //---------------------------------------------------------------------------------



            s->krylov_basis[0].apply_mask(s->defects_mask);
            s->krylov_basis[0].normalize();

            recorder rec;
            rec.set(s->krylov_basis[0].data(), s->defects_mask.data(), s->k, x0 - opt->frame_w/2, y0 - opt->frame_h/2, x0 + opt->frame_w/2, y0 + opt->frame_h/2);//nehme einf enster der groeße 100x100 auf

            for (int t=0;t<s->T;t++) {//propagate in time
                cout << "\nSim t " << t;

                //compute timestep
                //construct krylov basis and construct hamiltonian in krylov space
                K.process();//MPI -> hamilton ist fuer alle processe gleich, deshalb zeitentwicklung identisch, kein problem
                cplx* Vk = U.evolve(K.getHamilton());//evolve with hamiltonian from krylov space and given timestep
                K.convert(Vk);//write new state back into world space and into the krylov basis

                if (t%opt->frame_skip == 0) rec.take("test.vid");
            }
        }

    public:
        Simulator() {
            opt = 0;
            gr = 0;

            id = 0;
        }

        void start() {
            opt = options::get();
            gr = mpi_grid::get();
            if (gr) id = gr->id();

            switch (opt->job) {
                case 'c':
                    compute_correlation_function();//system s, number of steps
                    break;
                case 'd':
                    if (id == 0) compute_dos();
                    break;
                case 'w':
                    compute_diffusion_constant();
                    break;
                case 'p':
                    propagate_wavepacket();
                    break;
            }

            if (gr) gr->endMPI();
        }

};

#endif // SIMULATOR_H_INCLUDED
