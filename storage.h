#ifndef STORAGE_H_INCLUDED
#define STORAGE_H_INCLUDED

#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>
#include "Grid.h"

//implement addition, substraction, multiplication, etc..
class state {
    public:
        cplx* v;
        int N;

        mpi_grid* grid;

        state() {
            v = 0;
            N = 0;
            grid = mpi_grid::get();
        }

        void allocate(int n) {
            N = n;
            v = new cplx[n];
            for (int i=0;i<N;i++) v[i] = 0;
        }

        void kill() { delete[] v; }

        void save(string path, int offset) {
            ofstream ofile(path.c_str(), fstream::in | fstream::out | fstream::binary);
            ofile.seekp(ios_base::beg + offset);
            ofile.write((char*)v, N*sizeof(cplx));
            ofile.close();
        }

        void load(string path, int offset) {
            //exit if not present
            ifstream ifile(path.c_str(), fstream::in | fstream::binary);
            if (!ifile.is_open()) {
                cout << "\nError while loading state, no file " << path << "\n";
                exit(0);
            }

            ifile.seekg(offset);
            ifile.read((char*)v, N*sizeof(cplx));
            ifile.close();
        }

        void setRandom(int seed) {
            srand(seed);
            float rm = RAND_MAX;

            for (int i=0;i<N;i++) {
                float a = rand() - rm/2;
                float b = rand() - rm/2;
                cplx cab = cplx(a,b);
                v[i] = cab*(1./rm);
            }
        }

        void setDelta(int i) {
            for (int j=0;j<N;j++) {
                v[j] = 0;
            }
            v[i] = 1;
        }

        void apply_mask(state mask) {
            for (int i=0;i<N;i++) v[i] = v[i]*mask.v[i];
        }

        int size() {return N;}

        cplx* data() {return v;}

        cplx& operator[](int i) { return v[i]; }


        //--------------MATH-------------------

        cplx mult(state _v) {
            if (_v.N != N) cout << "\nWarning! in fkt mult, sizes of vectors don't match!\n";

            cplx c(0,0);

            for (int i=0; i<N; i++) c += std::conj(v[i]) * _v.v[i];

            if (grid != 0) return grid->gatherSum(c);
            else return c;
        }

        cplx norm() {
            return mult(*this);
        }

        double normalize() {
            double n = sqrt(real(norm()));
            if (n < 1e-15) cout << "\nWarning! norm of Vector nearly zero in fkt normalize!\n";
            for (int i=0;i<N;i++) v[i] /= n;
            return n;
        }

        cplx sqsum() {
            cplx c = cplx(0,0);
            for (int i=0;i<N;i++) c += std::norm(v[i]);

            if (grid != 0) return grid->gatherSum(c);
            else return c;
        }

        void copy(state _v) {
            if (_v.N != N) cout << "\nWarning! in fkt mult, sizes of vectors don't match!\n";
            for (int i=0;i<N;i++) v[i] = _v.v[i];
        }

        void conjugate() { for (int i=0;i<N;i++) v[i] = std::conj(v[i]); }

};

struct corrFkt {
    int N;
    cplx* corrfunc;

    int N_buffer;//buffer
    cplx* cf_buffer;


    void allocate(int n, int n_b) {
        N = n;
        N_buffer = n_b;

        corrfunc = new cplx[n];
        cf_buffer = new cplx[n_b];

        for (int i=0;i<N;i++) corrfunc[i] = 0;
        for (int i=0;i<N_buffer;i++) cf_buffer[i] = 0;
    }

    void kill() { delete[] corrfunc; delete[] cf_buffer; }

    void append(string path) {
        ofstream ofile(path.c_str(), fstream::out | fstream::binary | fstream::app);
        ofile.write((char*)cf_buffer, N_buffer*sizeof(cplx));
        ofile.close();
    }

    void load(string path, int offset) {
        //exit if not present
        ifstream ifile(path.c_str(), fstream::in | fstream::binary);
        if (!ifile.is_open()) {
            cout << "\nError while loading state, no file " << path << "\n";
            exit(0);
        }

        ifile.seekg(offset);
        ifile.read((char*)corrfunc, N*sizeof(cplx));
        ifile.close();
    }

    int size() {return N;}
    int bufferSize() {return N_buffer;}

    cplx* data() {return corrfunc;}
    cplx* buffer() {return cf_buffer;}

    cplx& operator[](int i) { return corrfunc[i]; }
};

//everything for the krylov basis
struct storage : public head {
    options* opt;
    bool saved;
    string path;

    state initial_state;
    state defects_mask;
    vector<state> krylov_basis;

    corrFkt dos;

    //mpi lattice bounds, 2*k
    cplx* westbound;
    cplx* eastbound;
    cplx* northbound;
    cplx* southbound;


    storage (options* op) {
        westbound = 0;
        eastbound = 0;
        northbound = 0;
        southbound = 0;

        saved = false;

        head* h = this;
        head* ho = op;
        (*h) = (*ho);

        opt = op;
    }

    void distributeRandomDefects(int seed, float def_A, float def_B) {
        int nA=def_A*k*k*0.5;
        int nB=def_B*k*k*0.5;
        int nA_t = 0;
        int nB_t = 0;
        char l=0;

        for (int i=0;i<defects_mask.size();i++) defects_mask[i] = cplx(1,0);

        srand(seed);
        int d_n = 0;
        for (int i=0;i<nA+nB;i++) {
            int dx = rand()%k;
            int dy = rand()%k;

            if ((dx%2 + dy%2)%2) l='A';
            else l='B';

            switch(l) {
                case 'A':
                    if (nA_t<nA) nA_t++;
                    else {dx++;nB_t++;}
                    break;
                case 'B':
                    if (nB_t<nB) nB_t++;
                    else {dx++;nA_t++;}
                    break;
            }

            if (dx >= k) dx -= k - k%2;

            cplx u = cplx(0,0);

            bool wrong=false;
            do {
                if(defects_mask[dx + dy*k] == u) {
                    wrong=true;
                    dx+=2;
                    if (dx >= k) {
                        dx -= k - k%2 -1;
                        if (dx == 2) dx = 0;
                        dy++;
                    }
                    if (dy >= k) dy -= k;
                    break;
                } else wrong=false;
            } while(wrong);

            defects_mask[dx + dy*k] = cplx(0,0); d_n++;
        }
        cout << "\nNumber of defects : " << d_n << " witch are " << 100.0*d_n/(k*k) << "%\n";
    }

    void allocate() {
        if (k>0) {
            initial_state.allocate(k*k);
            defects_mask.allocate(k*k);
            for (int i=0;i<m;i++) krylov_basis.push_back(state());
            for (int i=0;i<m;i++) krylov_basis[i].allocate(k*k);

            westbound = new cplx[2*k];
            eastbound = new cplx[2*k];
            northbound = new cplx[2*k];
            southbound = new cplx[2*k];
        }

        dos.allocate(N, opt->N_buffer);
    }

    void deallocate() {
        initial_state.kill();
        defects_mask.kill();
        for (int i=0;i<m;i++) krylov_basis[i].kill();
        krylov_basis = vector<state>();

        delete[] westbound;
        delete[] eastbound;
        delete[] northbound;
        delete[] southbound;

        dos.kill();
    }

    //speichert den aktuellen Zustand, grundzustand, und korr funktion, sowie alles weitere das für die DOS gebraucht wird
    void initial_write() {
        cout << "\nSave";

        path = "results/cf";
        path += getPath();

        //create file!
        ofstream of(path.c_str(), fstream::out | fstream::binary);
        of.close();

        writeHead(path);

        saved = true;
        int k2 = k*k;

        initial_state.save(path, ios_base::beg + sizeof(head));
        defects_mask.save(path, ios_base::beg + sizeof(head) + k2*sizeof(cplx));
        krylov_basis[0].save(path, ios_base::beg + sizeof(head) + 2*k2*sizeof(cplx));
    }

    void append() {
        int k2 = k*k;
        N += opt->N_buffer;

        if (!saved) initial_write();
        cout << "\nAppend to " << path;

        //rewrite head to change N
        writeHead(path);
        krylov_basis[0].save(path, ios_base::beg + sizeof(head) + 2*k2*sizeof(cplx));
        dos.append(path);
    }

    void load(string _path) {
        if (saved) deallocate();

        path = _path;
        cout << "\nload file " << path << "\n";
        readHead(path);
        //printHead();

        allocate();

        saved = true;
        int k2 = k*k;

        initial_state.load(path, ios_base::beg + sizeof(head));
        defects_mask.load(path, ios_base::beg + sizeof(head) + k2*sizeof(cplx));
        krylov_basis[0].load(path, ios_base::beg + sizeof(head) + 2*k2*sizeof(cplx));
        dos.load(path, ios_base::beg + sizeof(head) + 3*k2*sizeof(cplx));
    }
};

#endif // STORAGE_H_INCLUDED
