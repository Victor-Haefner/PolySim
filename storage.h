#ifndef STORAGE_H_INCLUDED
#define STORAGE_H_INCLUDED

#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>

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
};

struct options : public head {
    char job;
    string path;//datapath
    int N_buffer;//buffer
    float dos_crop;

    bool append;
    bool serial;
    bool graphene;
    bool debug;
};

struct storage : public head {
    options* opt;
    bool started;
    string path;

    cplx* v0;//root state
    cplx* vsys;//actual state

    //bounds, 2*k
    cplx* westbound;
    cplx* eastbound;
    cplx* northbound;
    cplx* southbound;

    vector<int> defects;

    cplx* corrfunc;

    int N_buffer;//buffer
    cplx* cf_buffer;

    storage () {
        v0 = 0;
        vsys = 0;
        corrfunc = 0;
        N_buffer = 0;
        cf_buffer = 0;
        westbound = 0;
        eastbound = 0;
        northbound = 0;
        southbound = 0;
        started = false;
    }

    void distributeRandomDefects(int gridid, float def_A, float def_B) {
        int nA=def_A*k*k*0.5;
        int nB=def_B*k*k*0.5;
        int nA_t = 0;
        int nB_t = 0;
        char l=0;

        srand(seed_disorder + gridid);
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

            bool wrong=false;
            do {
                for(int j=0;j<defects.size();j++) {
                    if(defects[j] == dx + dy*k) {
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
                }
            } while(wrong);

            defects.push_back(dx + dy*k);
        }
        cout << "\nNumber of defects : " << defects.size() << " witch are " << 100.0*defects.size()/(k*k) << "%\n";
    }

    void allocate() {
        if(v0 != 0) delete[] v0;
        if(vsys != 0) delete[] vsys;
        if(cf_buffer != 0) delete[] cf_buffer;
        if(corrfunc != 0) delete[] corrfunc;
        if(westbound != 0) delete[] westbound;
        if(eastbound != 0) delete[] eastbound;
        if(northbound != 0) delete[] northbound;
        if(southbound != 0) delete[] southbound;


        if (k>0) {
            v0 = new cplx[k*k];
            vsys = new cplx[k*k*m];
            westbound = new cplx[2*k];
            eastbound = new cplx[2*k];
            northbound = new cplx[2*k];
            southbound = new cplx[2*k];
        }
        if (N_buffer > 0) cf_buffer = new cplx[N_buffer];
        if (N > 0) corrfunc = new cplx[N];
    }

    void set(options* op) {
        head* h = this;
        head* ho = opt;
        (*h) = (*ho);

        opt = op;

        allocate();
    }

    //speichert den aktuellen Zustand, grundzustand, und alle projektionen auf diesen, sowie alles weitere das für die DOS gebraucht wird
    void save() {
        cout << "\nSAVE";

        path = "cf";
        path += getPath();
        writeHead(path);
        ofstream ofile(path.c_str(), fstream::out | fstream::in | fstream::binary);
        ofile.seekp(ios_base::beg + sizeof(head));

        started = true;
        int k2 = k*k;

        //defects
        int kd = defects.size();
        ofile.write((char*)&kd, sizeof(int));

        //v0
        ofile.write((char*)v0, k2*sizeof(cplx));
        //vt
        ofile.write((char*)vsys, k2*sizeof(cplx));
        //defekte
        ofile.write((char*)&defects[0], kd*sizeof(int));
        //korrelation fkt
        ofile.write((char*)cf_buffer, N_buffer*sizeof(cplx));


        ofile.close();

        cout << "\nFootprint : " << footprint() << endl;
    }

    void load(string path) {
        cout << "\nload file\n";
        readHead(path);

        //exit if not present
        ifstream ifile(path.c_str(), fstream::in | fstream::binary);
        if (!ifile.is_open()) {
            cout << "\nError while loading state, no file " << path << "\n";
            exit(0);
        }
        ifile.seekg(sizeof(head));

        started = true;

        //read defects
        int kd;
        ifile.read((char*)&kd, sizeof(int));


        cout << "\nfile contains " << k << " A" << dA << "% B" << dB << "% " << m << " " << N << " " << dt << " " << "\n";
        int k2 = k*k;

        allocate();

        //v0 vt
        if (k2>0) ifile.read((char*)v0, k2*sizeof(cplx));
        if (k2>0) ifile.read((char*)vsys, k2*sizeof(cplx));
        if (kd>0) {
            int def_tmp[kd];
            ifile.read((char*)def_tmp, kd*sizeof(int));
            for (int i=0;i<kd;i++) defects.push_back(def_tmp[i]);
        }
        if (N>0) ifile.read((char*)corrfunc, N*sizeof(cplx));

        ifile.close();

        //cout << "\nfile loaded successfully\n";
        cout << "\nFootprint : " << footprint() << endl;
    }

    void append() {
        int k2 = k*k;
        N += N_buffer;

        head tmp; tmp.readHead(path);
        if (tmp.N == 0 or started == false) {
            save();
            return;
        }
        cout << "\nAppend to " << path;

        //rewrite head to change head change N
        writeHead(path);

        fstream ofile;
        ofile.open(path.c_str(), fstream::in | fstream::out | fstream::binary);

        //change state
        ofile.seekp(ios_base::beg + sizeof(head) + k2*sizeof(cplx));
        ofile.write((char*)vsys, k2*sizeof(cplx));

        //korrelation fkt
        ofile.close();
        ofile.open(path.c_str(), fstream::out | fstream::binary | fstream::app);
        ofile.write((char*)cf_buffer, N_buffer*sizeof(cplx));
        ofile.close();
        cout << "\nFootprint : " << footprint() << endl;
    }

    //hint to help identifiing a system
    cplx footprint() {
        cplx c = 0;

        for (int i=0;i<k*k;i++) {
            c += v0[i];
            c += vsys[i];
        }

        for (int i=0;i<defects.size();i++) {
            c += defects[i];
        }

        return c;
    }
};

//sequence of states
struct sequence {
    int width;
    int height;

    int cx;
    int cy;

    int lenght;
    cplx* frames;

    int kd;
    int* defects;

    string path;

    storage* s;

    sequence() {
        width = 0;
        height = 0;
        cx = 0;//offset, corner
        cy = 0;
        lenght = 0;
        kd = 0;

        frames = 0;
        defects = 0;
    }

    void set(options* opt, storage* _s) {
        s = _s;
        width = s->k;
        height = s->k;
        lenght = s->T;

        if (width > s->k) width = s->k;
        if (height > s->k) height = s->k;

        cx = s->k/2-width/2;
        cy = s->k/2-height/2;

        frames = new cplx[lenght*width*height];

        kd = 0;
        for (int j=0;j<s->defects.size();j++)
            if (s->defects[j] >= cy*s->k+cx and s->defects[j] <= (cy+height)*s->k+cx+width)
                kd++;

        defects = new int[kd];

        int di = 0;
        int dj = 0;
        for (int j=0;j<s->defects.size();j++)
            if (s->defects[j] >= cy*s->k+cx and s->defects[j] <= (cy+height)*s->k+cx+width) {
                dj = s->defects[j];
                defects[di] = dj-(cy*s->k+cx)-(s->k-width)*(dj/s->k-cy);
                di++;
            }
    }

    //speichert einen Ausschnitt vom system
    void save() {
        char ctmp[100];
        sprintf(ctmp, "sys%i.%d.bin", s->k, (int)time(0));
        path = string(ctmp);

        ofstream ofile(path.c_str(), fstream::out | fstream::binary);

        ofile.write((char*)&width, sizeof(int));
        ofile.write((char*)&height, sizeof(int));
        ofile.write((char*)&cx, sizeof(int));
        ofile.write((char*)&cy, sizeof(int));
        ofile.write((char*)&lenght, sizeof(int));


        ofile.write((char*)&kd, sizeof(int));
        ofile.write((char*)defects, kd*sizeof(int));
    }

    void copyData() {
        int x;
        int y;
        int ii;

        for (int i=0; i<height; i++) {//geh durch frame
            for (int j=0; j<width; j++) {
                x = cx+j;
                y = cy+i;
                ii = x*s->k + y;

                frames[i*width+j] = s->vsys[ii];//speicher ausschnitt
            }
        }
    }

    void append() {
        fstream ofile(path.c_str(), fstream::in | fstream::out | fstream::binary | fstream::app);
        ofile.write((char*)frames, width*height*sizeof(cplx));
        ofile.close();
    }

    void load(string path) {
        //exit if not present
        ifstream ifile(path.c_str(), fstream::in | fstream::binary);
        if (!ifile.is_open()) {
            cout << "\nError while loading sequence, no file " << path << "\n";
            exit(0);
        }

        cout << "\nload file\n";

        //read overhead
        ifile.read((char*)&width, sizeof(int));
        ifile.read((char*)&height, sizeof(int));
        ifile.read((char*)&cx, sizeof(int));
        ifile.read((char*)&cy, sizeof(int));
        ifile.read((char*)&lenght, sizeof(int));

        ifile.read((char*)&kd, sizeof(int));
        ifile.read((char*)defects, kd*sizeof(int));

        ifile.read((char*)frames, width*height*lenght*sizeof(cplx));

        ifile.close();

        cout << "\nfile loaded successfully\n";

        cout << "\nLatice size : " << width << " x " << height << "\n";
        cout << "\nTime steps : " << lenght << "\n";
        cout << "\nDefects : " << kd << "\n";
    }
};

//im grid ne fkt die ein state laden kann
struct eigenvector {
    int N;
    cplx* v;
    string path;

    storage* s;

    void allocate() {
        if(v != 0) delete[] v;
        if (k>0) v = new cplx[k*k];
    }

    void writeData() {//baustelle
        s->writeHead(path);

        ofstream ofile(path.c_str(), fstream::in | fstream::out | fstream::binary);
        ofile.seekp(ios_base::beg + sizeof(head));

        //hier zeugs speichern!!!
        ofile.write((char*)&k, sizeof(int));

        ofile.close();
    }

    void readData() {//baustelle
        ;
    }

    void save() {
        path = "ev";
        path += s->getPath();
        writeData();
    }

    void append() {//baustelle
        ;
    }
};

#endif // STORAGE_H_INCLUDED
