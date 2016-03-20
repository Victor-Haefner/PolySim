#ifndef CURRENT_H_INCLUDED
#define CURRENT_H_INCLUDED

#include "storage.h"
#include "GenericBinaryIO.h"

using namespace std;

class P_rt : private binIO {
    private:
        int rN;
        int tN;
        float dt;
        float t_skip;
        int N;
        bool order;//the way the data is sorted, true means t wise

        double* data;//this pointer has to be the last member of this class for working IO operations
        int* rads;

    public:
        P_rt() : order(true) {data = 0; rads = 0; t_skip = 5; }

        void allocate(int _N, int _rN, int _tN, float _dt, float _t_skip) {
            //cout << "\nallocate : " << _t_skip << " " << _rN << " " << _tN << " " << _dt << " " << data << " " << rads << " " << endl;

            if (data != 0) delete[] data;
            if (rads != 0) delete[] rads;

            N = _N;
            rN = _rN;
            tN = _tN;
            dt = _dt;
            t_skip = _t_skip;
            data = new double[N];
            rads = new int[rN];
            for (int i=0;i<N;i++) data[i] = 0;
            for (int i=0;i<rN;i++) rads[i] = 0;
        }

        void setRads(vector<int> _r) { memcpy(rads, &(_r[0]), sizeof(int)*rN); }

        void save(string path) {
            write(path, sizeof(P_rt) - sizeof(double*) - sizeof(int*));
            append(path, N*sizeof(double), (char*)data);
        }

        void load(string path) {
            cout << "\n load " << path << endl;
            read(path, sizeof(P_rt) - sizeof(double*) - sizeof(int*));
            allocate(N, rN, tN, dt, options::get()->t_skip);
            extract(path, sizeof(P_rt) - sizeof(double*) - sizeof(int*), N*sizeof(double), (char*)data);
        }

        void sort(bool b) {//the way the data is sorted, true means t wise
            if (order == b) return;
            order = b;

            double* tmp = new double[N];

            for (int i=0;i<tN;i++) {
                for (int j=0;j<rN;j++) {
                    if(b) tmp[j*tN+i] = data[i*rN+j];
                    else tmp[i*rN+j] = data[j*tN+i];
                }
            }

            for (int i=0;i<N;i++) data[i] = tmp[i];
            delete[] tmp;
        }

        vector<double> getP_t(int r) {
            vector<double> vec;
            sort(true);
            for (int i=r*tN;i<(r+1)*tN;i++) vec.push_back(data[i]);

            return vec;
        }

        vector<double> getP_r(int t) {
            vector<double> vec;
            sort(false);
            for (int i=t*rN;i<(t+1)*rN;i++) vec.push_back(data[i]);

            return vec;
        }

        void writeAscii(string path) {
            cout << "\nsave to " << path << endl;
            ofstream file(path.c_str());

            for (int i=0;i<N;i++) {
                if(order) {// t wise
                    if (i%tN == 0 and i>0) file << "\n\n";
                    file << (i%tN)*dt << " " << data[i] << endl;
                } else {// r wise
                    if (i%rN == 0 and i>0) file << "\n\n";
                    int _t = round(t_skip/dt);
                    if ((i/rN)%_t == 0 ) {
                        //cout << "\n tskip " << t_skip << " dt " << dt << " tn " << i/rN << " _t " << _t;
                        file << sqrt(rads[i%rN]) << " " << data[i] << endl;
                    }
                }
            }

            file.close();
        }

        double& operator[](int i) { return data[i]; }

        int size() { return N; }

        void print() {
            cout << "\n Prt : " << N << " " << rN << " " << tN << " " << dt << " " << data << " " << rads << " " << endl;
        }
};

class diffusion {
    private:
        storage* s;
        options* opt;

        //double* dvec;
        P_rt prt;

        int x0, y0, rmax;

        //maps to all possible integer distances
        vector<int> Rmap;

        //maps to every radius a set of 2d coordinates
        vector< vector< int* > > Pmap;

        void fillRmap() {//all possible R values
            for (int i=-rmax; i<=rmax; i++) {
                for (int j=-rmax; j<=rmax; j++) {
                    int r = getIntegerDistance(i+x0,j+y0);
                    if(find(Rmap.begin(), Rmap.end(), r) == Rmap.end())
                        Rmap.push_back(r);
                }
            }

            sort(Rmap.begin(), Rmap.end());
        }

        void fillPmap() {
            fillRmap();

            Pmap.resize(Rmap.size(), vector<int*>());

            for (int i=-rmax; i<=rmax; i++) {
                for (int j=-rmax; j<=rmax; j++) {
                    int r = getIntegerDistance(i+x0,j+y0);
                    int id = int(find(Rmap.begin(), Rmap.end(), r) - Rmap.begin());
                    int* pos = new int[2];
                    pos[0] = i+x0; pos[1] = j+y0;
                    Pmap[id].push_back(pos);
                }
            }
        }

        bool isSameSublatice(int x, int y) { return (x+y-x0-y0)%2 == 0; }

        //returns r*r
        // r is the distance from the site xy to x0y0
        int getIntegerDistance(int x, int y) {
            int X = x-x0;
            int Y = y-y0;
            bool ssl = isSameSublatice(x,y);

            if (ssl) return (3*X*X + 9*Y*Y)/4;
            else {
                if ((y0 + x0)%2 == 0) {
                    return (3*X*X + 9*Y*Y + 1 - 6*Y)/4;
                } else {
                    return (3*X*X + 9*Y*Y + 1 + 6*Y)/4;
                }
            }
        }

        //returns the distance from the site xy to x0y0
        float getDistance(int x, int y) {
            int r = getIntegerDistance(x,y);

            return sqrt(r);
        }

        vector<string> getPaths(string subpath) {
            string dir = subpath;

            size_t pos_slash = dir.find_last_of('/');
            subpath = dir.substr(pos_slash+1, 1000);
            size_t beg = 0;
            size_t end = pos_slash;
            dir = dir.substr(beg, end);
            cout << "\nLooking in path : " << dir << " for " << subpath;

            DIR* dp = opendir(dir.c_str());
            dirent* ep;

            vector<string> vec;

            if (dp != 0) {
                while ((ep = readdir(dp))) {
                    string name = ep->d_name;
                    if (name.find(subpath) == 0) {
                        string tmp = dir;
                        tmp += '/';
                        tmp += name;
                        vec.push_back(tmp);
                    }
                }

                closedir(dp);
            }

            return vec;
        }

    public:
        diffusion() {
        }

        void set(storage* _s, options* _opt, int _x0, int _y0) {
            s = _s;
            opt = _opt;

            x0 = _x0;
            y0 = _y0;

            rmax = opt->R;

            fillPmap();
            prt.allocate(opt->T*Rmap.size(), Rmap.size(), opt->T, opt->dt, opt->t_skip);
            prt.setRads(Rmap);

            //debug output
            /*cout << "\nMAPS\n";
            for (int i=0;i<Rmap.size();i++) cout << " " << Rmap[i] << " ";
            cout << "\n";
            for (int i=0;i<Pmap.size();i++) {
                cout << "\n";
                for (int j=0;j<Pmap[i].size();j++)
                    cout << " " << Pmap[i][j][0] - x0 << "," << Pmap[i][j][1] - y0 << " ";
            }*/
        }

        void process(int t) {
            for (int i=0; i<Pmap.size(); i++) {
                for (int j=0; j<Pmap[i].size(); j++) {
                    int id = Pmap[i][j][0] + Pmap[i][j][1]*s->k;
                    cplx c = s->krylov_basis[0][id];
                    prt[t + i*s->T] += norm(c)/Pmap[i].size();
                }
            }
        }

        void average(string path) {
            cout << "\naverage over diffusion data\n";
            vector<string> paths = getPaths(path);
            int N = paths.size();
            if (N == 0) {
                cout << "\nNo files found\n";
                return;
            }

            prt.load(paths[0]);
            prt.setRads(Rmap);

            P_rt tmp;
            for (int i=1; i<N; i++) {
                tmp.load(paths[i]);
                for (int j=0; j<prt.size();j++) prt[j] += tmp[j];
            }

            //generate path string
            string ptmp = "avg_";
            string dir = paths[0];
            size_t pos_slash = dir.find_last_of('/');
            dir = dir.substr(pos_slash+1, 1000);
            ptmp.append(dir);

            //prt.save(ptmp);
            if (opt->order == 't') prt.sort(true);
            if (opt->order == 'r') prt.sort(false);
            prt.writeAscii(ptmp);
        }

        void save() {
            string path = "D_rt_";
            path += s->getPath();
            prt.save(path);
        }
};

#endif // CURRENT_H_INCLUDED
