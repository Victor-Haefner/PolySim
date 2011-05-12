#ifndef ZUSTANDSDICHTE_H_INCLUDED
#define ZUSTANDSDICHTE_H_INCLUDED

#include <fstream>
#include <sys/stat.h>

#include "storage.h"

class Zustandsdichte {
    private:

        int M;//stats
        systm* s;

        cplx* ffin;
        cplx* ffout;
        cplx* time;

        graph* gr;//visualisierung

        void writeData(string datapath, cplx* data, long n) {
            ofstream file(datapath.c_str(), ios::out | ios::binary);
            file.write((char*)data, n*sizeof(cplx));
            file.close();
        }

        void readData(string datapath, cplx* data, long n) {
            cout << "\n\nread : " << datapath << "\n\n";
            ifstream file(datapath.c_str(), ios::in | ios::binary);
            file.read((char*)data, n*sizeof(cplx));
            file.close();
            cout << "  done";
        }

        long getDataLenght(string datapath) {
            long n = 0;
            struct stat results;
            if (stat(datapath.c_str(), &results) == 0) n = results.st_size;
            return n;
        }

        void fourier(int N, cplx* in, cplx* out, int dir) {
            fftw_plan p;

            //p = fftw_plan_dft_1d(N, (double(*)[2])in, (double(*)[2])out, FFTW_FORWARD, FFTW_ESTIMATE);//FFTW_FORWARD?------------------------------------?
            p = fftw_plan_dft_1d(N, (double(*)[2])in, (double(*)[2])out, dir, FFTW_ESTIMATE);//FFTW_FORWARD?------------------------------------?

            fftw_execute(p);
            fftw_destroy_plan(p);
        }

        double gaus(double x, float s, double x0) {
            return exp(-(x-x0)*(x-x0)/(2*s*s));
        }

        double time2energy(int i, int N) {
            return ((float)i - N/2.0)*2*pi/(N*s->dt);
        }

        double hanningWindow(int i, int N) {
            return 0.5 + 0.5*cos(pi*i/N);//old, best
            //return 0.5 + 0.5*cos(2*pi*i/s->N+pi);
	    //return 1;
        }

    public:
        Zustandsdichte() {
            //allocate space for fftw
            ffin = 0;
            ffout = 0;
            time = 0;

            gr = 0;
        }

        ~Zustandsdichte() {
            if (ffin) fftw_free(ffin);
            if (ffout) fftw_free(ffout);
            if (time) fftw_free(time);
        }

        void set(systm* _s, int m) {
            s = _s;
            M = m;

            //allocate space for fftw
            ffin = (cplx*) fftw_malloc(sizeof(cplx) * s->N);
            ffout = (cplx*) fftw_malloc(sizeof(cplx) * s->N);
            time = (cplx*) fftw_malloc(sizeof(cplx) * s->N);
            reset();
        }

        void add(int i, cplx& c, int N) {
            cplx tmp = c;

            tmp *= s->dt/pi;//normiere den input
            tmp *= hanningWindow(i, N);

            if (i%2) tmp *= -1;//verschiebt den frequenz nullpunkt in die mitte!

            ffin[i] = tmp;
            time[i] = time2energy(i, N);
        }

        void process(float dos_crop) {
	    if (dos_crop<0) dos_crop=0;
	    if (dos_crop>1) dos_crop=1;
	    int N = s->N*dos_crop;
            for (int i=0;i<N;i++) add(i,s->v0vt[i], N);

            //output start vector
            saveData("DOS_in", ffin);

            cout << "\nStart FFTW\n";

            fourier(N,ffin, ffin, FFTW_BACKWARD);

            for (int i=0;i<N;i++) {
                cplx c = ffin[i];
                c *= 1./M;
                ffout[i] += c;
            }

            cout << "\nOffset Control : " << ffout[0] << " " << s->dt/pi/2 << "\n";
            for (int i=0;i<N;i++) ffout[i] -= ffout[N-1];//subtract offset!
        }

        vector<string> getPaths(string subpath) {
            string dir = "results";
            DIR* dp = opendir(dir.c_str());
            dirent* ep;

            vector<string> vec;

            if (dp != 0) {
                while ((ep = readdir(dp))) {
                    string name = ep->d_name;
                    if (name.find(subpath) == 0) {
                        cout << "\nfound data : " << name;
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

        void testfourier() {
            for (int i=0;i<s->N;i++) {
                cplx c = 0;
                double x = i;
                c += gaus(x,10,0);

                if (i%2) c *= -1;//verschiebt den frequenz nullpunkt in die mitte!

                ffin[i] = c;
                time[i] = time2energy(i, s->N);
            }

            cout << "\nStart FFTW\n";

            fourier(s->N, ffin, ffout, FFTW_BACKWARD);//shortcut when no statistics
            cout << "\nOffset Control : " << ffout[0] << "\n";

            for (int i=0;i<s->N;i++) {//subtract offset!
                ffout[i] -= ffout[s->N-1];
            }

            return;
        }

        void draw(float x1, float x2, float y1, float y2) {
            gr = new graph(800, 600);
            //gr->setVector(s->N, time, ffin, x1, x2, y1, y2, 255, 0, 0);
            gr->setVector(s->N, time, ffout, x1, x2, y1, y2, 255, 0, 0);
            gr->draw();
        }

        void reset() {//setzt in und output auf 0
            for (int i=0;i<s->N;i++) {
                time[i] = 0;
                ffin[i] = 0;
                ffout[i] = 0;
            }
        }

        //save the computed DOS in ascii
        void saveData(string path, cplx* vec = 0) {
            if (vec==0) vec=ffout;

            char ctmp[100];
            sprintf(ctmp, "%ix%i_A%f\%_B%f\%_seedS%i_seedD%i_%d.dat", s->k*s->grid_w, s->k*s->grid_h, s->dA*100, s->dB*100, s->seed_system, s->seed_disorder, (int)std::time(0));
            path += string(ctmp);

            ofstream file(path.c_str());
            //for (int i=0;i<s->N;i++) file << "\n" << real(time[i]) << " " << real(vec[i]);
            for (int i=0;i<s->N;i++) file << "\n" << real(time[i]) << " " << real(vec[i]) << " " << imag(vec[i]) << " " << norm(vec[i]);
        }

        void savePlot() {
            if (gr != 0) gr->saveCanvas("StateDensity.png");
        }
};

#endif // ZUSTANDSDICHTE_H_INCLUDED
