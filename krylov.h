#ifndef KRYLOV_H_INCLUDED
#define KRYLOV_H_INCLUDED

#include "storage.h"
#include "cimgVisual.h"

using namespace std;

class krylovRaum {
    private:
        systm* s;
        int k2;

        grid* gr;

        //Operatoren im krylov raum
        cplx* H_k;
        cplx* J_k;

        //zwischenspeicher der beim gram schmidt berechneten groesen!
        cplx* o_hh_o;//ok
        double* normC;//ok
        cplx* k_h_o;

        cplx mult(cplx* v1, cplx* v2, int N) {
            cplx c(0,0);

            for (int i=0; i<N; i++) c += std::conj(v1[i]) * v2[i];

            if (s->serial) return c;
            else return gr->gatherSum(c);
        }

        double normalize(cplx* v, int N) {
            long double norm;

            norm = real(mult(v,v,N));
            norm = sqrt(norm);

            if (norm < 1e-15) { cout << "\nError! norm of Vector nearly zero in fkt normalize! terminate\n"; exit(-1); }

            for (int i=0;i<N;i++) v[i] /= norm;

            return norm;
        }

        cplx khio(int k, int i) {
            cplx c = o_hh_o[k+i];
            for (int kk = 0; kk<k; kk++) c -= conj(k_h_o[k*s->m+kk]) * k_h_o[i*s->m+kk];
            c *= 1./normC[k];
            return c;
        }

        //copy boundaries in mpi buffer
        void getBoundsFrom(cplx* v) {
            for (int i=0;i<s->k;i++) {
                s->eastbound[s->k+i] = v[i*s->k +s->k -1];
                s->westbound[s->k+i] = v[i*s->k];
                s->northbound[s->k+i] = v[i];
                s->southbound[s->k+i] = v[k2 -s->k +i];
            }
        }

        //Current in direction n
        void JPsi_mpi(cplx* v1, cplx* v2, float nx, float ny, float nz) {
            for (int i=0;i<s->defects.size();i++) {//add defects on old vector
                v1[s->defects[i]] = 0;
            }

            getBoundsFrom(v1);
            gr->getBounds(s);

            float n = sqrt(nx*nx+ny*ny+nz*nz); nx/=n; ny/=n; nz/=n;//normiere n
            int x,y,g;
            cplx t_n,t_s,t_w,t_e;

            t_n = cplx(-ny,0);
            t_s = cplx(ny,0);

            //square lattice (later overwritten if graphene)
            t_w = cplx(nx,0);
            t_e = cplx(-nx,0);

            for (int i=0;i<k2;i++) {
                v2[i] = 0;
                x = i%s->k;
                y = i/s->k;
                g = x + y;

                #ifdef GRAPHENE
                if (g%2 == 1) {//sublatice A
                    t_w = cplx(nx*sqrt(3)/2 + ny/2,0);
                    t_e = cplx(-nx*sqrt(3)/2 + ny/2,0);
                } else {//sublatice B
                    t_w = cplx(nx*sqrt(3)/2 - ny/2,0);
                    t_e = cplx(-nx*sqrt(3)/2 - ny/2,0);
                }
                #endif

                //Innerhalb einer Zeile----------------------------------
                if (x != 0) v2[i] += v1[i-1]*t_w;
                if (x != s->k-1) v2[i] += v1[i+1]*t_e;

                //zwischen den Zeilen------------------------------------
                #ifdef GRAPHENE
                if (g%2 == 0)
                #endif
                if (i < k2-s->k) v2[i] += t_s*v1[i+s->k];

                #ifdef GRAPHENE
                if (g%2 == 1)
                #endif
                if (i >= s->k) v2[i] += t_n*v1[i-s->k];

                //periodizitaet horizontal------- mpi version------------------
                if (x == 0) v2[i] += t_w*s->westbound[y];
                if (x == s->k-1) v2[i] += t_e*s->eastbound[y];

                //periodizitaet vertikal-------
                #ifdef GRAPHENE
                if (g%2 == 0)
                #endif
                if (i >= k2-s->k) v2[i] += t_s*s->southbound[x];
                #ifdef GRAPHENE
                if (g%2 == 1)
                #endif
                if (i < s->k) v2[i] += t_n*s->northbound[x];
            }

            for (int i=0;i<s->defects.size();i++) {//add defects on new vector
                v2[s->defects[i]] = 0;
            }
        }

        void HPsi_mpi(cplx* v1, cplx* v2) {//MPI
            for (int i=0;i<s->defects.size();i++) {//add defects on old vector
                v1[s->defects[i]] = 0;
            }

            getBoundsFrom(v1);
            gr->getBounds(s);

            int g;
            for (int i=0;i<k2;i++) {
                v2[i] = 0;
                int x = i%s->k;
                int y = i/s->k;

                //Innerhalb einer Zeile----------------------------------
                if (x != 0) v2[i] += v1[i-1];
                if (x != s->k-1) v2[i] += v1[i+1];

                //zwischen den Zeilen------------------------------------
                #ifdef GRAPHENE
                g = x + y;
                if (g%2 == 0)
                #endif
                if (i < k2-s->k) v2[i] += v1[i+s->k];

                #ifdef GRAPHENE
                if (g%2 == 1)
                #endif
                if (i >= s->k) v2[i] += v1[i-s->k];

                //periodizitaet horizontal------- mpi version------------------
                if (x == 0) v2[i] += s->westbound[y];
                if (x == s->k-1) v2[i] += s->eastbound[y];

                //periodizitaet vertikal-------
                #ifdef GRAPHENE
                if (g%2 == 0)
                #endif
                if (i >= k2-s->k) v2[i] += s->southbound[x];
                #ifdef GRAPHENE
                if (g%2 == 1)
                #endif
                if (i < s->k) v2[i] += s->northbound[x];
            }

            for (int i=0;i<s->defects.size();i++) {//add defects on new vector
                v2[s->defects[i]] = 0;
            }
        }

        void HPsi_serial(cplx* v1, cplx* v2) {//calc basis vector j+1
            for (int i=0;i<s->defects.size();i++) {//add defects on old vector
                v1[s->defects[i]] = 0;
            }

            int g;
            for (int i=0;i<k2;i++) {
                v2[i] = 0;

                //Innerhalb einer Zeile----------------------------------
                if (i%s->k != 0) v2[i] += v1[i-1];
                if (i%s->k != s->k-1) v2[i] += v1[i+1];

                //periodizitaet horizontal-------
                if (i%s->k == 0) v2[i] += v1[i+s->k-1];
                if (i%s->k == s->k-1) v2[i] += v1[i-s->k+1];


                //zwischen den Zeilen------------------------------------
                #ifdef GRAPHENE
                g = i%s->k + i/s->k;
                if (g%2 == 0)
                #endif
                if (i < k2-s->k) v2[i] += v1[i+s->k];

                #ifdef GRAPHENE
                if (g%2 == 1)
                #endif
                if (i >= s->k) v2[i] += v1[i-s->k];

                //periodizitaet vertikal-------
                #ifdef GRAPHENE
                if (g%2 == 0)
                #endif
                if (i >= k2-s->k) v2[i] += v1[s->k-k2+i];
                #ifdef GRAPHENE
                if (g%2 == 1)
                #endif
                if (i < s->k) v2[i] += v1[k2-s->k+i];
            }

            for (int i=0;i<s->defects.size();i++) {//add defects on new vector
                v2[s->defects[i]] = 0;
            }
        }

        double computeBasis() {
            //normalize v0
            normalize(s->vsys, k2);
            normC[0] = 1;

            //calc v Hv HHv HHHv ...
            if (s->serial) for (int i=0;i<s->m-1;i++) HPsi_serial(s->vsys + i*k2, s->vsys + i*k2 + k2);
            else for (int i=0;i<s->m-1;i++) HPsi_mpi(s->vsys + i*k2, s->vsys + i*k2 + k2);

            //vor gram schmidt!
            for (int i=0;i<s->m;i++) {//OK
                o_hh_o[i] = mult(s->vsys + k2*i, s->vsys, k2);//i = 0...m-1
                //o_hh_o[i+m] = mult(i,m);//i = m...2m-1
                o_hh_o[i+s->m-1] = mult(s->vsys + k2*i, s->vsys + k2*(s->m-1), k2);//i = m...2m-1
            }

            //gram schmidt : b(i) = H i> - |j><j Hi 0> = w(i) - ...
            cplx c; int i1,i2,i3;
            {
                for (i1=0; i1<s->m; i1++) {// calc b(i)
                    for (i2=0; i2<i1; i2++) {
                        c = mult(s->vsys + k2*i2, s->vsys + k2*i1, k2);//<j Hi 0> hier vlt mit der function khio(int k, int i) ?!? könnte viel zeit sparen!!
                        for (i3=0; i3<k2; i3++) s->vsys[i1*k2+i3] -= s->vsys[i2*k2+i3] * c;//|i> -= c|j>
                    }

                    normC[i1] = normalize(s->vsys + k2*i1, k2);
                }
            }

            for (int i=0; i<=s->m; i++) {
                for (int j=0; j<s->m; j++) {//k
                    k_h_o[i*s->m+j] = khio(j, i);
                }
            }

            //teste die orthogonalität und norm
            if (s->debug) {
                cplx c;
                double q = 0;
                for (int i=1;i<s->m;i++) {
                    c = mult(s->vsys + k2*i, s->vsys + k2*(i-1), k2);
                    q+=norm(c);
                }
                q /= s->m;

                return abs(log10(q));
            }
            return 0;
        }

        void computeH_k() {
            for (int i=0;i<s->m;i++) {
                for (int j=0;j<s->m;j++) {
                    int ih = j*s->m+i;

                    H_k[ih] = o_hh_o[i+j+1];

                    for (int k=0;k<i;k++) {
                        H_k[ih] -= conj(k_h_o[i*s->m+k]) * k_h_o[j*s->m+s->m+k];

                        for (int _k=0;_k<j;_k++) {
                            H_k[ih] += conj(k_h_o[i*s->m+k]) * H_k[k*s->m+_k] * k_h_o[j*s->m+_k];
                        }
                    }

                    for (int _k=0;_k<j;_k++) {
                        H_k[ih] -= conj(k_h_o[i*s->m+s->m+_k]) * k_h_o[j*s->m+_k];
                    }

                    H_k[ih] /= normC[i]*normC[j];
                }
            }
        }

        void printStats() {
            long double mem = sizeof(cplx)*k2*s->m;
            long double mem_mb = mem/1048576.0;
            long double mem_gb = mem/1073741824.0;

            cout << "\n\n System : ";
            cout << "\n  size " << s->k << " x " << s->k << " sites";
            cout << "\n  krylov dim " << s->m ;

            cout << "\n\n Memory used : ";
            cout << "\n  " << mem << " bytes\n   = " << mem_mb << " mb\n   = " << mem_gb << " gb\n\n";
        }

        void copy(cplx* v1, cplx* v2) {
            for (int i=0;i<k2;i++) v2[i] = v1[i];
        }

        void conjugate(cplx* v) {
            for (int i=0;i<k2;i++) v[i] = std::conj(v[i]);
        }

    public:
        krylovRaum() {
            H_k = 0;
            J_k = 0;
            normC = 0;
            o_hh_o = 0;
            k_h_o = 0;

            s = 0;
            gr = 0;
        }

        ~krylovRaum() {
            printStats();

            if (normC) delete[] normC;
            if (o_hh_o) delete[] o_hh_o;
            if (k_h_o) delete[] k_h_o;
        }

        void set(systm* _s, grid* _gr) {
            s = _s;
            gr = _gr;
            int m = s->m;
            k2 = s->k*s->k;

            H_k = new cplx[m*m];
            J_k = new cplx[3*m*m];
            normC = new double[m];
            o_hh_o = new cplx[2*m];
            k_h_o = new cplx[m*m+m];
        }

        double process() {
            double d=0;
            d = computeBasis();
            computeH_k();
            return d;
        }

        //rechnet den Krylov vektor zurueck
        void convert(cplx* k_v) {
            int j;
            for (int i = 0; i<k2;  i++) {
                s->vsys[i] *= k_v[0];
                for (j=1 ; j<s->m ; j++) {
                    s->vsys[i] += s->vsys[j*k2+i] * k_v[j];
                }
            }
        }

        cplx* getHamilton() {
            return H_k;
        }

        cplx* getCurrent() {
            return J_k;
        }

        void setRandom(int gridid) {
            srand(s->seed_system + gridid);
            float rm = RAND_MAX;

            for (int i=0;i<k2;i++) {
                float a = rand() - rm/2;
                float b = rand() - rm/2;
                cplx cab = cplx(a,b);
                s->vsys[i] = cab*(1./rm);
            }

            for (int i=0;i<s->defects.size();i++) {//add defects on old vector
                s->vsys[s->defects[i]] = 0;
            }
        }

        void saveState() {
            if (s->v0 == 0) s->v0 = new cplx[k2];

            normalize(s->vsys, k2);
            copy(s->vsys, s->v0);
            conjugate(s->v0);
        }

        cplx mult0() {//mpi
            if (s->v0 == 0) saveState();

            cplx c = 0;
            for (int i=0; i<k2; i++) c += s->v0[i]*s->vsys[i];//v0 ist schon vorab conjugiert worden

            if (s->serial) return c;
            else return gr->gatherSum(c);
        }

        cplx multE() {//mpi multiply with EV
            JPsi_mpi(s->EV, s->vsys + k2);
            cplx c = mult(s->vsys, s->vsys + k2, k2);
            cplx cn = mult(s->v0, s->EV, k2);
            return c*(1/std::norm(cn));
        }

        void do_J_p0(float nx, float ny, float nz) {
            JPsi_mpi(s->vsys, s->vsys+k2, nx, ny, nz);
            copy(s->vsys+k2, s->vsys);
        }

        //gibt in der console den vektor j aus
        void print(cplx* v) {
            cout << "\nprint vector\n";
            for (int i=0;i<k2;i++) cout << v[i] << " ";
        }

        void print0() {
            cout << "\nprint base vector\n";
            for (int i=0;i<k2;i++) cout << s->v0[i] << " ";
        }

        void printH() {
            for (int i2=0; i2<s->m; i2++) {
                cout << "\n";
                for (int i3=0; i3<s->m; i3++) {
                    cout << " " << H_k[i2*s->m + i3];
                }
            }
        }
};

#endif // KRYLOV_H_INCLUDED
