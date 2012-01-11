#ifndef KRYLOV_H_INCLUDED
#define KRYLOV_H_INCLUDED

#include "storage.h"

using namespace std;

class krylovRaum {
    private:
        storage* s;
        int k2;

        mpi_grid* gr;

        //Operatoren im krylov raum
        cplx* H_k;
        //cplx* J_k;

        //zwischenspeicher der beim gram schmidt berechneten groesen!
        cplx* o_hh_o;//ok
        double* normC;//ok
        cplx* k_h_o;

        cplx khio(int k, int i) {
            cplx c = o_hh_o[k+i];
            for (int j = 0; j<k; j++) c -= conj(k_h_o[k*s->m+j]) * k_h_o[i*s->m+j];
            c *= 1./normC[k];
            return c;
        }

        //copy boundaries in mpi buffer
        void getBoundsFrom(state v) {
            for (int i=0;i<s->k;i++) {
                s->eastbound[s->k+i] = v[i*s->k +s->k -1];
                s->westbound[s->k+i] = v[i*s->k];
                s->northbound[s->k+i] = v[i];
                s->southbound[s->k+i] = v[k2 -s->k +i];
            }
        }

        //Current in direction n, never used/tested
        void JPsi_mpi(state& src, state& res, float nx, float ny, float nz) {
            src.apply_mask(s->defects_mask);

            getBoundsFrom(src);
            gr->getBounds(s->k, s->northbound, s->southbound, s->westbound, s->eastbound);

            float n = sqrt(nx*nx+ny*ny+nz*nz); nx/=n; ny/=n; nz/=n;//normiere n
            int x,y,g;
            cplx t_n,t_s,t_w,t_e;

            t_n = cplx(-ny,0);
            t_s = cplx(ny,0);

            //square lattice (later overwritten if graphene)
            t_w = cplx(nx,0);
            t_e = cplx(-nx,0);

            for (int i=0;i<k2;i++) {
                res[i] = 0;
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
                if (x != 0) res[i] += src[i-1]*t_w;
                if (x != s->k-1) res[i] += src[i+1]*t_e;

                //zwischen den Zeilen------------------------------------
                #ifdef GRAPHENE
                if (g%2 == 0)
                #endif
                if (i < k2-s->k) res[i] += t_s*src[i+s->k];

                #ifdef GRAPHENE
                if (g%2 == 1)
                #endif
                if (i >= s->k) res[i] += t_n*src[i-s->k];

                //periodizitaet horizontal------- mpi version------------------
                if (x == 0) res[i] += t_w*s->westbound[y];
                if (x == s->k-1) res[i] += t_e*s->eastbound[y];

                //periodizitaet vertikal-------
                #ifdef GRAPHENE
                if (g%2 == 0)
                #endif
                if (i >= k2-s->k) res[i] += t_s*s->southbound[x];
                #ifdef GRAPHENE
                if (g%2 == 1)
                #endif
                if (i < s->k) res[i] += t_n*s->northbound[x];
            }

            res.apply_mask(s->defects_mask);
        }

        void HPsi_mpi_graphene(state& src, state& res) {//MPI
            src.apply_mask(s->defects_mask);
            float c_r, c_i;

            getBoundsFrom(src);
            gr->getBounds(s->k, s->northbound, s->southbound, s->westbound, s->eastbound);

            int g;
            for (int i=0;i<k2;i++) {
                res[i] = 0;

                //c_r = real(src[i]);
                //c_i = imag(src[i]);
                //if (c_r == 0 and c_i == 0) continue;

                int x = i%s->k;
                int y = i/s->k;

                //Innerhalb einer Zeile----------------------------------
                if (x != 0) res[i] += src[i-1];
                if (x != s->k-1) res[i] += src[i+1];

                //zwischen den Zeilen------------------------------------
                g = x + y;
                if (g%2 == 0)
                if (i < k2-s->k) res[i] += src[i+s->k];

                if (g%2 == 1)
                if (i >= s->k) res[i] += src[i-s->k];

                //periodizitaet horizontal------- mpi version------------------
                if (x == 0) res[i] += s->westbound[y];
                if (x == s->k-1) res[i] += s->eastbound[y];

                //periodizitaet vertikal-------
                if (g%2 == 0)
                if (i >= k2-s->k) res[i] += s->southbound[x];

                if (g%2 == 1)
                if (i < s->k) res[i] += s->northbound[x];
            }

            res.apply_mask(s->defects_mask);
        }

        void HPsi_mpi_square(state& src, state& res) {//MPI
            src.apply_mask(s->defects_mask);
            float c_r, c_i;

            getBoundsFrom(src);
            gr->getBounds(s->k, s->northbound, s->southbound, s->westbound, s->eastbound);

            int g;
            for (int i=0;i<k2;i++) {
                res[i] = 0;

                //c_r = real(src[i]);
                //c_i = imag(src[i]);
                //if (c_r == 0 and c_i == 0) continue;

                int x = i%s->k;
                int y = i/s->k;

                //Innerhalb einer Zeile----------------------------------
                if (x != 0) res[i] += src[i-1];
                if (x != s->k-1) res[i] += src[i+1];

                //zwischen den Zeilen------------------------------------
                if (i < k2-s->k) res[i] += src[i+s->k];
                if (i >= s->k) res[i] += src[i-s->k];

                //periodizitaet horizontal------- mpi version------------------
                if (x == 0) res[i] += s->westbound[y];
                if (x == s->k-1) res[i] += s->eastbound[y];

                //periodizitaet vertikal-------
                if (i >= k2-s->k) res[i] += s->southbound[x];
                if (i < s->k) res[i] += s->northbound[x];
            }

            res.apply_mask(s->defects_mask);
        }

        void HPsi_serial_graphene(state& src, state& res) {//calc basis vector j+1
            src.apply_mask(s->defects_mask);
            float c_r, c_i;

            int g;
            int x,y;
            for (int i=0;i<k2;i++) {
                res[i] = 0;

                x = i%s->k;
                y = i/s->k;

                //Innerhalb einer Zeile----------------------------------
                if (x != 0) res[i] += src[i-1];
                if (x != s->k-1) res[i] += src[i+1];

                //periodizitaet horizontal-------
                if (x == 0) res[i] += src[i+s->k-1];
                if (x == s->k-1) res[i] += src[i-s->k+1];

                //zwischen den Zeilen------------------------------------
                g = x + y;
                if (g%2 == 0)
                if (i < k2-s->k) res[i] += src[i+s->k];

                if (g%2 == 1)
                if (i >= s->k) res[i] += src[i-s->k];

                //periodizitaet vertikal-------
                if (g%2 == 0)
                if (i >= k2-s->k) res[i] += src[s->k-k2+i];

                if (g%2 == 1)
                if (i < s->k) res[i] += src[k2-s->k+i];
            }

            res.apply_mask(s->defects_mask);
        }

        void HPsi_serial_square(state& src, state& res) {//calc basis vector j+1
            src.apply_mask(s->defects_mask);
            float c_r, c_i;

            int g;
            int x,y;
            for (int i=0;i<k2;i++) {
                res[i] = 0;

                x = i%s->k;
                y = i/s->k;

                //Innerhalb einer Zeile----------------------------------
                if (x != 0) res[i] += src[i-1];
                if (x != s->k-1) res[i] += src[i+1];

                //periodizitaet horizontal-------
                if (x == 0) res[i] += src[i+s->k-1];
                if (x == s->k-1) res[i] += src[i-s->k+1];

                //zwischen den Zeilen------------------------------------
                if (i < k2-s->k) res[i] += src[i+s->k];
                if (i >= s->k) res[i] += src[i-s->k];

                //periodizitaet vertikal-------
                if (i >= k2-s->k) res[i] += src[s->k-k2+i];
                if (i < s->k) res[i] += src[k2-s->k+i];
            }

            res.apply_mask(s->defects_mask);
        }

        double computeBasis() {

            //normalize(s->krylov_basis[0]);
            normC[0]=1;


            //calc v Hv HHv HHHv ...
            for (int i=0;i<s->m-1;i++) {
                if ( s->opt->serial and  s->opt->graphene) HPsi_serial_graphene(s->krylov_basis[i], s->krylov_basis[i+1]);
                if (!s->opt->serial and  s->opt->graphene) HPsi_mpi_graphene(s->krylov_basis[i], s->krylov_basis[i+1]);
                if ( s->opt->serial and !s->opt->graphene) HPsi_serial_square(s->krylov_basis[i], s->krylov_basis[i+1]);
                if (!s->opt->serial and !s->opt->graphene) HPsi_mpi_square(s->krylov_basis[i], s->krylov_basis[i+1]);
            }

            //vor gram schmidt!
            for (int i=0;i<s->m;i++) {//OK
                o_hh_o[i] = s->krylov_basis[i].mult(s->krylov_basis[0]);
                o_hh_o[i+s->m-1] = s->krylov_basis[i].mult(s->krylov_basis[s->m-1]);
            }

            //gram schmidt : b(i) = H i> - |j><j Hi 0> = w(i) - ...
            cplx c; int i1,i2,i3;
            {
                for (i1=1; i1<s->m; i1++) {// calc b(i)
                    for (i2=0; i2<i1; i2++) {
                        //c = s->krylov_basis[i2].mult(s->krylov_basis[i1]);//<j Hi 0> //to be removed
                        c = khio(i2, i1);
                        for (i3=0; i3<k2; i3++) s->krylov_basis[i1][i3] -= s->krylov_basis[i2][i3] * c;//|i> -= c|j>
                    }

                    normC[i1] = s->krylov_basis[i1].normalize();
                }
            }

            for (int i=0; i<=s->m; i++) {
                for (int j=0; j<s->m; j++) {//k
                    k_h_o[i*s->m+j] = khio(j, i);
                }
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
            if (!s) return;

            long double mem = sizeof(cplx)*k2*s->m;
            long double mem_mb = mem/1048576.0;
            long double mem_gb = mem/1073741824.0;

            cout << "\n\n System : ";
            cout << "\n  size " << s->k << " x " << s->k << " sites";
            cout << "\n  krylov dim " << s->m ;

            cout << "\n\n Memory used : ";
            cout << "\n  " << mem << " bytes\n   = " << mem_mb << " mb\n   = " << mem_gb << " gb\n\n";
        }

    public:
        krylovRaum() {
            H_k = 0;
            //J_k = 0;
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

        void set(storage* _s, mpi_grid* _gr) {
            s = _s;
            gr = _gr;
            int m = s->m;
            k2 = s->k*s->k;

            H_k = new cplx[m*m];
            normC = new double[m];
            o_hh_o = new cplx[2*m];
            k_h_o = new cplx[m*m+m];

            normC[0] = 1;//norm of state allways one!
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
                s->krylov_basis[0][i] *= k_v[0];
                for (j=1 ; j<s->m ; j++) {
                    s->krylov_basis[0][i] += s->krylov_basis[j][i] * k_v[j];
                }
            }
        }

        cplx* getHamilton() { return H_k; }

        /*cplx* getCurrent() {
            return J_k;
        }*/

        void do_J_p0(float nx, float ny, float nz) {
            JPsi_mpi(s->krylov_basis[0], s->krylov_basis[1], nx, ny, nz);
            s->krylov_basis[0].copy(s->krylov_basis[1]);
        }

        //gibt in der console den vektor j aus
        void print(cplx* v) {
            cout << "\nprint vector\n";
            for (int i=0;i<k2;i++) cout << v[i] << " ";
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
