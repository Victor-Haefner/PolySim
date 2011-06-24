#ifndef CURRENT_H_INCLUDED
#define CURRENT_H_INCLUDED

#include "storage.h"

using namespace std;

class diffusion {
    private:
        storage* s;
        options* opt;

        double* dvec;

        int x0, y0, t_min, dr, rmax;

        float getDistance(int x, int y) {
            int Dxi = x-x0;
            int Dyi = y-y0;

            //lattice cst
            float la = 1;

            //Same SubLattice?
            bool ssl;
            if ((Dxi%2 + Dyi%2)%2 == 0) ssl = true;
            else ssl = false;

            //site distance in space
            float Dx = Dxi*sqrt(3)/2.*la;
            float Dy = Dyi*1.5*la;
            if (!ssl) {
                if (y0%2 + x0%2 == 1) Dy += 0.5*la; // +/- ..je nachdem auf welchem untergitter ich mich in x0y0 befinde!!
                else Dy -= 0.5*la;
            }

            return sqrt(Dx*Dx + Dy*Dy);
        }

    public:
        diffusion() {
        }

        void set(storage* _s, options* _opt, int _x0, int _y0) {
            s = _s;
            opt = _opt;

            //build D(r,t) matrix
            t_min = opt->t_min;
            dr = opt->dr;
            //s->N = opt->R*(s->T/t_min);//shape
            t_min *= round(1./s->dt);
            s->N = opt->R*(s->T-t_min);
            dvec = new double[s->N];
            for (int i=0;i<s->N;i++) dvec[i] = 0;

            x0 = _x0;
            y0 = _y0;

            rmax = opt->R*dr;
        }

        void process(int t) {
            if (t >= t_min) {//let it propagate for a given total time before taking data

                for (int i=-rmax; i<rmax; i++) {
                    for (int j=-rmax; j<rmax; j++) {//gehe durch das gebiet das den kreis einschileßt
                        float _r = getDistance(i+x0,j+y0);
                        int r = round(_r/dr - 0.5);
                        if (r>=opt->R) continue;
                        if (r<0) r = 0;

                        //float A;
                        //if(r == 0) A = 1;//sitzt in der mitte
                        //else A = 2*pi*r*dr*dr;//area of a circle with border 1, (no error)

                        cplx c = s->krylov_basis[0][(i+x0)*s->k + (j+y0)];
                        //c *= 1./N[r];
                        //c *= 1./A;

                        dvec[t-t_min + r*(s->T-t_min)] += norm(c);
                    }
                }
            }
        }

        void getShape(int t) {
            int tt = t_min*s->dt;
            if (t%tt != 0) return;

            //sum over all sites in the area r_i r_i+1
            int k = t/tt;
            for (int r=0; r<opt->R; r++) {
                cplx c = s->krylov_basis[0][x0*s->k+y0+r*2];//r*2 damit die punkte geometrisch auf einer reihe liegen!
                dvec[k*opt->R + r] = norm(c);
            }
        }

        void save() {

            /*---------------------------EXPORT-------------------------------*/

            string path = "D_rt_";
            path += s->getPath();

            ofstream file(path.c_str());

            for (int i=0;i<s->N;i++) {
                if (i%(s->T-t_min) == 0 and i>0) file << "\n\n";
                double _t = (i%(s->T-t_min) + t_min)*s->dt;
                double _v = dvec[i];
                file << "\n" << _t << " " << _v;
            }
            file << "\n";

            file.close();
        }

        void saveShape() {

            /*---------------------------EXPORT-------------------------------*/

            string path = "D_rt_";
            path += s->getPath();

            ofstream file(path.c_str());

            for (int i=0;i<s->N;i++) {
                if (i%opt->R == 0 and i>0) file << "\n\n";
                double _t = i%opt->R;
                double _v = norm(s->dos[i]);
                //file << _t*_t << " " << log(_v) << "\n";
                file << _t*_t << " " << _v << "\n";
            }



            file.close();
        }
};

#endif // CURRENT_H_INCLUDED
