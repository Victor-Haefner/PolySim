#ifndef WAVEPACKET_H_INCLUDED
#define WAVEPACKET_H_INCLUDED

#include "storage.h"

using namespace std;

class wavepacket {
    private:
        int k, x0, y0, x, y;
        float kx0, ky0, s;

        bool isSameSublatice(int x, int y) { return (x+y-x0-y0)%2 == 0; }

        cplx getGausSquare(int i) {
            x = i%k;
            y = i/k;

            cplx a = cplx(0, (kx0*(x-x0) + ky0*(y-y0)) );

            cplx A = (x - x0)*(x - x0) + (y - y0)*(y - y0);

            A *= -1/(2*s*s);

            return exp(A + a);
        }

        cplx getGausHoneycomb(int i, float la = 1) {
            x = i%k;
            y = i/k;

            int Dxi = x-x0;
            int Dyi = y-y0;

            //Same SubLattice?
            bool ssl = isSameSublatice(x,y);

            //site distance in space
            float Dx = Dxi*sqrt(3)/2.*la;
            float Dy = Dyi*1.5*la;
            if (!ssl) {
                if ((y0 + x0)%2 == 1) Dy += 0.5*la; // +/- ..je nachdem auf welchem untergitter ich mich in x0y0 befinde!!
                else Dy -= 0.5*la;
            }

            //gaus
            cplx a = cplx(0, (kx0*Dx + ky0*Dy) );
            cplx A = Dx*Dx + Dy*Dy;

            A *= -1/(2*s*s);

            return exp(A + a);
        }

    public:
        wavepacket() {
            ;
        }

        void set(cplx* v, int k_, int x0_, int y0_, float p, float p_phi, float s_) {
            k = k_;
            x0 = x0_;
            y0 = y0_;

            s = s_;


            kx0 = p*cos(p_phi);
            ky0 = p*sin(p_phi);

            for (int i=0;i<k*k;i++) {
                if (options::get()->graphene) v[i] = getGausHoneycomb(i);
                else v[i] = getGausSquare(i);
            }
        }

};

#endif // WAVEPACKET_H_INCLUDED
