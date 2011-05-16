#ifndef TIMEEVOLUTION_H_INCLUDED
#define TIMEEVOLUTION_H_INCLUDED


//LAPACK Methode zur diagonalisierung von Matrizen
extern "C" {
    extern void zheev_(char *jobz, char *uplo, int *n, std::complex<double> *a, int *lda, double *w, std::complex<double> *work, int *lwork, double *rwork, int *info);
}

using namespace std;

class timeEvolution {
    private:
        cplx* S_k;//unitaere matrix vom diagonalisieren H_k S = S D
        cplx* H_k;
        cplx* V_k;
        double* EW;

        storage* s;

        double qlt[3];//qualitaet
        double Q, Qd;

        void copy(int l, cplx* target, cplx* source) {
            for (int i=0; i<l; i++) target[i] = source[i];
        }

        void transpose(int l, cplx* v) {
            for (int i=0;i<l;i++)
                for (int j=0;j<i;j++)
                    swap(v[i*l + j], v[j*l + i]);
        }

        void conjugate(int l, cplx* v) {
            for (int i=0;i<l;i++) v[i] = conj(v[i]);
        }

        void multMatrix(int l, cplx* m1, cplx* m2, cplx* result) {
            cplx* tmp = new cplx[s->m*s->m];

            for (int i=0;i<l;i++)
                for (int j=0;j<l;j++) {
                    int ii = i*l+j;
                    tmp[ii] = 0;
                    for (int k=0;k<l;k++)
                        tmp[ii] += m1[i*l+k]*m2[k*l+j];
                }

            copy(s->m*s->m, result, tmp);
	    delete[] tmp;
        }

        double normMatrix(int l, cplx* m1) {
            double sum = 0;
            cplx tmp;
            for (int i=0;i<l;i++) {
                for (int j=0;j<l;j++) {
                    tmp = std::conj(m1[i*l+j])*m1[i*l+j];
                    sum += real(tmp);
                }
            }
            return sum;
        }

        void printMatrix(cplx* m1) {
            for (int i2=0; i2<s->m; i2++) {
                cout << "\n";
                for (int i3=0; i3<s->m; i3++) {
                    cout << " " << m1[i2*s->m + i3];
                }
            }
        }


        void checkDiagQ() {
            int m2 = s->m*s->m;

            //S_k_1
            cplx* S_k_1 = new cplx[m2];
            copy(m2, S_k_1, S_k);//ok
            conjugate(m2, S_k_1);//ok
            transpose(s->m, S_k_1);//ok

            multMatrix(s->m, S_k_1, H_k, S_k_1);
            multMatrix(s->m, S_k_1, S_k, S_k_1);

            for (int i=0;i<s->m;i++) S_k_1[i*s->m+i] -= EW[i];// St H S - D

            Qd = normMatrix(s->m, S_k_1) / normMatrix(s->m, H_k);
            Qd = abs(log10(Qd));

	    delete[] S_k_1;
        }

        void checkQ() {
            //merke mir die letzten drei
            qlt[0] = qlt[1];
            qlt[1] = qlt[2];
            qlt[2] = sqrt(norm(V_k[s->m-1]));

            //berechne den mittelwert
            Q = qlt[0]+qlt[1]+qlt[2];
            Q /= 3;
            Q = abs(log10(Q));
        }

        void diagonalize(int m, cplx* mat, double* EWvec) {
            cplx* workSpace = new cplx;//ok
            double* dWorkSpace = new double[3*m-2];//ok

            int info;
            char V = 'V';
            char L = 'L';

            int wsSize = -1;
            zheev_(&V, &L, &m, mat, &m, EWvec, workSpace, &wsSize, dWorkSpace, &info);
            if (info != 0) return;

            wsSize = workSpace[0].real();
            delete workSpace;
            workSpace = new cplx[wsSize];

            zheev_(&V, &L, &m, mat, &m, EWvec, workSpace, &wsSize, dWorkSpace, &info);
            if (info != 0) return;

            delete[] dWorkSpace;
            delete[] workSpace;
        }

        void computeU_k(double dt) {//baustelle
            //berechne U = SDS-1
            cplx* tmp = new cplx[s->m];

            //berechne die diagonale exponential matrix D
            cplx a;
            for (int i=0; i<s->m; i++) {
                a = cplx(0,1);
                a *= dt;
                a *= EW[i];
                V_k[i] = exp(a)*S_k[i];//koennte statt i*m auch nur i sein
            }

            for (int i=0;i<s->m;i++) {
                tmp[i] = 0;
                for (int j=0;j<s->m;j++) {
                    tmp[i] += S_k[i*s->m+j]*V_k[j];
                }
            }

            for (int i=0;i<s->m;i++) V_k[i] = tmp[i];

	    delete[] tmp;

            //teste die hermitesche invertierung
            /*if (s->debug) {
                matrix tmp = S_k_1 * S_k;


                cplx c = 0;

                for (int i=0;i<m;i++) {
                    tmp[i][i] -= 1;// Ut U - 1
                    c += tmp[i][i];
                }

                if (norm(c) > 0.000001)
                    cout << "\nKrylov Zeitentwicklung, test failed : " << norm(c);
            }*/
        }

    public:
        timeEvolution() {
            EW = 0;
            V_k = 0;
            S_k = 0;
            H_k = 0;

            qlt[0] = qlt[1] = qlt[2] =0;
            Q = Qd = 0;
        }

        ~timeEvolution() {
            if (EW) delete[] EW;
            if (V_k) delete[] V_k;
            if (S_k) delete[] S_k;
            if (H_k) delete[] H_k;
        }

        void set(storage* _s) {
            s = _s;

            EW = new double[s->m];
            V_k = new cplx[s->m];
            S_k = new cplx[s->m*s->m];
            H_k = new cplx[s->m*s->m];
        }

        cplx* evolve(cplx* H) {
            copy(s->m*s->m, H_k, H);
            copy(s->m*s->m, S_k, H);

            diagonalize(s->m, S_k, EW);
            transpose(s->m,S_k);

            if (s->opt->debug) {
                checkDiagQ();
                checkQ();
            }

            computeU_k(s->dt);

            //cout << "\nVk : "; for (int i=0;i<m;i++) cout << V_k[i];
            return V_k;
        }


        void printS() {
            printMatrix(S_k);
        }

        double getQ() {
            return Q;
        }

        double getDiagQ() {
            return Qd;
        }
};

#endif // TIMEEVOLUTION_H_INCLUDED
