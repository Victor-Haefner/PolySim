#ifndef RECORDER_H_INCLUDED
#define RECORDER_H_INCLUDED

#include <algorithm>

using namespace std;

class recorder {
    private:
        vector<cplx* > rows;//rows of the 2D state to save
        vector<char* > defects;//rows of the 2D state to save

        vector<string> files;

        int rN, cN;//number of rows, collums

        void init(string path) {
            files.push_back(path);
            ofstream file(path.c_str(), fstream::out | fstream::binary);
            file.write((char*)&rN, sizeof(int));
            file.write((char*)&cN, sizeof(int));
            for (int i=0;i<rN;i++)  file.write((char*)defects[i], cN*sizeof(char));
            file.close();
        }

        void append(string path) {
            ofstream file(path.c_str(), fstream::out | fstream::binary | fstream::app);
            for (int i=0;i<rN;i++)  file.write((char*)rows[i], cN*sizeof(cplx));
            file.close();
        }

    public:
        recorder() { rN = cN = 0; }

        void set(cplx* v, cplx* def, int k, int x0, int y0, int xm, int ym) {
            cout << "\n set recorder at : " << x0 << " " << xm << " x " << y0 << " " << ym;

            rN = ym-y0;
            cN = xm-x0;

            for (int i=y0;i<ym;i++) {
                char* _d = new char[cN];
                defects.push_back(_d);

                //cplx* _v = new cplx[cN];
                //rows.push_back(_v);

                rows.push_back(v + x0 + i*k);

                for (int j=x0;j<xm;j++) {
                    //_v[j-x0] = v[j + i*k];
                    if(def[j + i*k] == cplx(0,0)) _d[j-x0] = 'D';
                    else if((i+j)%2 == 0) _d[j-x0] = 'B';
                    else _d[j-x0] = 'A';
                }
            }
        }

        void take(string path) {
            if(find(files.begin(), files.end(), path) != files.end()) append(path);
            else init(path);
        }
};

#endif // RECORDER_H_INCLUDED
