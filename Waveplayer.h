#ifndef WAVEPLAYER_H_INCLUDED
#define WAVEPLAYER_H_INCLUDED

using namespace std;

class player {
    private:
        //visualization
        topView<cplx>* top;

        //data
        sequence* seq;

    public:
        player() {
            top = 0;

            top = new topView<cplx>(100, 50, 2, 1000);
        }

        ~player() {
            if (top) delete top;
        }

        void set(sequence* _seq) {
            seq = _seq;
            resize(seq->width, seq->height, 2);
        }

        void resize(int w, int h, double s) {
            if (top == 0) return;
            top->resize(w,h,s);
        }

        void draw(int key) {
            if (top == 0) return;
            if (seq->frames == 0) return;

            for (int i=0; i<seq->height; i++) {
                for (int j=0; j<seq->width; j++) {
                    top->drawCplx(seq->frames[i*seq->width + j + key*seq->width*seq->height], i, j);
                }
            }

            for (int i=0;i<seq->kd;i++) top->drawDefect(seq->defects[i]%seq->width, seq->defects[i]/seq->height);

            top->update();
        }

        void play() {
            cout << "\nStart visualising\n";

            for (int i=0;i<seq->lenght;i++) {
                cout << "\nDraw frame : " << i;
                draw(i);
            }

            cout << "\n\nEnd of visualisation.\n";
        }
};


#endif // WAVEPLAYER_H_INCLUDED
