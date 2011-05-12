#ifndef CAIROVISUAL_H_INCLUDED
#define CAIROVISUAL_H_INCLUDED

#include "CImg.h"
#include <math.h>
#include <string>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

class cimgVisual_base {
    protected:
        typedef unsigned char color;

        cimg_library::CImg<unsigned char> canvas;
        cimg_library::CImgDisplay main_disp;
        int H;
        int W;

        int scale;
        int colScale;
        int mode;

        bool save; //if true the frames will be saved for a vid
        string fileTag;

        void drawPixel(int x, int y, int r, int g, int b) {
            color col[] = {r,g,b};
            canvas.draw_point( x, y, col) ;
        }

        long frameCounter;
        void saveFrame() {
            string path = "vid/frames/frame";
            path += fileTag;

            int id = frameCounter;

            char a0 = id%10 + '0';
            id /= 10;
            char a1 = id%10 + '0';
            id /= 10;
            char a2 = id%10 + '0';
            id /= 10;
            char a3 = id%10 + '0';
            id /= 10;
            char a4 = id%10 + '0';


            path.push_back(a4);
            path.push_back(a3);
            path.push_back(a2);
            path.push_back(a1);
            path.push_back(a0);

            path.append(".png");

            canvas.save_png(path.c_str());
            frameCounter++;
        }

        void resizeDisp() {
            main_disp.resize(W*scale, H*scale);
            canvas.resize(W*scale, H*scale);
        }

    public:
        cimgVisual_base(int w, int h, int s = 1, int c = 1) {
            W = w;
            H = h;

            scale = s;
            colScale = c;
            frameCounter = 0;
            mode = 0;

            save = 0;

            canvas = cimg_library::CImg<unsigned char>(w*scale, h*scale, 1, 3, 0);

            main_disp = cimg_library::CImgDisplay(canvas, "A CImg-Example");
        }

        void setMode(int m) {
            mode = m;
        }

        void logQuary(string tag, char c = ' ') {
            while(true) {
                if (c == 'y') {
                    save = true;
                    fileTag = tag;
                    return;
                } else if (c == 'n') {
                    save = false;
                    return;
                }
                cout << "\nDo you want to log the visual output? (y/n) : ";
                cin >> c;
            }
        }

        void saveCanvas(string path) {
            canvas.save_png(path.c_str());
        }

        void resize(int w, int h, double s = 0) {
            W = w;
            H = h;
            if (s != 0) scale = s;
            resizeDisp();
        }

        void update() {
            if (save) saveFrame();
            main_disp.display(canvas);
        }
};

class graph : public cimgVisual_base {
    private:
        cplx* Xvalues;
        cplx* Yvalues;

        float xmin, xmax, ymin, ymax;

        int zero[2];

        int N;

        color r;
        color g;
        color b;

        void drawYAxis() {
            string sp = "";
            string sn = "";
            char str[100];
            color bg[] = {255,255,255};
            color tmp[] = {0,0,0};

            float sca = (ymax-ymin)/H;
            int delta = 60;

            for (int i=0;i<H;i++) drawPixel(zero[0],i,150,150,150);//vertikale graue linie
            for (int i=0;i<H; i+= delta) {
                for (int j=-3; j<4; j++) {
                    drawPixel(zero[0] +j, zero[1]+i, 150,150,150);
                    drawPixel(zero[0] +j, zero[1]-i, 150,150,150);
                }
            }

            for (int i = 0; i<H; i+= delta) {
                float e = i*sca;
                sp = sn = "";
                sp.append(str, sprintf(str,"%.*f", 3, -e));
                sn.append(str, sprintf(str,"%.*f", 3,  e));
                sp += " "; sn += " ";
                canvas.draw_text(sp.c_str(), 0, zero[1] + i-3, tmp, bg);//beschrifte positive axe
                canvas.draw_text(sn.c_str(), 0, zero[1] - i-3, tmp, bg);//beschrifte negative axe
            }

            canvas.draw_text("D(E)", 20, zero[1] + 20, tmp, bg);
        }

        void drawXAxis() {
            float sca = (float)(xmax-xmin)/W;
            int delta = 60;

            for (int i=0;i<W;i++) drawPixel(i,zero[1],150,150,150);
            for (int i=0;i<W; i+= delta) {
                for (int j=-3; j<4; j++) {
                    drawPixel(zero[0] +i, zero[1]+j, 150,150,150);
                    drawPixel(zero[0] -i, zero[1]+j, 150,150,150);
                }
            }

            //numbers
            string sp = "";
            string sn = "";
            char str[100];
            color bg[] = {255,255,255};
            color tmp[] = {0,0,0};

            for (int i = 0; i<W; i+= delta) {
                float e = i*sca;
                sp = sn = "";
                sp.append(str, sprintf(str,"%.*f", 3,  e));
                sn.append(str, sprintf(str,"%.*f", 3, -e));
                sp += " "; sn += " ";
                canvas.draw_text(sp.c_str(), zero[0] +i-3, H-15, tmp, bg);//beschrifte positive axe
                canvas.draw_text(sn.c_str(), zero[0] -i-3, H-15, tmp, bg);//beschrifte negative axe
            }

            canvas.draw_text("E/t (t = 1)", zero[0] + 20, H-40, tmp, bg);
        }

        void _draw() {
            reset();

            float scax = W/(xmax-xmin);
            float scay = H/(ymax-ymin);

            for (int i=0;i<N;i++) {
                double x = real(Xvalues[i])*scax;
                double y_r = real(Yvalues[i]);

                drawPixel(zero[0] + x, zero[1] - y_r*scay, 0, 0, 0);
            }

            update();
        }

        void zoom(float sx_r, float sx_l, float sy_d, float sy_u) {
            xmin *= sx_l;
            xmax *= sx_r;
            ymin *= sy_u;
            ymax *= sy_d;

            zero[0] = W*xmin/(xmin-xmax);
            zero[1] = H - H*ymin/(ymin-ymax);

            _draw();
        }

        void startInteraction() {
            while (!main_disp.is_closed) {
                main_disp.wait();
                if (main_disp.button == 2) return;
                if (main_disp.button == 1) {
                    float x = main_disp.mouse_x - W/2;
                    float y = main_disp.mouse_y - H/2;

                    x /= W/2;
                    y /= H/2;

                    float d = 0.95;

                    if (abs(x) > abs(y)) {
                        if (x > 0) zoom(d + 2*(1-d)*x, 1, 1, 1);
                        else zoom(1, d - 2*(1-d)*x, 1, 1);
                    }

                    if (abs(x) < abs(y)) {
                        if (y < 0) zoom(1, 1, d - 2*(1-d)*y, 1);
                        else zoom(1, 1, 1, d + 2*(1-d)*y);
                    }
                }
            }
        }

    public:
        graph(int w, int h) : cimgVisual_base(w,h) {
            Xvalues = 0;
            Yvalues = 0;
            zero[0] = w/2;
            zero[1] = h/2;

            r = 255; g = 255; b = 0;
        }

        void setVector(int n, cplx* xv, cplx* yv, float x1, float x2, float y1, float y2, int r_, int g_, int b_) {
            r = r_;
            g = g_;
            b = b_;
            Xvalues = xv;
            Yvalues = yv;
            N = n;

            xmin = x1;
            xmax = x2;
            ymin = y1;
            ymax = y2;

            cout << "\nvector SET\n";

            zoom(1,1,1,1);
            reset();
        }

        void draw() {
            _draw();
            startInteraction();
        }

        void reset() {
            canvas.fill(255);
            drawYAxis();
            drawXAxis();
        }
};

class timeline : public cimgVisual_base {
    private:
        vector<double>* values;
        int t;
        int legend;
        float scale;

        string legendStr[5];

        color r[5];
        color g[5];
        color b[5];

        void init() {
            //horizontal lines
            canvas.fill(0);
            for (int i=legend;i<W;i++) {
                for (int j=0;j<H;j+=5*scale) drawPixel(i,j,50,50,50);
            }

            //axis
            string s0 = "10e-";
            string s;
            char str[5];
            for (int i = 0; i<H - 15; i+= 5*scale) {
                int e = i/scale;

                color bg[] = {0,0,0};
                color tmp[] = {100,100,100};
                s = s0;
                s.append(str, sprintf(str,"%d",e));
                s += " ";
                canvas.draw_text(s.c_str(), legend - 50, 2+i, tmp, bg);
            }
        }

        void drawLegend() {
            color bg[] = {0,0,0};
            for (int i=0;i<5;i++) {
                color tmp[] = {r[i],g[i],b[i]};
                canvas.draw_text(legendStr[i].c_str(), 3, 2+i*19, tmp, bg);
            }
        }

    public:
        timeline(int w, int h, float s, vector<double>* v = 0) : cimgVisual_base(w+100,h) {
            values = v;
            legend = 100;
            t = legend;
            scale = s;

            r[0] = 255;
            g[0] = b[0] = 0;

            g[1] = 255;
            r[1] = b[1] = 0;

            b[2] = 255;
            g[2] = r[2] = 0;

            g[3] = r[3] = 255;
            b[3] = 0;

            r[4] = b[4] = 255;
            g[4] = 0;


            legendStr[0] = " Qual ";
            legendStr[1] = " Diag ";
            legendStr[2] = " Norm ";
            legendStr[3] = " Orth ";
            legendStr[4] = " dt ";

            init();
            drawLegend();
        }

        void setVector(vector<double>* v) {
            values = v;
        }

        void draw() {
            t++;
            if (t > W) t = 101;
            for (unsigned int i=0;i<values->size();i++) {
                drawPixel(t,(*values)[i]*scale, r[i], g[i], b[i]);
            }
        }

        void setLegend(int i, string s) {
            legendStr[i] = s;
        }

        void reset() {
            t = 0;
        }
};

template<typename type>
class histogram : public cimgVisual_base {
        int j;

    public:
        histogram(int w, int h, int s = 1, int c = 1) : cimgVisual_base(w,h,s,c) {
            ;
        }

        /*void draw() {
            canvas.fill(255);
            double w = 0;
            color col[] = {0,0,0};

            for (int i=0;i<W;i++) {
                switch(mode) {
                    case 0:
                        w = norm((*vec)[W*j+i]);
                        break;
                    case 1:
                        w = (*vec)[W*j+i].real();
                        break;
                    case 2:
                        w = (*vec)[W*j+i].imag();
                        break;
                }
                w = std::log(w + 1);
                w *= colScale;
                w = ceil(w);
                for (int u=0;u<scale;u++) canvas.draw_point(i*scale + u, (H-1)*scale - w, col);
            }
        }*/


};

template<typename type>
class topView : public cimgVisual_base {

        void drawSquare(int x, int y, int r, int g, int b) {
            for (int i=0;i<scale;i++)
                for (int j=0;j<scale;j++)
                    drawPixel( x*scale+i, y*scale+j, r,g,b) ;
        }

        void getColor(double a, int& r, int& g, int& b) {

            a = std::log(a + 1);
            a *= colScale;

            int i = ceil(a);
            if (i<0) i = -i;

            int c = i%255;
            if (i>=255) c = 255;

            r = c;
            g = c;
            b = c;
        }


    public:
        topView(int w, int h, int s = 1, int c = 1) : cimgVisual_base(w,h,s,c) {
            ;
        }

        void drawDefect(int i, int j) {
            drawSquare(i, j, 255, 0, 0);
        }

        void drawCplx(type& c, int i, int j) {
            //cout << "\ndraw : " << i << " " << j << " " << c;
            int r,g,b;
            double w;
            switch(mode) {
                case 0:
                    w = sqrt(norm(c));
                    break;
                case 1:
                    w = c.real();
                    break;
                case 2:
                    w = c.imag();
                    break;
            }
            getColor(w, r, g, b);
            drawSquare(i, j, r, g, b);
        }
};

#endif // CAIROVISUAL_H_INCLUDED
