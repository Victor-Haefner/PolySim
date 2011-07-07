#ifndef GENERICBINARYIO_H_INCLUDED
#define GENERICBINARYIO_H_INCLUDED

using namespace std;

class binIO {
    private:
    public:
        binIO() {;}

        //write the content of this class to file
        void write(string path, int size) {
            //ofstream file(path.c_str(), fstream::in | fstream::out | fstream::binary);
            ofstream file(path.c_str(), fstream::out | fstream::binary);
            file.seekp(ios_base::beg);
            file.write((char*)this, size);
            file.close();
        }

        //reads a chunk of data and puts it in destination
        void extract(string path, int begin, int size, char* destination) {
            ifstream file(path.c_str(), fstream::in | fstream::binary);
            if (!file.is_open()) {
                cout << "\nError while loading binary file, no file " << path << "\n";
                return;
            }

            file.seekg(ios_base::beg + begin);
            file.read(destination, size);
            file.close();
        }

        //read a file and put the data in this
        void read(string path, int size) {
            extract(path, 0, size, (char*)this);
        }

        //append data to a file
        void append(string path, int size, char* buffer) {
            ofstream ofile(path.c_str(), fstream::out | fstream::binary | fstream::app);
            ofile.write(buffer, size);
            ofile.close();
        }
};

#endif // GENERICBINARYIO_H_INCLUDED
