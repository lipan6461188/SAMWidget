#include "mainwindow.h"
#include <QApplication>
#include "calcrpkm.h"
#include <sstream>




int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();


    /*
    SamAlign samAlign("/Users/lee/test.sam");
    samAlign.calcBD();

    ofstream OUT("/Users/lee/test.bd", ofstream::out);
    if(OUT)
        samAlign.writeBD(OUT);
    OUT.close();

    cout << "finish..." << endl;
    */

    /*
    shared_ptr<ifstream> pSAM;
    if(pSAM)
        cout << "Hello" << endl;
    else
        cout << "Yes" << endl;
    */

   // cout << "Start to run..." << endl;
   // SamAlign samAlign("/Users/lee/DG_python_new.sam");

    //samAlign.printChrLen();

    return a.exec();
}
