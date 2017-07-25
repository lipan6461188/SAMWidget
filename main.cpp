#include "mainwindow.h"
#include <QApplication>
#include "calcrpkm.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();


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
