#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QString>
#include <QFileDialog>
#include <memory>
#include <QTreeWidget>
#include <QDir>
#include <fstream>
#include <QStringList>

#include "include_heder.h"
#include "calcrpkm.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    //get status bar
    const QStatusBar *getStatusBar(){ return this->statusBar(); }

private:
    Ui::MainWindow *ui;

    QString samFileName;
    //SAM File information
    SamAlign samAlign;


private slots:
    void chooseFile(QString file_prefix="sam");
    void loadFile();
    void showLen();
    void readFewReads();
    void calcRPKM();
};

#endif // MAINWINDOW_H
