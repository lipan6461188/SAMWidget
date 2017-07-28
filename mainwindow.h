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
#include <QDialog>
#include <QListView>
#include <QtCharts/QBarSet>
#include <QtCharts/QBarSeries>
#include <QtCharts/QChart>
#include <QChartView>
#include <QtCharts>
#include <QPainter>

#include "include_heder.h"
#include "calcrpkm.h"

namespace Ui {
class MainWindow;
class BDWindow;
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
    Ui::BDWindow *bd_ui;

    shared_ptr<QDialog> bd_dialog;

    QString samFileName;
    //SAM File information
    SamAlign samAlign;


private slots:
    void chooseFile(QString file_prefix="sam");
    void loadFile();
    void showLen();
    void readFewReads();
    void calcRPKM();
    void calcBDRT();

    void showBDRT();
    void plotBD();
};

#endif // MAINWINDOW_H
