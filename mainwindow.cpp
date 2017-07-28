#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "ui_bdwindow.h"

#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    bd_ui(new Ui::BDWindow)
{
    ui->setupUi(this);
    ui->chrLenList->setColumnCount(2);

    bd_dialog.reset(new QDialog(this));
    bd_ui->setupUi(bd_dialog.get());
    bd_dialog->setVisible(false);

    statusBar()->showMessage(tr("Ready"));
    //p_SAM = nullptr;
    connect(ui->choseFileButton, SIGNAL(clicked(bool)), this, SLOT(chooseFile()));
    connect(ui->loadFileButton, SIGNAL(clicked(bool)), this, SLOT(loadFile()));
    connect(ui->showLenButton, SIGNAL(clicked(bool)), this, SLOT(showLen()));
    connect(ui->fewReadsButton, SIGNAL(clicked(bool)), this, SLOT(readFewReads()));
    connect(ui->rpkmButton, SIGNAL(clicked(bool)), this, SLOT(calcRPKM()));
    connect(ui->bdButton, SIGNAL(clicked(bool)), this, SLOT(calcBDRT()));
    connect(bd_ui->showBDButton, SIGNAL(clicked(bool)), this, SLOT(plotBD()));
    connect(bd_ui->chrNameListWidget, SIGNAL(currentRowChanged(int)), this, SLOT(showBDRT()));
}

void MainWindow::chooseFile(QString file_prefix)
{
    const char *fileFormat = QString("%1 Files (*.%1)").arg(file_prefix).toStdString().c_str();
    //qDebug() << "fileFormat: " << fileFormat << "\n";
    samFileName = QFileDialog::getOpenFileName(this, tr("Choose a sam file"), QDir::homePath(), tr(fileFormat));

    ui->samFIleLineEdit->setText(samFileName);
}

void MainWindow::loadFile()
{

    samFileName = ui->samFIleLineEdit->text();

    std::ifstream *p_SAM = new std::ifstream(samFileName.toStdString(), std::ifstream::in);

    if(not p_SAM->is_open())
    {
        statusBar()->showMessage( QString("FATAL Error: open sam file ") + samFileName + " Failed"  );
        //std::cerr << "FATAL Error: open sam file " << samFileName.toStdString() << " Failed" << std::endl;
        //p_SAM = nullptr;
    }else{
        //std::cout << "Success: open sam file " << samFileName.toStdString() << std::endl;
        statusBar()->showMessage( QString("Success: open sam file ") + samFileName  );
        samAlign.setIn(p_SAM);
    }
}

void MainWindow::showLen()
{
    ui->chrLenList->clear();
    ui->chrLenList->setColumnCount(2);
    if(samAlign.is_good())
    {
        QList<QTreeWidgetItem *> items;

        auto chrLen = samAlign.getChrLen();
        for(auto a: chrLen)
        {
            QStringList chrLenLine = QStringList() << QString(a.first.c_str()) << QString(to_string(a.second).c_str());
            items.append(new QTreeWidgetItem( ui->chrLenList, chrLenLine ));
        }
        ui->chrLenList->insertTopLevelItems(0, items);
    }
}

void MainWindow::readFewReads()
{
    ui->chrLenList->clear();

    size_t num = ui->showLinesEdit->text().toInt();
    CStringMatrix stringMatrix = samAlign.readFewReads(num);
    if(stringMatrix.size() == 0)
    {
        statusBar()->showMessage("Sam file has a bad status or has no reads to read");
        return;
    }

    size_t readWidth = stringMatrix.at(0).size();
    ui->chrLenList->setColumnCount(readWidth);

    QList<QTreeWidgetItem *> items;
    for(auto read: stringMatrix)
    {
        QStringList readSample;
        for(auto readItem: read)
            readSample << QString(readItem.c_str());
        items.append(new QTreeWidgetItem( ui->chrLenList,  readSample));
    }
    ui->chrLenList->insertTopLevelItems(0, items);
}

void MainWindow::calcRPKM()
{
    statusBar()->showMessage("Start to calculate RPKM...");
    //void calcRPKM(bool removeMultiMap=true, bool removeGappedRead=true, bool removeReverseRead=true);
    samAlign.calcRPKM(not ui->multiBox->isChecked(), not ui->gappedBox->isChecked(), not ui->reverseBox->isChecked());
    statusBar()->showMessage("Finish calculating RPKM...");

    // show RPKM information
    const StringToUL &chrLen = samAlign.getChrLen();
    const StringToD &chrRPKM = samAlign.getRPKM();
    const StringToUL &chrSingleRead = samAlign.getChrSingleRead();
    const StringToDA &chrMultiRead = samAlign.getChrMultiRead();
    uLONG total_mapped_reads = samAlign.getTotalReadNum();

    ui->chrLenList->clear();
    ui->chrLenList->setColumnCount(6);

    QStringList ColumnNames = QStringList() << "chrName" << "chrLen" << "uniqMap" << "multi-map" << "multi-map-power" << "RPKM";
    ui->chrLenList->setHeaderLabels(ColumnNames);

    QList<QTreeWidgetItem *> items;
    for(auto chrLenPair: chrLen)
    {
        string chrName = chrLenPair.first;
        QStringList chrRPKMInfo;


        chrLen.at(chrName);

        double averaged_mutiple_mapped_reads = accumulate(chrMultiRead.at(chrName).begin(), chrMultiRead.at(chrName).end(), 0.0);
      //  std::cout << chrName << ":  " << averaged_mutiple_mapped_reads << endl;
        uLONG multiMapReads = chrMultiRead.at(chrName).size();
        chrRPKMInfo << QString(chrName.c_str()) << QString( std::to_string(chrLen.at(chrName)).c_str() ) << QString( std::to_string(chrSingleRead.at(chrName)).c_str() )
                    << QString( std::to_string(multiMapReads).c_str() ) << QString::number(averaged_mutiple_mapped_reads)
                    << QString( std::to_string(chrRPKM.at(chrName)).c_str() );
        items.append(new QTreeWidgetItem(ui->chrLenList, chrRPKMInfo));
    }
    std::cout << "Mapped Reads: " << total_mapped_reads << endl;
    ui->chrLenList->insertTopLevelItems(0, items);
}

/*
 *  //get chrLen
    const StringToUL& getChrLen(){ return chrLen; }
    //get RPKM
    const StringToD & getRPKM(){ return chrRPKM; }
    //get chrSingleRead
    const StringToUL & getChrSingleRead(){ return chrSingleRead; }
    //get RPKM
    const StringToDA & getChrMultiRead(){ return chrMultiRead; }
    //get total num
    uLONG getTotalReadNum(){ return total_mapped_reads; }
*/

void MainWindow::calcBDRT()
{
    std::clog << "Start to calculate BD..." << endl;
    statusBar()->showMessage("Start to calculate BD...");
    //void calcRPKM(bool removeMultiMap=true, bool removeGappedRead=true, bool removeReverseRead=true);
    samAlign.calcBDRT(not ui->multiBox->isChecked(), not ui->gappedBox->isChecked(), not ui->reverseBox->isChecked());
    statusBar()->showMessage("Finish calculating BD...");

    std::clog << "Finish calculating BD..." << endl;

    bd_ui->chrNameListWidget->clear();

    // get chr info
    const StringToUL &chrLen = samAlign.getChrLen();
    for(const auto &chr_length: chrLen)
        bd_ui->chrNameListWidget->addItem( QString::fromStdString(chr_length.first) );

    bd_ui->chrNameListWidget->setCurrentRow(0);
    bd_dialog->setVisible(true);

    std::clog << "Finish Load All BD/RT..." << endl;
}

void MainWindow::showBDRT()
{
    string chrName = bd_ui->chrNameListWidget->currentItem()->text().toStdString();
    auto BD_RT = samAlign.getChrBDRT(chrName);
    //auto chrBD = BD_RT.first;
    //auto chrRT = BD_RT.second;

    bd_ui->bdTreeWidget->clear();
    bd_ui->bdTreeWidget->setColumnCount(2);
    QStringList ColumnNames = QStringList() << "base" << "Base Density" << "RT";
    bd_ui->bdTreeWidget->setHeaderLabels(ColumnNames);

    QList<QTreeWidgetItem *> items;
  //  size_t baseIndx = 1;
    for(size_t baseIndx=0; baseIndx<BD_RT.first.size(); baseIndx++)
    {
        QStringList chrBDInfo;

        chrBDInfo << QString::number(baseIndx+1) << QString::number(BD_RT.first.at(baseIndx)) << QString::number(BD_RT.second.at(baseIndx));

     //   baseIndx++;

        //double averaged_mutiple_mapped_reads = accumulate(chrMultiRead.at(chrName).begin(), chrMultiRead.at(chrName).end(), 0.0);
      //  std::cout << chrName << ":  " << averaged_mutiple_mapped_reads << endl;
       // uLONG multiMapReads = chrMultiRead.at(chrName).size();


        //chrRPKMInfo << QString(chrName.c_str()) << QString( std::to_string(chrLen.at(chrName)).c_str() ) << QString( std::to_string(chrSingleRead.at(chrName)).c_str() )
          //          << QString( std::to_string(multiMapReads).c_str() ) << QString::number(averaged_mutiple_mapped_reads)
            //        << QString( std::to_string(chrRPKM.at(chrName)).c_str() );
        items.append(new QTreeWidgetItem(bd_ui->bdTreeWidget, chrBDInfo));
    }


}



void MainWindow::plotBD()
{

   /*
    string chrName = bd_ui->chrNameListWidget->currentItem()->text().toStdString();
    vector<uLONG> chrBD = samAlign.getChrBD(chrName);

    QtCharts::QBarSet *set0 = new QtCharts::QBarSet( QString::fromStdString(chrName) );
    for(auto bd: chrBD)
        *set0 << bd;

    QtCharts::QBarSeries *series = new QtCharts::QBarSeries();
    series->append(set0);

    QtCharts::QChart *chart = new QtCharts::QChart();
    chart->addSeries(series);
    chart->setTitle("Simple barchart example");
    chart->setAnimationOptions(QtCharts::QChart::SeriesAnimations);


    QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QMainWindow *window = new QMainWindow();
    window->setCentralWidget(chartView);
    window->resize(420, 300);
    window->show();
    */


}

MainWindow::~MainWindow()
{
    delete ui;
}
