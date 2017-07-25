#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)//:p_SAM(nullptr)
{
    ui->setupUi(this);
    ui->chrLenList->setColumnCount(2);


    statusBar()->showMessage(tr("Ready"));
    //p_SAM = nullptr;
    connect(ui->choseFileButton, SIGNAL(clicked(bool)), this, SLOT(chooseFile()));
    connect(ui->loadFileButton, SIGNAL(clicked(bool)), this, SLOT(loadFile()));
    connect(ui->showLenButton, SIGNAL(clicked(bool)), this, SLOT(showLen()));
    connect(ui->fewReadsButton, SIGNAL(clicked(bool)), this, SLOT(readFewReads()));
    connect(ui->rpkmButton, SIGNAL(clicked(bool)), this, SLOT(calcRPKM()));
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
        std::cout << chrName << ":  " << averaged_mutiple_mapped_reads << endl;
        uLONG multiMapReads = chrMultiRead.at(chrName).size();
        chrRPKMInfo << QString(chrName.c_str()) << QString( std::to_string(chrLen.at(chrName)).c_str() ) << QString( std::to_string(chrSingleRead.at(chrName)).c_str() )
                    << QString( std::to_string(multiMapReads).c_str() ) << QString::number(averaged_mutiple_mapped_reads)
                    << QString( std::to_string(chrRPKM.at(chrName)).c_str() );
        items.append(new QTreeWidgetItem(ui->chrLenList, chrRPKMInfo));
    }
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

MainWindow::~MainWindow()
{
    delete ui;
}
