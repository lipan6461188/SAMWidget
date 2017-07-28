#ifndef CALCRPKM_H
#define CALCRPKM_H

#include "include_heder.h"
using namespace std;

// Const String Matrix
using CStringMatrix = const vector< vector<string> >;
// String Matrix
using StringMatrix = vector< vector<string> >;
// String to uLong
using StringToUL = map<string, uLONG>;
// String to double array
using StringToDA = map< string, vector<double> >;
// String to double
using StringToD = map<string, double>;
// String to uLONG Array
using StringToULMatrix = map<string, vector<uLONG>>;

// print stdout read information
void printRead(const vector<string> &read);

// Read Attribute
bool read_is_mapped(const vector<string> &read);
bool read_is_gapped(const vector<string> &read);
bool read_is_reverse(const vector<string> &read);

// add one read to a map to calc RPKM
void addOneRead_RPKM(CStringMatrix &sameReadArray, StringToUL &chrSingleRead,
                     StringToDA &chrMultiRead, uLONG &total_mapped_reads,
                     bool removeMultiMap=false, bool removeGappedRead=false, bool removeReverseRead=false);

// fomula to calculate RPKM
double RPKM(double read_in_gene, uLONG gene_len, uLONG total_reads);

// read chrLen from sam file
void readChrLen(ifstream &SAM, StringToUL &chrLen);

// calculate RPKM from sam file
void samRPKM( string samFileName, bool removeMultiMap=false, bool removeGappedRead=false, bool removeReverseRead=false);



class SamAlign
{
public:
    SamAlign() = default;

    SamAlign(string samFileName)
    {
        pSAM.reset(new ifstream(samFileName, ifstream::in));
        if(not *pSAM)
        {
            cerr << "FATAL Error: open sam file " << samFileName << " failed\n" << endl;
            pSAM.reset();
        }else{
            if( readChrLen() != NORMAL )
                pSAM.reset();
        }
    }

    SamAlign(ifstream *SAM)
    {
        this->pSAM.reset(SAM);
        if( readChrLen() != NORMAL )
            pSAM.reset();
    }

    ~SamAlign()
    {
        if(pSAM && *pSAM)
            pSAM->close();
    }

    void setIn(ifstream *SAM)
    {
        cout << "Set IN " << endl;
        if(pSAM && *pSAM)
            pSAM->close();
        pSAM.reset(SAM);
        if( readChrLen() != NORMAL )
            pSAM.reset();
    }

    void printChrLen()
    {
        cout << "# ************Chromosome Length**************" << endl;
        for(auto p: chrLen)
            cout << p.first << ":   " << p.second << endl;
    }

    // if the current pSAM is good
    bool is_good(){ if(pSAM and *pSAM) return true; else return false; }
    // read next few reads
    CStringMatrix readFewReads(size_t num=5);

    //calcRPKM
    void calcRPKM(bool removeMultiMap=true, bool removeGappedRead=true, bool removeReverseRead=true);

    //get chrLen
    const StringToUL& getChrLen(){ return chrLen; }
    //get RPKM
    const StringToD & getRPKM(){ return chrRPKM; }
    //get chrSingleRead
    const StringToUL & getChrSingleRead(){ return chrSingleRead; }
    //get RPKM
    const StringToDA & getChrMultiRead(){ return chrMultiRead; }
    //get total num
    uLONG getTotalReadNum(){ return total_mapped_reads; }

    //calcBD/RT
    void calcBDRT(bool removeMultiMap=true, bool removeGappedRead=true, bool removeReverseRead=true);

    //output BD/RT information
    void writeBDRT(ofstream &RT_OUT, ofstream &BD_OUT);
    //get chr BD/RT
    pair<vector<uLONG>,vector<uLONG>>  getChrBDRT( const string &chrName );


private:
    shared_ptr<ifstream> pSAM;      // pointer to sam iftream
    StringToUL chrLen;              // map<string, usigned long> to record chr length
    StringToUL chrSingleRead;       // uniq mapped reads
    StringToDA chrMultiRead;        // mutiple mapped reads
    uLONG total_mapped_reads = 0UL; // total mapped reads
    StringToD chrRPKM;              // RPKM of chromosome
    StringToULMatrix BD;           // BD of chromosome
    StringToULMatrix RT;           // RT of chromosome

    // read chrLen from sam file
    int readChrLen();

};

#endif // CALCRPKM_H
