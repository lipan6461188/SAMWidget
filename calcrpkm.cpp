
#include "calcrpkm.h"
#include "parsestring.h"
#include "mainwindow.h"

void printRead(const vector<string> &read)
{
    for_each(read.begin(), read.end(), [](string readItem){ cout << readItem << " "; });
    cout << "\n";
}

bool read_is_mapped(const vector<string> &read)
{
    if( stoi(read.at(1)) & 4 )
    {
        //printRead(read);
        return false;
    }
    else
        return true;
}

bool read_is_gapped(const vector<string> &read)
{
    if( read.at(5).find('N') != string::npos )
    {
        //printRead(read);
        return true;
    }
    else
        return false;
}

bool read_is_reverse(const vector<string> &read)
{
    if ( stoi(read.at(1)) & 16 )
    {
        //printRead(read);
        return true;
    }
    else
        return false;
}

void addOneRead_RPKM(CStringMatrix &sameReadArray, StringToUL &chrSingleRead, StringToDA &chrMultiRead,
                     uLONG &total_mapped_reads, bool removeMultiMap, bool removeGappedRead, bool removeReverseRead)
{
    // collect the valid map
    vector< vector<string> > filteredVector;
    for(auto readInfo: sameReadArray)
    {
        if( not read_is_mapped(readInfo) )
            continue;
        if( removeGappedRead and read_is_gapped(readInfo))
            continue;
        if( removeReverseRead and read_is_reverse(readInfo) )
            continue;
        filteredVector.push_back(readInfo);
    }

    if( (removeMultiMap and filteredVector.size() > 1) or filteredVector.size() == 0 )
    {
        //skip this read
    }
    else{
        if( removeMultiMap )
            // it must map to one place
            ++chrSingleRead[ filteredVector.at(0).at(2) ];
        else{
            // it may map to multiple places
            size_t mapped_places = filteredVector.size();
            if( mapped_places == 1 )
                ++chrSingleRead[ filteredVector.at(0).at(2) ];
            else
                for(auto readInfo: filteredVector)
                    chrMultiRead[ readInfo.at(2) ].push_back( 1.0 / mapped_places );
        }
        ++total_mapped_reads;
    }
}

double RPKM(uLONG read_in_gene, uLONG gene_len, uLONG total_reads)
{
    return 1.0 * read_in_gene / ( (1.0 * gene_len / 1000) * (1.0 * total_reads/1000000) + 0.01 );
}

void readChrLen(ifstream &SAM, StringToUL &chrLen)
{
    string thisLine;
    vector<string> samItems;

    auto lastLinePos = SAM.tellg();
    while( getline(SAM, thisLine) )
    {
        if(thisLine.find("@SQ") == 0)
            // collect read length and name information
        {
            split(thisLine, '\t', samItems);
            vector<string> namePair, lenPair;
            split(samItems.at(1), ':', namePair);
            split(samItems.at(2), ':', lenPair);

            chrLen[ namePair[1] ] = stoul(lenPair[1]);

            continue;
        }else if(thisLine.find("@") == 0)
        {
            continue;
        }else{
            // line is an alignment
            SAM.seekg(lastLinePos);
        }
        // keep last line
        lastLinePos = SAM.tellg();
    }
}

int SamAlign::readChrLen()
{
    if(not pSAM)
    {
        cerr << "FATAL Error: pSam is NULL" << endl;
        return NULL_PTR;
    }
    if(not *pSAM)
    {
        cerr << "FATAL Error: Sam file is not open" << endl;
        return IO_ERROR;
    }

    chrLen.clear();
    chrMultiRead.clear();
    chrRPKM.clear();
    chrSingleRead.clear();
    total_mapped_reads = 0;

    ifstream &SAM = *pSAM;
    chrLen.clear();

    string thisLine;
    vector<string> samItems;

    auto lastLinePos = SAM.tellg();
    while( getline(SAM, thisLine) )
    {
        if(thisLine.find("@SQ") == 0)
            // collect read length and name information
        {
            split(thisLine, '\t', samItems);
            vector<string> namePair, lenPair;
            split(samItems.at(1), ':', namePair);
            split(samItems.at(2), ':', lenPair);

            chrLen[ namePair[1] ] = stoul(lenPair[1]);
            chrMultiRead[ namePair[1] ];
            chrRPKM[ namePair[1] ];
            chrSingleRead[ namePair[1] ];

        }else if(thisLine.find("@") == 0)
        {
            // do nothing
        }else{
            // line is an alignment
            SAM.seekg(lastLinePos);
            break;
        }
        // keep last line
        lastLinePos = SAM.tellg();
    }

    if(chrLen.size()==0)
    {
        cerr << "Unexpected Error: sam file has no header" << endl;
        return FORMAT_ERROR;
    }

    return NORMAL;
}

CStringMatrix SamAlign::readFewReads(size_t num)
{
    StringMatrix stringMatrix;
    string thisLine;
    size_t idx = 0;
    vector<string> samItems;

    if(not is_good())
        return stringMatrix;

    ifstream &SAM = *pSAM;

    // preserve location
    auto lastLinePos = SAM.tellg();

    while(getline(SAM, thisLine) and idx < num)
    {
        split(thisLine, '\t', samItems);
        stringMatrix.push_back(samItems);
        idx++;
    }

    // return back to location
    SAM.seekg(lastLinePos);

    return stringMatrix;
}

void samRPKM( string samFileName, bool removeMultiMap, bool removeGappedRead, bool removeReverseRead)
{
    ifstream SAM(samFileName, ifstream::in);
    if(not SAM)
    {
        cerr << "Fatal Error: open sam file " << samFileName << "failed\n";
        exit(IO_ERROR);
    }else{
        clog << "Start to read sam file " << samFileName << endl;
    }

    string thisLine;
    string currentReadName;
    vector< vector<string> > sameReadArray;

    map<string, uLONG> chrLen;
    map<string, uLONG> chrSingleRead;
    map< string, vector<double> > chrMultiRead;
    uLONG total_mapped_reads = 0;

    vector<string> samItems;
    uLONG count = 0;
    while( getline(SAM, thisLine) )
    {
        ++count;
        if(count % 100000 == 0)
            clog << "Have read " << count << " lines..." << endl;
        if(thisLine.find("@SQ") == 0)
            //collect read length and name information
        {
            split(thisLine, '\t', samItems);
            vector<string> namePair, lenPair;
            split(samItems.at(1), ':', namePair);
            split(samItems.at(2), ':', lenPair);

            chrLen[ namePair[1] ] = stoul(lenPair[1]);
            chrSingleRead[ namePair[1] ] = 0UL;
            chrMultiRead[ namePair[1] ];

         //   cout << namePair[1] << "\t" << lenPair[1] << endl;
            continue;
        }
        if(thisLine.find("@") == 0)
            continue;
        split(thisLine, '\t', samItems);
        if( samItems[0] == currentReadName )
        {   // it's a multi-map
            sameReadArray.push_back(samItems);
        }else{
            // process those reads
            addOneRead_RPKM(sameReadArray, chrSingleRead, chrMultiRead, total_mapped_reads, removeMultiMap, removeGappedRead, removeReverseRead);
            sameReadArray.clear();
            sameReadArray.push_back(samItems);
            currentReadName = samItems[0];
        }
        //break;
    }

    SAM.close();

    /* output statistic information */

    map<string, double> chrRPKM;

    cout << "#TR:" << total_mapped_reads << endl;


    for(auto singleMap: chrSingleRead)
    {
        string chrName = singleMap.first;
        uLONG single_mapped_reads = singleMap.second;
        uLONG mutiple_mapped_reads = chrMultiRead[chrName].size();
        double averaged_mutiple_mapped_reads = accumulate(chrMultiRead[chrName].begin(), chrMultiRead[chrName].end(), 0.0);

        chrRPKM[chrName] = RPKM(single_mapped_reads+mutiple_mapped_reads, chrLen[chrName], total_mapped_reads);
        cout << chrName << ":   " << chrLen[chrName] << "\t" << single_mapped_reads << "   " << mutiple_mapped_reads << "   " << averaged_mutiple_mapped_reads << "   " << chrRPKM[chrName] << endl;
    }
}

void SamAlign::calcRPKM(bool removeMultiMap, bool removeGappedRead, bool removeReverseRead)
{
    if(not is_good())
        return;

    clog << "Start to calculate RPKM..." << endl;
    // initinize
    for(auto &pair: chrMultiRead)
        pair.second.clear();
    for(auto &pair: chrSingleRead)
        pair.second = 0UL;
    for(auto &pair: chrRPKM)
        pair.second = 0;
    total_mapped_reads = 0UL;

    // record the start of sam
    ifstream &SAM = *pSAM;
    auto lastLinePos = SAM.tellg();

    string thisLine;
    string currentReadName;
    vector< vector<string> > sameReadArray;
    vector<string> samItems;
    uLONG count = 0;

    while( getline(SAM, thisLine) )
    {
        ++count;
        if(count % 100000 == 0)
            clog << "Have read " << count << " lines..." << endl;

        split(thisLine, '\t', samItems);
        if( samItems[0] == currentReadName )
        {   // it's a multi-map
            sameReadArray.push_back(samItems);
        }else{
            // process those reads
            addOneRead_RPKM(sameReadArray, chrSingleRead, chrMultiRead, total_mapped_reads, removeMultiMap, removeGappedRead, removeReverseRead);
            sameReadArray.clear();
            sameReadArray.push_back(samItems);
            currentReadName = samItems[0];
        }
    }

    // return back to start position
    //SAM.seekg(lastLinePos);
    SAM.close();

    /* output statistic information */
    for(auto singleMap: chrSingleRead)
    {
        string chrName = singleMap.first;
        uLONG single_mapped_reads = singleMap.second;
        uLONG mutiple_mapped_reads = chrMultiRead[chrName].size();
       // double averaged_mutiple_mapped_reads = accumulate(chrMultiRead[chrName].begin(), chrMultiRead[chrName].end(), 0.0);

        chrRPKM[chrName] = RPKM(single_mapped_reads+mutiple_mapped_reads, chrLen[chrName], total_mapped_reads);
       // cout << chrName << ":   " << chrLen[chrName] << "\t" << single_mapped_reads << "   " << mutiple_mapped_reads << "   " << averaged_mutiple_mapped_reads << "   " << chrRPKM[chrName] << endl;
    }
}










