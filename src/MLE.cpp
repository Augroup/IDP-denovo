//estimate exp for each tx, with input of sam output of sr aligned to pseudo loci seq & pseudo loci GPD file & single-ended read length
#include <string>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <map>
#include <tr1/unordered_map>
#include <vector>
#include <set>
#include <numeric>
#include <cmath>
#include <unistd.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace boost;

#define MINJUNCTIONOVERLAP 10

int nThreads= 24;

int runCmd(string cmd) {
cerr<< cmd << endl;
system(cmd.c_str() );
}

int extractColumn(string fileName) {
cout << "beginning awk" << endl;
ifstream gpdFile(fileName.c_str() );
ofstream listFile("positive_candidate_list1");
assert( gpdFile.is_open() && listFile.is_open() );

string line= "";
getline(gpdFile, line);
while( !gpdFile.eof() ) {
     if(line != "" && line[0] != '#') {
        vector<string> contentVec;
        if(line.find("\t")== string::npos)
           split(contentVec, line, is_any_of(" ") );
        else
           split(contentVec, line, is_any_of("\t") );
        listFile << contentVec[2] << "\t" << contentVec[1] << endl;
       }
     getline(gpdFile, line);
     }


gpdFile.close();
listFile.close();
}

int main(int argc, char* argv[]) {

cerr<< "./MLE <SR.sam> <report.gpd> python_path code_path SR_length(single)" << endl;

assert(argc== 6);

string samFileName= argv[1];
string gpdFileName= argv[2];
string pythonPath= argv[3];
string codePath= argv[4];
string srLength= argv[5];

string gpdPrefix= gpdFileName.substr(0, gpdFileName.find_last_of(".") );

//run parseRef.py
string cmd=  pythonPath + " "+  codePath +  "/parseRef.py " + gpdFileName +" "+ srLength + " " + lexical_cast<string>(MINJUNCTIONOVERLAP);
runCmd(cmd);

//run parseSAM
cmd=  pythonPath + " "+ codePath + "/parseSAM_MT.py " + gpdPrefix + "_regions.gpd " + samFileName +  " " + lexical_cast<string>(nThreads) +  " " + pythonPath + " " + srLength + " " + lexical_cast<string>(MINJUNCTIONOVERLAP) + " >parseSAM_MT1.log";  
runCmd(cmd);

//run awk
extractColumn(gpdFileName);

cmd= "mv refSeq_MLE_input.txt refSeq_MLE_input1.txt";
runCmd(cmd);

//run markknownTranscripts.py
cmd= pythonPath + " " + codePath + "/markKnownTranscripts.py refSeq_MLE_input1.txt positive_candidate_list1 refSeq_MLE_input_marked1.txt";
runCmd(cmd);

//run MLE_MT.py
cmd= pythonPath + " " + codePath + "/MLE_MT.py refSeq_MLE_input_marked1.txt refSeq_MLE_output1.txt " + lexical_cast<string>(nThreads) +  " " + pythonPath;
runCmd(cmd);

//run reformat.py
cmd= pythonPath + " " + codePath + "/reformat.py refSeq_MLE_output1.txt > refSeq_MLE_output1.txt_";
runCmd(cmd);

//run MLEout2tab.py
cmd= pythonPath + " " + codePath + "/MLEout2tab.py refSeq_MLE_output1.txt_ > refSeq_MLE_output1.tab";
runCmd(cmd);


}
