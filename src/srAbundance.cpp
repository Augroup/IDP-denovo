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

int runHisat(string srFileName_1, string srFileName_2, string pseudoRefFileName){

string cmd= "mkdir pseudoHisatIndex";
runCmd(cmd);
cmd= "hisat-build " + pseudoRefFileName + " pseudoHisatIndex >hisatLog_sr 2>hisatlog_sr";
runCmd(cmd);

cmd= "hisat -f -N1 --no-head --no-sq -k1 --no-unal -x pseudoHisatIndex -1 " + srFileName_1 + " -2 "+ srFileName_2 +" -S hisatOutput.sam -p " + lexical_cast<string>(nThreads) +" >hisatLog_sr 2>hisatlog_sr";
runCmd(cmd);
}

string getSrLength(string srFileName_1) {

ifstream srFile(srFileName_1.c_str() );
assert(srFile.is_open() );

string line= "";
getline(srFile, line);
getline(srFile, line);

srFile.close();
int srLength= line.length();

return lexical_cast<string>(srLength);
}

int getPythonPath(string& pythonPath) {

string cmd= "which python >pythonLog";
runCmd(cmd);
ifstream pythonFile("pythonLog");
assert(pythonFile.is_open() );

string line= "";
getline(pythonFile, line);
pythonPath= line.substr(line.find('\'')+1 );
pythonPath= pythonPath.substr(0, pythonPath.find('\'') );

pythonFile.close();
}

int main(int argc, char* argv[]) {

cerr<< "./srQuantify sr_1.fa sr_2.fa <report.gpd> <report_seq> python_path code_path SR_length(single) nThreads" << endl;

assert(argc== 6);

string srFileName_1= argv[1];
string srFileName_2= argv[2];
string gpdFileName= argv[3];
string pseudoRefFileName= argv[4];
//string pythonPath= argv[5];
//string codePath= argv[6];
nThreads= atoi(argv[5]);

runHisat(srFileName_1, srFileName_2, pseudoRefFileName);

string cmd= ("awk '{if($6!~/N/ && $1!~/@/) print}' hisatOutput.sam > hisatOutput1.sam");
runCmd(cmd);
string samFileName= "hisatOutput1.sam";
string srLength= getSrLength(srFileName_1);
cout <<"length of SRs is " << srLength << endl;

string codePath= "";
string pythonPath= "" ;
getPythonPath(pythonPath);
cout <<"python path " << pythonPath << endl;

string gpdPrefix= gpdFileName.substr(0, gpdFileName.find_last_of(".") );

//run parseRef.py
cmd=  codePath +  "parseRef.py " + gpdFileName +" "+ srLength + " " + lexical_cast<string>(MINJUNCTIONOVERLAP);
runCmd(cmd);

//run parseSAM
cmd=   codePath + "parseSAM_MT.py " + gpdPrefix + "_regions.gpd " + samFileName +  " " + lexical_cast<string>(nThreads) +  " " + pythonPath + " " + srLength + " " + lexical_cast<string>(MINJUNCTIONOVERLAP) + " >parseSAM_MT1.log";  
runCmd(cmd);

//exit(1);

//run awk
extractColumn(gpdFileName);

cmd= "mv refSeq_MLE_input.txt refSeq_MLE_input1.txt";
runCmd(cmd);

//run markknownTranscripts.py
cmd= codePath + "markKnownTranscripts.py refSeq_MLE_input1.txt positive_candidate_list1 refSeq_MLE_input_marked1.txt";
runCmd(cmd);

//run MLE_MT.py
cmd= codePath + "MLE_MT.py refSeq_MLE_input_marked1.txt refSeq_MLE_output1.txt " + lexical_cast<string>(nThreads) +  " " + pythonPath;
runCmd(cmd);

//run reformat.py
cmd= codePath + "reformat.py refSeq_MLE_output1.txt > refSeq_MLE_output1.txt_";
runCmd(cmd);

//run MLEout2tab.py
cmd= codePath + "MLEout2tab.py refSeq_MLE_output1.txt_ > refSeq_MLE_output1.tab";
runCmd(cmd);


}
