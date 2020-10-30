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
#include <bitset>
#include <deque>
#include <omp.h>
#include <unordered_set>
#include <time.h>

#include <boost/bind.hpp>

using namespace std;
using namespace boost;

string nThreads= "1", kmer= "15", kCutoff= "0.05";
string scaffoldFileName= "test_data/scaffold.fa", lrFileName= "test_data/lr.fa", leftSrFileName= "left.fa", rightSrFileName= "right.fa";

int detectParameter(string str) {

if(str.find("scaffold=") != string::npos) 
   scaffoldFileName= str.substr(str.find("=")+1);

if(str.find("lr=") != string::npos)
   lrFileName= str.substr(str.find("=")+1);

if(str.find("leftSr=") != string::npos)
   leftSrFileName= str.substr(str.find("=")+1);

if(str.find("rightSr=") != string::npos)
   rightSrFileName= str.substr(str.find("=")+1);

if(str.find("nThreads=") != string::npos)
   nThreads= str.substr(str.find("=")+1);

if(str.find("kmer=") != string::npos)
   kmer= str.substr(str.find("=")+1);

if(str.find("kcutoff=") != string::npos)
   kCutoff= str.substr(str.find("=")+1);
}

int runCmd(string cmd) {
cerr<< cmd << endl;
system(cmd.c_str() );
}

int getSrLength(string srFileName) {

ifstream srFile(srFileName.c_str() );
assert(srFile.is_open() );

string line= "";
getline(srFile, line);
getline(srFile, line);
srFile.close();

return line.length();
}

int runIdpdenovo( ) {

string cmd= "";

cmd= "lrTag " + lrFileName + " lr > lr_input";
runCmd(cmd);

lrFileName= "lr_input";

system("mkdir scaffoldIndex");
//--------------scaffold extension
cmd= "gmap_build -D scaffoldIndex -d scaffoldIndex " + scaffoldFileName + " >gmapLog 2>gmaplog";
runCmd(cmd);

cmd= "gmap -D scaffoldIndex -d scaffoldIndex -B5 -A -f samse --no-sam-header --failed-input=unalignedLr -n1 -t" + nThreads + " " + lrFileName + " >lrGmap2Scaffold.sam 2>>gmaplog";
runCmd(cmd);

cmd= "scaffoldExtend lrGmap2Scaffold.sam " + scaffoldFileName + " " + lrFileName + " " + nThreads + " >extendLog 2>extendlog";
runCmd(cmd);

cmd= "combineSrAssemblyLr longerSeq "+ scaffoldFileName + " unusedLr >combineAssembly 2>>extendlog";
runCmd(cmd);
//cluster 
cmd= "extension2Cluster longerSeq "+ scaffoldFileName + " lrGmap2Scaffold.sam unusedLr > seq_cluster 2>Log"; 
runCmd(cmd);

cmd= "kCluster unalignedLr " + kmer + " " + kCutoff + " " + nThreads + " >>seq_cluster 2>clusterLog ";
runCmd(cmd);
//annotate
system("cat combineAssembly unalignedLr >combine_seq");

cmd= "annotate seq_cluster combine_seq " + nThreads + " >annotateResult 2>annotatelog";
runCmd(cmd);

//----------------SR confirm
cmd= "hisat-build report_seq report_seq >hisatLog 2>hisatlog "; 
runCmd(cmd);

cmd= "hisat -x report_seq -f -N1 --no-unal --no-hd --no-sq -S srHisat2Ref.sam -k1 -1 " + leftSrFileName + " -2 "+ rightSrFileName + " -p " + nThreads + " 2>>hisatlog";
runCmd(cmd);

cmd= "srConfirm report_seq gapInfo srHisat2Ref.sam > confirmed_gap 2>confirmLog";
runCmd(cmd);


cmd= "updatedExtendGpd extend.gpd > updated_extend.gpd";
runCmd(cmd);

cmd= "uniq_spliceTx report.gpd > uniq_report.gpd";
runCmd(cmd);

cmd= "getConsensus uniq_report.gpd updated_extend.gpd combine_seq > extendTx";
runCmd(cmd);

cmd= "cat extendTx largeCluster.fa unannotated.fa > idpdenovoSeq";
runCmd(cmd);
//-----------quantify

//-----------LR quantify
system("mkdir extendTx_gmapIndex");
cmd= "gmap_build -D extendTx_gmapIndex -d extendTx_gmapIndex extendTx >> gmapLog 2>>gmaplog";
runCmd(cmd);

cmd= "gmap -D extendTx_gmapIndex -d extendTx_gmapIndex -B5 -A -f samse --no-sam-header --nofails -n1 -t" + nThreads + " " + lrFileName + " >lrGmap2ExtendTx.sam 2>>gmaplog"; 
runCmd(cmd);

cmd="lrCount lrGmap2ExtendTx.sam > lr_quantify";
runCmd(cmd);

//----------SR quantify

cmd= "srAbundance " +  leftSrFileName + " " + rightSrFileName + " report.gpd report_seq " + nThreads + " >srAbundanceLog 2>srAbundancelog";
runCmd(cmd);

return 0;


}

int readConfigFile(string configFileName) {

ifstream configFile(configFileName.c_str() );
assert(configFile.is_open() );

string line= "";
getline(configFile, line);
while( !configFile.eof() ) {
      if(line.find("scaffold=") != string::npos) 
         scaffoldFileName= line.substr(line.find("=")+1);
      if(line.find("lr=") != string::npos) 
         lrFileName= line.substr(line.find("=")+1);
      if(line.find("leftSr=") != string::npos) 
         leftSrFileName= line.substr(line.find("=")+1);
      if(line.find("rightSr=") != string::npos) 
         rightSrFileName= line.substr(line.find("=")+1);
      if(line.find("nThreads=") != string::npos) 
         nThreads= line.substr(line.find("=")+1);
      if(line.find("kmer=") != string::npos) 
         kmer= line.substr(line.find("=")+1);
      if(line.find("kcutoff=") != string::npos) 
         kCutoff= line.substr(line.find("=")+1);
      getline(configFile, line);
     }

configFile.close();
}

inline bool exists_test2 (const std::string& fileName) {
return ( access( fileName.c_str(), F_OK ) != -1 );
}
/*
int changeEnvPath(  ) {

char cwd[1024];
if(getcwd(cwd, sizeof(cwd) ) !=NULL) 
   fprintf(stdout, "Current working dir: %s\n", cwd);
else
   perror("getcwd() error");

char* pPath;
pPath = getenv ("PATH");
if(pPath!=NULL)
   printf ("environment variable PATH is: %s\n",pPath);

string envPath= pPath, cwdPath= cwd;

string addedPath_total= "";
string addedPath= cwd + "/bin";
if(envPath.find(addedPath) == string::npos)
   addedPath_total += ":" + addedPath;

addedPath= cwd + "/plugins/BIN";
if(envPath.find(addedPath) == string::npos)
   addedPath_total += ":" + addedPath;

addedPath= cwd + "/IDP_0.1.9/bin";
if(envPath.find(addedPath) == string::npos)
   addedPath_total += ":" + cmd;
}
*/

int main(int argc, char* argv[]) {

//cerr<< "usage 1: ./idpdenvo scaffold=scaffold.fa lr=lr.fa leftSr=left.fa rightSr=right.fa " << endl;
cerr<< "usage 1: ./idpdenvo config_file " << endl;
//cerr<<" optional: ncpus=1 kmer=15 kcutoff=0.05" << endl;

if(argc <5 && argc != 2) {
   cerr<< "wrong argument\nPlease input scaffold.fa & lr.fa & leftSr.fa & rightSr.fa" << endl;
   exit(1);
  }

if(argc >= 5) {
   cerr<< "./idpdenovo ";
   for(int i= 1; i< argc; i++) {
       cerr<< argv[i] <<" ";
       detectParameter(lexical_cast<string>(argv[i]) );
      }
   cerr<< endl;
   }

if(argc== 2) {
   string configFileName= argv[1];
   readConfigFile(configFileName);
  }
cout <<"***************************" << endl;
cout <<"number of CPUs used: " << nThreads << endl;
//--------check existence--------
if(exists_test2(scaffoldFileName)== 0) {
   cout<<"Please check path to scaffold file "<< scaffoldFileName << endl;
   exit(1);
  }
else
   cout<<"scaffold file exists. " << scaffoldFileName << endl;

if(exists_test2(lrFileName)== 0) {
   cout<<"Please check path to lr file "<< lrFileName << endl;
   exit(1);
  }
else
   cout<<"LR file exists. " << lrFileName << endl;

if(exists_test2(leftSrFileName)== 0) {
   cout<<"Please check path to SR mate 1s file "<< leftSrFileName << endl;
   exit(1);
  }
else
  cout<<"SR mate 1s file exists. " << leftSrFileName << endl;

if(exists_test2(rightSrFileName)== 0) {
   cout<<"Please check path to SR amte 2s file "<< rightSrFileName << endl;
   exit(1);
   }
else
   cout<<"SR mate 2s file exists. " << rightSrFileName << endl;

cout <<"***************************" << endl;

system("mkdir tmp_output");
system("mkdir idpdenovo_output");
chdir("./tmp_output");
runIdpdenovo( );
chdir("../");
system("pwd");
system("cp tmp_output/confirmed_gap idpdenovo_output/");
system("cp tmp_output/lr_input idpdenovo_output/");
system("cp tmp_output/report_seq idpdenovo_output/");
system("cp tmp_output/report.gpd idpdenovo_output/");
system("cp tmp_output/combine_seq idpdenovo_output/");
system("cp tmp_output/seq_cluster idpdenovo_output/");
system("cp tmp_output/lr_quantify idpdenovo_output/");
system("cp tmp_output/refSeq_MLE_output1.tab idpdenovo_output/sr_quantify");


}
