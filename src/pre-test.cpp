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

int runCmd(string cmd) {
cerr<< cmd << endl;
system(cmd.c_str() );
}

inline bool exists_test2 (const std::string& fileName) {
return ( access( fileName.c_str(), F_OK ) != -1 );
}

int removeFile(string fileName) {
if(remove(fileName.c_str() ) !=0)
   cerr<<"fail to delete file " << fileName << endl;
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

bool checkGmap( ) {

system("which gmap" );
system("gmap > gmapTest 2>gmapTestLog");

ifstream testFile("gmapTest");
assert(testFile.is_open() );

string line= "";
getline(testFile, line);
while( !testFile.eof() ) {
      if(line.find("Usage:") != string::npos)
         break;
      getline(testFile, line);
      }
testFile.close();
removeFile("gmapTest");
removeFile("gmapTestLog");

if(line.find("Usage:") != string::npos)
   return 1;
else 
   return 0;
}

bool checkCap3( ) {

system("which cap3" );
system("cap3 2>cap3TestLog");

ifstream testFile("cap3TestLog");
assert(testFile.is_open() );

string line= "";
getline(testFile, line);
while( !testFile.eof() ) {
      if(line.find("Usage:") != string::npos)
         break;
      getline(testFile, line);
      }
testFile.close();
removeFile("cap3TestLog");    

if(line.find("Usage:") != string::npos)
   return 1;
else
   return 0;
}

bool checkHisat( ) {

system("which hisat" );
system("hisat 2>hisatTestLog");

ifstream testFile("hisatTestLog");
assert(testFile.is_open() );

string line= "";
getline(testFile, line);
while( !testFile.eof() ) {
      if(line.find("Usage:") != string::npos)
          break;
       getline(testFile, line);
       }
testFile.close();
removeFile("hisatTestLog");

if(line.find("Usage:") != string::npos)
   return 1;
else
   return 0;
}

bool checkClustalo( ) {

system("which clustalo-1.2.0-Ubuntu-x86_64");
system("clustalo-1.2.0-Ubuntu-x86_64 -h  >clustaloTest");

ifstream testFile("clustaloTest");
assert(testFile.is_open() );

string line= "";
getline(testFile, line);
while( !testFile.eof() ) {
      if(line.find("Usage:") != string::npos)
         break;
      getline(testFile, line);
     }
testFile.close();
removeFile("clustaloTest");

if(line.find("Usage:") != string::npos)
   return 1;
else
   return 0;
}

bool checkIDP( ) {

system("which runIDP.py");
system("runIDP.py  >IDPTest");

ifstream testFile("IDPTest");
assert(testFile.is_open() );

string line= "";
getline(testFile, line);
while( !testFile.eof() ) {
      if(line.find("usage:") != string::npos)
         break;
      getline(testFile, line);
     }
testFile.close();
removeFile("IDPTest");

if(line.find("usage:") != string::npos)
   return 1;
else
   return 0;
}

bool checkPythonVersion( ) {

system("python --version 2>pythonVersionTest");
ifstream testFile("pythonVersionTest");
assert(testFile.is_open() );

string line= "";
getline(testFile, line);

testFile.close();
removeFile("pythonVersionTest");
cout << line << endl;
if(line.find("2.7") != string::npos)
   return 1;
else
   return 0;
}

int main(int argc, char* argv[]) {

cout <<"usage: pre-test config_file" << endl;
assert(argc== 2);

string configFileName= argv[1];

cout<<"=====TEST before Running IDP-denovo====" << endl;

char *pPath;
pPath = getenv ("PATH");
if (pPath!=NULL)
    printf ("environment variable PATH is: %s\n",pPath);

cout << "\n******TEST for Python Version ********* " << endl;
bool pythonResult= checkPythonVersion();
if(pythonResult== 0) {
   system("python --version");
   cout << "......................FAIL, not version 2.7"<< endl;
  }
else
   cout << "......................SUCCESS"<< endl;

cout<< "\n*******TEST for GMAP installation **********" << endl;
bool gmapResult= checkGmap();
if(gmapResult== 0) {
   system("gmap");
   cout << "......................FAIL"<< endl;
  }
else
   cout << "......................SUCCESS"<< endl;

cout<< "\n*******TEST for CAP3 **********" << endl;
bool cap3Result= checkCap3();
if(cap3Result== 0) {
   system("cap3");
   cout << "......................FAIL"<< endl;
  }
else
   cout << "......................SUCCESS"<< endl;

cout<< "\n*******TEST for HISAT **********" << endl;
bool hisatResult= checkHisat();
if(hisatResult== 0) {
   system("hisat");
   cout << "......................FAIL"<< endl;
  }
else
   cout << "......................SUCCESS"<< endl;

cout<< "\n*******TEST for Clustal Omega **********" << endl;
bool clustaloResult= checkClustalo();
if(clustaloResult== 0) {
   system("clustalo-1.2.0-Ubuntu-x86_64 -h ");
   cout << "......................FAIL"<< endl;
   }
else
  cout << "......................SUCCESS"<< endl;

cout<<"\n*******TEST for path to IDP executables *****" << endl;
bool IDPResult= checkIDP();
if(IDPResult== 0){
   system("runIDP.py");
   cout << "......................FAIL"<< endl;
  }
else
   cout << "......................SUCCESS"<< endl;

cout<<"\n*******TEST for path to IDP-denovo executables *****" << endl;  
system("pathtest");

cout<<"\n*******TEST for path to input files *****" << endl;  

readConfigFile(configFileName);

if(exists_test2(scaffoldFileName)== 0) {
   cout<<"Please check path to scaffold file "<< endl;
   exit(1);
   }
else
  cout<<"scaffold file exists. " << scaffoldFileName << endl;

if(exists_test2(lrFileName)== 0) {
   cout<<"Please check path to lr file "<< endl;
   exit(1);
   }
else
   cout<<"LR file exists. " << lrFileName << endl;

if(exists_test2(leftSrFileName)== 0) {
   cout<<"Please check path to SR mate 1s file "<< endl;
   exit(1);
   }
else
   cout<<"SR mate 1s file exists. " << leftSrFileName << endl;

if(exists_test2(rightSrFileName)== 0) {
   cout<<"Please check path to SR amte 2s file "<< endl;
   exit(1);
   }
else
  cout<<"SR mate 2s file exists. " << rightSrFileName << endl;


}
