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
#include <boost/bind.hpp>
//#include <omp.h>
#include <boost/multi_array.hpp>

#include <sys/time.h>
#include <time.h>

using namespace std;
using namespace boost;

unordered_set<string> coveredTagSet;

int readSeqFile(string seqFileName, unordered_map<string, string> &seqSet) {

ifstream seqFile(seqFileName.c_str() );
assert( seqFile.is_open() );

string line= "", tag= "", seq= "";
getline( seqFile, line);
while( !seqFile.eof() ) {
      if(line[0]== '>') {
         if(tag != "") {
            if(coveredTagSet.find(tag) != coveredTagSet.end() )
               seqSet.insert(make_pair(tag, seq) );
           }
        tag= line.substr(1, line.find(" ")-1);
        seq= "";
       }
    else
      seq+= line;
    getline( seqFile, line);
    }
if(tag != "") {
   if(coveredTagSet.find(tag) != coveredTagSet.end() )
      seqSet.insert(make_pair(tag, seq) );
  }
seqFile.close();
}

int getUniqTag(string reportGpdFileName) {

ifstream gpdFile(reportGpdFileName.c_str() );
assert(gpdFile.is_open() );

string line= "";
getline(gpdFile, line);
while( !gpdFile.eof() ) {
      if(line != "") {
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         coveredTagSet.insert(contentVec[1]);
        }
      getline(gpdFile, line);
     }

gpdFile.close();
}  

int main(int argc, char* argv[]) {

cerr<< "./getConsensus <uniq_report.gpd> <updated_extend.gpd> combine_seq" <<endl;
assert(argc== 4);


string reportGpdFileName= argv[1];
string extendGpdFileName= argv[2];
string seqFileName= argv[3];

getUniqTag(reportGpdFileName);

unordered_map<string, string> seqSet;
readSeqFile(seqFileName, seqSet);

for(unordered_map<string,string>::iterator i= seqSet.begin(); i!= seqSet.end(); i++) {
    cout << ">" << i->first << "\n" << i->second << endl;
   }



}
