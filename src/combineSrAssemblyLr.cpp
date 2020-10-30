//combine extension seq with SR-assembly seq
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

int readSeqFile(string seqFileName, unordered_map<string, string> &seqSet) {

ifstream seqFile(seqFileName.c_str() );
assert( seqFile.is_open() );

string line= "", tag= "", seq= "";
getline( seqFile, line);
while( !seqFile.eof() ) {
      if(line[0]== '>') {
         if(tag != "") {
            seqSet.insert(make_pair(tag, seq) );
           }
         tag= line.substr(1, line.find(" ")-1);
         seq= "";
         }
       else
          seq+= line;
       getline( seqFile, line);
      }
if(tag != "")
  seqSet.insert(make_pair(tag, seq) );
    seqFile.close();
}

int main(int argc, char* argv[]) {

cerr<< "./combineSrAssembly <extendedSeq.fa> <SR-scaffold.fa> <unusedLr.fa> " << endl;
assert(argc== 4);

string extendedSeqFileName= argv[1], assemblyFileName= argv[2], lrFileName= argv[3];
unordered_map<string, string> extendedSeqSet, assemblySeqSet, lrSeqSet;
readSeqFile(extendedSeqFileName, extendedSeqSet);
readSeqFile(assemblyFileName, assemblySeqSet);
readSeqFile(lrFileName, lrSeqSet);

unordered_set<string> assemblyTagSet, usedLrTagSet;
for(unordered_map<string, string>::iterator i= extendedSeqSet.begin(); i!= extendedSeqSet.end(); i++) {
    cout << ">" << i->first << "\n" << i->second << endl;
    string tag= i->first.substr(0, i->first.find("lr") );
    assemblyTagSet.insert(tag);
/*   
    string lrTag= i->first.substr(i->first.find("lr") );
    string lrUnitTag= lrTag.substr(lrTag.find_last_of("lr") );
    usedLrTagSet.insert(lrUnitTag);
 
    if(lrUnitTag!= lrTag) {
       lrUnitTag= lrTag.substr(0, lrTag.find("lr") );
       usedLrTagSet.insert(lrUnitTag);
      }   
*/
   }

for(unordered_map<string, string>::iterator i= assemblySeqSet.begin(); i!= assemblySeqSet.end(); i++) {
    if(assemblyTagSet.find(i->first) == assemblyTagSet.end() )
       cout << ">" << i->first << "\n" << i->second << endl;
   }



for(unordered_map<string,string>::iterator i= lrSeqSet.begin(); i!= lrSeqSet.end(); i++) {
    cout << ">" << i->first << "\n" << i->second << endl;
   }
}    
