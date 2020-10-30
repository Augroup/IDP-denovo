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

int main(int argc, char* argv[]) {

cerr<< "./lrCount <lr_gmap2_assembled_tx.sam>" << endl;

ifstream samFile(argv[1]);
assert( samFile.is_open() );

unordered_map<string,int> lrCountSet;

string line= "";
getline(samFile, line);
while( !samFile.eof() ) {
      if(line != "") {
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         string refTag= contentVec[2];
         if(refTag != "*") {
            unordered_map<string,int>::iterator it;
            if( (it= lrCountSet.find(refTag) )== lrCountSet.end() )
               lrCountSet.insert(make_pair(refTag, 1) );
            else
               it->second ++;
           }
        }
      getline(samFile, line);
     }

cout<< "#transcript\tcount"<< endl;
for(unordered_map<string,int>::iterator i= lrCountSet.begin(); i!= lrCountSet.end(); i++) 
    cout << i->first << "\t" << i->second << endl;







samFile.close();
}
