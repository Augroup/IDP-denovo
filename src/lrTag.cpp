#include <stdio.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <numeric>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <unordered_set>

using namespace std;

int main(int argc, char* argv[]) {
cerr << "./renameHeader <FataFile> <tag> " << endl;

ifstream seqFile(argv[1]);
assert(seqFile.is_open()  && argc==3);


string line= "", seq= "", tag= "";
int i=1;
getline(seqFile, line);
while( !seqFile.eof() ) {
      if(line[0]=='>')  {
         if(tag != "") {
            cout<< ">" << argv[2] << i << " " << tag << "\n" << seq << endl;
            i++;
            seq= "";
           }
         tag= line.substr(1);
        }
      else{
          seq += line;
         }
      getline(seqFile, line);
     }

if(tag != "") {
   cout << ">" << argv[2] << i << "\n" << seq << endl;
  }

seqFile.close();
}
