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

using namespace std;
using namespace boost;

int updatePos(vector<vector<string> > txVec, int leftMostPos) {

int addPos= (-1)*leftMostPos;
for(vector<vector<string> >::iterator i= txVec.begin(); i!= txVec.end(); i++) {
    cout<< (*i).at(0)<< "\t" << (*i).at(1) << "\t" << (*i).at(0) << "\t" << (*i).at(3) << "\t";
    int beginPos= lexical_cast<int>( (*i).at(4) )+ addPos;
    int endPos= lexical_cast<int>( (*i).at(5) )+ addPos;
    cout<< beginPos<< "\t"<< endPos<< "\t"<< beginPos << "\t" << endPos << "\t" << (*i).at(8)<< "\t";
    vector<string> startVec, endVec;
    split( startVec, (*i).at(9), is_any_of(",") );
    split( endVec, (*i).at(10), is_any_of(",") );
    for(vector<string>::iterator j= startVec.begin(); j!= startVec.end(); j++) {
        if(*j== "")
           continue;
        int pos= lexical_cast<int>(*j)+ addPos;
        cout << pos << ",";
       }
    cout << "\t";
    for(vector<string>::iterator j= endVec.begin(); j!= endVec.end(); j++) {
        if(*j== "")
          continue;
       int pos= lexical_cast<int>(*j)+ addPos;
       cout << pos << ",";
      }
    cout << endl;
   }

}

int main(int argc, char* argv[]) {

cerr<< "./updatedExtendGpd <extend.gpd>" << endl;

ifstream gpdFile(argv[1]);
assert( gpdFile.is_open() );

vector<vector<string> > txVec;

int leftMostPos= 0;
string prevGeneTag= "";
string line= "";
getline(gpdFile, line);
while( !gpdFile.eof() ) {
      if(line != "") {
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         string geneTag= contentVec[0];
         int startPos= lexical_cast<int>(contentVec[4]);
         if(prevGeneTag== "")
            prevGeneTag= geneTag;
         if(geneTag == prevGeneTag) {
            txVec.push_back(contentVec);          
            if(leftMostPos > startPos)
               leftMostPos= startPos;
           }
         else{
//cout <<":" << line << "\t" << leftMostPos << endl;
            updatePos(txVec, leftMostPos);
            vector<vector<string> >().swap(txVec);
            leftMostPos= startPos;
            txVec.push_back(contentVec);
           }
         prevGeneTag= geneTag;

        }
      getline(gpdFile, line);
     }
updatePos(txVec, leftMostPos);

gpdFile.close();
}
