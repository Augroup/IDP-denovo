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

//#include <omp.h>

using namespace std;
using namespace boost;

#define ALLOWANCE 10
#define MININTRONSIZE 68
#define MAXEXONSIZE 1700

struct txUnit{
       string content;
       int exonNum;
       vector<string> startSiteVec;      
       vector<string> endSiteVec;      
      };

unordered_set<string> startSiteSet_total, endSiteSet_total;

int formSiteSet(string siteList, string chrm, vector<string>& siteVec) {

vector<string> siteVec_tmp;
split(siteVec_tmp, siteList, is_any_of(",") );
for(vector<string>::iterator i= siteVec_tmp.begin(); i!= siteVec_tmp.end(); i++) {
    if(*i== "")
       continue;
    string site= chrm + ":" + *i;
    siteVec.push_back(site);
   }
}

bool compareByExonNum(const txUnit& t1, const txUnit& t2) {
return t1.exonNum > t2.exonNum;
}

int compareSite(string site, unordered_set<string>& siteSet_total, bool& novelSite) {
//cout <<"for site " << site << endl;
string chrm= site.substr(0, site.find(":") );
int pos= lexical_cast<int>(site.substr(site.find(":")+1) );
int deviation= 0;
for(deviation= 0; deviation<= ALLOWANCE; deviation++) {
    string querySite= chrm + ":" + lexical_cast<string>(pos+ deviation);

    if(siteSet_total.find(querySite)!= siteSet_total.end() ) {
       break;
      }
    if(deviation== 0)
       continue;
    querySite= chrm + ":" + lexical_cast<string>(pos- deviation);
    if(siteSet_total.find(querySite)!= siteSet_total.end() ) {
       break;
      }
   }

siteSet_total.insert(site);

if(deviation > ALLOWANCE) {
   novelSite= 1;
//   cout <<"cannot find "<< site << endl;
  }
else
   novelSite= 0;
}

int checkIntronSize(vector<string>& startSiteVec, vector<string>& endSiteVec, int& exonNum) {

vector<string> startSiteVec_tmp, endSiteVec_tmp;
startSiteVec_tmp.push_back( *(startSiteVec.begin())  );


vector<string>::iterator j= endSiteVec.begin();
for(vector<string>::iterator i= startSiteVec.begin()+1; i!= startSiteVec.end(); i++) {
    int intronEnd= lexical_cast<int>( (*i).substr( (*i).find(":")+1) ); 
    int intronBegin= lexical_cast<int>( (*j).substr( (*j).find(":")+1) ); 
//cout << "size " << intronEnd - intronBegin << endl;
    if(intronEnd - intronBegin >= MININTRONSIZE) {
       startSiteVec_tmp.push_back(*i);
       endSiteVec_tmp.push_back(*j);
      }
    else{
       vector<string>().swap(startSiteVec);
       vector<string>().swap(endSiteVec);
       return 1;
      }
    j++;
   }
endSiteVec_tmp.push_back( *(endSiteVec.rbegin() ) );

vector<string>().swap(startSiteVec);
vector<string>().swap(endSiteVec);

vector<string>::iterator n= endSiteVec_tmp.begin();
for(vector<string>::iterator m= startSiteVec_tmp.begin(); m!= startSiteVec_tmp.end(); m++) {
    int exonBegin= lexical_cast<int>( (*m).substr( (*m).find(":")+1) );
    int exonEnd= lexical_cast<int>( (*n).substr( (*n).find(":")+1) );
    if(exonEnd- exonBegin <= MAXEXONSIZE) {
       startSiteVec.push_back(*m);
       endSiteVec.push_back(*n);
      }
    else{
      vector<string>().swap(startSiteVec);
       vector<string>().swap(endSiteVec);

       break;
      }
    n++;
   }

//startSiteVec= startSiteVec_tmp;
//endSiteVec= endSiteVec_tmp;

exonNum= startSiteVec.size();
//cout << exonNum << endl;
//exit(1);
}

int main(int argc, char* argv[]) {

cerr<< "./uniq_spliceTx <splice.gpd>" << endl;

ifstream gpdFile(argv[1]);
assert( gpdFile.is_open() );

vector<txUnit> txUnitVec;

string line= "";
getline(gpdFile, line);
while( !gpdFile.eof() ) {
      if(line != "") {
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         string chrm= contentVec[2];
         int exonNum= lexical_cast<int>(contentVec[8]);
         string startSiteList= contentVec[9];
         string endSiteList= contentVec[10];

/*
         startSiteList= startSiteList.substr(startSiteList.find(",")+1);
         endSiteList= endSiteList.substr(0, endSiteList.length()-1);
         endSiteList= endSiteList.substr(0, endSiteList.find_last_of(",")+1);
*/
         vector<string> startSiteVec, endSiteVec;
         formSiteSet(startSiteList, chrm, startSiteVec);
         formSiteSet(endSiteList, chrm, endSiteVec);
//         checkIntronSize(startSiteVec, endSiteVec, exonNum);

         if(exonNum >1) {
            txUnit txUnitTmp= {line, exonNum, startSiteVec, endSiteVec};
            txUnitVec.push_back(txUnitTmp);
           }
         else
           cout << line << endl;
        }
      getline(gpdFile, line);
     }

sort(txUnitVec.begin(), txUnitVec.end(), compareByExonNum);
//cerr<<" sort "  << endl;

for(vector<txUnit>::iterator i= txUnitVec.begin(); i!= txUnitVec.end(); i++) {
    bool novelSite= 0, novelTx= 0;
    vector<string> startSiteVec= i->startSiteVec;
    vector<string> endSiteVec= i->endSiteVec;

     vector<string> startSiteVec_trimmed= startSiteVec, endSiteVec_trimmed= endSiteVec;
    int exonNum= 0;
//    checkIntronSize(startSiteVec_trimmed, endSiteVec_trimmed, exonNum);
    set<string> startSiteSet(startSiteVec_trimmed.begin(), startSiteVec_trimmed.end() );
    set<string> endSiteSet(endSiteVec_trimmed.begin(), endSiteVec_trimmed.end() );

    endSiteVec.pop_back();
    startSiteVec.erase(startSiteVec.begin() );
//cerr <<" size " << startSiteVec.size() << " " << endSiteVec.size() << endl;
    
    for(vector<string>::iterator j= startSiteVec.begin(); j!= startSiteVec.end(); j++) {
        compareSite(*j, startSiteSet_total, novelSite);
        if(novelSite== 1 && startSiteSet.find(*j)!= startSiteSet.end() )
           novelTx= 1;
       }
    for(vector<string>::iterator j= endSiteVec.begin(); j!= endSiteVec.end(); j++) {
        compareSite(*j, endSiteSet_total, novelSite);
        if(novelSite== 1 && endSiteSet.find(*j)!= endSiteSet.end() )
           novelTx= 1;
       }

    if(novelTx== 1)
       cout << i->content << endl;
   }


gpdFile.close();
}
