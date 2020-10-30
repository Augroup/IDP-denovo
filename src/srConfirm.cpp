//hisat align SR to pseudo-loci to get sam file 
//input is "gapInfo" to check if gaps are supported by SR alignment
//updated on 6/5/2016, vote for splice sites with SR alignment
//updated on 6/28/2016, for uniq mapping "NH:i:1"
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

#define AVGCOVCUTOFF 10
#define FLANKCUTOFF 5
#define DEVIATIONCUTOFF 10

struct regionUnit{
      vector<int> skipVec;
      vector<int> mapVec;
      };


unordered_map<string, regionUnit>alignmentSet;
unordered_map<string, string> seqSet;
unordered_map<string, unordered_map<int, int> >prevSiteSet, postSiteSet;

int readSeqFile(string seqFileName, unordered_map<string, string> &seqSet) {

ifstream seqFile(seqFileName.c_str() );
assert( seqFile.is_open() );

string line= "", tag= "", seq= "";
getline( seqFile, line);
while( !seqFile.eof() ) {
      if(line[0]== '>') {
         if(tag != "") {
            seqSet.insert(make_pair(tag, seq) );
            vector<int> skipVec(seq.length(), 0);
            vector<int> mapVec(seq.length(), 0);
            regionUnit regionUnitTmp= {skipVec, mapVec};
            alignmentSet.insert(make_pair(tag, regionUnitTmp) );
           }
         tag= line.substr(1, line.find(" ")-1);
         seq= "";
         }
      else
          seq+= line;
      getline( seqFile, line);
     }
if(tag != ""){
   seqSet.insert(make_pair(tag, seq) );
   vector<int> skipVec(seq.length(), 0);
   vector<int> mapVec(seq.length(), 0);
   regionUnit regionUnitTmp= {skipVec, mapVec};
   alignmentSet.insert(make_pair(tag, regionUnitTmp) );
  }
seqFile.close();
cerr<< "size:" << alignmentSet.size() << endl;
}

int add2Skip(vector<int>& skipVec, int beginPos, int length) {
//cout << skipVec.size() << "\t" << beginPos+ length << endl;
assert(skipVec.size() >= beginPos + length);
for(vector<int>::iterator i= skipVec.begin()+ beginPos; i!= skipVec.begin() + beginPos + length && i!= skipVec.end();  i++){
    (*i)++;
   }
}    

int add2Map(vector<int>& mapVec, int beginPos, int length) {
//cout << mapVec.size() << "\t" << beginPos+ length << endl;
//assert(mapVec.size() >= beginPos + length);
for(vector<int>::iterator i= mapVec.begin()+ beginPos; i!= mapVec.begin() + beginPos + length && i!= mapVec.end(); i++){
    (*i)++;
    }
}

int saveSite(string refTag, int position, int gapLength) {

unordered_map<string, unordered_map<int, int> >::iterator it;
if( (it= prevSiteSet.find(refTag) )== prevSiteSet.end() ) {
   unordered_map<int, int> site;
   site.insert(make_pair(position, 1) );
   prevSiteSet.insert(make_pair(refTag, site) );
  }
else {
   unordered_map<int,int>::iterator jt;
   if( (jt= it->second.find(position) )== it->second.end() ) {
      it->second.insert(make_pair(position, 1) );
     }
   else
     jt->second ++;
  }

int postPosition= position + gapLength;
if( (it= postSiteSet.find(refTag) )== postSiteSet.end() ) {
   unordered_map<int, int> site;
   site.insert(make_pair(postPosition, 1) );
   postSiteSet.insert(make_pair(refTag, site) );
   }
else {
     unordered_map<int,int>::iterator jt;
      if( (jt= it->second.find(postPosition) )== it->second.end() ) {
         it->second.insert(make_pair(postPosition, 1) );
        }
      else
         jt->second ++;
     }
}

int readSamFile(string samFileName) {

ifstream samFile(samFileName.c_str() );
assert( samFile.is_open() );

string line= "";
getline( samFile, line);
while( !samFile.eof() ) {
      if(line != "" && line.find("NH:i:1")!= string::npos) {
//cout<< line << endl; 
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );

         if(lexical_cast<int>(contentVec[4]) ==0) {
            getline( samFile, line);
            continue;
           }
         string cigar= contentVec[5];
/*
         if(cigar.find("N")== string::npos) {
            getline( samFile, line);
            continue;
           }
*/
         string refTag= contentVec[2];
//cout << refTag << "\t" << contentVec[0] << endl;
         unordered_map<string,  regionUnit>::iterator it;
         if( (it= alignmentSet.find(refTag) )== alignmentSet.end() ) {
            cerr<< "cannot find " << refTag << endl;
            exit(1);
           }
//cout << "size: " << it->second.mapVec.size() << endl;
         int beginPos= lexical_cast<int>(contentVec[3])-1;
         vector<string> lengthVec;
         split(lengthVec, cigar, is_any_of("DMIHNS") );
         int pLength= 0, pPosition= beginPos;
         for(vector<string>::iterator i= lengthVec.begin(); i!= lengthVec.end(); i++) {
             if(*i== "")
                continue;
             int length= lexical_cast<int>(*i);
             pLength += (*i).length();
             char op= lexical_cast<char>(cigar.substr(pLength, 1) );
             switch(op) {
                    case 'N':{add2Skip(it->second.skipVec, pPosition, length);
                              saveSite(refTag, pPosition, length);
                             }
                              break;
                    case 'M':
                    case 'D': add2Map(it->second.mapVec, pPosition, length);
                              break;
                    default: break;
                   }
             pLength += 1;
             if(op== 'M' || op== 'D' || op== 'N')
                pPosition += length;
            }
        }
      getline( samFile, line);
     }
samFile.close();
}

double getAvgCov(vector<int> alignmentVec, int startPos, int endPos) {
//cerr << "region:" << startPos << "\t" << endPos << "\t" << alignmentVec.size() << endl;
double avgCov= 0;
for(vector<int>::iterator i= alignmentVec.begin()+ startPos; i!= alignmentVec.begin()+ endPos+1; i++) 
    avgCov += *i;
avgCov= avgCov/ (double)(endPos- startPos +1);
return avgCov;
}

int readGpdFile(string gpdFileName) {

ifstream gpdFile(gpdFileName.c_str() );
assert( gpdFile.is_open() );

string line= "";
getline(gpdFile, line);
while( !gpdFile.eof() ) {
      if(line != "") {
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         string geneTag= contentVec[0];
         string startList= contentVec[9];
         if( count(startList.begin(), startList.end(), ',')<=1) {
cout << line << endl;
            getline(gpdFile, line);
            continue;
           }
//cout << line << endl;
         unordered_map<string, regionUnit >::iterator it;
         if( (it= alignmentSet.find(geneTag) )== alignmentSet.end() ){
            cerr << "cannot find alignment for " << geneTag << endl;
cout << line << endl;
            getline(gpdFile, line);
            continue;
           }
//cout << line << endl;
         string endList= contentVec[10];

         string startList_confirmed= startList.substr(0, startList.find(",")+1);
         startList= startList.substr(startList.find(",")+1);

         endList= endList.substr(0, endList.length()-1);
         string endList_tail= endList.substr(endList.find_last_of(",")+1) +","; 
         string endList_confirmed= "";
         endList= endList.substr(0, endList.find_last_of(",")+1);
//remove 1st in startList & last in endList
         vector<string> startVec, endVec;
         split(startVec, startList, is_any_of(",") );
         split(endVec, endList, is_any_of(",") );
         vector<string>::iterator j= endVec.begin();
         for(vector<string>::iterator i= startVec.begin(); i!= startVec.end(); i++) {
             if(*i== "")
                continue;
             int startPos= lexical_cast<int>(*i);
             int endPos= lexical_cast<int>(*j)-1;
//cout << it->second.size() << "\t" << endPos << "\t" << startPos << endl;
             double avgSkipCov= getAvgCov(it->second.skipVec, endPos, startPos);
             double avgMapCov= getAvgCov(it->second.mapVec, endPos, startPos);
             if(avgSkipCov >= AVGCOVCUTOFF && avgMapCov>= AVGCOVCUTOFF) {
                startList_confirmed += *i + ",";
                endList_confirmed += *j + ","; 
               }
//cout <<"avg: " <<  avgCov << endl;
//exit(1);
             j++;
            }
         endList_confirmed += endList_tail;
//cout << startList_confirmed << "\t" << endList_confirmed << endl;
         int nBlocks= count(endList_confirmed.begin(), endList_confirmed.end(),',');
         cout<< geneTag << "\t" << contentVec[1] << "\t" << geneTag << "\t" << contentVec[3] << "\t" << contentVec[4] << "\t" << contentVec[5] << "\t" << contentVec[4] << "\t" << contentVec[5] << "\t" << nBlocks << "\t" << startList_confirmed << "\t" << endList_confirmed << endl;
        }
      getline(gpdFile, line);
     }

gpdFile.close();
}

string compressSeq(string seq) {

string conciseSeq= "";
char prevNt= 'Z';
for(string::iterator i= seq.begin(); i!= seq.end(); i++) {
    char ch= boost::lexical_cast<char>(*i);
    ch= toupper(ch);
    if(prevNt == 'Z' || prevNt != ch) {
       conciseSeq += boost::lexical_cast<string>(ch);
       }
    prevNt= ch;
    }
if(seq.length() >0 && conciseSeq== "")
   conciseSeq= seq.substr(0,1);
return conciseSeq;
}

int tuneSite(string geneTag, int& beginPos, int& endPos) { 

unordered_map<string, unordered_map<int, int> >::iterator it;
if( (it= prevSiteSet.find(geneTag) )== prevSiteSet.end() ) {
   cerr<< "wrong cannot find prev site for " << geneTag << endl;
   exit(1);
  }
int cov= 0, prevSite= 0, postSite= 0;
for(unordered_map<int, int>::iterator i= it->second.begin(); i!= it->second.end(); i++) {
   if( abs(beginPos- i->first) <= DEVIATIONCUTOFF && i->second > cov) {
       cov= i->second;
       prevSite= i->first;
     }
  }
if(prevSite != 0)
   beginPos= prevSite;
cov= 0;

if( (it= postSiteSet.find(geneTag) )== postSiteSet.end() ) {
   cerr<< "wrong cannot find post site for " << geneTag << endl;
   exit(1);
  }
for(unordered_map<int, int>::iterator i= it->second.begin(); i!= it->second.end(); i++) {
   if( abs(endPos- i->first) <= DEVIATIONCUTOFF && i->second > cov) {
      cov= i->second;
      postSite= i->first;
      }
   }
if(postSite != 0) 
   endPos= postSite;
}

int readGapFile(string gapFileName) {

ifstream gapFile(gapFileName.c_str() );
assert( gapFile.is_open()  );
string line= "";
getline(gapFile, line);
while( !gapFile.eof() ) {
      if(line != "") {
         string geneTag= line.substr(0, line.find(":") );
//cerr<< line << "\t" << geneTag << endl;
         unordered_map<string, regionUnit >::iterator it;
         if( (it= alignmentSet.find(geneTag) )== alignmentSet.end() ){
            cerr<< line << "\ncannot find alignment for " << geneTag << endl;
            exit(1);
           }
         line= line.substr(line.find(":")+1);
         string lrTag= line.substr(0, line.find("\t") );
         line= line.substr(line.find("\t")+1);
         vector<string> contentVec;
         vector<string> confirmedVec;
         split(contentVec, line, is_any_of("\t") );
         for(vector<string>::iterator i= contentVec.begin(); i!= contentVec.end(); i++) {
             if(*i== "")
                continue;
             string range= (*i).substr((*i).find(":")+1); 
             int length= lexical_cast<int>(range.substr(0,range.find("-") ) );
             int beginPos= lexical_cast<int>(range.substr(range.find("-")+1) );
//cout << "for " << geneTag << "\t" << beginPos << " " << length << "\t" << it->second.skipVec.size() << "\t" << it->second.mapVec.size() << endl;
             double avgSkipCov= getAvgCov(it->second.skipVec, beginPos, beginPos + length);
//             double avgSkipCov= 10;
//             double avgMapCov= 10;
             double avgMapCov= getAvgCov(it->second.mapVec, beginPos, beginPos + length);
             double avgFlankMapCov_prev= getAvgCov(it->second.mapVec, beginPos- FLANKCUTOFF, beginPos+ FLANKCUTOFF);
             double avgFlankMapCov_post= getAvgCov(it->second.mapVec, beginPos +length- FLANKCUTOFF, beginPos+length + FLANKCUTOFF);
//             if(avgSkipCov >= AVGCOVCUTOFF && avgMapCov>= AVGCOVCUTOFF && avgFlankMapCov_prev >= AVGCOVCUTOFF && avgFlankMapCov_post >= AVGCOVCUTOFF)
             if(avgSkipCov >= AVGCOVCUTOFF && avgMapCov>= AVGCOVCUTOFF)
                confirmedVec.push_back(*i);
            }
         if(confirmedVec.size() >0) {
            cout << geneTag << ":" << lrTag << "\t";

            unordered_map<string, string>::iterator it;
            if( (it= seqSet.find(geneTag) )== seqSet.end() ) {
               cerr<< "cannot find " << geneTag << endl;
               exit(1);
              }

            for(vector<string>::iterator j= confirmedVec.begin(); j!= confirmedVec.end(); j++) {
//                cout << *j << "\t";
                string range= (*j).substr((*j).find(":")+1);
                int length=lexical_cast<int>(range.substr(0,range.find("-") ) );
                int beginPos= lexical_cast<int>(range.substr(range.find("-")+1) );
                int endPos= beginPos +length;
                tuneSite(geneTag, beginPos, endPos);

                string gap= (*j).substr(0, (*j).find(":")+1) + lexical_cast<string>(endPos- beginPos) + "-" + lexical_cast<string>(beginPos);

                string gapSeq= it->second.substr(beginPos, length);
                string compressedSeq= compressSeq(gapSeq);
                if(length<=970&& compressedSeq.length() *2 > gapSeq.length() ){
//                   cout << *j << "\t";
                   cout << gap << "\t";
                   cerr<<">" << geneTag << ":" << beginPos << ":" << length <<":" << compressedSeq.length() << "\n" << gapSeq << endl;
                  }

               }
            cout << endl;
           }
        }
      getline(gapFile, line);
     }

gapFile.close();
} 

int main(int argc, char* argv[]) {

cerr<< "./srConfirmGap <report_seq.fa> gapInfo <report.sam> " << endl;

assert(argc== 4);

string seqFileName= argv[1];
string gapFileName= argv[2];
string samFileName= argv[3];

readSeqFile(seqFileName, seqSet);
cerr << "finish reading seqFile" << endl;

readSamFile(samFileName);
cerr << "finish reading samFile" << endl;
readGapFile(gapFileName);
cerr << "finish reading gapFile" << endl;
}
