//form cluster from extensionOutput
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

//#include <omp.h>
#include "bloomfilter.h"

#define CLUSTERMAXSIZE 100

using namespace std;
using namespace boost;

struct tagUnit{
       string tag;
       bool revComp;
       };

struct groupUnit{
       vector<tagUnit> tagVec;
       BloomFilter* kmerSet_filter_p;
       BloomFilter* kmerSet_filter_twin_p;
      };

unordered_map<string, groupUnit> groupSet;
int kmer= 17;

string getComplementSeq(string sequence) {

string complementSeq="";
for(string::iterator i= sequence.begin(); i!= sequence.end(); i++) {
    switch(*i) {
           case 'A':
           case 'a': complementSeq += 'T';
                     break;
           case 'T':
           case 't': complementSeq += 'A';
                     break;
           case 'G':
           case 'g': complementSeq += 'C';
                     break;
           case 'C':
           case 'c': complementSeq += 'G';
                     break;
            default: complementSeq += 'N';
           }
   }
return complementSeq;
}

int decomposeIntoKmer(string seq, BloomFilter& kmerSet_filter, BloomFilter& kmerSet_filter_twin, string seqTag) {
if(seq.length() < kmer)
   return 0;
for(int i= 0; i< seq.length()- kmer +1; i++) {
    string kmerSeq= seq.substr(i, kmer);
    kmerSet_filter.add(kmerSeq);
    }
string seqTwin= getComplementSeq(seq);
reverse(seqTwin.begin(), seqTwin.end() );
for(int i= 0; i< seqTwin.length()- kmer +1; i++) {
    string kmerSeq= seqTwin.substr(i, kmer);
    kmerSet_filter_twin.add(kmerSeq);
    }
}

int decomposeIntoKmerCompare(string seq, string representativeTag, BloomFilter *bfP, BloomFilter *bfTwinP, bool &revComp) {

int similarCount= 0, similarCount_twin= 0;
if(seq.length() < kmer)
   return 0;

for(int i= 0; i< seq.length()- kmer +1; i++) {
    string kmerSeq= seq.substr(i, kmer);
    if( (*bfP).contains(kmerSeq ) ) {
       similarCount ++;
      }
    if( (*bfTwinP).contains(kmerSeq) )
       similarCount_twin ++;
   }
if(similarCount_twin > similarCount) {
   similarCount= similarCount_twin;
   revComp= 1;
  }
return similarCount;
}

bool compareByLength(const pair<string,string> &p1, const pair<string,string> &p2) {
return p1.second.length() > p2.second.length();
}

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
if(tag != ""){
  seqSet.insert(make_pair(tag, seq) );
  }
seqFile.close();
}

int main(int argc, char* argv[]) {

cerr<< "./extension2Cluster <longerSeq> <scaffold.fa> <lrGmap2Scaffold.sam> <lr.fa>" << endl;
assert(argc== 5);

string extensionSeqFileName= argv[1];
string scaffoldSeqFileName= argv[2];
string samFileName= argv[3];
string lrSeqFileName= argv[4];

unordered_map<string, string> scaffoldSeqSet, extensionSeqSet, lrSeqSet;
readSeqFile(scaffoldSeqFileName, scaffoldSeqSet);
readSeqFile(extensionSeqFileName, extensionSeqSet);
readSeqFile(lrSeqFileName, lrSeqSet);

unordered_map<string, set<string> >clusterSet;
unordered_map<string, int> strandSet;
unordered_set<string> tagSet;

ifstream samFile(argv[3]);
assert(samFile.is_open() );

set<string> visitedTagSet;
string line= "";
getline(samFile, line);
while( !samFile.eof() ) {
      if(line != "") {
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         string lrTag= contentVec[0];
         if(visitedTagSet.find(lrTag) != visitedTagSet.end() ) {
            getline(samFile, line);
            continue;
           }
         visitedTagSet.insert(lrTag);
         int strand= lexical_cast<int>(contentVec[1]);
         if(strand >=16)
            strand= 1;
         string scaffold= contentVec[2];
         string locus= scaffold;
         locus= locus.substr(locus.find("_")+1);
         locus= locus.substr(0, locus.find("_") );
         strandSet.insert(make_pair(lrTag, strand) );        
         unordered_map<string, set<string> >::iterator it;
         if( (it= clusterSet.find(locus) )== clusterSet.end() ) {
            set<string> lrTagSet;
            lrTagSet.insert(scaffold);
            lrTagSet.insert(lrTag);
            clusterSet.insert(make_pair(locus, lrTagSet) );
           }
         else
           it->second.insert(lrTag);
        }
      getline(samFile, line);
     }

unordered_set<string> extendedScaffoldTagSet;

for(unordered_map<string, string>::iterator i= extensionSeqSet.begin(); i!= extensionSeqSet.end(); i++) {
    string locus= i->first;
    string scaffold= locus.substr(0, locus.find("lr") );
    extendedScaffoldTagSet.insert(scaffold);
    locus= locus.substr(locus.find("_")+1);
    locus= locus.substr(0, locus.find("_") );
    unordered_map<string, set<string> >::iterator it;
    if( (it= clusterSet.find(locus) )== clusterSet.end() ) {
       set<string> lrTagSet;
       lrTagSet.insert(i->first);
       clusterSet.insert(make_pair(locus, lrTagSet) );
      }
    else
       it->second.insert(i->first);
   }


for(unordered_map<string, string>::iterator i= scaffoldSeqSet.begin(); i!= scaffoldSeqSet.end(); i++) {
    string scaffold= i->first;
    if(extendedScaffoldTagSet.find(scaffold) != extendedScaffoldTagSet.end() )
       continue;
    string locus= scaffold;
    locus= locus.substr(locus.find("_")+1);
    locus= locus.substr(0, locus.find("_") );
    tagSet.insert(scaffold);
    unordered_map<string, set<string> >::iterator it;
    if( (it= clusterSet.find(locus) )== clusterSet.end() ) {
       set<string> lrTagSet;
       lrTagSet.insert(scaffold);
       clusterSet.insert(make_pair(locus, lrTagSet) );
      }
    else
       it->second.insert(scaffold);
    }


cerr<< "finish input with size " << clusterSet.size() << endl;
for(unordered_map<string, set<string> >::iterator i= clusterSet.begin(); i!= clusterSet.end(); i++) {
//cout << i->first << endl;
//continue;
    vector<pair<string, string> >seqVec;
    for(set<string>::iterator j= i->second.begin(); j!= i->second.end(); j++) {
        string seq= "";
        unordered_map<string, string>::iterator it;
        if( (*j).find("Locus") != string::npos) {
           if( (*j).find("lr") != string::npos) {
              if( (it= extensionSeqSet.find(*j) )== extensionSeqSet.end() ) {
                 cerr<< "cannot find extension "<< *j << endl;
                 exit(1);
                }
             }
           else{
              if(extendedScaffoldTagSet.find(*j)!= extendedScaffoldTagSet.end())
                continue;
              if( (it= scaffoldSeqSet.find(*j) )== scaffoldSeqSet.end() ) {
                 cerr<< "cannot find scaffold "<< *j << endl;
                 exit(1);
                } 
             }
          }
        else{
            if( (it= lrSeqSet.find(*j) )== lrSeqSet.end() ) {
               continue;
               cerr<< "cannot find lr "<< *j << endl;
               exit(1);
              }
          }
        seqVec.push_back(make_pair(*j, it->second) );
       }
     if(seqVec.size()== 1){
        cout << (*(seqVec.begin() ) ).first << " head:\n" << endl; 
        continue;
       }
    if(seqVec.size() >CLUSTERMAXSIZE) {
       for(vector<pair<string, string> >::iterator j= seqVec.begin(); j!= seqVec.end(); j++) {
           if(j== seqVec.begin() ) {
              cout << j->first << " head:" << endl;
             }
           else
             cout << j->first << ":0" << "\t";
          }
       cout << endl;
       continue;
      }
   sort(seqVec.begin(), seqVec.end(), compareByLength);

   for(vector<pair<string,string> >::iterator j= seqVec.begin(); j!= seqVec.end(); j++) {
       groupUnit groupUnitTmp;
       string representativeTag= "";
       if(j== seqVec.begin() ) {
         representativeTag= j->first;
         cout << j->first << " head:" << endl;
         string representativeSeq= j->second;
         vector<tagUnit> tagVec;
         int insertSize= representativeSeq.length()- kmer +1 ;
         BloomFilter *bfP= (new BloomFilter(insertSize, 0.001) );
         BloomFilter *bfTwinP= (new BloomFilter(insertSize, 0.001) );
         decomposeIntoKmer(j->second, *bfP, *bfTwinP, j->first);
         groupUnitTmp= {tagVec, bfP, bfTwinP};
         }
      else{
         bool revComp= 0;
         decomposeIntoKmerCompare(j->second, representativeTag, groupUnitTmp.kmerSet_filter_p, groupUnitTmp.kmerSet_filter_twin_p, revComp);
         cout << j->first << ":" << revComp << "\t";
         }
     }
//    if(seqVec.size() >0)
       cout << endl;
continue;


    string head= "";
    for(set<string>::iterator j= i->second.begin(); j!= i->second.end(); j++) {
        if( (*j).find("Locus") != string::npos) {
           head= *j;
           break;
          }
       }
    cout << head << " head:" << endl;
    for(set<string>::iterator j= i->second.begin(); j!= i->second.end(); j++) {
        if(*j== head)
           continue;
        int direction= 0;


        if( (*j).find("Locus") == string::npos) {
           unordered_map<string, int>::iterator it;
           if( (it= strandSet.find(*j) )== strandSet.end() ) {
              cerr<< "cannot detect strand " << *j << endl;
              exit(1);
             }
           direction= it->second;
          }
        cout << *j << ":" << direction << "\t";
       }
    if(i->second.size() >0)
       cout << endl;
   }


samFile.close();
}
