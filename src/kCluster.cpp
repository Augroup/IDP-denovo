//updated groupUnit structure just has vector<string>, use vector only before
//updated on 3/25/2016 multiple bloomfilter for each representativeSeq
//pointer store in structure groupUnit
//output if reverse complement(0/1 for has) is better for kmer counting 3/31/2016 
//updated on 4/7/2016 reverse complement and itself two bloom filters
//updated on 4/7/2016 add gettimeofday to report clock time
//updated on 5/19/2016 compression before clustering
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
#include <sys/time.h>

#include <omp.h>
#include "bloomfilter.h"

using namespace std;
using namespace boost;

#define KMERCOUNTCUTOFF 50
#define PVALCUTOFF 0.05

struct tagUnit{
        string tag;
        bool revComp;
        };

struct groupUnit{
        vector<tagUnit> tagVec;
        BloomFilter* kmerSet_filter_p;
        BloomFilter* kmerSet_filter_twin_p;
       };

double percentCutoff= 0.05;
unordered_map<string, groupUnit> groupSet; 
int mpNum= 32, kmer= 7;
//BloomFilter *kmerSet_filter_p;
//unordered_map<string, BloomFilter*> bfPointerSet;

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

int readSeqFile(string seqFileName, vector<pair<string, string> > &seqVec) {

ifstream seqFile(seqFileName.c_str() );
assert( seqFile.is_open() );

string line= "", tag= "", seq= "";
getline( seqFile, line);
while( !seqFile.eof() ) {
      if(line[0]== '>') {
          if(tag != "") {
             seq= compressSeq(seq);
             if(seq.length() >= kmer)
                seqVec.push_back(make_pair(tag, seq) );
             }
          tag= line.substr(1, line.find(" ")-1);
          seq= "";
          }
      else
         seq+= line;
      getline( seqFile, line);
      }
if(tag != "") {
   seq= compressSeq(seq);
   if(seq.length() >= kmer)
      seqVec.push_back(make_pair(tag, seq) );
   }
seqFile.close();
}

bool compareByLength(const pair<string,string> &p1, const pair<string,string> &p2) {
return p1.second.length() > p2.second.length();
}

int decomposeIntoKmer(string seq, BloomFilter& kmerSet_filter, BloomFilter& kmerSet_filter_twin, string seqTag) {
if(seq.length() < kmer)
   return 0;
//cout << "begin decompose" << endl;
for(int i= 0; i< seq.length()- kmer +1; i++) {
    string kmerSeq= seq.substr(i, kmer);
//    if(! kmerSet_filter.contains(kmerSeq ) )
       kmerSet_filter.add(kmerSeq);
    }
string seqTwin= getComplementSeq(seq);
reverse(seqTwin.begin(), seqTwin.end() );
for(int i= 0; i< seqTwin.length()- kmer +1; i++) {
    string kmerSeq= seqTwin.substr(i, kmer);
    kmerSet_filter_twin.add(kmerSeq);
   }
//cout << "finish decompose" << endl;
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
/*
string seqTwin= getComplementSeq(seq);
reverse(seqTwin.begin(), seqTwin.end() );
for(int i= 0; i< seqTwin.length()- kmer +1; i++) {
    string kmerSeq= seqTwin.substr(i, kmer);
    if( (*bfP).contains(kmerSeq ) )
       similarCount_twin ++;
   }
*/
if(similarCount_twin > similarCount) {
   similarCount= similarCount_twin;
   revComp= 1;
   }

//cout << "similarCount " << representativeTag <<" " << similarCount  <<" " << seq.length()-kmer+1 << endl;
return similarCount;
}    

int clusterKmerSeq(string seqTag,  string seq) {

omp_lock_t writeLock;
omp_init_lock(&writeLock);

double maxValue= INT_MIN;
string representativeTag= "";

vector<pair<string, groupUnit> >groupVec(groupSet.begin(), groupSet.end());
unordered_map<string, bool> revCompSet;

#pragma omp parallel
{
#pragma omp for
for(int j= 0; j< groupVec.size(); j++) {
    vector<pair<string, groupUnit > >::iterator i;
    if(j<= groupVec.size()/2)
       i= groupVec.begin() +j;
    else
       i= groupVec.end() - (groupVec.size()-j);

    bool revComp= 0;
    int similarKmerCount= decomposeIntoKmerCompare(seq,  i->first, i->second.kmerSet_filter_p, i->second.kmerSet_filter_twin_p, revComp);
    double percent= (double)similarKmerCount/(double)(seq.length()- kmer+1);
    if(percent >= percentCutoff) {
omp_set_lock(&writeLock);
       if(percent > maxValue) {
          maxValue= percent;
          representativeTag= i->first;
          revCompSet.insert(make_pair(representativeTag, revComp) );
          }
omp_unset_lock(&writeLock);
      }
   }
}
omp_destroy_lock(&writeLock);

cerr<< "representativeTag: " << seqTag << "\t" << representativeTag << "\t" << maxValue << endl;
if(representativeTag == "") {
   vector<tagUnit> tagVec;
   int insertSize= seq.length()- kmer+1;
   BloomFilter *bfP= (new BloomFilter(insertSize, 0.001) );
   BloomFilter *bfTwinP= (new BloomFilter(insertSize, 0.001) );
   decomposeIntoKmer(seq, *bfP, *bfTwinP, seqTag);
   groupUnit groupUnitTmp= {tagVec, bfP, bfTwinP};
   groupSet.insert(make_pair(seqTag, groupUnitTmp) );
   }
else{
     unordered_map<string, bool>::iterator jt;
     if( (jt= revCompSet.find(representativeTag) )== revCompSet.end() )
        cerr<< "cannot found revComp for " << representativeTag << endl;
     bool revComp= jt->second;
     tagUnit tagUnitTmp= {seqTag, revComp};
    groupSet.find(representativeTag)->second.tagVec.push_back(tagUnitTmp);
    }
}

int main(int argc, char* argv[]) {

cerr<< "./kmerCluster <seq.fa> kmer percentCutoff nThreads" << endl;
assert(argc== 5);

string seqFileName= argv[1];
kmer= atoi(argv[2]);
percentCutoff= atof(argv[3]);
mpNum= atoi(argv[4]);

struct timeval startClock, endClock;
gettimeofday(&startClock, NULL);

clock_t startT, endT;
startT= clock();

cerr << "mp available: " << omp_get_num_procs() << endl;
mpNum= (mpNum > omp_get_num_procs()? omp_get_num_procs(): mpNum);
omp_set_num_threads(mpNum);

//BloomFilter kmerSet_filter(pow(4,10), 0.00001);
cerr << "finish optimizaing " << endl;

//kmerSet_filter_p= &(kmerSet_filter);

vector<pair<string, string> > seqVec;
readSeqFile(seqFileName, seqVec);
sort(seqVec.begin(), seqVec.end(), compareByLength);

for(vector<pair<string,string> >::iterator i= seqVec.begin(); i!= seqVec.end(); i++) {
    if(i== seqVec.begin() ) {
       string representativeSeq= i->second;
       vector<tagUnit> tagVec;
//cout << "aaa" << endl;
       int insertSize= representativeSeq.length()- kmer +1 ;
       BloomFilter *bfP= (new BloomFilter(insertSize, 0.001) );
       BloomFilter *bfTwinP= (new BloomFilter(insertSize, 0.001) );
//       BloomFilter kmerSet_filter_lr(insertSize, 0.001);
//cout << "xxx" << endl;
       decomposeIntoKmer(i->second, *bfP, *bfTwinP, i->first);
//cout << "yyy" << endl;
//       BloomFilter* bfP= &(kmerSet_filter_lr);
//       bfPointerSet.insert(make_pair(i->first, bfP) );
//cout << "zzz" << endl;       
       groupUnit groupUnitTmp= {tagVec, bfP, bfTwinP};
       groupSet.insert(make_pair(i->first,groupUnitTmp) );
      }
    else{
       clusterKmerSeq(i->first, i->second );
      }
   }
cerr << "fnish clustering" << endl;
endT= clock();
float diff ((float)endT- (float)startT) ;
float seconds= diff/CLOCKS_PER_SEC;

gettimeofday(&endClock, NULL);

for(unordered_map<string, groupUnit>::iterator i= groupSet.begin(); i!= groupSet.end(); i++) {
    cout << i->first <<" head:" << endl;
    for(vector<tagUnit>::iterator j= i->second.tagVec.begin(); j!= i->second.tagVec.end(); j++) {
       cout << j->tag << ":" << j->revComp <<"\t";
       }
    cout << endl;
   }
cerr << "total cluster for k=" << kmer << " with cutoff " << percentCutoff <<  " : " << groupSet.size() << endl;
cerr << "running time(s) :" << seconds << endl;

long double clockSecond= endClock.tv_sec - startClock.tv_sec; 
long double clockUsecond= (endClock.tv_usec- startClock.tv_usec)/1000000;
long double mtime= clockSecond + clockUsecond;
cerr << "clock time(s) is "<< mtime << endl; 
}
