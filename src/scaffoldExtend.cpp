//openMP
//updated on 2/21/2016 faster, use GMAP output as input rather than BLAST output
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

#define CONFIDENCECUTOFF 1
#define EVALCUTOFF  1e-20
#define ALIGNCUTOFF 50
#define EXTCUTOFF 1
#define COVCUTOFF 3
#define BITSET boost::dynamic_bitset<>
#define SOLIDCUTOFF 1
#define MAXSEQNUMFOREXTEND 300
#define MAPQCUTOFF 0
#define FPKMCUTOFF 1

namespace std {
template<> struct hash<dynamic_bitset<> > {
       size_t operator()(const dynamic_bitset<>& t) const {
                         string s;
                         to_string(t,s);
                         return hash<string>()(s);
                        }
      };
}

struct addedPartUnit{
        string rawSeq;
        unordered_map<string, string> prevAddedSet;
        unordered_map<string, string> prevCap3Set;
        unordered_map<string, string> postAddedSet;
        unordered_map<string, string> postCap3Set;
       };

struct contigUnit{
        string contigTag;
        string contigSeq;
        string representativeTag;
        int lrNum;
       };

int kmer= 15, mpNum= 16;
unordered_map<BITSET, int> kmerFreqSet, kmerPercentSet;
unordered_map<string, set<string> >blastPairSet;
unordered_map<string,string> lrSeqSet, oasesSeqSet;
unordered_map<string, addedPartUnit>contigAddedPartSet;

vector<pair<string, set<string> > >blastPairVec;

//unordered_map<BITSET, bool> visitedSet;

int process_mem_usage(double & vm_usage, double& resident_set) {

vm_usage= 0.0;
resident_set= 0.0;

ifstream stat_stream("/proc/self/stat", ios_base::in);

string pid, comm, state, ppid, pgrp, session, tty_nr;
string tpgid, flags, minflt, cminflt, majflt, cmajflt;
string utime, stime, cutime, cstime, priority, nice;
string O, itrealvalue, starttime;

unsigned long vsize;
long rss;

stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
>> utime >> stime >> cutime >> cstime >> priority >> nice
>> O >> itrealvalue >> starttime >> vsize >> rss;

stat_stream.close();

long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
vm_usage     = vsize / 1024.0;
resident_set = rss * page_size_kb;
}

int runCmd(string cmd) {
cerr<< cmd << endl;
system(cmd.c_str() );
}

int createNt(int value, BITSET& head) {

switch(value) {
       case 0: head= BITSET(string("00") );
               break;
       case 1: head= BITSET(string("01") );
               break;
       case 2: head= BITSET(string("10") );
               break;
       case 3: head= BITSET(string("11") );
               break;
      }
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

string convertBit2seq(BITSET kmerBitset) {
string seq= "";
assert(kmerBitset.size() %2== 0);
bitset<2> nt;
for(int i=  kmerBitset.size()-1; i>= 0; i -=2) {
    nt[0]= kmerBitset[i-1];
    nt[1]= kmerBitset[i];
    if(nt== 00)
       seq += "A";
    if(nt== 01)
       seq += "C";
    if(nt== 10)
       seq += "G";
    if(nt== 11)
       seq += "T";
    }
return seq;
}

BITSET combineBitset(BITSET prevB, BITSET postB) {

boost::dynamic_bitset<> prevCopy(prevB), postCopy(postB);
int totalSize= prevB.size() + postB.size();
prevCopy.resize(totalSize);
postCopy.resize(totalSize);
prevCopy <<= postB.size();
prevCopy |= postCopy;
return prevCopy;
}

BITSET convertBitSet(string kmerSeq) {

BITSET bitSeq;
BITSET ba(string("00") ), bc(string("01")), bg(string("10")), bt(string("11"));

for(string::iterator i= kmerSeq.begin(); i!= kmerSeq.end(); i++) {
    switch(*i) {
           case 'A':
           case 'a': bitSeq= combineBitset(bitSeq, ba);
                     break;
           case 'T':
           case 't': bitSeq= combineBitset(bitSeq, bt);
                     break;
           case 'G':
           case 'g': bitSeq= combineBitset(bitSeq, bg);
                     break;
           case 'C':
           case 'c': bitSeq= combineBitset(bitSeq, bc);
                     break;
           case 'N':
           case 'n': bitSeq= combineBitset(bitSeq, ba);
                     break;
         }     
   }
return bitSeq;
}

BITSET reverseBit(BITSET original) {

BITSET reversed(original.size() );
int j= 0;
for(int i= original.size()-1; i>= 0; i-=2) {
    reversed[j]= original[i-1];
    reversed[j+1]= original[i];
    j+= 2;
   }
return reversed;
}

BITSET findStoreKmerBitset(BITSET kmerBitset) {

BITSET twinKmerBitset= reverseBit(~kmerBitset);
BITSET storeKmerBitset= (kmerBitset< twinKmerBitset? kmerBitset:twinKmerBitset);
return storeKmerBitset;
}

bool detectRepeat(BITSET kmerBitset) {

set<int> valueSet;
for(int i= 0; i< kmerBitset.size(); i+=2) {
    int value= 2* kmerBitset[i+1] + kmerBitset[i];
    valueSet.insert(value);
    }
if(valueSet.size()< 2)
   return 1;
return 0;
}

int getPercent(BITSET prefix) {

int value= 2*prefix[1]+ prefix[0];
prefix.resize(prefix.size()+2);
prefix= prefix << 2;
BITSET head;
int totalFreq_unit= 0;
vector<pair<BITSET, int> >frKmerPercentVec;
for(int j= 0; j< 4; j++) {
    createNt(j, head);
    head.resize(prefix.size() );
    BITSET newKmer= prefix | head;
    BITSET storeKmerBitset= findStoreKmerBitset(newKmer);
    unordered_map<BITSET,int >::iterator it;
    int freq= 0;
    if( (it= kmerFreqSet.find(storeKmerBitset) ) != kmerFreqSet.end() )
       freq= it->second;
       if(freq !=0) {
          frKmerPercentVec.push_back(pair<BITSET, int>(newKmer,freq));
          totalFreq_unit += freq;
          }
       }
for(vector<pair<BITSET, int> >::iterator k= frKmerPercentVec.begin(); k!= frKmerPercentVec.end(); k++) {
    int percent= (int)(100*k->second/totalFreq_unit);
    kmerPercentSet.insert(pair<BITSET, int>(k->first, percent) );
//    visitedSet.insert(make_pair(k->first,0) );
    }
}

int reportPercent( ) {

unordered_set<BITSET> prefixSet;

for(unordered_map<BITSET,int>::iterator i= kmerFreqSet.begin(); i!= kmerFreqSet.end(); i++) {
    BITSET prefix= (i->first) >>2;
    prefix.resize(prefix.size()-2);
    if(prefixSet.find(prefix) == prefixSet.end() ) {
       prefixSet.insert(prefix);
       getPercent(prefix);
       }
    prefix= reverseBit( (~(i->first) ) )>>2;
    prefix.resize(prefix.size()-2);
    if(prefixSet.find(prefix) == prefixSet.end() ) {
       prefixSet.insert(prefix);
       getPercent(prefix);
      }
   }

ofstream percentFile("kmerPercent");
assert(percentFile.is_open()  );
for(unordered_map<BITSET, int>::iterator i= kmerPercentSet.begin(); i!= kmerPercentSet.end(); i++) {
    percentFile << convertBit2seq(i->first) << "\t" << i->second << endl;
    }

percentFile.close();
}

int readPercent( ) {

ifstream percentFile("kmerPercent");

string line= "";
getline(percentFile, line);
while( !percentFile.eof() ) {
      if(line != "") {
         string kmerSeq= line.substr(0, line.find("\t") );
         int percent= boost::lexical_cast<int>(line.substr(line.find("\t")+1) );
         BITSET kmerBitset= convertBitSet(kmerSeq);
         kmerPercentSet.insert(make_pair(kmerBitset, percent) );
//         visitedSet.insert(make_pair(kmerBitset,0) );
         }
     getline(percentFile, line);
    }
percentFile.close();
}

inline bool exist_test(const string& fileName) {
return (access(fileName.c_str(), F_OK) != -1);
}

int readJellyFile(string jellyFileName) {

ifstream jellyFile(jellyFileName.c_str() );
assert(jellyFile.is_open() );

vector<pair<string, int> >kmerFreqVec;
string line= "";
getline(jellyFile, line);
while( !jellyFile.eof() ) {
      if(line != "") {
         int freq= boost::lexical_cast<int>(line.substr(line.find("\t")+1) );
         string kmerSeq= line.substr(0, line.find("\t") );
         if(freq > COVCUTOFF)
            kmerFreqVec.push_back(make_pair(kmerSeq, freq) );
        }
      getline(jellyFile, line);
     }

omp_lock_t writelock;
omp_init_lock(&writelock);

#pragma omp parallel
{
#pragma omp for
for(vector<pair<string,int> >::iterator i= kmerFreqVec.begin(); i< kmerFreqVec.end(); i++) {
    BITSET kmerBitset= convertBitSet(i->first);
    BITSET storeBitset=findStoreKmerBitset(kmerBitset);
    int freq= i->second;
    if(detectRepeat(kmerBitset) ==1)
       freq /= 10;
       omp_set_lock(&writelock);
       kmerFreqSet.insert(pair<BITSET, int>(storeBitset, freq) );
       omp_unset_lock(&writelock);
      }
}
omp_destroy_lock(&writelock);

jellyFile.close();
}

int removeFile(string fileName) {
if(remove(fileName.c_str() ) !=0)
   cerr << "error deleting file " << fileName << endl;
else
   cerr << "successful deleting file " << fileName << endl;
}

int copyFile(string inputFileName, string outputFileName) {

ifstream inputFile(inputFileName.c_str() ); 
ofstream outputFile(outputFileName.c_str() );
assert( inputFile.is_open() && outputFile.is_open() );

string line= "";
getline(inputFile, line);
while( !inputFile.eof() ) {
      if(line != "") {
         outputFile << line << endl;
        }
      getline(inputFile, line);
     }

inputFile.close();
outputFile.close();
}    

int getBlastPair(string blastOutputFileName) {

ifstream blastOutputFile(blastOutputFileName.c_str() );
assert( blastOutputFile.is_open() );

set<string> queryTagSet;
set<string> multiAlignQueryTagSet;
unordered_map<string, set<string> >blastPairSetTmp;

string line= "";
getline(blastOutputFile, line);
while( !blastOutputFile.eof() ) {
      if(line != "") {
         string queryTag= line.substr(0, line.find("\t") );
         if(queryTagSet.find(queryTag) != queryTagSet.end() ) {
            getline(blastOutputFile, line);
            multiAlignQueryTagSet.insert(queryTag);
            continue;
           }
         queryTagSet.insert(queryTag);
         line= line.substr(line.find("\t")+1);
         string hitTag= line.substr(0, line.find("\t") );
         string confidence= hitTag.substr(hitTag.find("Confidence")+11);
         double confidenceNum= boost::lexical_cast<double>(confidence.substr(0, confidence.find("_") ) );
         if(confidenceNum < CONFIDENCECUTOFF) {
            getline(blastOutputFile, line);
            continue;
           }
         line= line.substr(0, line.find_last_of("\t") );
         double eval= boost::lexical_cast<double>(line.substr(line.find_last_of("\t")+1) );
         if(eval > EVALCUTOFF) {
            getline(blastOutputFile, line);
            continue;
           }
         unordered_map<string, set<string> >::iterator it;
         if( (it= blastPairSetTmp.find(hitTag) )== blastPairSetTmp.end() ) {
            set<string> queryTagSet;
            queryTagSet.insert(queryTag);
            blastPairSetTmp.insert(make_pair(hitTag, queryTagSet) );
           }
         else
           it->second.insert(queryTag);
        }
      getline(blastOutputFile, line);
     }

for(unordered_map<string, set<string> >::iterator i= blastPairSetTmp.begin(); i!= blastPairSetTmp.end(); i++) {
   set<string> tagSet;
   for(set<string>::iterator j= i->second.begin(); j!= i->second.end(); j++) {
       if(multiAlignQueryTagSet.find(*j)== multiAlignQueryTagSet.end() )
          tagSet.insert(*j);
      }
   if(tagSet.size() >0) {
//      blastPairSet.insert(make_pair(i->first, tagSet) );
      blastPairVec.push_back(make_pair(i->first, tagSet) );
     }
  }

blastOutputFile.close();
}

int getGmapPair(string gmapFileName) {

ifstream gmapFile(gmapFileName.c_str() );
assert( gmapFile.is_open() );

unordered_map<string,string> unusedLrSeqSet;
unordered_set<string> visitedQueryTagSet;
string line= "";
getline(gmapFile, line);
while( !gmapFile.eof() ) {
//cout << line << endl;    
      if(line != "") {
         vector<string>contentVec;
         boost::split(contentVec, line, boost::is_any_of("\t") );
         int posCount= 0;
         string queryTag= "", sbjctTag= "", cigar= "", seq= "";
         int position= 0, mapQ= 0;
         for(vector<string>::iterator i= contentVec.begin(); i!= contentVec.end(); i++) {
//cout << *i << endl;
             switch(posCount) {
                    case 0: queryTag= *i;
                            break;
                    case 2: sbjctTag= *i;
                            break;
                    case 3: position= boost::lexical_cast<int>(*i);
                            break;
                    case 4: mapQ= boost::lexical_cast<int>(*i);
                            break;
                    case 5: cigar= *i;
                            break;
                    case 9: seq= *i;
                            break;
                    default: break;        
                   }
             posCount ++;
            }
         if(sbjctTag== "*") {
            getline(gmapFile, line);
            continue;
           }
         double fpkm= boost::lexical_cast<double>(sbjctTag.substr(sbjctTag.find_last_of("_")+1 ) );
         if(visitedQueryTagSet.find(queryTag) != visitedQueryTagSet.end() ) {
            getline(gmapFile, line);
            continue;
           }
         visitedQueryTagSet.insert(queryTag);
         if(mapQ < MAPQCUTOFF || cigar.find("N")!= string::npos || cigar.find("S")== string::npos || fpkm < FPKMCUTOFF) {
            unusedLrSeqSet.insert(make_pair(queryTag, seq) );
            getline(gmapFile, line);
            continue;
           }
//         sbjctTag= sbjctTag.substr(0, sbjctTag.find_last_of("_") );
         int txLength;
         unordered_map<string, string>::iterator kt;
         if( (kt= oasesSeqSet.find(sbjctTag) ) == oasesSeqSet.end() ) {
            cerr<< sbjctTag << "  not exist" << endl;
            exit(1);
           }
         else
           txLength= kt->second.length();
         string rawSeq= kt->second;
         vector<string> numFragVec;
         boost::split(numFragVec, cigar, boost::is_any_of("MIDNSH") );
         int alignLength= 0, refAlignLength= 0, mutualAlignLength= 0, prevClip= 0, postClip= 0;
         string::iterator it= cigar.begin();
         for(vector<string>::iterator i= numFragVec.begin(); i!= numFragVec.end(); i++) {
             if(*i== "")
                continue;
             it += (*i).length();
             if(*it== 'M')
                mutualAlignLength += boost::lexical_cast<int>(*i);
             if(*it== 'M' || *it== 'D')
                refAlignLength += boost::lexical_cast<int>(*i);
             if(*it != 'D')
                alignLength += boost::lexical_cast<int>(*i);
            if(*it== 'S') {
               if(i== numFragVec.begin() )
                  prevClip= boost::lexical_cast<int>(*i);
               else
                  postClip= boost::lexical_cast<int>(*i);
              }
           it++;
          }
        if(mutualAlignLength < ALIGNCUTOFF) {
           unusedLrSeqSet.insert(make_pair(queryTag, seq) );
           getline(gmapFile, line);
           continue;
          }
        unordered_map<string, addedPartUnit>::iterator jt;
        if(postClip - (txLength- position- refAlignLength) >=30) {
           string postCap3Part= seq.substr(seq.length()- postClip + (txLength -position - refAlignLength) );
           string postAddedPart= seq.substr(prevClip);
           if( (jt= contigAddedPartSet.find(sbjctTag) )== contigAddedPartSet.end() ) {
              unordered_map<string, string> prevAddedSet, postAddedSet, prevCap3Set, postCap3Set;
              postAddedSet.insert(make_pair(queryTag,postAddedPart) );
              addedPartUnit addedPartUnitTmp= {rawSeq,prevAddedSet, prevCap3Set, postAddedSet, postCap3Set};
              contigAddedPartSet.insert(make_pair(sbjctTag, addedPartUnitTmp));
             }
           else {
              jt->second.postAddedSet.insert(make_pair(queryTag,postAddedPart) );
              jt->second.postCap3Set.insert(make_pair(queryTag,postCap3Part) );
              }
          }
       if(prevClip - position >=30) {
          string prevCap3Part= seq.substr(0, prevClip- position);
          string prevAddedPart= seq.substr(0, seq.length()- postClip);
          if( (jt= contigAddedPartSet.find(sbjctTag) )== contigAddedPartSet.end() ) {
             unordered_map<string, string> prevAddedSet, postAddedSet, prevCap3Set, postCap3Set;
             postAddedSet.insert(make_pair(queryTag,prevAddedPart) );
             postCap3Set.insert(make_pair(queryTag,prevCap3Part) );
             addedPartUnit addedPartUnitTmp= {rawSeq,prevAddedSet, prevCap3Set, postAddedSet, postCap3Set};
             contigAddedPartSet.insert(make_pair(sbjctTag, addedPartUnitTmp));
            }
         else {
            jt->second.prevAddedSet.insert(make_pair(queryTag,prevAddedPart) );
            jt->second.prevCap3Set.insert(make_pair(queryTag,prevCap3Part) );
            }
         }

       }
      getline(gmapFile, line);
     }
gmapFile.close();

ofstream unusedLrFile("unusedLr");
assert(unusedLrFile.is_open() );
for(unordered_map<string,string>::iterator i= unusedLrSeqSet.begin(); i!= unusedLrSeqSet.end(); i++) {
    unusedLrFile<< ">" << i->first << "\n" << i->second << endl;
   }
unusedLrFile.close();
}    


int readSeqFile(string seqFileName, unordered_map<string, string> &seqSet) {

ifstream seqFile(seqFileName.c_str() );
if(! (seqFile.is_open() ) )
   cerr<< "not such file " << seqFileName << endl;
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

bool testContentExist(string fileName) {

bool contentExist= 0;
ifstream file(fileName.c_str() );
assert( file.is_open() );

string line= "";
getline(file, line);
while( !file.eof() ) {
      if(line != "")
         contentExist= 1;
      getline(file, line);
     }
file.close();
return contentExist;
}

string readClustaloFile(string contigTag, string partSeqTag,string fileNameTag){

unordered_map<string, string> alignmentSet;
readSeqFile("clustaloOutput_"+ fileNameTag, alignmentSet);

int prevSkip= 0, postSkip= 0;
string contigSeq= alignmentSet.find(contigTag)->second, partSeq= alignmentSet.find(partSeqTag)->second;
if(contigSeq.substr(0, 10) == string(10, '-') ) {

   string alignment= contigSeq.substr(contigSeq.find_first_not_of("-") );
   prevSkip= contigSeq.length()- alignment.length();
  }
else if(contigSeq.substr(contigSeq.length()-10)== string(10, '-') ) {
        string alignment= contigSeq.substr(0, contigSeq.find_last_not_of("-") );
        postSkip= contigSeq.length()- alignment.length();
       }
assert(contigSeq.length()== partSeq.length() );

string correctedSeq= "";
string::iterator j= partSeq.begin() + prevSkip;
for(string::iterator i= contigSeq.begin()+ prevSkip; i!= contigSeq.end()- postSkip; i++) {
    if(*i != '-')   
       correctedSeq += *i;
    j++;
   }
if(prevSkip != 0) 
   correctedSeq = partSeq.substr(0, prevSkip)+ correctedSeq;
if(postSkip != 0) 
   correctedSeq += partSeq.substr(partSeq.length()- postSkip);
//cout << correctedSeq << endl;
return correctedSeq;
}

int getPrevAttach(string logFileName, string contigFileName, unordered_map<string, string> partSeqSet) {

ifstream logFile(logFileName.c_str() );
assert( logFile.is_open() );

unordered_map<string, string> representativeSet;

string line= "";
getline(logFile, line);
while( !logFile.eof() ) {
      if(line.find("DETAILED") != string::npos)
         break;
      if(line.find("Contig") != string::npos) {
         line= line.substr(line.find("Contig")+ 7);
         string contigTag= "Contig"+ line.substr(0, line.find(" ") );
         set<string> lrTagSet;
         getline(logFile, line);
         string representativeTag= line.substr(0, line.length()-1);
         representativeSet.insert(make_pair(contigTag, representativeTag) );
cout << contigTag << " " << representativeTag << endl;
        }
   getline(logFile, line);
   }

logFile.close();

unordered_map<string, string> contigSeqSet;
readSeqFile(contigFileName, contigSeqSet);

for(unordered_map<string,string>::iterator i= representativeSet.begin(); i!= representativeSet.end(); i++) {
    string contigSeq= contigSeqSet.find(i->first)->second;
    string lrPartSeq= partSeqSet.find(i->second)->second;
    transform(lrPartSeq.begin(), lrPartSeq.end(), lrPartSeq.begin(), ::toupper);

    ofstream clustaloInputFile("clustaloSeq");
    assert(clustaloInputFile.is_open() );
    clustaloInputFile << ">" << i->first << "\n" << contigSeq << endl;
    clustaloInputFile << ">" << i->second << "\n" << lrPartSeq << endl;
    clustaloInputFile.close();
    string cmd= "clustalo-1.2.0-Ubuntu-x86_64  -i clustaloSeq  -o clustaloOutput --force --threads=24";
    runCmd(cmd);
    string correctedSeq= readClustaloFile(i->first, i->second, "1");
cout << ">" << i->second << "\n" << correctedSeq << endl;    
//exit(1);
   }
}

int readGmapFile(int txLength, unordered_map<string, string>& prevPartSeqSet, unordered_map<string, string>& postPartSeqSet, string fileNameTag) {

ifstream gmapFile("gmapOutput" + fileNameTag);
ofstream prevCap3InputFile("prevCap3" + fileNameTag), postCap3InputFile("postCap3" + fileNameTag);
assert(gmapFile.is_open() && prevCap3InputFile.is_open() && postCap3InputFile.is_open() );

unordered_map<string,string> prevCap3SeqSet, postCap3SeqSet;
//unordered_map<string, string> prevPartSeqSet, postPartSeqSet;

set<string> queryTagSet, multiAlignTagSet;
string line= "";
getline(gmapFile, line);
while( !gmapFile.eof() ) {
      if(line != "") {
//cerr << line << endl;
         vector<string> contentVec;
         boost::split(contentVec, line, boost::is_any_of("\t") );
         string queryTag= contentVec[0];
         if(queryTagSet.find(queryTag) != queryTagSet.end() ) {
            getline(gmapFile, line);
            multiAlignTagSet.insert(queryTag);
            continue;
           }
         queryTagSet.insert(queryTag);
         int flag= boost::lexical_cast<int>(contentVec[1]);
         int position= boost::lexical_cast<int>(contentVec[3])-1;
         string cigar= contentVec[5];
//cout << cigar << endl;
         if(cigar.find("N") != string::npos || cigar.find("S")== string::npos) {
            getline(gmapFile, line);
            continue;
           }
         string seq= contentVec[9];
         vector<string> numFragVec;
         boost::split(numFragVec, cigar, boost::is_any_of("MIDNSH") );
         int alignLength= 0, refAlignLength= 0, mutualAlignLength= 0, prevClip= 0, postClip= 0;
         string::iterator it= cigar.begin();
         for(vector<string>::iterator i= numFragVec.begin(); i!= numFragVec.end(); i++) {
             if(*i== "")
                continue;
             it += (*i).length();
             if(*it== 'M')
                mutualAlignLength += boost::lexical_cast<int>(*i);
             if(*it== 'M' || *it== 'D')
                refAlignLength += boost::lexical_cast<int>(*i);
             if(*it != 'D')
                alignLength += boost::lexical_cast<int>(*i);
             if(*it== 'S') {
                if(i== numFragVec.begin() )
                   prevClip= boost::lexical_cast<int>(*i);
                else
                   postClip= boost::lexical_cast<int>(*i);
               }
//             cout << *it << endl;
             it ++;
            }
         if(mutualAlignLength < ALIGNCUTOFF) {
            getline(gmapFile, line);
            continue;
           }
         if(postClip - (txLength- position- refAlignLength) >30) { 
cout<< cigar<< "\t" << seq.length() << "\t" << postClip << "\t" << alignLength- prevClip << endl;
cout<<">"<< queryTag<<"post\n"<< seq.substr(prevClip)<< endl;    
postPartSeqSet.insert(make_pair(queryTag, seq.substr(prevClip) ) );
            postCap3SeqSet.insert(make_pair(queryTag, seq.substr(seq.length()- postClip + (txLength - position- refAlignLength) ) ) );
//       postCap3InputFile<< ">" << queryTag << "\n" << seq.substr(seq.length()-postClip+ (txLength-position- refAlignLength) ) << endl;
           }
         if(prevClip - position > 30) {
cout<<">" << queryTag<<"prev\n"<<seq.substr(0, seq.length() -postClip)<< endl;
prevPartSeqSet.insert(make_pair(queryTag, seq.substr(0, seq.length()- postClip)));
            prevCap3SeqSet.insert(make_pair(queryTag, seq.substr(0, prevClip- position) ) );
//       prevCap3InputFile<< ">" << queryTag << "\n" << seq.substr(0, prevClip-position) << endl;
           }
        }
    getline(gmapFile, line);
   }

int prevCount= 0, postCount= 0;
for(unordered_map<string,string>::iterator i= prevCap3SeqSet.begin(); i!= prevCap3SeqSet.end(); i++) {
    if(multiAlignTagSet.find(i->first) == multiAlignTagSet.end() ) {
       prevCount ++;
       prevCap3InputFile << ">" << i->first << "\n" << i->second << endl;
      }
   }
for(unordered_map<string,string>::iterator i= postCap3SeqSet.begin(); i!= postCap3SeqSet.end(); i++) {
    if(multiAlignTagSet.find(i->first) == multiAlignTagSet.end() ) {
       postCount ++;
       postCap3InputFile << ">" << i->first << "\n" << i->second << endl;
      }
   }
/*
if(prevCount >1 || postCount >1) {
   cout << "achieve " << endl;
   exit(1);
  }
*/

gmapFile.close();  
prevCap3InputFile.close();
postCap3InputFile.close();
/*
if(prevCount == 1) {
   copyFile("prevCap3", "prevCap3.prev.contigs");
  }
*/
if(prevCount >1) {
   string cmd= "cap3 prevCap" + fileNameTag+" > prevLog" + fileNameTag;
   system(cmd.c_str() );
//   getPrevAttach("prevLog", "prevCap3.prev.contigs", prevPartSeqSet);
  }
/*
if(postCount == 1) {
   copyFile("postCap3", "postCap3.post.contigs");
  }
*/
if(postCount >1) {
   string cmd= "cap3 postCap" + fileNameTag+ " > postLog" + fileNameTag;
   system(cmd.c_str() );
//   getPrevAttach("postLog", "postCap3.post.contigs", postPartSeqSet);
//   exit(1);
  }
removeFile( "prevCap" + fileNameTag);
removeFile( "postCap" + fileNameTag);
}

bool compareByContig(const contigUnit &p1, const contigUnit &p2) {

if(p1.lrNum > p2.lrNum)
   return 1;
if(p1.lrNum < p2.lrNum)
   return 0;
if(p1.lrNum == p2.lrNum)
   return p1.contigSeq.length() > p2.contigSeq.length();
}

int selectExt(unordered_map<string, string>& seqSet, string logFileName, unordered_map<string, string>& partSeqSet, string fileNameTag) {

ifstream logFile(logFileName.c_str() );
assert( logFile.is_open() );

vector<contigUnit> contigUnitVec;

string line= "";
getline(logFile, line);
while( !logFile.eof() ) {
      if(line.find("DETAILED") != string::npos) 
         break;
      if(line.find("Contig") != string::npos) {
         line= line.substr(line.find("Contig")+ 7);
         string contigTag= "Contig"+ line.substr(0, line.find(" ") );
         set<string> lrTagSet;
         getline(logFile, line);
         string representativeTag= "";
         while(line.find("**")== string::npos && line != "") {
               line= line.substr(line.find_first_not_of(" ") );
               string lrTag= line.substr(0, line.find(" ") );
               lrTag= lrTag.substr(0, lrTag.length()-1);
//cout << "lr " << lrTag << endl;
               if(representativeTag == "")
                  representativeTag= lrTag;
               lrTagSet.insert(lrTag);
               getline(logFile, line);
              }
         string contigSeq= seqSet.find(contigTag)->second;
         int lrNum= lrTagSet.size();
         contigUnit contigUnitTmp={contigTag,contigSeq,representativeTag,lrNum};
//cout << contigTag << " " << representativeTag << " " << lrNum << endl;
         contigUnitVec.push_back(contigUnitTmp);
         continue;
        }
      getline(logFile, line);
     }

logFile.close();

if(contigUnitVec.size() >1)
   sort(contigUnitVec.begin(), contigUnitVec.end(), compareByContig);
//cout << "finish sorting " << endl;
/*
vector<pair<string,string> >v(seqSet.begin(), seqSet.end() );
sort(v.begin(), v.end(), compareByLength);
*/

unordered_map<string, string> correctedPartSeqSet;

int seqCount= 0;
//for(vector<pair<string, string> >::iterator i= v.begin(); i!= v.end(); i++) {
for(vector<contigUnit>::iterator i= contigUnitVec.begin(); i!= contigUnitVec.end(); i++) {
    if(seqCount >= EXTCUTOFF)
       break;
    seqCount ++;
//cout << i->representativeTag << endl;
    string lrPartSeq= partSeqSet.find(i->representativeTag)->second;
    transform(lrPartSeq.begin(), lrPartSeq.end(), lrPartSeq.begin(),::toupper);
    
    ofstream clustaloInputFile("clustaloSeq_" + fileNameTag);
    assert(clustaloInputFile.is_open() );
    clustaloInputFile<< ">" << i->contigTag << "\n" << i->contigSeq<< endl;
    clustaloInputFile<< ">" << i->representativeTag << "\n" << lrPartSeq<< endl;
    clustaloInputFile.close();
    string cmd= "clustalo-1.2.0-Ubuntu-x86_64  -i clustaloSeq_" + fileNameTag + " -o clustaloOutput_" + fileNameTag +" --force --threads=24";
    runCmd(cmd);
    string correctedSeq= readClustaloFile(i->contigTag, i->representativeTag, fileNameTag);
cout << ">" << i->representativeTag << "\n" <<  correctedSeq << endl;    
    correctedPartSeqSet.insert(make_pair(i->representativeTag, correctedSeq) );

    removeFile("clustaloSeq_" + fileNameTag);
    removeFile("clustaloOutput_" + fileNameTag);
   }



unordered_map<string,string>().swap(partSeqSet);
partSeqSet= correctedPartSeqSet;
}    

int outputPercent(string seq) { 

for(string::iterator i= seq.begin(); i!= seq.end()- kmer+1; i++) {
    string kmerSeq= string(i, i+kmer);
    BITSET kmerBitset= convertBitSet(kmerSeq);
    int percent= 0;
    unordered_map<BITSET, int>::iterator it;
    if( (it= kmerPercentSet.find(kmerBitset) )!= kmerPercentSet.end() )
        percent= it->second;
     cout << percent << "\t";
   }
cout << endl;
}

int printAllPathsUtil(BITSET srcKmerBitset, set<BITSET> destKmerBitsetSet, unordered_map<BITSET, bool>& visitedSet, BITSET path[], int &path_index, vector<int> pathFreq, int oriLength) {

//cerr<< "to " << convertBit2seq(srcKmerBitset) << endl;
if(pathFreq.size() > oriLength + 50)
   return 0;

BITSET storeKmerBitset= findStoreKmerBitset(srcKmerBitset);
unordered_map<BITSET, bool>::iterator it;
if( (it= visitedSet.find(storeKmerBitset) )== visitedSet.end() ) { 
   return 0;
  }

else if(pathFreq.size() >=10) {
   cout << "freq \t";
   int totalFreq= 0;
   for(vector<int>::iterator i= pathFreq.begin(); i!= pathFreq.end(); i++) {
       cout << *i << "\t";
       totalFreq += *i;
      }
   cout << endl;
   if( (double)totalFreq/(double)pathFreq.size() <20)
      return 0;
  }

it->second= 1;

path[path_index]= srcKmerBitset;
path_index ++;

//if(srcKmerBitset== destKmerBitset) {
if( destKmerBitsetSet.find(srcKmerBitset)!= destKmerBitsetSet.end() ) {
cout<< " PATH \t";
   for(int i=0; i< path_index; i++) {
       cout << convertBit2seq(path[i]) << " ";
      }
   cout << endl;
  }
else{
  BITSET prefix= srcKmerBitset;
//  prefix.resize(prefix.size()+2);
  prefix= prefix << 2;
  BITSET head;
  for(int j= 0; j< 4; j++) {
      createNt(j, head);
      head.resize(prefix.size() );
      BITSET newKmer= prefix | head;
//cerr<< "look for " << convertBit2seq(newKmer) << endl;
//      BITSET storeKmerBitset= findStoreKmerBitset(newKmer);
      unordered_map<BITSET, int>::iterator it;
      if( (it= kmerPercentSet.find(newKmer) )!= kmerPercentSet.end() ) { 
         if(visitedSet.find(storeKmerBitset)->second== 0) {
            int freq= it->second;
            pathFreq.push_back(freq);
            printAllPathsUtil(newKmer, destKmerBitsetSet, visitedSet, path, path_index, pathFreq, oriLength);
           }
        }
     }
  }
path_index --;
it->second= 0;
pathFreq.pop_back();
}
/*
int keepSmallSize(vector<string> &seqVec, string seqRemained) {
cerr << "exceed 300" << endl;

vector<pair<string, int> >seqScoreVec;
for(vector<string>::iterator i= seqVec.begin(); i!= seqVec.end(); i++) {
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    string query= *i;
    bool alignSucceed= 0;
    int seqRemained_length= seqRemained.size();
    alignSucceed= aligner.Align(query.c_str(), seqRemained.c_str(), seqRemained_length, filter, &alignment);
    cerr<< " align " << alignSucceed <<" " << alignment.sw_score << endl;
    seqScoreVec.push_back(make_pair(query, alignment.sw_score) );
    alignment.Clear();
   }
sort(seqScoreVec.begin(), seqScoreVec.end(), 
     boost::bind(&std::pair<string, int>::second, _1) >
     boost::bind(&std::pair<string, int>::second, _2));

vector<string> seqVec_small;
int seqCount= MAXSEQNUMFOREXTEND;
for(vector<pair<string, int> >::iterator i= seqScoreVec.begin(); i!= seqScoreVec.end(); i++) {
    seqCount --;
    seqVec_small.push_back(i->first);
    if(seqCount <=0)
       break;
   }
vector<string>().swap(seqVec);
seqVec= seqVec_small;
cerr<< "finish selection" << endl;
}

int dfs(vector<string>& seqVec, set<BITSET>& destKmerBitsetSet, set<BITSET> visitedKmerSet, int& oriLength, string& seqRemained) {

string sampleSeq= *(seqVec.begin() );
int seqAvgLength= sampleSeq.length();
cerr<< "begin dfs " << endl;
if( seqAvgLength - oriLength >50) {
   cerr<< "too long " << seqAvgLength << " " << oriLength << endl;
   return 0;
  }
bool detectDest= 0;
vector<string> seqVec_extend;
for(vector<string>::iterator i= seqVec.begin(); i!= seqVec.end(); i++) {
    string oriSeq= *i;
cout << oriSeq << endl;    
cerr << "ori:" << oriSeq << endl;    
    string srcKmerSeq= (*i).substr( (*i).length()- kmer); 
    BITSET srcKmerBitset= convertBitSet(srcKmerSeq);
    BITSET prefix= srcKmerBitset;
    prefix = prefix<<2;
    BITSET head;
    for(int j= 0; j< 4; j++) {
        createNt(j, head);
        head.resize(prefix.size() );
        BITSET newKmer= prefix | head;
        BITSET storeKmerBitset= findStoreKmerBitset(newKmer);

        if(destKmerBitsetSet.find(newKmer) != destKmerBitsetSet.end() ) {
           cout << "found it " << *i << j << endl;
           detectDest= 1;
           break;
          }
        unordered_map<BITSET, int>::iterator it;
        if( (it=kmerFreqSet.find(storeKmerBitset) ) == kmerFreqSet.end() )
           continue;
        else if(it->second==0){
                continue;
           }

        if(visitedKmerSet.find(newKmer) == visitedKmerSet.end() ) {
           visitedKmerSet.insert(newKmer);
           string seq_extend= oriSeq;
           switch(j) {
                  case 0: seq_extend += "A";
                          break;
                  case 1: seq_extend += "C";
                          break;
                  case 2: seq_extend += "G";
                          break;
                  case 3: seq_extend += "T";
                          break;
                  default: break;
                 }
           seqVec_extend.push_back(seq_extend);
          }
       }
    if(detectDest== 1)
       break;
   }
vector<string>().swap(seqVec);
seqVec= seqVec_extend;
vector<string>().swap(seqVec_extend);

if(seqVec.size() > MAXSEQNUMFOREXTEND) {
   keepSmallSize(seqVec, seqRemained);
  }
if(detectDest== 0 && seqVec.size() >0) {
   dfs(seqVec, destKmerBitsetSet, visitedKmerSet, oriLength, seqRemained);
  }
}
*/

/*
int attachPrev(string oasesTx, int rightMostPosition, string attachedSeq) {

string srcKmerSeq= attachedSeq.substr(rightMostPosition - kmer+1, kmer);
string seqRemained= attachedSeq.substr(rightMostPosition -kmer+1);
string destKmerSeq= oasesTx.substr(0, kmer);
set<BITSET> destKmerBitsetSet;
for(int i= 0; i< 15; i++) {
    string destKmerSeq= oasesTx.substr(i, kmer);
    BITSET destKmerBitset= convertBitSet(destKmerSeq);
    destKmerBitsetSet.insert(destKmerBitset);
   }

BITSET srcKmerBitset= convertBitSet(srcKmerSeq);
BITSET destKmerBitset= convertBitSet(destKmerSeq);
cerr<< srcKmerSeq << "\t" << srcKmerBitset << endl;
cerr<< destKmerSeq << "\t" << destKmerBitset << endl;

//BITSET *path= new BITSET[kmerPercentSet.size() ];
//vector<int> pathFreq;
//int path_index= 0;
int oriLength= attachedSeq.length() - rightMostPosition;
set<BITSET> visitedKmerSet;
vector<string> seqVec;
seqVec.push_back(srcKmerSeq);
cerr<<"remain: " << seqRemained << endl;

double vm, rss;
process_mem_usage(vm,rss);
cerr<< "VM: " << vm <<";RSS: " << rss << endl;
dfs(seqVec, destKmerBitsetSet, visitedKmerSet, oriLength, seqRemained);

}    
*/

/*
int errorCorrect(unordered_map<string, string> attachedSeqSet) {

string cmd= " /Shared/Au/jason/Source/bowtie2-2.2.1/bin/bowtie2-build attachedPart attachedSeqBowtie2Index >bowtie2Log";
runCmd(cmd);

cmd= "/Shared/Au/jason/Source/bowtie2-2.2.1/bin/bowtie2 -x attachedSeqBowtie2Index -f -U ~/data/oxford/SR_uniq.fa -S attached.sam -p "+ boost::lexical_cast<string>(mpNum)+ " --no-unal --no-hd --no-sq 2>> bowtie2Log";
runCmd(cmd);

readSeqFile("attachedPart", attachedSeqSet);

unordered_map<string, vector<int> >readCoveredSet;
for(unordered_map<string, string>::iterator i= attachedSeqSet.begin(); i!= attachedSeqSet.end(); i++) {
    vector<int> seqVec(i->second.length(), 0);
    readCoveredSet.insert(make_pair(i->first, seqVec) );
   }
cerr<< "finish allocate" << endl;


ifstream samFile("attached.sam");
assert(samFile.is_open() );
string line= "";
getline(samFile, line);
while( !samFile.eof() ) {
      if(line != "") {
         vector<string> contentVec;
         boost::split(contentVec, line, boost::is_any_of("\t") );
         string srTag= contentVec[0];
         int freq= boost::lexical_cast<int>(srTag.substr(srTag.find("_")+1) );
         string cigar= contentVec[5];
         if(cigar.find("D") != string::npos || cigar.find("S") != string::npos || cigar.find("I") != string::npos || cigar.find("H") !=string::npos || cigar.find("N") != string::npos) {
             getline(samFile, line);
             continue;
            }
        string match= "MD:Z:"+ cigar.substr(0, cigar.length()-1);
//cout << match << endl;        
        if(contentVec[17].find(match)== string::npos) {
            getline(samFile, line);
            continue;
           }
        int position= boost::lexical_cast<int>(contentVec[3]) -1;
        int coveredLength= boost::lexical_cast<int>(cigar.substr(0, cigar.length()-1) );
        string refTag= contentVec[2];
        vector<int>* p= &(readCoveredSet.find(refTag)->second);
        for(vector<int>::iterator i= (*p).begin()+ position; i!= (*p).begin()+ position+ coveredLength; i++) {
            (*i)+= freq;
           }
        }
      getline(samFile, line);
     }

samFile.close();
cerr<< "finish reading sam File" << endl;

for(unordered_map<string, vector<int> >::iterator i= readCoveredSet.begin(); i!= readCoveredSet.end(); i++) {
//    cout<< i->first << endl;
    string oasesTxTag= (i->first).substr(0,i->first.find("prev") );
    oasesTxTag= oasesTxTag.substr(0,oasesTxTag.find("post") );
cerr<<  ">" << oasesTxTag << endl;
cout<<  ">" << oasesTxTag << endl;
    string oasesTx= oasesSeqSet.find(oasesTxTag)->second;
cerr<< oasesTx << endl;    
    int position= 0, leftMostPosition= INT_MAX, rightMostPosition= INT_MIN;   
    for(vector<int>::iterator j= i->second.begin(); j!= i->second.end(); j++) {
        if(*j >= SOLIDCUTOFF) {
           leftMostPosition= (leftMostPosition< position? leftMostPosition: position); 
           rightMostPosition= (rightMostPosition> position? rightMostPosition: position); 
          }
        position ++;
       }
cerr<<"region: " << leftMostPosition << "\t" << rightMostPosition << endl;
    string attachedSeq= attachedSeqSet.find(i->first)->second;
    if(i->first.find("prev")!= string::npos && rightMostPosition != INT_MIN) {  
       attachPrev(oasesTx, rightMostPosition, attachedSeq);
      }

    for(vector<int>::iterator j= i->second.begin(); j!= i->second.end(); j++) {
        cout <<*j << "\t";
       }
    cout << endl;
   }
exit(1);
}
*/
int runClustalo(unordered_map<string,string>& rawSeqSet, unordered_map<string,string> partSeqSet, string fileNameTag) {

unordered_map<string,string> extendedSeqSet;
for(unordered_map<string,string>::iterator i= partSeqSet.begin(); i!= partSeqSet.end(); i++) {
    for(unordered_map<string,string>::iterator j= rawSeqSet.begin(); j!= rawSeqSet.end(); j++) {
       ofstream clustaloInputFile("clustaloSeq_" + fileNameTag );
       assert(clustaloInputFile.is_open() );
       clustaloInputFile<< ">" << j->first << "\n" << j->second << endl;
       clustaloInputFile<< ">" << i->first << "\n" << i->second << endl;
       clustaloInputFile.close();

       string cmd= "clustalo-1.2.0-Ubuntu-x86_64  -i clustaloSeq_" + fileNameTag + " -o clustaloOutput_" + fileNameTag+" --force --threads=24";
       runCmd(cmd);
       string mergedSeq= readClustaloFile(j->first, i->first, fileNameTag);
       extendedSeqSet.insert(make_pair(j->first + i->first, mergedSeq) );
       cout << ">merged" << "\n" <<  mergedSeq << endl;

       removeFile("clustaloOutput_"+ fileNameTag);
       removeFile("clustaloSeq_"+ fileNameTag);
      }
//exit(1);
   }
unordered_map<string,string>().swap(rawSeqSet);
rawSeqSet= extendedSeqSet;
}

int extendSeq(string rawSeq, string rawSeqTag, unordered_map<string, string>prevPartSeqSet, unordered_map<string, string> postPartSeqSet, string fileNameTag) {

unordered_map<string, string> prevContigSeqSet, postContigSeqSet, rawSeqSet;
rawSeqSet.insert(make_pair(rawSeqTag, rawSeq) );

if(prevPartSeqSet.size() == 1) { 
   runClustalo(rawSeqSet, prevPartSeqSet, "1");
  }
if(postPartSeqSet.size() == 1) { 
   runClustalo(rawSeqSet, postPartSeqSet, "1");
  }

if(exist_test("prevCap3" +fileNameTag+ ".cap.contigs") ==1 && exist_test("prevLog" + fileNameTag)==1) {
   readSeqFile("prevCap3" +fileNameTag+ ".cap.contigs", prevContigSeqSet);
   removeFile("prevCap3" +fileNameTag+ ".cap.contigs");
   selectExt(prevContigSeqSet, "log_cap" + fileNameTag, prevPartSeqSet, fileNameTag);
//   removeFile("prevLog" + fileNameTag);
   runClustalo(rawSeqSet, prevPartSeqSet, "1");
   }
if(exist_test("postCap3" +fileNameTag+ ".cap.contigs") ==1 && exist_test("postLog" + fileNameTag)==1) {
   readSeqFile("postCap3" +fileNameTag+ ".cap.contigs", postContigSeqSet);
   removeFile("postCap3" +fileNameTag+ ".cap.contigs");
   selectExt(postContigSeqSet, "log_cap" + fileNameTag, postPartSeqSet, fileNameTag);
//   removeFile("postLog" + fileNameTag);
   runClustalo(rawSeqSet, postPartSeqSet, "1");
   }

system( ("rm prevCap3" +fileNameTag+ ".cap.*").c_str() );
system( ("rm postCap3" +fileNameTag+ ".cap.*").c_str() );
//removeFile("prevLog"+ fileNameTag);
//removeFile("postLog"+ fileNameTag);

ofstream outputSeqFile("longerSeq" , ios::out|ios::app);
for(unordered_map<string,string>::iterator i= rawSeqSet.begin(); i!= rawSeqSet.end(); i++) {
    outputSeqFile << ">" << i->first << "\n" << i->second << endl;
   }
outputSeqFile.close();

return 1;




/*
for(unordered_map<string,string>::iterator i= prevPartSeqSet.begin(); i!= prevPartSeqSet.end(); i++) {
    cout << ">" << i->first << "\n" << i->second << endl;
   }
for(unordered_map<string,string>::iterator i= postPartSeqSet.begin(); i!= postPartSeqSet.end(); i++)
    cout << ">" << i->first << "\n" << i->second << endl;
*/


exit(1);







/*
if(prevContigSeqSet.size() > EXTCUTOFF) {
   selectExt(prevContigSeqSet, "prevLog");
//   exit(1);
   }
if(postContigSeqSet.size() > EXTCUTOFF) 
   selectExt(postContigSeqSet, "postLog");

//if(prevSeqSet.size() >0) {
//   errorCorrect(prevSeqSet);
  }
*/

/*
ofstream attachedSeqFile("attachedPart", ios::out|ios::app);
assert( attachedSeqFile.is_open() );
int attachedCount= 0;
for(unordered_map<string,string>::iterator i= prevSeqSet.begin(); i!= prevSeqSet.end(); i++) {
    attachedSeqFile << ">" << rawSeqTag << "prev-" << attachedCount << "\n" << i->second << endl;
    attachedCount ++;
   }
for(unordered_map<string,string>::iterator i= postSeqSet.begin(); i!= postSeqSet.end(); i++) {
    attachedSeqFile << ">" << rawSeqTag << "post-" << attachedCount << "\n" << i->second << endl;
    attachedCount ++;
    }

attachedSeqFile.close();
*/

return 1;











exit(1);



/*
int seqCount= 0;
if(prevSeqSet.size() >0) {
   for(unordered_map<string, string>::iterator i= prevSeqSet.begin(); i!= prevSeqSet.end(); i++) {
       if(postSeqSet.size() >0) {
          for(unordered_map<string,string>::iterator j= postSeqSet.begin(); j!= postSeqSet.end(); j++) {
              cout<< ">" << rawSeqTag << "ext-" <<  seqCount << "\n" << i->second << rawSeq << j->second << endl;
              seqCount ++;
              outputPercent(i->second);
              outputPercent(rawSeq);
              outputPercent(j->second);
             }
         }
       else{
           cout<< ">" << rawSeqTag<< "prev-" <<seqCount << "\n" << i->second << rawSeq<< endl;
           seqCount ++;
           outputPercent(i->second);
           outputPercent(rawSeq);
         }
      }
  }
else if(postSeqSet.size() >0) {
        for(unordered_map<string,string>::iterator j= postSeqSet.begin(); j!= postSeqSet.end(); j++) {
            cout << ">" << rawSeqTag<< "post-" << seqCount << "\n" << rawSeq << j->second << endl;
           seqCount ++;
           outputPercent(rawSeq);
           outputPercent(j->second);
           }
       }
*/

}

int runGmap4Pair( ) {

omp_lock_t writeLock;
omp_init_lock(&writeLock);
#pragma omp parallel
{
#pragma omp for
for(int j= 0; j< blastPairVec.size(); j++) {
//for(unordered_map<string, set<string> >::iterator i= blastPairSet.begin(); i != blastPairSet.end(); i++) {

    vector<pair<string, set<string> > >::iterator i= blastPairVec.begin()+ j;
    string j_str= boost::lexical_cast<string>(j);
   
    string txFileName= "oasesTx"+ j_str, lrSeqFileName= "lrSeq"+ j_str;
    ofstream txFile(txFileName.c_str() ), lrSeqFile(lrSeqFileName.c_str() );
    assert(txFile.is_open()  && lrSeqFile.is_open() );
//cerr<< i->first << endl;
    string txSeq= oasesSeqSet.find(i->first)->second;
    int txLength= txSeq.length();
    txFile << ">" << i->first << "\n" << txSeq << endl;
    txFile.close();
//    cout << i->first << endl;
    
    for(set<string>::iterator j= i->second.begin(); j!= i->second.end(); j++) {
        string lrSeq= lrSeqSet.find(*j)->second; 
        lrSeqFile<< ">" << *j << "\n" << lrSeq << endl;
//        cout << *j << "\t";
       }
    lrSeqFile.close();
    string cmd= "gmap_build  -D locusDir -d " + txFileName + " "+ txFileName + " >> gmapLog 2>>gmapLog2";
    system(cmd.c_str() );
    cmd= "gmap -D locusDir -t 1 -B 5 -A -f samse --no-sam-header --nofails  -d " + txFileName+ " " +lrSeqFileName+ " > gmapOutput" + j_str;
    system(cmd.c_str() );
    system( ("rm locusDir/oasesTx"+ j_str+ ".*").c_str() );
    system( ("rm  -r locusDir/oasesTx"+ j_str).c_str() );

    unordered_map<string,string> prevPartSeqSet, postPartSeqSet;
    readGmapFile(txLength, prevPartSeqSet, postPartSeqSet, j_str);
omp_set_lock(&writeLock);
//    extendSeq(txSeq, i->first, prevPartSeqSet, postPartSeqSet, j_str);
omp_unset_lock(&writeLock);

     removeFile("gmapOutput" + j_str);
     removeFile("oasesTx"+ j_str);
     removeFile("lrSeq"+ j_str);
//    exit(1);
//    cout << endl;
   }
}
omp_destroy_lock(&writeLock);
}    

int formCap3InputFile(string inputFileName, unordered_map<string, string> cap3Set) {
ofstream cap3InputFile(inputFileName.c_str() );
assert( cap3InputFile.is_open() );

for(unordered_map<string, string>::iterator i= cap3Set.begin(); i!= cap3Set.end(); i++) {
    cap3InputFile<< ">" << i->first << "\n" << i->second << endl;
   }
cap3InputFile.close();

string cmd= "cap3 " + inputFileName + " > log_" + inputFileName;
runCmd(cmd);
removeFile(inputFileName);
}

int processConsensus(unordered_map<string,string> prevCap3Set, unordered_map<string, string> prevAddedSet, unordered_map<string, string> postCap3Set, unordered_map<string, string> postAddedSet, string rawSeqTag, string rawSeq, string fileNameTag, unordered_map<string, string>& rawSeqSet) {
/*
ofstream cap3InputFile( ("prevCap" + fileNameTag).c_str() );
assert(cap3InputFile.is_open() );
for(unordered_map<int,string>::iterator i= prevCap3Set.begin(); i!= prevCap3Set.end(); i++) {
    cap3InputFile<< ">" << i->first << "\n" << i->second << endl;
    }
cap3InputFile.close();

string cmd= "~/app/cap3/CAP3/cap3 prevCap" + fileNameTag +" >prevLog" + fileNameTag;
runCmd(cmd);
removeFile("prevCap3"+ fileNameTag);
*/

unordered_map<string, string> prevContigSeqSet, postContigSeqSet;
rawSeqSet.insert(make_pair(rawSeqTag, rawSeq) );

if(prevAddedSet.size() ==1) {
   runClustalo(rawSeqSet, prevAddedSet, fileNameTag);
  }
if(postAddedSet.size() ==1) {
   runClustalo(rawSeqSet, postAddedSet, fileNameTag);
  }

if(prevCap3Set.size() >1)
   formCap3InputFile("prevCap" + fileNameTag, prevCap3Set);
if(postCap3Set.size() >1)
   formCap3InputFile("postCap" + fileNameTag, postCap3Set);

if(exist_test("prevCap" +fileNameTag+ ".cap.contigs") ==1 && exist_test("log_prevCap" + fileNameTag)==1) {
   readSeqFile("prevCap" +fileNameTag+ ".cap.contigs", prevContigSeqSet);
   removeFile("prevCap" +fileNameTag+ ".cap.contigs");
   selectExt(prevContigSeqSet, "log_prevCap" + fileNameTag, prevAddedSet, fileNameTag);
   removeFile("log_prevCap" + fileNameTag);
   runClustalo(rawSeqSet, prevAddedSet, fileNameTag);
   }

if(exist_test("postCap" +fileNameTag+ ".cap.contigs") ==1 && exist_test("log_postCap" + fileNameTag)==1) {
   readSeqFile("postCap" +fileNameTag+ ".cap.contigs", postContigSeqSet);
   removeFile("postCap" +fileNameTag+ ".cap.contigs");
   selectExt(postContigSeqSet, "log_postCap" + fileNameTag, postAddedSet, fileNameTag);
   removeFile("log_postCap" + fileNameTag);
   runClustalo(rawSeqSet, postAddedSet, fileNameTag);
  }

system( ("rm prevCap" +fileNameTag+ ".cap.*").c_str() );
system( ("rm postCap" +fileNameTag+ ".cap.*").c_str() );
//removeFile("prevLog"+ fileNameTag);
//removeFile("postLog"+ fileNameTag);
/*
ofstream outputSeqFile("longerSeq" , ios::out|ios::app);
for(unordered_map<string,string>::iterator i= rawSeqSet.begin(); i!= rawSeqSet.end(); i++) {
    outputSeqFile << ">" << i->first << "\n" << i->second << endl;
    }
outputSeqFile.close();
*/
}

int writeResultIntoFile(unordered_map<string, string> rawSeqSet) {


ofstream outputSeqFile("longerSeq" , ios::out|ios::app);
for(unordered_map<string,string>::iterator i= rawSeqSet.begin(); i!= rawSeqSet.end(); i++) {
if(i->first.find("lr") == string::npos) {
   continue;
  }
    outputSeqFile << ">" << i->first << "\n" << i->second << endl;
   }
outputSeqFile.close();
}

int getConsensus( ) {

vector<pair<string, addedPartUnit> >contigAddedPartVec(contigAddedPartSet.begin(), contigAddedPartSet.end() );

omp_lock_t writelock;
omp_init_lock(&writelock);

#pragma omp parallel
{
#pragma omp for
for(int j= 0; j< contigAddedPartVec.size(); j++) {
    string fileNameTag= boost::lexical_cast<string>(j);
   
    unordered_map<string, string> rawSeqSet;
    vector<pair<string, addedPartUnit> >::iterator i= contigAddedPartVec.begin()+j;
    processConsensus(i->second.prevCap3Set, i->second.prevAddedSet, i->second.postCap3Set, i->second.postAddedSet, i->first, i->second.rawSeq, fileNameTag, rawSeqSet);
omp_set_lock(&writelock);
    writeResultIntoFile(rawSeqSet);
omp_unset_lock(&writelock);
   }
}
omp_destroy_lock(&writelock);
} 

int main(int argc, char* argv[]) {

clock_t startT, endT;
startT= clock();

cerr<< "./scaffoldExtend <LR_gmap2_SrTx.sam> <Oases_output:Transcripts.fa/soap Output> <lrseq.fa> nThreads" << endl;
assert(argc== 5);
string gmapSamFileName= argv[1], oasesFileName= argv[2], lrSeqFileName= argv[3];
mpNum= atoi(argv[4]);

ofstream outputSeqFile("longerSeq");
assert(outputSeqFile.is_open() );
outputSeqFile.close();

mpNum= 24;
cerr << "mp available: " << omp_get_num_procs() << endl;
mpNum= (mpNum > omp_get_num_procs()? omp_get_num_procs(): mpNum);
omp_set_num_threads(mpNum);

#pragma omp parallel
{
if(omp_get_thread_num()== 0)
   cerr << "use threads num: " << omp_get_num_threads() << endl;
}

//readJellyFile(jellyFileName);
cerr << "finish scanning with size " << kmerFreqSet.size()  << endl;

#pragma omp parallel
{
#pragma omp for
for(int i=0 ; i< 4; i++) {
    if(i== 0) {
//       getBlastPair(blastOutputFileName);
//       getGmapPair(gmapSamFileName);      
       cerr<< "finish reading blast Output file" << endl;
      }
/*
    if(i== 1) {
       if( !exist_test("kmerPercent") )
          reportPercent();
       else
          readPercent();
      }
*/
    if(i== 2) {
       readSeqFile(oasesFileName, oasesSeqSet);
       cerr<< "finish Oases file" << endl;
      }
    if(i== 3) {
       readSeqFile(lrSeqFileName, lrSeqSet);
       cerr<< "finish LR seq file" << endl;
      }
   }
}
getGmapPair(gmapSamFileName);      
getConsensus();

endT= clock();
float diff ((float)endT- (float)startT) ;
float seconds= diff/CLOCKS_PER_SEC;
cout << "running time(s) :" << seconds << endl;

/*
getBlastPair(blastOutputFileName);
cerr<< "finish reading blast Output file" << endl;
readSeqFile(oasesFileName, oasesSeqSet);
cerr<< "finish Oases file" << endl;
readSeqFile(lrSeqFileName, lrSeqSet);
cerr<< "finish LR seq file" << endl;
*/
//runGmap4Pair( );

unordered_map<string, string> seqSetTmp;
//errorCorrect(seqSetTmp);



}
