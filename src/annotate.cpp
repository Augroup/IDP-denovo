//updated on 6/4/2016, add reverse complement from cluster file to get seq for msa
//updated on 6/6/2016, annot from extension/lr/scaffolds set
//updated on 6/15/2016, with extended scaffold as 1st seq in MSA is possible
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
#include <omp.h>
#include <boost/multi_array.hpp>

#include <sys/time.h>
#include <time.h>

#define CLUSTERSIZECUTOFF 30
#define PERCENTCUTOFF 30
#define MAXEXONLENGTH 480
#define MINEXONLENGTH 43

using namespace std;
using namespace boost;

struct ntUnit{
       int Acount;
       int Tcount;
       int Gcount;
       int Ccount;
       double confidence;
       int ntCount;
       char nt;
      };

int mpNum= 32;
unordered_map<string, string> seqSet;
unordered_set<string> flTagSet;
unordered_map<string, bool> strandSet;
vector<vector<string> > largeClusterVec;

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
                     break;
           }
     }
return complementSeq;
}
/*
int readAnnotFile(string annotFileName) {

ifstream annotFile(annotFileName.c_str() );
assert(annotFile.is_open() );

string line= "";
getline(annotFile, line);
while( !annotFile.eof() ) {
      if(line != "" && line[0] != '#') {
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         if(contentVec[4]== "1")
            flTagSet.insert(contentVec[0]);
           }
         getline(annotFile, line);
       }
annotFile.close();
cout << "finish reading annot file" << endl;
}
*/
int runCmd(string cmd) {
cerr<< cmd << endl;
system(cmd.c_str() );
}
   
int removeFile(string fileName) {
if(remove(fileName.c_str() ) !=0)
   cerr << "error deleting file " << fileName << endl;
else
   cerr << "successful deleting file " << fileName << endl;
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
         tag= line.substr(1);
         tag= tag.substr(0, tag.find(" ") );
//         tag= line.substr(1, line.find("_")-1);
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

int readClusterFile(string clusterFileName, vector<vector<string> >& tagVec_total) {

ifstream clusterFile(clusterFileName.c_str() );
assert( clusterFile.is_open() );
string line= "";
getline(clusterFile, line);
while( !clusterFile.eof() ) {
      if(line.find("head:") != string::npos) {
         vector<string> tagVec;
         string tag= line.substr(0, line.find(" ") );
         string lrTag= tag.substr(0, tag.find("_") );
         if(lrTag[0] != 'l') {
            lrTag= tag;
           }
         strandSet.insert(make_pair(tag, 0) );
//         if(flTagSet.find(lrTag) != flTagSet.end() )
            tagVec.push_back(tag);
         getline(clusterFile, line);
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         for(vector<string>::iterator i= contentVec.begin(); i!= contentVec.end(); i++) {
             if(*i== "")
                continue;
             string tag= (*i).substr(0, (*i).find(":") );
             string strand= (*i).substr(0, (*i).find_last_of(":") );
             bool direction= 0;
             if(strand== "1") 
                direction= 1;
             strandSet.insert(make_pair(tag, direction) );
             string lrTag= tag.substr(0, tag.find("_") );
             if(lrTag[0] != 'l')
                lrTag= tag;
//             if(flTagSet.find(lrTag) != flTagSet.end() )
                tagVec.push_back(tag);
             }
            if(tagVec.size() <= CLUSTERSIZECUTOFF)
               tagVec_total.push_back(tagVec);
            else
               largeClusterVec.push_back(tagVec);
         }
      getline(clusterFile, line);
      }
clusterFile.close();
}
   
string getSeq(string tag) {
//tag= tag.substr(0, tag.find("_") );
string seq= "";
unordered_map<string, string>::iterator it;
if( (it= seqSet.find(tag) )== seqSet.end() ) {
   cerr<< "cannot find seq " << tag << endl;
   exit(1);
   }
unordered_map<string, bool>::iterator jt;
if( (jt= strandSet.find(tag) )== strandSet.end() ) {
   cerr<< "cannot find direction " << tag << endl;
   exit(1);
  }
if(jt->second== 1) {
   reverse(seq.begin(), seq.end() );
   seq= getComplementSeq(seq);
  }

seq= it->second;
return seq;
}

bool compareByLength(const pair<string,string> &p1, const pair<string,string> &p2) {
return p1.second.length() > p2.second.length();
}

int sortByLength(vector<string>& tagVec) {

vector<pair<string, string> >seqVec;
for(vector<string>::iterator i= tagVec.begin(); i!= tagVec.end(); i++) {
    string seq= getSeq(*i);
    seqVec.push_back(make_pair(*i, seq) );
   }
sort(seqVec.begin(), seqVec.end(), compareByLength);
vector<string>().swap(tagVec);
for(vector<pair<string,string> >::iterator i= seqVec.begin(); i!= seqVec.end(); i++) {
    tagVec.push_back(i->first);
   }
}

int fillNt(char nt, ntUnit &ntUnitTmp) {
switch(nt) {
      case 'A':
      case 'a': ntUnitTmp.Acount++;
                break;
      case 'T':
      case 't': ntUnitTmp.Tcount++;
                break;
      case 'C':
      case 'c': ntUnitTmp.Ccount++;
                break;
      case 'G':
      case 'g': ntUnitTmp.Gcount++;
                break;
      default:  break;
     }
}

int getNt(vector<ntUnit>& ntUnitVec) {

for(vector<ntUnit>::iterator i= ntUnitVec.begin(); i!= ntUnitVec.end(); i++) {
    (*i).ntCount= (*i).Acount;
    (*i).nt= 'A';
    if( (*i).ntCount < (*i).Tcount) {
       (*i).nt= 'T';
       (*i).ntCount= (*i).Tcount;
       }
    if( (*i).ntCount < (*i).Ccount) {
       (*i).nt= 'C';
       (*i).ntCount= (*i).Ccount;
       }
    if( (*i).ntCount < (*i).Gcount) {
       (*i).nt= 'G';
       (*i).ntCount= (*i).Gcount;
       }
    (*i).confidence= (double)(*i).ntCount/(double)( (*i).Acount+ (*i).Tcount + (*i).Gcount+ (*i).Ccount);
    }
}

int readSeqFileVec(string seqFileName, vector<pair<string, string> > &seqVec) {

ifstream seqFile(seqFileName.c_str() );
assert( seqFile.is_open() );

string line= "", tag= "", seq= "";
getline( seqFile, line);
while( !seqFile.eof() ) {
      if(line[0]== '>') {
         if(tag != "") {
            seqVec.push_back(make_pair(tag, seq) );
            }
         tag= line.substr(1, line.find(" ")-1);
         seq= "";
         }
      else
         seq+= line;
      getline( seqFile, line);
     }
if(tag != "")
   seqVec.push_back(make_pair(tag, seq) );
seqFile.close();
}


string getAlignConsensus(string clustaloOutputFileName) {

vector<pair<string, string> >seqVec;
readSeqFileVec(clustaloOutputFileName, seqVec);

string representativeSeq= (*(seqVec.begin() ) ).second;

unordered_map<int, int> gapPosSet;
int gapBeginPos= -1, pos= 0;
for(string::iterator i= representativeSeq.begin(); i!= representativeSeq.end(); i++) {
    if(*i== '-' && gapBeginPos== -1) {
       gapBeginPos= pos;
      }
    if(*i != '-' && gapBeginPos != -1) {
       int gapLength= pos- gapBeginPos;
       if(gapLength >= MINEXONLENGTH)
          gapPosSet.insert(make_pair(gapBeginPos, pos) );
       gapBeginPos= -1;
      }
    pos ++;
   }

for(unordered_map<int, int>::iterator i= gapPosSet.begin(); i!= gapPosSet.end(); i++) {
    cout << i->first << "\t" << i->second << endl;
    int regionLength= i->second - i->first; 
    ntUnit ntUnitTmp= {0,0,0,0,0,0,'N'};
    vector<ntUnit> ntUnitVec(regionLength, ntUnitTmp);
    for(vector<pair<string, string> >::iterator j= seqVec.begin()+1; j!= seqVec.end(); j++) {
       vector<ntUnit>::iterator it= ntUnitVec.begin();
       for(string::iterator k= j->second.begin()+ i->first; k!= j->second.begin()+ i->second; k++) {
           fillNt(*k, *it);
           it++;
           }
       }
    getNt(ntUnitVec);
    string gapSeq= "";
    for(vector<ntUnit>::iterator m= ntUnitVec.begin();m!= ntUnitVec.end(); m++){
        gapSeq += (*m).nt;
        cout << (*m).nt;
       }
    cout << endl;
    representativeSeq= representativeSeq.substr(0, i->first) + gapSeq + representativeSeq.substr(i->second); 
   }

representativeSeq.erase(remove(representativeSeq.begin(), representativeSeq.end(), '-'), representativeSeq.end() );

cout << representativeSeq << endl;
return representativeSeq;




string consensusSeq= "";
unordered_map<string, string> alignedSeqSet;
readSeqFile(clustaloOutputFileName, alignedSeqSet);

int regionLength= (*(alignedSeqSet.begin() ) ).second.length();
ntUnit ntUnitTmp= {0,0,0,0,0,0,'N'};
vector<ntUnit> ntUnitVec(regionLength, ntUnitTmp);

for(unordered_map<string, string>::iterator i= alignedSeqSet.begin(); i!= alignedSeqSet.end(); i++) {
    vector<ntUnit>::iterator it= ntUnitVec.begin();
    for(string::iterator j= i->second.begin(); j!= i->second.end(); j++) {
        fillNt(*j, *it);
        it++;
       }
    }
getNt(ntUnitVec);

for(vector<ntUnit>::iterator i= ntUnitVec.begin(); i!= ntUnitVec.end(); i++){
    cout << (*i).nt;
    consensusSeq += (*i).nt;
   }
cout << endl;
return consensusSeq;
}

int readMatrixFile(string fileTag, vector<pair<string,string> >& inputSeqVec, vector<pair<string,string> >& unannotatedSeqVec){

ifstream matrixFile("matrix_" +fileTag);
assert(matrixFile.is_open() );

string line= "";
getline(matrixFile, line);
getline(matrixFile, line);
vector<string> contentVec, tmpVec;
split(tmpVec, line, is_any_of(" ") );
for(vector<string>::iterator i= tmpVec.begin(); i!= tmpVec.end(); i++) {
    if(*i != "")
      contentVec.push_back(*i);
   }
cout<< "%:"<<contentVec[1]<< "\t"<<contentVec[2]<< "\t"<< contentVec[3] << endl;
if(lexical_cast<double>(contentVec[3]) < PERCENTCUTOFF) {
   string tag= inputSeqVec[2].first + ":" + contentVec[3];
   unannotatedSeqVec.push_back(make_pair(tag, inputSeqVec[2].second) );
   inputSeqVec.pop_back();
   cout << "3rd seq bad" << endl;
  }
if(lexical_cast<double>(contentVec[2]) < PERCENTCUTOFF) {
   string tag= inputSeqVec[1].first + ":" + contentVec[2];
   unannotatedSeqVec.push_back(make_pair(tag, inputSeqVec[1].second) );
   inputSeqVec.erase(inputSeqVec.begin()+1);
   cout << "2nd seq bad" << endl;
  }

matrixFile.close();
removeFile("matrix_" + fileTag);
if(inputSeqVec.size() ==1)
   return 1;
//cout << "size:" << inputSeqVec.size() << endl;
if(inputSeqVec.size() != 3) {
   ofstream clustalInputFile( ("clustaloInput_"+ fileTag).c_str() );
   assert( clustalInputFile.is_open() );
   for(vector<pair<string,string> >::iterator i= inputSeqVec.begin(); i!= inputSeqVec.end(); i++) {
       clustalInputFile << ">" << i->first << "\n" << i->second << endl;
      }
   clustalInputFile.close();
   string cmd= "clustalo-1.2.0-Ubuntu-x86_64 --percent-id --force --threads=1 --iter=1 -i clustaloInput_" + fileTag + " -o clustaloOutput_" + fileTag;
    runCmd(cmd);
  }
return 0;
}

int getGapInfo(string templateSeq, unordered_map<string, string>seqSet_subset, string fileTag, unordered_map<string, vector<string> >& gapSet, unordered_map<string, string>& querySeqSet, string &extendedSeq, vector<string> &reportGpdContentVec, vector<string>& extendGpdContentVec) {

string indexFileName= "gmapIndex_"+ fileTag;
ofstream indexFile(indexFileName.c_str() );
assert(indexFile.is_open() );
indexFile<< ">consensus" << "\n" << templateSeq << endl;
indexFile.close();

string cmd= "gmap_build -D indexDir -d " + indexFileName  + " " + indexFileName + " >> gmapLog1 2>>gmapLog2";
runCmd(cmd);

string queryFileName= "queryFile_" + fileTag;
ofstream queryFile(queryFileName.c_str());
assert(queryFile.is_open() );
for(unordered_map<string, string>::iterator i= seqSet_subset.begin(); i!= seqSet_subset.end(); i++) {
    queryFile<< ">" << i->first << "\n" << i->second << endl;
   }
queryFile.close();


cmd= "gmap  -D indexDir -d " + indexFileName + " -t 16 -B 5 -A -f samse --no-sam-header --nofails "+ queryFileName + " >samOutput_" + fileTag + " 2>> gmapLog2";
runCmd(cmd);

cmd= "rm -r indexDir/gmapIndex_"+ fileTag;
runCmd(cmd);
removeFile(queryFileName);

ifstream samFile( ("samOutput_"+ fileTag).c_str() );
assert(samFile.is_open() );

int prevAttachmentLength= 0;
string reportSeq= templateSeq;
extendedSeq= templateSeq;
unordered_set<string> visitedTagSet;
string prevTag= "" ;
string line= "";
getline(samFile, line);
while( !samFile.eof() ) {
      if(line != "") {
         vector<string> contentVec;
         split(contentVec, line, is_any_of("\t") );
         int mapQ= lexical_cast<int>(contentVec[4]);
         if(mapQ== 0) {
            getline(samFile, line);
            continue;
           }
        string queryTag= contentVec[0];
        string strand= contentVec[1];
        if(prevTag == queryTag) {
           getline(samFile, line);
           continue;
           }
        prevTag= queryTag;
        int beginPos= lexical_cast<int>(contentVec[3]) -1;
        string cigar= contentVec[5];
        string querySeq= contentVec[9];
        string oriQuerySeq= querySeq;
        int oriQuerySeqLength= querySeq.length();
        string startList= "", endList= "",reportStartList="", reportEndList="";
        vector<string> lengthVec;
        split(lengthVec, cigar, is_any_of("DMINHS") );
        int pLength= 0, prevSoft= 0, postSoft= 0, pPosition= beginPos;
        int matchBegin= -1;
        bool prevSoftFound= 0, postSoftFound= 0;
        int prevSoftLength= 0, postSoftLength= 0;
        int alignedLength= 0, unalignedLength= 0;
        for(vector<string>::iterator i= lengthVec.begin(); i!= lengthVec.end(); i++) {
            if(*i== "")
               continue;
            int length= lexical_cast<int>(*i);
            pLength +=  (*i).length();
            char op= lexical_cast<char>(cigar.substr(pLength, 1) );
             switch(op) {
                    case 'S': if(pLength== cigar.length()-1) {
                                 oriQuerySeq= oriQuerySeq.substr(0, oriQuerySeq.length()- length);
                                 postSoftLength= length;
                                 postSoftFound= 1;
                                 assert(querySeq.length() >= length);
                                 string exonStr= querySeq.substr(querySeq.size()-length);
                                 endList += lexical_cast<string>(pPosition) +",";
                                 startList += lexical_cast<string>(extendedSeq.length()- prevAttachmentLength ) +",";
                                 extendedSeq += exonStr;
                                 pPosition= extendedSeq.length()- prevAttachmentLength -1;
                                 }
                              else{
                                  prevSoftLength= length;
                                  oriQuerySeq= oriQuerySeq.substr(length);
                                  prevSoftFound= 1;
                                  endList += lexical_cast<string>(prevAttachmentLength*(-1) )+",";
                                  prevAttachmentLength += length;
                                  startList += lexical_cast<string>(prevAttachmentLength*(-1) )+",";
                                  assert(querySeq.length() >= length);
                                  string exonStr= querySeq.substr(0, length);
                                  extendedSeq= exonStr + extendedSeq;
                                  }
                               break;
                      case 'M': if(matchBegin== -1) {
                                   matchBegin= pPosition;
                                   startList += lexical_cast<string>(matchBegin)+ ",";
                                   }
                                 break;
                    case 'N':{endList += lexical_cast<string>(pPosition)+",";
                              matchBegin= -1;
                              int queryPos= oriQuerySeq.length()- querySeq.length();
                              string gap= lexical_cast<string>(queryPos) + ":" + lexical_cast<string>(length) + "-" + lexical_cast<string>(pPosition);
                              unordered_map<string, vector<string> >::iterator it;
                              if(length <= MAXEXONLENGTH && length >= MINEXONLENGTH) {
                                 if( (it= gapSet.find(queryTag))==gapSet.end()){
                                     vector<string> gapVec;
                                     gapVec.push_back(gap);
                                     gapSet.insert(make_pair(queryTag, gapVec));
                                    }
                                 else
                                  it->second.push_back(gap);
                                }
                              }
                             break;
                    default: break;
                   }
             pLength += 1;
             if(op== 'M' || op== 'D' || op== 'N')
                pPosition += length;
             if(op== 'M' || op== 'I' || op== 'S')
                querySeq = querySeq.substr(length);

             if(op== 'M')
                alignedLength += length;
             if(op== 'S')
                unalignedLength += length;
             }
          querySeqSet.insert(make_pair(queryTag, oriQuerySeq) );
          if(unalignedLength > alignedLength){
             gapSet.erase(queryTag);
          cout << "erase:" << queryTag << " " << unalignedLength << " " << alignedLength << endl;           
          }
         endList += lexical_cast<string>(pPosition) +",";
         string end= lexical_cast<string>(pPosition);
         string reportEnd= end;
                                                                                         int nBlocks= count(startList.begin(), startList.end(), ',');
         int reportNBlocks= nBlocks;
         reportStartList= startList;
         reportEndList= endList;
        if(prevSoftFound== 1) {
           reportStartList= reportStartList.substr(reportStartList.find(",")+1);
           reportEndList= reportEndList.substr(reportEndList.find(",")+1);
           reportNBlocks --; 
          } 
        if(postSoftFound== 1) {
           reportStartList= reportStartList.substr(0, reportStartList.size()-1);
           reportStartList= reportStartList.substr(0, reportStartList.find_last_of(",")+1);
           reportEndList= reportEndList.substr(0, reportEndList.size()-1);
           reportEndList= reportEndList.substr(0, reportEndList.find_last_of(",")+1 );
           reportEnd= reportEndList.substr(0, reportEndList.size()-1);
           if(reportEnd.find(",") != string::npos)
              reportEnd= reportEnd.substr(reportEnd.find_last_of(",")+1);
              reportNBlocks --;
             }
          char strand_char= '+';
         if(strand== "16")
            strand_char= '-';
         string start= startList.substr(0, startList.find(",") );
         string reportStart= reportStartList.substr(0, reportStartList.find(",") );

cout << startList << "\n" << endList << "\n" << endl;
cout << reportStartList << "\n" << reportEndList << "\n" << endl;
//exit(1);
string extendGpdContent= "Gene"+ fileTag + "\t" + queryTag + "\tGene" + fileTag+ "\t" + lexical_cast<string>(strand_char) + "\t" + start + "\t" + end + "\t" + start + "\t" + end + "\t" + lexical_cast<string>(nBlocks) + "\t" + startList + "\t" + endList;
extendGpdContentVec.push_back(extendGpdContent);
string reportGpdContent= "Gene"+ fileTag + "\t" + queryTag + "\tGene" + fileTag+ "\t" + lexical_cast<string>(strand_char) + "\t" + reportStart + "\t" + reportEnd + "\t" + reportStart + "\t" + reportEnd + "\t" + lexical_cast<string>(reportNBlocks) + "\t" + reportStartList + "\t" + reportEndList;
reportGpdContentVec.push_back(reportGpdContent);
        }
      getline(samFile, line);
     }
samFile.close();

removeFile("samOutput_"+ fileTag);
removeFile("gmapIndex_"+ fileTag);
}

int findExtendedTemplate(vector<string>& tagVec) {

string templateTag= "";
for(vector<string>::iterator i= tagVec.begin(); i!= tagVec.end(); i++) {
    if( (*i).find("Locus") != string::npos && (*i).find("lr") != string::npos) {
       if( (*i).length() > templateTag.length() )
          templateTag= *i;
      }
   }
if(templateTag != "") {
   tagVec.erase(remove(tagVec.begin(), tagVec.end(), templateTag), tagVec.end() );
   vector<string>::iterator it= tagVec.begin();
   tagVec.insert(it, templateTag);
  }
}

int msa(vector<string> tagVec, string fileTag, string& templateSeq, string& extendedSeq, vector<string>& reportGpdContentVec, vector<string>& extendGpdContentVec, unordered_map<string, vector<string> >& gapSet, unordered_map<string, string>& querySeqSet, vector<pair<string,string> >& unannotatedSeqVec) {

unordered_map<string, string> seqSet_subset;

sortByLength(tagVec);
findExtendedTemplate(tagVec);

string templateTag= *(tagVec.begin() );
templateSeq= getSeq(templateTag);
seqSet_subset.insert(make_pair(templateTag, templateSeq) );
for(vector<string>::iterator i= tagVec.begin()+1; i!= tagVec.end(); i++) {
    string cmd= "";
    vector<pair<string, string> > inputSeqVec;
    inputSeqVec.push_back(make_pair(templateTag, templateSeq) );
//cout << templateTag << endl;
    cmd= "clustalo-1.2.0-Ubuntu-x86_64 --percent-id --force --threads=1 --iter=1 -i clustaloInput_" + fileTag + " -o clustaloOutput_" + fileTag;
    ofstream clustalInputFile( ("clustaloInput_"+ fileTag).c_str() );
    assert( clustalInputFile.is_open() );
    clustalInputFile<< ">" << templateTag << "\n" << templateSeq << endl;
    string lrTag= *i;
    string lrSeq= getSeq(lrTag);
    clustalInputFile<< ">" << lrTag << "\n" << lrSeq << endl;
//cout << lrTag << endl;
    inputSeqVec.push_back(make_pair(lrTag, lrSeq) );
    seqSet_subset.insert(make_pair(lrTag, lrSeq) );
    if(i+1 != tagVec.end() ) {
       string lrTag_next= *(i+1);
       string lrSeq_next= getSeq(lrTag_next);
       clustalInputFile<< ">" << lrTag_next << "\n" << lrSeq_next << endl;
//cout << lrTag_next << endl;
       inputSeqVec.push_back(make_pair(lrTag_next, lrSeq_next) );
       seqSet_subset.insert(make_pair(lrTag_next, lrSeq_next) );
       cmd= "clustalo-1.2.0-Ubuntu-x86_64 --percent-id --force --full --threads=1 --iter=1 --distmat-out=matrix_" + fileTag + " -i clustaloInput_" + fileTag + " -o clustaloOutput_" + fileTag;
      }
    clustalInputFile.close();

    bool errorCluster= 0;
    runCmd(cmd);
    if(i+1 != tagVec.end() ) {
       errorCluster= readMatrixFile(fileTag, inputSeqVec, unannotatedSeqVec);
       i++;
      }
    if(errorCluster== 0)
       templateSeq= getAlignConsensus( "clustaloOutput_"+ fileTag);
    templateTag= "consensus";
//cout << "pos:" << distance(tagVec.begin(), i) << endl;
    removeFile("clustaloInput_"+ fileTag);
    removeFile("clustaloOutput_"+ fileTag);
//exit(1);
    if(i== tagVec.end() )
       break;
   }

getGapInfo(templateSeq, seqSet_subset, fileTag, gapSet, querySeqSet, extendedSeq, reportGpdContentVec, extendGpdContentVec);
}

int main(int argc, char* argv[]) {

cerr<< "./annotate <cluster> <seq.fa> nThreads" << endl;

assert(argc== 4);

string clusterFileName= argv[1];
string seqFileName= argv[2];
mpNum= atoi(argv[3]);

//readAnnotFile(annotFileName);

ofstream extendSeqFile("extendedSeq"), reportSeqFile("report_seq"), extendGpdFile("extend.gpd"), reportGpdFile("report.gpd"), gapFile("gapInfo"), queryFile("querySeq");
assert(extendSeqFile.is_open()  && reportSeqFile.is_open()  && extendGpdFile.is_open()  && reportGpdFile.is_open()  && gapFile.is_open() && queryFile.is_open() );
extendSeqFile.close();
reportSeqFile.close();
extendGpdFile.close();
reportGpdFile.close();
gapFile.close();
queryFile.close();

vector<vector<string> >tagVec_total;

cerr << "mp available: " << omp_get_num_procs() << endl;
mpNum= (mpNum > omp_get_num_procs()? omp_get_num_procs(): mpNum);
omp_set_num_threads(mpNum);

#pragma omp parallel
{ 
if(omp_get_thread_num()== 0)
   cerr << "use threads num: " << omp_get_num_threads() << endl;
}

struct timeval startClock, endClock;
clock_t startT, endT;
timespec time1, time2;

clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time1);
gettimeofday(&startClock, NULL);
startT= clock();

#pragma omp parallel
{
#pragma omp for
for(int i= 0; i< 2; i++) {
    if(i== 0)
       readSeqFile(seqFileName, seqSet);
    if(i== 1)
       readClusterFile(clusterFileName, tagVec_total);
    }
}

ofstream largeClusterFile("largeCluster.fa");
assert(largeClusterFile.is_open() );
for(vector<vector<string> >::iterator i= largeClusterVec.begin(); i!= largeClusterVec.end(); i++) {
   for(vector<string>::iterator j= (*i).begin(); j!= (*i).end(); j++) {  
       string seq= getSeq(*j);
       largeClusterFile<< ">" <<  *j << "\n" << seq << endl;
      }
  }
largeClusterFile.close();

cout << "finish reading input files" << endl;
string cmd= "mkdir indexDir";
runCmd(cmd);

omp_lock_t writelock;
omp_init_lock(&writelock);

ofstream unannotatedFile("unannotated.fa");
assert(unannotatedFile.is_open() );

cerr<< "num of clusters " << tagVec_total.size() << endl;
#pragma omp parallel
{
#pragma omp for
for(int j= 0; j< tagVec_total.size(); j++) {
    vector<vector<string> >::iterator i= tagVec_total.begin() +j;
    string fileTag= lexical_cast<string>(j);
    if( (*i).size() >=2 ) {
       string reportSeq= "", extendedSeq= "";
       vector<string> reportGpdContentVec, extendGpdContentVec;
       unordered_map<string, vector<string> >gapSet;
       unordered_map<string, string> querySeqSet;
       vector<pair<string,string> >unannotatedSeqVec;
       msa(*i, fileTag, reportSeq, extendedSeq, reportGpdContentVec, extendGpdContentVec, gapSet, querySeqSet, unannotatedSeqVec); 

omp_set_lock(&writelock);
       ofstream extendedFile("extendedSeq", ios::out|ios::app), reportSeqFile("report_seq", ios::out|ios::app), reportGpdFile("report.gpd", ios::out|ios::app), extendGpdFile("extend.gpd", ios::out|ios::app), gapFile("gapInfo", ios::out|ios::app), queryFile("querySeq", ios::out|ios::app);
        assert( extendedFile.is_open() && reportSeqFile.is_open() && reportGpdFile.is_open() && extendGpdFile.is_open()  && gapFile.is_open()  && queryFile.is_open() );
        extendedFile<< ">Gene"  << fileTag << "\n" << extendedSeq << endl;
        reportSeqFile<< ">Gene"  << fileTag << "\n" << reportSeq << endl;
        for(vector<string>::iterator m= reportGpdContentVec.begin(); m!= reportGpdContentVec.end(); m++)
            reportGpdFile << *m << endl;
            for(vector<string>::iterator m= extendGpdContentVec.begin(); m!= extendGpdContentVec.end(); m++)
                 extendGpdFile << *m << endl;

                 for(unordered_map<string, vector<string> >::iterator m= gapSet.begin(); m!= gapSet.end(); m++) {
                     gapFile << "Gene" << fileTag << ":" <<  m->first << "\t" ;
                     for(vector<string>::iterator n= m->second.begin(); n!= m->second.end(); n++) {
                         gapFile << *n << "\t";
                        }
                     gapFile << endl;
                     }

                 for(unordered_map<string, string>::iterator m= querySeqSet.begin(); m!= querySeqSet.end(); m++)
                     queryFile << ">" << m->first << "\n" << m->second << endl;

                 extendedFile.close();
                 reportSeqFile.close();
                 reportGpdFile.close();
                 extendGpdFile.close();
                 gapFile.close();
                 queryFile.close();

        for(vector<pair<string,string> >::iterator m= unannotatedSeqVec.begin(); m!= unannotatedSeqVec.end(); m++) {
            unannotatedFile<< ">" << (*m).first << "\n" << (*m).second << endl;
           }
omp_unset_lock(&writelock);

      }
   else{
        string tag= *((*i).begin() );
        string seq= getSeq(tag);
omp_set_lock(&writelock);
        ofstream reportSeqFile("report_seq", ios::out|ios::app), reportGpdFile("report.gpd", ios::out|ios::app);
        reportSeqFile<< ">Gene" << fileTag << "\n" << seq << endl;
        reportGpdFile<< "Gene" << fileTag << "\t" << tag << "\tGene" << fileTag << "\t+\t0\t" << seq.length() << "\t0\t" << seq.length() << "\t1\t0,\t" << seq.length() << "," << endl;
        reportSeqFile.close();
        reportGpdFile.close();
omp_unset_lock(&writelock);
       }
   }
}

omp_destroy_lock(&writelock);

unannotatedFile.close();
system("rmdir indexDir");
}
