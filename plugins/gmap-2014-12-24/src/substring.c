static char rcsid[] = "$Id: substring.c 154591 2014-12-04 02:00:32Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "substring.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>		/* For log and exp */

#include "assert.h"
#include "mem.h"
#include "maxent_hr.h"
#include "listdef.h"
#include "list.h"
#include "complement.h"
#include "genome128_hr.h"
#include "mapq.h"
#include "pair.h"		/* For Pair_print_gsnap */
#include "pairdef.h"		/* For State_T */
#include "comp.h"


#define TRIM_MATCH_SCORE 1
#define KNOWN_SPLICESITE_EDGE 1	/* Cannot be 0; otherwise will hit found splice sites */

#define CONSISTENT_TEXT "consistent"
#define TRANSLOCATION_TEXT "translocation"
#define INVERSION_TEXT "inversion"
#define SCRAMBLE_TEXT "scramble"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* mark_mismatches */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Substring_new */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Substring_overlap_p and Substring_insert_length */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


/* splice site probs */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


/* Substring_convert_to_pairs */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif


/* contains known splicesite */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif


/* trimming.  may also want to turn on DEBUG8 in pair.c */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* bad_stretch_p */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


#define LOG_99 -0.01005033585
#define LOG_01 -4.605170186

#if 0
/* Switches on 5 consecutive mismatches */
#define LOG_99_9999 -0.01015034085
#define LOG_99_0001 -9.220390708
#define LOG_25_0001 -10.59663473
#define LOG_25_9999 -1.386394366
#define LOG_01_9999 -4.605270191
#define LOG_01_0001 -13.81551056
#define LOG_75_0001 -9.498022444
#define LOG_75_9999 -0.2877820775
#endif

#if 0
#define LOG_99_999 -0.01105083619
#define LOG_99_001 -6.917805615
#define LOG_25_001 -8.29404964
#define LOG_25_999 -1.387294861
#define LOG_01_999 -4.606170686
#define LOG_01_001 -11.51292546
#define LOG_75_001 -7.195437351
#define LOG_75_999 -0.2886825728
#endif

/* Switches on 4 consecutive mismatches */
#define LOG_99_99 -0.02010067171
#define LOG_99_01 -4.615220522
#define LOG_25_01 -5.991464547
#define LOG_25_99 -1.396344697
#define LOG_01_99 -4.615220522
#define LOG_01_01 -9.210340372
#define LOG_75_01 -4.892852258
#define LOG_75_99 -0.2977324083


static bool print_nsnpdiffs_p;
static bool print_snplabels_p;
static bool show_refdiff_p;

static IIT_T snps_iit;
static int *snps_divint_crosstable;

static IIT_T genes_iit;
static int *genes_divint_crosstable;

static IIT_T splicesites_iit;
static int *splicesites_divint_crosstable;

static int donor_typeint;
static int acceptor_typeint;

static int trim_mismatch_score;
static bool novelsplicingp;
static bool knownsplicingp;
static bool output_sam_p;

static Mode_T mode;

static double genomelength;	/* For BLAST E-value */
static int reject_trimlength;


char *
Endtype_string (Endtype_T endtype) {
  switch (endtype) {
  case END: return "end";
  case INS: return "ins";
  case DEL: return "del";
  case DON: return "don";
  case ACC: return "acc";
  case AMB_DON: return "amb_don";
  case AMB_ACC: return "amb_acc";
  case TERM: return "term";
  default:
    fprintf(stderr,"Unexpected endtype %d\n",endtype);
    abort();
  }
  return "";
}


static char complCode[128] = COMPLEMENT_LC;

#if 0
static char *
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC_OUT(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}
#endif

static char *
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return sequence;
}



#define T Substring_T

struct T {
  bool exactp;
  int nmismatches_whole;        /* Over total substring.  Comes from calling procedure */
  int nmismatches_bothdiff;	/* Over region left after trimming */
  int nmismatches_refdiff;	/* Over region left after trimming */
  /* nsnpdiffs = nmismatches_bothdiff - nmismatches_refdiff */
  int nmatches;			/* Over region left after trimming */

  int trim_left;
  int trim_right;
  bool trim_left_splicep;
  bool trim_right_splicep;

  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;

  Univcoord_T left_genomicseg;	/* left needed to retrieve genomicseg */
  Univcoord_T left;	 /* adjusted by Substring_new for aligndiff */
  Univcoord_T genomicstart;	/* For region corresponding to entire querylength (if extrapolated) */
  Univcoord_T genomicend;

  Endtype_T start_endtype;
  Endtype_T end_endtype;

  int querystart_orig;		/* For part that aligns to genome, pre-trim */
  int queryend_orig;
  int querystart;		/* For part that aligns to genome, post-trim */
  int queryend;
  int querylength;

  Univcoord_T alignstart;	/* For part that aligns to genome, including part that is trimmed (pre-trim) */
  Univcoord_T alignend;

  Univcoord_T alignstart_trim;	/* For part that aligns to genome, excluding part that is trimmed (post-trim) */
  Univcoord_T alignend_trim;

  int alignoffset;
  int extraleft;
  int extraright;

  int genomiclength;
  bool plusp;
  int genestrand;
  bool first_read_p;

  char *genomic_bothdiff; /* In same direction as query.  NULL if same
  			     as query.  Has dashes outside of
  			     querystart..(queryend-1) for indels and
  			     splices.  Has lowercase to indicate
  			     trimmed regions, terminal regions,
  			     deletions, splice dinucleotides, and
  			     mismatches from ref and alt.  Use for
  			     MAPQ computations and scoring. */

  char *genomic_refdiff;  /* Same as above, but lowercase for
			     mismatches from ref only.  For
			     non-SNP-tolerant alignment, this is just
			     a pointer to genomic_bothdiff.  Use for
			     NM and MD computations.  */

  float mapq_loglik;

  /* for splices */
  bool chimera_sensep;

  Univcoord_T splicecoord;
  int splicesites_knowni;	/* Needed for intragenic_splice_p in stage1hr.c */

  bool chimera_knownp;
  bool chimera_novelp;
  Univcoord_T chimera_modelpos;
  int chimera_pos;
  double chimera_prob;

  /* for shortexon (always use *_1 for acceptor and *_2 for donor) */
  /* for donor/acceptor: the ambiguous position */
  Univcoord_T splicecoord_2;
  int splicesites_knowni_2;

  bool chimera_knownp_2;
  bool chimera_novelp_2;
  Univcoord_T chimera_modelpos_2;
  int chimera_pos_2;
  double chimera_prob_2;
};


void
Substring_alias_circular (T this) {
  Chrpos_T chrlength;

  if (this != NULL) {
    chrlength = this->chrlength;

    this->left_genomicseg += chrlength;
    this->left += chrlength;
    this->genomicstart += chrlength;
    this->genomicend += chrlength;
    this->alignstart += chrlength;
    this->alignend += chrlength;
    this->alignstart_trim += chrlength;
    this->alignend_trim += chrlength;
    this->chimera_modelpos += chrlength;
    this->chimera_modelpos_2 += chrlength;
  }

  return;
}


void
Substring_unalias_circular (T this) {
  Chrpos_T chrlength;

  if (this != NULL) {
    chrlength = this->chrlength;

    this->left_genomicseg -= chrlength;
    this->left -= chrlength;
    this->genomicstart -= chrlength;
    this->genomicend -= chrlength;
    this->alignstart -= chrlength;
    this->alignend -= chrlength;
    this->alignstart_trim -= chrlength;
    this->alignend_trim -= chrlength;
    this->chimera_modelpos -= chrlength;
    this->chimera_modelpos_2 -= chrlength;
  }

  return;
}



static void
fill_w_dashes (char *string, int start, int end) {
  int i;

  for (i = start; i < end; i++) {
    string[i] = '-';
  }
  return;
}


#if 0
static bool
dibase_mismatch_p (char leftnt, char leftcolor, char nt, char rightcolor, char rightnt) {
  char left_predict, right_predict;

  switch (toupper(leftnt)) {
  case 'A':
    switch (leftcolor) {
    case '0': left_predict = 'A'; break;
    case '1': left_predict = 'C'; break;
    case '2': left_predict = 'G'; break;
    case '3': left_predict = 'T'; break;
    }
    break;
  case 'C':
    switch (leftcolor) {
    case '0': left_predict = 'C'; break;
    case '1': left_predict = 'A'; break;
    case '2': left_predict = 'T'; break;
    case '3': left_predict = 'G'; break;
    }
    break;
  case 'G':
    switch (leftcolor) {
    case '0': left_predict = 'G'; break;
    case '1': left_predict = 'T'; break;
    case '2': left_predict = 'A'; break;
    case '3': left_predict = 'C'; break;
    }
    break;
  case 'T':
    switch (leftcolor) {
    case '0': left_predict = 'T'; break;
    case '1': left_predict = 'G'; break;
    case '2': left_predict = 'C'; break;
    case '3': left_predict = 'A'; break;
    }
    break;
  }

  switch (toupper(rightnt)) {
  case 'A':
    switch (rightcolor) {
    case '0': right_predict = 'A'; break;
    case '1': right_predict = 'C'; break;
    case '2': right_predict = 'G'; break;
    case '3': right_predict = 'T'; break;
    }
    break;
  case 'C':
    switch (rightcolor) {
    case '0': right_predict = 'C'; break;
    case '1': right_predict = 'A'; break;
    case '2': right_predict = 'T'; break;
    case '3': right_predict = 'G'; break;
    }
    break;
  case 'G':
    switch (rightcolor) {
    case '0': right_predict = 'G'; break;
    case '1': right_predict = 'T'; break;
    case '2': right_predict = 'A'; break;
    case '3': right_predict = 'C'; break;
    }
    break;
  case 'T':
    switch (rightcolor) {
    case '0': right_predict = 'T'; break;
    case '1': right_predict = 'G'; break;
    case '2': right_predict = 'C'; break;
    case '3': right_predict = 'A'; break;
    }
    break;
  }

  if (left_predict != right_predict) {
    /* printf("left_predict %c != right_predict %c\n",left_predict,right_predict); */
    return false;
  } else if (left_predict == toupper(nt)) {
    return false;
  } else {
    return true;
  }
}
#endif


#if 0
/* query is numbers */
static void
mark_mismatches_dibase (char *genome, char *query, int start, int end, bool plusp) {
  int i;

  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",genome));
  debug1(printf("mark:   "));

  if (plusp) {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start],genome[start+1])) {
      genome[start] = (char) tolower(genome[start]);
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i-1],genome[i],query[i],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i-1],genome[i],query[i],genome[i+1]); */
	genome[i] = (char) tolower(genome[i]);
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-2],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      genome[end-1] = (char) tolower(genome[end-1]);
    }

  } else {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start+1],genome[start+1])) {
      genome[start] = (char) tolower(genome[start]);
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i],genome[i],query[i+1],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i],genome[i],query[i+1],genome[i+1]); */
	genome[i] = (char) tolower(genome[i]);
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-1],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      genome[end-1] = (char) tolower(genome[end-1]);
    }
  }

  debug1(printf("\n"));

  return;
}
#endif


#if 0
static int
count_mismatches_dibase (char *genome, char *query, int start, int end, bool plusp) {
  int nmismatches = 0, i;

  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",genome));
  debug1(printf("mark:   "));

  if (plusp) {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start],genome[start+1])) {
      nmismatches++;
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i-1],genome[i],query[i],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i-1],genome[i],query[i],genome[i+1]); */
	nmismatches++;
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-2],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      nmismatches++;
    }

  } else {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start+1],genome[start+1])) {
      nmismatches++;
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i],genome[i],query[i+1],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i],genome[i],query[i+1],genome[i+1]); */
	nmismatches++;
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-1],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      nmismatches++;
    }
  }

  debug1(printf("\n"));
  debug1(printf("%d mismatches\n\n",nmismatches));

  return nmismatches;
}
#endif


static int
trim_left_end (Compress_T query_compress, Univcoord_T left, int querystart, int queryend,
	       int querylength, bool plusp, int genestrand, bool first_read_p, int trim_mismatch_score) {
  int bestscore, score;
  int trim5, alignlength, pos, prevpos, i;
  int nmismatches;

#ifdef HAVE_ALLOCA
  int *mismatch_positions = (int *) ALLOCA(querylength*sizeof(int));
#else
  int mismatch_positions[MAX_READLENGTH];
#endif


  debug8(printf("Entered trim_left_end with querystart %d, queryend %d\n",querystart,queryend));

  alignlength = queryend - querystart;
  bestscore = 0;
  score = 0;
  trim5 = 0;

  if (plusp == true) {
    nmismatches = Genome_mismatches_right_trim(mismatch_positions,/*max_mismatches*/alignlength,
					       query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
					       plusp,genestrand,first_read_p);
    debug8(printf("%d mismatches:",nmismatches));
    debug8(
	   for (i = 0; i < nmismatches; i++) {
	     printf(" %d",mismatch_positions[i]);
	   }
	   printf("\n");
	   );

    prevpos = queryend;
    for (i = 0; i < nmismatches; i++) {
      pos = mismatch_positions[i] + 1; /* the position just after the mismatch */
      score += (prevpos - pos)*TRIM_MATCH_SCORE;
      if (score >= bestscore) {	/* want >= and not >, so ties go to end */
	bestscore = score;
	trim5 = pos;
      }
      debug8(printf("Trim left pos %d, score %d, trim5 %d\n",pos,score,trim5));

      score += trim_mismatch_score; /* For the mismatch */
      if (score < 0) {
	score = 0;
      }
      prevpos = pos - 1;		/* On the mismatch */
    }

  } else {
    nmismatches = Genome_mismatches_left_trim(mismatch_positions,/*max_mismatches*/alignlength,
					      query_compress,left,/*pos5*/querylength - queryend,
					      /*pos3*/querylength - querystart,plusp,genestrand,first_read_p);

    debug8(printf("%d mismatches:",nmismatches));
    debug8(
	   for (i = 0; i < nmismatches; i++) {
	     printf(" %d",querylength - 1 - mismatch_positions[i]);
	   }
	   printf("\n");
	   );

    prevpos = queryend;
    for (i = 0; i < nmismatches; i++) {
      pos = (querylength - 1 - mismatch_positions[i]) + 1; /* the position just after the mismatch */
      score += (prevpos - pos)*TRIM_MATCH_SCORE;
      if (score >= bestscore) {	/* want >= and not >, so ties go to end */
	bestscore = score;
	trim5 = pos;
      }
      debug8(printf("Trim left pos %d, score %d, trim5 %d\n",pos,score,trim5));

      score += trim_mismatch_score; /* For the mismatch */
      if (score < 0) {
	score = 0;
      }
      prevpos = pos - 1;		/* On the mismatch */
    }
  }

  score += prevpos*TRIM_MATCH_SCORE;
  if (score >= bestscore) {	/* want >= and not >, so ties go to end */
    bestscore = score;
    trim5 = 0;
  }

  debug8(printf("Trim left pos 0, score %d, trim5 %d\n",score,trim5));
  debug8(printf("\n"));

  return trim5;
}



/* querystart (typically 0) and queryend (typically querylength) are exclusive */
/* sequences may have had lower case characters marked */
static int
trim_right_end (Compress_T query_compress, Univcoord_T left, int querystart, int queryend,
		int querylength, bool plusp, int genestrand, bool first_read_p, int trim_mismatch_score) {
  int bestscore, score;
  int trim3, alignlength, pos, prevpos, i;
  int nmismatches;

#ifdef HAVE_ALLOCA
  int *mismatch_positions = (int *) ALLOCA(querylength*sizeof(int));
#else
  int mismatch_positions[MAX_READLENGTH];
#endif

  debug8(printf("Entered trim_right_end with querystart %d, queryend %d\n",querystart,queryend));

  alignlength = queryend - querystart;
  bestscore = 0;
  score = 0;
  trim3 = 0;

  if (plusp == true) {
    nmismatches = Genome_mismatches_left_trim(mismatch_positions,/*max_mismatches*/alignlength,
					      query_compress,left,/*pos5*/querystart,/*pos3*/queryend,
					      plusp,genestrand,first_read_p);

    debug8(printf("%d mismatches:",nmismatches));
    debug8(
	   for (i = 0; i < nmismatches; i++) {
	     printf(" %d",mismatch_positions[i]);
	   }
	   printf("\n");
	   );

    prevpos = querystart - 1;
    for (i = 0; i < nmismatches; i++) {
      pos = mismatch_positions[i] - 1; /* the position just before the mismatch */
      score += (pos - prevpos)*TRIM_MATCH_SCORE;
      if (score >= bestscore) {	/* want >= and not >, so ties go to end */
	bestscore = score;
	trim3 = querylength - pos - 1;
      }
      debug8(printf("Trim right pos %d, score %d, trim3 %d\n",pos,score,trim3));
      
      score += trim_mismatch_score; /* For the mismatch */
      if (score < 0) {
	score = 0;
      }
      prevpos = pos + 1;		/* On the mismatch */
    }

  } else {
    nmismatches = Genome_mismatches_right_trim(mismatch_positions,/*max_mismatches*/alignlength,
					       query_compress,left,/*pos5*/querylength - queryend,
					       /*pos3*/querylength - querystart,plusp,genestrand,first_read_p);

    debug8(printf("%d mismatches:",nmismatches));
    debug8(
	   for (i = 0; i < nmismatches; i++) {
	     printf(" %d",querylength - 1 - mismatch_positions[i]);
	   }
	   printf("\n");
	   );

    prevpos = querystart - 1;
    for (i = 0; i < nmismatches; i++) {
      pos = (querylength - 1 - mismatch_positions[i]) - 1; /* the position just before the mismatch */
      score += (pos - prevpos)*TRIM_MATCH_SCORE;
      if (score >= bestscore) {	/* want >= and not >, so ties go to end */
	bestscore = score;
	trim3 = querylength - pos - 1;
      }
      debug8(printf("Trim right pos %d, score %d, trim3 %d\n",pos,score,trim3));
      
      score += trim_mismatch_score; /* For the mismatch */
      if (score < 0) {
	score = 0;
      }
      prevpos = pos + 1;		/* On the mismatch */
    }
  }

  score += (queryend - 1 - prevpos)*TRIM_MATCH_SCORE;
  if (score >= bestscore) { /* want >= and not >, so ties go to end */
    bestscore = score;
    trim3 = 0;
  }

  debug8(printf("Trim right pos %d, score %d, trim3 %d\n",queryend-1,score,trim3));
  debug8(printf("\n"));

  return trim3;
}


/* Uses hmm */
/* Modified for Stage3_bad_stretch_p */
bool
Substring_bad_stretch_p (T this, Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  int alignlength, startpos, endpos, pos, i;
  float vprob_good, vprob_bad, prev_vprob_good, prev_vprob_bad, good_incr_prob, bad_incr_prob;
#ifdef DEBUG9
  bool result;
  State_T vstate_good[MAX_READLENGTH], vstate_bad[MAX_READLENGTH], state;
#endif

#ifdef HAVE_ALLOCA
  int *mismatch_positions;
#else
  int mismatch_positions[MAX_READLENGTH];
#endif


  debug9(printf("Entered bad_stretch_p with querystart_orig %d, queryend_orig %d, plusp %d\n",
		this->querystart_orig,this->queryend_orig,this->plusp));

  alignlength = this->queryend_orig - this->querystart_orig;
#ifdef HAVE_ALLOCA
  mismatch_positions = (int *) ALLOCA(alignlength*sizeof(int));
#endif

  if (this->plusp == true) {
    startpos = this->querystart_orig;
    endpos = this->queryend_orig;
    debug9(printf("Calling Genome_mismatches_left from pos5 %d to pos3 %d\n",startpos,endpos));
#ifdef DEBUG9
    nmismatches = Genome_mismatches_left(mismatch_positions,/*max_mismatches*/alignlength,
					 query_compress_fwd,this->left,/*pos5*/startpos,/*pos3*/endpos,
					 /*plusp*/true,this->genestrand,this->first_read_p);
#endif
  } else {
    startpos = this->querylength - this->queryend_orig;
    endpos = this->querylength - this->querystart_orig;
    debug9(printf("Calling Genome_mismatches_left from pos5 %d to pos3 %d\n",startpos,endpos));
#ifdef DEBUG9
    nmismatches = Genome_mismatches_left(mismatch_positions,/*max_mismatches*/alignlength,
					 query_compress_rev,this->left,/*pos5*/startpos,/*pos3*/endpos,
					 /*plusp*/false,this->genestrand,this->first_read_p);
#endif
  }

  debug9(printf("%d mismatches:",nmismatches));
  debug9(
	 for (i = 0; i < nmismatches; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");
	 );

  /* Initialize priors */
#ifdef COMPUTE_LOG
  prev_vprob_good = log(0.99);
  prev_vprob_bad = log(0.01);
#else
  prev_vprob_good = LOG_99;
  prev_vprob_bad = LOG_01;
#endif

  i = 0;
  for (pos = startpos; pos < endpos; pos++) {
    if (pos != mismatch_positions[i]) {
      /* match */
#ifdef DEBUG9
      printf("querypos %d (match): ",pos);
#endif

      /* state: GOOD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.99) + log(/*transition_prob*/0.99);
      bad_incr_prob = log(/*emission_prob*/0.99) + log(/*transition_prob*/0.01);
#else
      good_incr_prob = LOG_99_99;
      bad_incr_prob = LOG_99_01;
#endif

#ifdef DEBUG9
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	printf("%.2f + %.2f (GOOD->GOOD*) or %.2f + %.2f (BAD->GOOD ).  ",
	       prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob);
	vprob_good = prev_vprob_good + good_incr_prob;
	vstate_good[pos] = GOOD;
      } else {
	printf("%.2f + %.2f (GOOD->GOOD ) or %.2f + %.2f (BAD->GOOD!).  ",
	       prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob);
	vprob_good = prev_vprob_bad + bad_incr_prob;
	vstate_good[pos] = BAD;
      }
#else
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_good = prev_vprob_good + good_incr_prob;
      } else {
	/* vprob_good = prev_vprob_bad + bad_incr_prob; */
	return true;
      }
#endif

      /* state: BAD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.25) + log(/*transition_prob*/0.01);
      bad_incr_prob = log(/*emission_prob*/0.25) + log(/*transition_prob*/0.99);
#else
      good_incr_prob = LOG_25_01;
      bad_incr_prob = LOG_25_99;
#endif
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_bad = prev_vprob_good + good_incr_prob;
#ifdef DEBUG9
	printf("%.2f + %.2f (GOOD->BAD*) or %.2f + %.2f (BAD->BAD ).\n",
	       prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob);
	vstate_bad[pos] = GOOD;
#endif
      } else {
	vprob_bad = prev_vprob_bad + bad_incr_prob;
#ifdef DEBUG9
	printf("%.2f + %.2f (GOOD->BAD ) or %.2f + %.2f (BAD->BAD*).\n",
	       prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob);
	vstate_bad[pos] = BAD;
#endif
      }

    } else {
      /* mismatch */
#ifdef DEBUG9
      printf("querypos %d (mismatch): ",pos);
#endif

      /* state: GOOD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.01) + log(/*transition_prob*/0.99);
      bad_incr_prob = log(/*emission_prob*/0.01) + log(/*transition_prob*/0.01);
#else
      good_incr_prob = LOG_01_99;
      bad_incr_prob = LOG_01_01;
#endif
#ifdef DEBUG9
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	printf("%.2f + %.2f (GOOD->GOOD*) or %.2f + %.2f (BAD->GOOD ).  ",
	       prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob);
	vprob_good = prev_vprob_good + good_incr_prob;
	vstate_good[pos] = GOOD;
      } else {
	printf("%.2f + %.2f (GOOD->GOOD ) or %.2f + %.2f (BAD->GOOD!).  ",
	       prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob);
	vprob_good = prev_vprob_bad + bad_incr_prob;
	vstate_good[pos] = BAD;
      }
#else
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_good = prev_vprob_good + good_incr_prob;
      } else {
	/* vprob_good = prev_vprob_bad + bad_incr_prob; */
	return true;
      }
#endif

      /* state: BAD */
#ifdef COMPUTE_LOG
      good_incr_prob = log(/*emission_prob*/0.75) + log(/*transition_prob*/0.01);
      bad_incr_prob = log(/*emission_prob*/0.75) + log(/*transition_prob*/0.99);
#else
      good_incr_prob = LOG_75_01;
      bad_incr_prob = LOG_75_99;
#endif
      if (prev_vprob_good + good_incr_prob > prev_vprob_bad + bad_incr_prob) {
	vprob_bad = prev_vprob_good + good_incr_prob;
#ifdef DEBUG9
	printf("%.2f + %.2f (GOOD->BAD*) or %.2f + %.2f (BAD->BAD ).\n",
	       prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob);
	vstate_bad[pos] = GOOD;
#endif
      } else {
	vprob_bad = prev_vprob_bad + bad_incr_prob;
#ifdef DEBUG9
	printf("%.2f + %.2f (GOOD->BAD ) or %.2f + %.2f (BAD->BAD*).\n",
	       prev_vprob_good,good_incr_prob,prev_vprob_bad,bad_incr_prob);
	vstate_bad[pos] = BAD;
#endif
      }

      i++;
    }
    
    prev_vprob_good = vprob_good;
    prev_vprob_bad = vprob_bad;
  }

#ifdef DEBUG9
  result = false;
  if (prev_vprob_good > prev_vprob_bad) {
    printf("Final state: GOOD\n");
    state = GOOD;
  } else {
    printf("Final state: BAD\n");
    state = BAD;
    result = true;
  }

  for (pos = endpos-1; pos >= startpos; pos--) {
    if (state == GOOD) {
      printf("querypos %d: GOOD\n",pos);
      state = vstate_good[pos];
    } else {
      printf("querypos %d: BAD\n",pos);
      state = vstate_bad[pos];
      result = true;
    }
  }
  return result;

#else
  return false;
#endif

}





void
Substring_free (T *old) {
  if ((*old)->genomic_bothdiff != NULL) {
    if ((*old)->genomic_refdiff != (*old)->genomic_bothdiff) {
      FREE_OUT((*old)->genomic_refdiff);
    }
    FREE_OUT((*old)->genomic_bothdiff);
  }

  FREE_OUT(*old);
  return;
}


bool
Substring_contains_p (T this, int querypos) {
  if (querypos >= this->querystart && querypos < this->queryend) {
    return true;
  } else {
    return false;
  }
}


void
Substring_print_ends (T this, int chrnum) {
  if (this == NULL) {
    printf("NA ");
  } else {
    printf("#%d:%u..%u ",chrnum,this->alignstart,this->alignend);
  }
  return;
}


int
Substring_compare (T substring1, T substring2, int alias1, int alias2, Chrpos_T chrlength1, Chrpos_T chrlength2) {
  Univcoord_T alignstart1, alignend1, alignstart2, alignend2;

  if (substring1 == NULL && substring2 == NULL) {
    return 0;
  } else if (substring1 == NULL) {
    return -1;
  } else if (substring2 == NULL) {
    return +1;
  } else {
    alignstart1 = substring1->alignstart;
    alignend1 = substring1->alignend;
    if (alias1 < 0) {
      alignstart1 += chrlength1;
      alignend1 += chrlength1;
    }

    alignstart2 = substring2->alignstart;
    alignend2 = substring2->alignend;
    if (alias2 < 0) {
      alignstart2 += chrlength2;
      alignend2 += chrlength2;
    }

    if (alignstart1 < alignstart2) {
      return -1;
    } else if (alignstart1 > alignstart2) {
      return +1;
    } else if (alignend1 < alignend2) {
      return -1;
    } else if (alignend1 > alignend2) {
      return +1;
    } else {
      return 0;
    }
  }
}


#if 0
bool
Substring_equal_p (T substring1, T substring2) {

  if (substring1 == NULL && substring2 == NULL) {
    return true;
  } else if (substring1 == NULL) {
    return false;
  } else if (substring2 == NULL) {
    return false;
  } else if (substring1->alignstart != substring2->alignstart) {
    return false;
  } else if (substring1->alignend != substring2->alignend) {
    return false;
  } else {
    return true;
  }
}
#endif


bool
Substring_overlap_p (T substring1, T substring2) {
  Univcoord_T low1, high1, low2, high2;

  if (substring1->plusp == true) {
    low1 = substring1->alignstart;
    high1 = substring1->alignend;
  } else {
    low1 = substring1->alignend;
    high1 = substring1->alignstart;
  }

  if (substring2->plusp == true) {
    low2 = substring2->alignstart;
    high2 = substring2->alignend;
  } else {
    low2 = substring2->alignend;
    high2 = substring2->alignstart;
  }

  debug3(printf("Checking overlap between %u..%u and %u..%u",low1,high1,low2,high2));

  if (high2 < low1) {
    debug3(printf(" => no because %u < %u\n",high2,low1));
    return false;
  } else if (low2 > high1) {
    debug3(printf(" => no because %u > %u\n",low2,high1));
    return false;
  } else {
    debug3(printf(" => yes\n"));
    return true;
  }
}


Chrpos_T
Substring_insert_length (T substring5, T substring3) {
  Univcoord_T pos5, pos3;

  pos5 = substring5->genomicstart;
  debug3(printf("pos5 %u\n",substring5->genomicstart));

  pos3 = substring3->genomicend;
  debug3(printf("pos3 %u\n",substring3->genomicend));

  if (pos5 > pos3) {
    return pos5 - pos3;
  } else {
    return pos3 - pos5;
  }
}


bool
Substring_overlap_point_trimmed_p (T substring, Univcoord_T endpos) {
  Univcoord_T low, high;

  if (substring == NULL) {
    return false;
  }

  if (substring->plusp == true) {
    low = substring->alignstart_trim;
    high = substring->alignend_trim;
    debug3(printf("Checking overlap between plus %u..%u and %u",low,high,endpos));
  } else {
    low = substring->alignend_trim;
    high = substring->alignstart_trim;
    debug3(printf("Checking overlap between minus %u..%u and %u",low,high,endpos));
  }


  if (endpos < low) {
    debug3(printf(" => no because %u < %u\n",endpos,low));
    return false;
  } else if (endpos > high) {
    debug3(printf(" => no because %u > %u\n",endpos,high));
    return false;
  } else {
    debug3(printf(" => yes\n"));
    return true;
  }
}


Univcoord_T
Substring_overlap_segment_trimmed (T substring1, T substring2) {
  Univcoord_T maxlow, minhigh;
  Univcoord_T low1, high1, low2, high2;

  if (substring1->plusp == true) {
    low1 = substring1->alignstart_trim;
    high1 = substring1->alignend_trim;
  } else {
    low1 = substring1->alignend_trim;
    high1 = substring1->alignstart_trim;
  }

  if (substring2->plusp == true) {
    low2 = substring2->alignstart_trim;
    high2 = substring2->alignend_trim;
  } else {
    low2 = substring2->alignend_trim;
    high2 = substring2->alignstart_trim;
  }

  debug3(printf("Checking overlap between %u..%u and %u..%u",low1,high1,low2,high2));

  if (high2 < low1) {
    debug3(printf(" => no because %u < %u\n",high2,low1));
    return 0;
  } else if (low2 > high1) {
    debug3(printf(" => no because %u > %u\n",low2,high1));
    return 0;
  } else {
    maxlow = (low1 > low2) ? low1 : low2;
    minhigh = (high1 < high2) ? high1 : high2;
    debug3(printf(" => yes.  maxlow %llu, minhigh %llu.  returning %llu\n",
		  maxlow,minhigh,maxlow + (minhigh - maxlow)/2));
    return maxlow + (minhigh - maxlow)/2;
  }
}



#if 0
Chrpos_T
Substring_insert_length_trimmed (T substring5, T substring3, int hit5_trim_left, int hit5_trim_right, int hit3_trim_left, int hit3_trim_right) {
  Univcoord_T rightpos, leftpos, pos;
  int right_extent, left_extent;
  
  debug3(printf("substring5 %d..%d out of %d.  aligned to %u..%u\n",
		substring5->querystart,substring5->queryend,substring5->querylength,
		substring5->alignstart_trim - substring5->chroffset,substring5->alignend_trim - substring5->chroffset));
  debug3(printf("substring3 %d..%d out of %d.  aligned to %u..%u\n",
		substring3->querystart,substring3->queryend,substring3->querylength,
		substring3->alignstart_trim - substring3->chroffset,substring3->alignend_trim - substring3->chroffset));
  if (substring5->plusp == true) {
    if (substring5->alignend_trim > substring3->alignend_trim) {
      rightpos = substring5->alignend_trim;
      right_extent = substring5->querylength - substring5->queryend - hit5_trim_right;
    } else {
      rightpos = substring3->alignend_trim;
      right_extent = substring3->querylength - substring3->queryend - hit3_trim_right;
    }
    debug3(printf("rightpos = %u, right_extent = %d\n",rightpos,right_extent));

    if (substring5->alignstart_trim < substring3->alignstart_trim) {
      leftpos = substring5->alignstart_trim;
      left_extent = substring5->querystart - hit5_trim_left;
    } else {
      leftpos = substring3->alignstart_trim;
      left_extent = substring3->querystart - hit3_trim_left;
    }
    debug3(printf("leftpos = %u, left_extent = %d\n",leftpos,left_extent));

    return (rightpos - leftpos) + right_extent + left_extent;

  } else {
    if (substring3->alignstart_trim > (pos = substring5->alignend_trim + (substring5->queryend - substring5->querystart /* - hit5_trim_left */) +
				       (substring3->querylength - substring3->querystart - hit3_trim_right))) {
      return 0;
    } else {
      debug3(printf("minus: substring5->alignend_trim %u - substring3->alignstart_trim %u +\n",
		    substring5->alignend_trim - substring5->chroffset,substring3->alignstart_trim - substring3->chroffset));
      debug3(printf("      substring5->queryend %d - substring5->querystart %d +\n",substring5->queryend,substring5->querystart));
      debug3(printf("      substring3->querylength %d - substring3->querystart %d - hit3_trim_right %d\n",
		    substring3->querylength,substring3->querystart,hit3_trim_right));
      return (pos - substring3->alignstart_trim);
    }
  }
}
#endif



static void
mark_mismatches_cmet_gsnap (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("\n"));
  debug1(printf("first_read_p %d\n",first_read_p));
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'G' && query[i] == 'A') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'C' && query[i] == 'T') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }

  return;
}


static void
mark_mismatches_cmet_sam (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'G' && query[i] == 'A') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'A';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'C' && query[i] == 'T') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'T';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }
  
  return;
}



static void
mark_mismatches_atoi_gsnap (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'T' && query[i] == 'C') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'A' && query[i] == 'G') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }
  
  return;
}



static void
mark_mismatches_atoi_sam (char *gbuffer, char *query, int start, int end, int genestrand) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  if (genestrand == +2) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'T' && query[i] == 'C') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'C';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'A' && query[i] == 'G') {
	debug1(printf("."));
#if 0
	/* Want to show mismatches */
	gbuffer[i] = 'G';		/* Avoids showing mismatches in MD and NM strings */
#endif
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }

  return;
}



void
Substring_setup (bool print_nsnpdiffs_p_in, bool print_snplabels_p_in,
		 bool show_refdiff_p_in, IIT_T snps_iit_in, int *snps_divint_crosstable_in,
		 IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		 IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
		 int donor_typeint_in, int acceptor_typeint_in, int trim_mismatch_score_in,
		 bool novelsplicingp_in, bool knownsplicingp_in,
		 bool output_sam_p_in, Mode_T mode_in, Univcoord_T genomelength_in,
		 int reject_trimlength_in) {
  print_nsnpdiffs_p = print_nsnpdiffs_p_in;
  print_snplabels_p = print_snplabels_p_in;
  show_refdiff_p = show_refdiff_p_in;

  snps_iit = snps_iit_in;
  snps_divint_crosstable = snps_divint_crosstable_in;

  genes_iit = genes_iit_in;
  genes_divint_crosstable = genes_divint_crosstable_in;

  splicesites_iit = splicesites_iit_in;
  splicesites_divint_crosstable = splicesites_divint_crosstable_in;

  donor_typeint = donor_typeint_in;
  acceptor_typeint = acceptor_typeint_in;

  trim_mismatch_score = trim_mismatch_score_in;
  
  novelsplicingp = novelsplicingp_in;
  knownsplicingp = knownsplicingp_in;

  output_sam_p = output_sam_p_in;
  mode = mode_in;

  genomelength = (double) genomelength_in;
  reject_trimlength = reject_trimlength_in;

  return;
}



static char *
embellish_genomic (char *genomic_diff, char *query, int querystart, int queryend, int querylength,
		   int alignoffset, int extraleft, int extraright, int genestrand) {
  char *result;
  int i, j, k;

  result = (char *) MALLOC_OUT((querylength+1) * sizeof(char));
  result[querylength] = '\0';

  /* Add aligned region with lower-case diffs, surrounded by dashes */
  fill_w_dashes(result,0,querystart);
  debug2(printf("g1: %s\n",result));

  strncpy(&(result[querystart]),&(genomic_diff[alignoffset]),queryend-querystart);
  debug2(printf("g1: %s\n",result));

  if (mode == STANDARD) {
    /* Skip */
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    mark_mismatches_cmet_gsnap(result,query,querystart,queryend,genestrand);
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    mark_mismatches_atoi_gsnap(result,query,querystart,queryend,genestrand);
  } else {
    abort();
  }

  fill_w_dashes(result,queryend,querylength);
  debug2(printf("g1: %s\n",result));

  /* Add terminal ends as lower-case */
  for (k = 0, i = querystart-1, j = alignoffset-1; k < extraleft && i >= 0; k++, i--, j--) {
    result[i] = (char) tolower(genomic_diff[j]);
#if 0
    printf("k=%d i=%d result[i]=%c\n",k,i,result[i]);
#endif
    assert(result[i] == 'a' || result[i] == 'c' || result[i] == 'g' || result[i] == 't' || result[i] == 'n');
    }
  for (k = 0, i = queryend, j = alignoffset+queryend-querystart; k < extraright && i < querylength; k++, i++, j++) {
    result[i] = (char) tolower(genomic_diff[j]);
#if 0
    printf("k=%d i=%d result[i]=%c\n",k,i,result[i]);
#endif
    assert(result[i] == 'a' || result[i] == 'c' || result[i] == 'g' || result[i] == 't' || result[i] == 'n');
  }
  debug2(printf("g1: %s\n",result));

  return result;
}


static char *
embellish_genomic_sam (char *genomic_diff, char *query, int querystart, int queryend, int querylength,
		       int genomiclength, int alignoffset, int genestrand) {
  char *result;
  int i, j, k;

  result = (char *) MALLOC_OUT((querylength+1) * sizeof(char));
  result[querylength] = '\0';

  strncpy(&(result[querystart]),&(genomic_diff[alignoffset]),queryend-querystart);

  if (mode == STANDARD) {
    /* Skip */
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    mark_mismatches_cmet_sam(result,query,querystart,queryend,genestrand);
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    mark_mismatches_atoi_sam(result,query,querystart,queryend,genestrand);
  } else {
    abort();
  }

  /* Add terminal ends as lower-case */
  for (k = 0, i = querystart-1, j = alignoffset-1; i >= 0 && j >= 0; k++, i--, j--) {
    if (query[i] == genomic_diff[j]) {
      result[i] = genomic_diff[j];
    } else {
      result[i] = (char) tolower(genomic_diff[j]);
    }
#if 0
    printf("k=%d i=%d j=%d result[i]=%c\n",k,i,j,result[i]);
#endif
  }

  for (k = 0, i = queryend, j = alignoffset+queryend-querystart; i < querylength && j < genomiclength; k++, i++, j++) {
    if (query[i] == genomic_diff[j]) {
      result[i] = genomic_diff[j];
    } else {
      result[i] = (char) tolower(genomic_diff[j]);
    }
#if 0
    printf("k=%d i=%d j=%d result[i]=%c\n",k,i,j,result[i]);
#endif
  }

  return result;
}





/************************************************************************/


T
Substring_new (int nmismatches_whole, Chrnum_T chrnum, Univcoord_T chroffset,
	       Univcoord_T chrhigh, Chrpos_T chrlength, Univcoord_T left,
	       Univcoord_T genomicstart, Univcoord_T genomicend,
	       Compress_T query_compress, Endtype_T start_endtype, Endtype_T end_endtype,
	       int querystart, int queryend, int querylength,
	       Univcoord_T alignstart, Univcoord_T alignend, int genomiclength,
	       int extraleft, int extraright, bool exactp,
	       bool plusp, int genestrand, bool first_read_p,
	       bool trim_left_p, bool trim_right_p, int minlength) {
  T new;
  int aligndiff;
  int nmatches;
  double prob1, prob2;
  int nonterminal_trim = 0;


  /* General test for goodness over original region */
  nmatches = queryend - querystart - nmismatches_whole;
  if (nmatches - 3*nmismatches_whole < 0) {
    debug2(printf("Substring fails general test for goodness with %d matches and %d mismatches\n",
		  nmatches,nmismatches_whole));
    return (T) NULL;
  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
  }

  new->exactp = exactp;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;

  new->left_genomicseg = left;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  new->start_endtype = start_endtype;
  new->end_endtype = end_endtype;

  new->querystart_orig = new->querystart = querystart;
  new->queryend_orig = new->queryend = queryend;
  new->querylength = querylength;

  new->alignstart = new->alignstart_trim = alignstart;
  new->alignend = new->alignend_trim = alignend;

  new->extraleft = extraleft;
  new->extraright = extraright;

  new->genomiclength = genomiclength;
  new->plusp = plusp;
  new->genestrand = genestrand;
  new->first_read_p = first_read_p;

  new->chimera_knownp = false;
  new->chimera_knownp_2 = false;
  new->chimera_novelp = false;
  new->chimera_novelp_2 = false;

  debug2(printf("Entered Substring_new with chrnum %d (chroffset %u, chrhigh %u)\n",
		chrnum,chroffset,chrhigh));

  /* Compute coordinates */
  if (plusp == true) {
    new->alignoffset = alignstart - genomicstart;
    aligndiff = /* (alignstart - genomicstart) - querystart = */ new->alignoffset - querystart;
    left += aligndiff;
    new->left = left;

    debug2(printf("\n"));
    debug2(printf("querylength is %d, genomiclength is %d, alignstart is %u, alignend is %u, genomicstart is %u, genomicend is %u, alignoffset is %d\n",
		  querylength,genomiclength,alignstart,alignend,genomicstart,genomicend,new->alignoffset));
    
  } else {
    new->alignoffset = genomicstart - alignstart;
    aligndiff = (alignend - genomicend) - (querylength - queryend);
    left += aligndiff;
    new->left = left;

    debug2(printf("\n"));
    debug2(printf("querylength is %d, genomiclength is %d, alignstart is %u, alignend is %u, genomicstart is %u, genomicend is %u, alignoffset is %d\n",
		  querylength,genomiclength,alignstart,alignend,genomicstart,genomicend,new->alignoffset));
  }
  assert(alignstart <= chrhigh);
  assert(alignend <= chrhigh);

  /* Assign new->nmismatches_whole */
  new->nmismatches_whole = nmismatches_whole;

  /* Initialize these so Substring_free knows what to do */
  new->genomic_bothdiff = (char *) NULL;
  new->genomic_refdiff = (char *) NULL;

  /* Do trimming */
  debug8(printf("trim_left_p %d, trim_right_p %d\n",trim_left_p,trim_right_p));

  if (trim_left_p == false) {
    new->trim_left = 0;

  } else if (new->start_endtype == TERM) {
    /* Accept true terminals generated by GSNAP procedure */
    new->trim_left = trim_left_end(query_compress,left,querystart,queryend,querylength,plusp,genestrand,first_read_p,
				   /*trim_mismatch_score*/-3);

  } else {
    new->trim_left = trim_left_end(query_compress,left,querystart,queryend,querylength,plusp,genestrand,first_read_p,
				   trim_mismatch_score);
    nonterminal_trim += new->trim_left;
  }

  if (trim_right_p == false) {
    new->trim_right = 0;

  } else if (new->end_endtype == TERM) {
    /* Accept true terminals generated by GSNAP procedure */
    new->trim_right = trim_right_end(query_compress,left,querystart,queryend,querylength,plusp,genestrand,first_read_p,
				     /*trim_mismatch_score*/-3);

  } else {
    new->trim_right = trim_right_end(query_compress,left,querystart,queryend,querylength,plusp,genestrand,first_read_p,
				     trim_mismatch_score);
    nonterminal_trim += new->trim_right;
  }

  debug8(printf("Nonterminal trim %d vs reject_trimlength %d\n",nonterminal_trim,reject_trimlength));
  if (nonterminal_trim >= reject_trimlength) {
    /* Reject non-terminal alignments (including those by sarray search) with excessive trim */
    /* Keep true terminals for now in case they help lead to GMAP alignments */
    debug8(printf("Nonterminal trims %d exceeds reject_trimlength %d, so returning NULL\n",
		  nonterminal_trim,reject_trimlength));
    Substring_free(&new);
    return (T) NULL;
  } else {
    new->querystart += new->trim_left;
    new->queryend -= new->trim_right;
    /* Check for minlength.  Needed to avoid nonsensical terminal alignments */
    if (new->queryend - new->querystart <= minlength) {
      debug8(printf("queryend %d - querystart %d <= minlength %d, so returning NULL\n",
		    new->queryend,new->querystart,minlength));
      Substring_free(&new);
      return (T) NULL;
    }
  }

  /* ? Should we spend the time to determine trim_left_splicep and
     trim_right_splicep, especially since trimming may not be perfect */
  if (plusp == true) {
    /* nmatches: Counts matches over whole region including trims */
    new->nmatches = (new->alignend - new->alignstart) - new->nmismatches_whole;

    new->alignstart_trim += new->trim_left;
    new->alignend_trim -= new->trim_right;
    
    if (novelsplicingp == false) {
      new->trim_left_splicep = new->trim_right_splicep = false;
    } else {
      if (new->trim_left == 0) {
	new->trim_left_splicep = false;
      } else {
	prob1 = Maxent_hr_acceptor_prob(left + new->trim_left,chroffset);
	prob2 = Maxent_hr_antidonor_prob(left + new->trim_left,chroffset);
	/* fprintf(stderr,"At %u, acceptor prob %f, antidonor prob %f\n",left+new->trim_left,prob1,prob2); */
	if (prob1 > 0.90 || prob2 > 0.90) {
	  new->trim_left_splicep = true;
	} else {
	  new->trim_left_splicep = false;
	}
      }

      if (new->trim_right == 0) {
	new->trim_right_splicep = false;
      } else {
	prob1 = Maxent_hr_donor_prob(left + querylength - new->trim_right,chroffset);
	prob2 = Maxent_hr_antiacceptor_prob(left + querylength - new->trim_right,chroffset);
	/* fprintf(stderr,"At %u, donor prob %f, antiacceptor prob %f\n",left + querylength - new->trim_right,prob1,prob2); */
	if (prob1 > 0.90 || prob2 > 0.90) {
	  new->trim_right_splicep = true;
	} else {
	  new->trim_right_splicep = false;
	}
      }
    }

  } else {
    /* nmatches: Counts matches over whole region including trims */
    new->nmatches = (new->alignstart - new->alignend) - new->nmismatches_whole;

    new->alignstart_trim -= new->trim_left;
    new->alignend_trim += new->trim_right;
    
    if (novelsplicingp == false) {
      new->trim_left_splicep = new->trim_right_splicep = false;
    } else {
      if (new->trim_left == 0) {
	new->trim_left_splicep = false;
      } else {
	prob1 = Maxent_hr_donor_prob(left + querylength - new->trim_left,chroffset);
	prob2 = Maxent_hr_antiacceptor_prob(left + querylength - new->trim_left,chroffset);
	/* fprintf(stderr,"At %u, donor prob %f, antiacceptor prob %f\n",left + querylength - new->trim_left,prob1,prob2); */
	if (prob1 > 0.90 || prob2 > 0.90) {
	  new->trim_left_splicep = true;
	} else {
	  new->trim_left_splicep = false;
	}
      }

      if (new->trim_right == 0) {
	new->trim_right_splicep = false;
      } else {
	prob1 = Maxent_hr_acceptor_prob(left + new->trim_right,chroffset);
	prob2 = Maxent_hr_antidonor_prob(left + new->trim_right,chroffset);
	/* fprintf(stderr,"At %u, acceptor prob %f, antidonor prob %f\n",left+new->trim_right,prob1,prob2); */
	if (prob1 > 0.90 || prob2 > 0.90) {
	  new->trim_right_splicep = true;
	} else {
	  new->trim_right_splicep = false;
	}
      }
    }
  }


  /* nmismatches_bothdiff: Counts matches of trimmed region */
  if (new->trim_left == 0 && new->trim_right == 0) {
    new->nmismatches_bothdiff = nmismatches_whole;

  } else if (plusp == true) {
    new->nmismatches_bothdiff = 
      Genome_count_mismatches_substring(query_compress,left,/*pos5*/new->alignstart_trim-left,
					/*pos3*/new->alignend_trim-left,plusp,genestrand,first_read_p);
  } else {
    new->nmismatches_bothdiff = 
      Genome_count_mismatches_substring(query_compress,left,/*pos5*/new->alignend_trim-left,
					/*pos3*/new->alignstart_trim-left,plusp,genestrand,first_read_p);
  }

  /* General test for goodness over trimmed region */
  nmatches = new->queryend - new->querystart - new->nmismatches_bothdiff;
  if (nmatches - 3*new->nmismatches_bothdiff < 0) {
    debug2(printf("Substring fails general test for goodness with %d matches and %d mismatches\n",
		  nmatches,new->nmismatches_bothdiff));
    Substring_free(&new);
    return (T) NULL;
  }

  return new;
}


/* Look also at Pair_compute_mapq */
float
Substring_compute_mapq (T this, Compress_T query_compress, char *quality_string, bool trim_terminals_p) {
  int mapq_start, mapq_end;

  /* mapq */
  mapq_start = this->querystart_orig;
  mapq_end = this->queryend_orig;

  /* It appears from simulated reads that it is better not to trim in
     computing MAPQ.  The correct mapping then tends to be selected
     with a higher MAPQ score. */
  /* But if all ends are terminals, then terminal parts should not be
     included in MAPQ scoring */
  if (trim_terminals_p == true) {
    if (this->start_endtype == TERM) {
      mapq_start += this->trim_left;
    }
    if (this->end_endtype == TERM) {
      mapq_end -= this->trim_right;
    }
  }

  if (this->exactp == true) {
    /* this->mapq_loglik = MAPQ_loglik_exact(quality_string,0,querylength); */
    this->mapq_loglik = 0.0;
  } else {
    debug2(printf("trim_left %d, trim_right %d, mapq_start = %d, mapq_end = %d\n",
		  this->trim_left,this->trim_right,mapq_start,mapq_end));
    this->mapq_loglik = MAPQ_loglik(query_compress,this->left,mapq_start,mapq_end,
				    this->querylength,quality_string,this->plusp,this->genestrand,this->first_read_p);
  }

  return this->mapq_loglik;
}


/* Note: query needed only for dibase */
int
Substring_display_prep (char **deletion, T this, char *query, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Genome_T genome, int deletion_pos, int deletion_length) {
  char *genomic_diff;
  unsigned char *gbuffer;
  int mismatch_offset;
#ifndef HAVE_ALLOCA
  unsigned char gbuffer_alloc[MAX_READLENGTH/*+MAX_END_DELETIONS*/+1];
  bool allocp;
#endif

  mismatch_offset = this->alignoffset - this->querystart_orig;

  /* genomic_bothdiff, genomic_refdiff, and nmismatches_refdiff */
  if (this->exactp == true) {
    this->genomic_bothdiff = (char *) NULL;
    this->genomic_refdiff = (char *) NULL;
    this->nmismatches_refdiff = this->nmismatches_whole;

  } else if (this->plusp == true) {
#ifdef HAVE_ALLOCA
    gbuffer = (unsigned char *) ALLOCA((this->genomiclength+1) * sizeof(unsigned char));
#else
    if (this->genomiclength < MAX_READLENGTH) {
      gbuffer = gbuffer_alloc;
      allocp = false;
    } else {
      gbuffer = (unsigned char *) MALLOC((this->genomiclength+1) * sizeof(unsigned char));
      allocp = true;
    }
#endif

    Genome_fill_buffer_simple(genome,this->left_genomicseg,this->genomiclength,gbuffer);
    if (deletion_pos >= 0) {
      *deletion = (char *) CALLOC_OUT(deletion_length+1,sizeof(char));
      strncpy(*deletion,&(gbuffer[deletion_pos]),deletion_length);
    }
    genomic_diff = gbuffer;

    Genome_mark_mismatches(genomic_diff,this->querylength,query_compress_fwd,
			   this->left,/*pos5*/this->querystart_orig,/*pos3*/this->queryend_orig,
			   mismatch_offset,/*plusp*/true,this->genestrand,this->first_read_p);

    this->genomic_bothdiff = embellish_genomic(genomic_diff,query,this->querystart_orig,this->queryend_orig,
					       this->querylength,this->alignoffset,this->extraleft,this->extraright,
					       this->genestrand);

    if (snps_iit == NULL) {
      this->genomic_refdiff = this->genomic_bothdiff;
      this->nmismatches_refdiff = this->nmismatches_bothdiff;

    } else {
      this->nmismatches_refdiff =
	Genome_count_mismatches_substring_ref(query_compress_fwd,this->left,
					      /*pos5*/this->alignstart_trim - this->left,
					      /*pos3*/this->alignend_trim - this->left,
					      /*plusp*/true,this->genestrand,this->first_read_p);

      Genome_mark_mismatches_ref(genomic_diff,this->querylength,query_compress_fwd,this->left,
				 /*pos5*/this->querystart_orig,/*pos3*/this->queryend_orig,
				 mismatch_offset,/*plusp*/true,this->genestrand,this->first_read_p);
      if (output_sam_p == false) {
	this->genomic_refdiff = embellish_genomic(genomic_diff,query,this->querystart_orig,this->queryend_orig,
						  this->querylength,this->alignoffset,this->extraleft,this->extraright,
						  this->genestrand);
      }
    }

    if (output_sam_p == true) {
      this->genomic_refdiff = embellish_genomic_sam(genomic_diff,query,this->querystart_orig,this->queryend_orig,
						    this->querylength,this->genomiclength,this->alignoffset,
						    this->genestrand);
    }

#ifdef HAVE_ALLOCA
    FREEA(gbuffer);
#else
    if (allocp == true) {
      FREE(gbuffer);
    }
#endif

  } else {
#ifdef HAVE_ALLOCA
    gbuffer = (unsigned char *) ALLOCA((this->genomiclength+1) * sizeof(unsigned char));
#else
    if (this->genomiclength < MAX_READLENGTH) {
      gbuffer = gbuffer_alloc;
      allocp = false;
    } else {
      gbuffer = (unsigned char *) MALLOC((this->genomiclength+1) * sizeof(unsigned char));
      allocp = true;
    }
#endif

    Genome_fill_buffer_simple(genome,this->left_genomicseg,this->genomiclength,gbuffer);
    genomic_diff = make_complement_inplace(gbuffer,this->genomiclength);
    if (deletion_pos >= 0) {
      *deletion = (char *) CALLOC_OUT(deletion_length+1,sizeof(char));
      strncpy(*deletion,&(gbuffer[deletion_pos]),deletion_length);
    }

    Genome_mark_mismatches(genomic_diff,this->querylength,query_compress_rev,
			   this->left,/*pos5*/this->querylength - this->queryend_orig,
			   /*pos3*/this->querylength - this->querystart_orig,
			   mismatch_offset,/*plusp*/false,this->genestrand,this->first_read_p);

    this->genomic_bothdiff = embellish_genomic(genomic_diff,query,this->querystart_orig,this->queryend_orig,
					       this->querylength,this->alignoffset,this->extraleft,this->extraright,
					       this->genestrand);

    if (snps_iit == NULL) {
      this->genomic_refdiff = this->genomic_bothdiff;
      this->nmismatches_refdiff = this->nmismatches_bothdiff;

    } else {
      this->nmismatches_refdiff = 
	Genome_count_mismatches_substring_ref(query_compress_rev,this->left,
					      /*pos5*/this->alignend_trim - this->left,
					      /*pos3*/this->alignstart_trim - this->left,/*plusp*/false,
					      this->genestrand,this->first_read_p);

      Genome_mark_mismatches_ref(genomic_diff,this->querylength,query_compress_rev,this->left,
				 /*pos5*/this->querylength - this->queryend_orig,
				 /*pos3*/this->querylength - this->querystart_orig,
				 mismatch_offset,/*plusp*/false,this->genestrand,this->first_read_p);

      if (output_sam_p == false) {
	this->genomic_refdiff = embellish_genomic(genomic_diff,query,this->querystart_orig,this->queryend_orig,
						  this->querylength,this->alignoffset,this->extraleft,this->extraright,
						  this->genestrand);
      }
    }

    if (output_sam_p == true) {
      this->genomic_refdiff = embellish_genomic_sam(genomic_diff,query,this->querystart_orig,this->queryend_orig,
						    this->querylength,this->genomiclength,this->alignoffset,
						    this->genestrand);
    }

#ifdef HAVE_ALLOCA
    FREEA(gbuffer);
#else
    if (allocp == true) {
      FREE(gbuffer);
    }
#endif
  }

  return this->nmismatches_refdiff;
}


Univcoord_T
Substring_splicecoord (T this) {
  return this->splicecoord;
}

int
Substring_splicesites_knowni (T this) {
  return this->splicesites_knowni;
}

Univcoord_T
Substring_splicecoord_A (T this) {
  return this->splicecoord;
}

Univcoord_T
Substring_splicecoord_D (T this) {
  return this->splicecoord_2;
}

bool
Substring_plusp (T this) {
  return this->plusp;
}

int
Substring_genestrand (T this) {
  return this->genestrand;
}

bool
Substring_first_read_p (T this) {
  return this->first_read_p;
}

char *
Substring_genomic_bothdiff (T this) {
  return this->genomic_bothdiff;
}

char *
Substring_genomic_refdiff (T this) {
  return this->genomic_refdiff;
}


int
Substring_nmismatches_whole (T this) {
  return this->nmismatches_whole;
}


int
Substring_nmismatches_bothdiff (T this) {
  return this->nmismatches_bothdiff;
}


int
Substring_nmismatches_refdiff (T this) {
  return this->nmismatches_refdiff;
}

int
Substring_nmatches (T this) {
  return this->nmatches;
}

int
Substring_nmatches_posttrim (T this) {
  return this->queryend - this->querystart - this->nmismatches_bothdiff;
}


void
Substring_set_nmismatches_terminal (T this, int nmismatches_whole, int nmismatches_bothdiff) {
  this->nmismatches_whole = nmismatches_whole;
  this->nmismatches_bothdiff = nmismatches_bothdiff;
  if (this->plusp == true) {
    this->nmatches = (this->alignend_trim - this->alignstart_trim) - nmismatches_whole;
  } else {
    this->nmatches = (this->alignstart_trim - this->alignend_trim) - nmismatches_whole;
  }
  return;
}

Endtype_T
Substring_start_endtype (T this) {
  return this->start_endtype;
}

Endtype_T
Substring_end_endtype (T this) {
  return this->end_endtype;
}

void
Substring_set_endtypes (T this, Endtype_T start_endtype, Endtype_T end_endtype) {
  this->start_endtype = start_endtype;
  this->end_endtype = end_endtype;
  return;
}



float
Substring_mapq_loglik (T this) {
  return this->mapq_loglik;
}

int
Substring_trim_left (T this) {
  return this->trim_left;
}

int
Substring_trim_right (T this) {
  return this->trim_right;
}

bool
Substring_trim_left_splicep (T this) {
  return this->trim_left_splicep;
}

bool
Substring_trim_right_splicep (T this) {
  return this->trim_right_splicep;
}


int
Substring_querystart (T this) {
  return this->querystart;
}

int
Substring_querystart_orig (T this) {
  return this->querystart_orig;
}

int
Substring_queryend (T this) {
  return this->queryend;
}

int
Substring_queryend_orig (T this) {
  return this->queryend_orig;
}

int
Substring_querylength (T this) {
  return this->querylength;
}

int
Substring_match_length (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return this->queryend - this->querystart;
  }
}

/* Before trimming */
int
Substring_match_length_orig (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return this->queryend_orig - this->querystart_orig;
  }
}

/* Used only by Goby */
Chrpos_T
Substring_genomic_alignment_length (T this) {
  if (this == NULL) {
    return 0U;
  } else if (this->alignend_trim > this->alignstart_trim) {
    return this->alignend_trim - this->alignstart_trim;
  } else {
    return this->alignstart_trim - this->alignend_trim;
  }
}


Chrnum_T
Substring_chrnum (T this) {
  return this->chrnum;
}

Univcoord_T
Substring_chroffset (T this) {
  return this->chroffset;
}

Univcoord_T
Substring_chrhigh (T this) {
  return this->chrhigh;
}

Chrpos_T
Substring_chrlength (T this) {
  return this->chrlength;
}

Univcoord_T
Substring_alignstart (T this) {
  return this->alignstart;
}

Univcoord_T
Substring_alignend (T this) {
  return this->alignend;
}

Univcoord_T
Substring_alignstart_trim (T this) {
  return this->alignstart_trim;
}

Univcoord_T
Substring_alignend_trim (T this) {
  return this->alignend_trim;
}


#if 0
Univcoord_T
Substring_alignmid_trim (T this) {
  if (this->alignend_trim > this->alignstart_trim) {
    return this->alignstart_trim + (this->alignend_trim - this->alignstart_trim)/2;
  } else {
    return this->alignend_trim + (this->alignstart_trim - this->alignend_trim)/2;
  }
}
#endif


Univcoord_T
Substring_left_genomicseg (T this) {
  return this->left_genomicseg;
}

Univcoord_T
Substring_genomicstart (T this) {
  return this->genomicstart;
}

Univcoord_T
Substring_genomicend (T this) {
  return this->genomicend;
}

Chrpos_T
Substring_genomiclength (T this) {
  return this->genomiclength;
}

Chrpos_T
Substring_chrstart (T this) {
  return this->genomicstart - this->chroffset;
}

Chrpos_T
Substring_chrend (T this) {
  return this->genomicend - this->chroffset;
}


double
Substring_chimera_prob (T this) {
  return this->chimera_prob;
}

double
Substring_chimera_prob_2 (T this) {
  return this->chimera_prob_2;
}

int
Substring_chimera_pos (T this) {
  return this->chimera_pos;
}

/* For shortexon */
int
Substring_chimera_pos_A (T this) {
  return this->chimera_pos;
}

/* For shortexon */
int
Substring_chimera_pos_D (T this) {
  return this->chimera_pos_2;
}

bool
Substring_chimera_knownp (T this) {
  return this->chimera_knownp;
}

int
Substring_nchimera_known (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return (int) this->chimera_knownp + (int) this->chimera_knownp_2;
  }
}

int
Substring_nchimera_novel (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return (int) this->chimera_novelp + (int) this->chimera_novelp_2;
  }
}


int
Substring_chimera_sensedir (T this) {
  if (this->chimera_sensep == true) {
    return SENSE_FORWARD;
  } else {
    return SENSE_ANTI;
  }
}


bool
Substring_chimera_sensep (T this) {
  return this->chimera_sensep;
}

/* circularpos measures query distance from SAM chrlow to origin */
int
Substring_circularpos (T this) {
  if (this == NULL) {
    return -1;

  } else if (this->plusp == true) {
    /* printf("substring plus: looking at %u (trim %u) vs %u+%u\n",
       this->alignend,this->alignend_trim,this->chroffset,this->chrlength); */
    if (this->alignend_trim > this->chroffset + this->chrlength) {
      /* return (this->querystart - this->trim_left) + (this->chroffset + this->chrlength) - this->alignstart; */
      return this->querystart + (this->chroffset + this->chrlength) - this->alignstart_trim;
    } else {
      return -1;
    }

  } else {
    /* printf("substring minus: looking at %u (trim %u) vs %u+%u\n",
       this->alignstart,this->alignstart_trim,this->chroffset,this->chrlength); */
    if (this->alignstart_trim > this->chroffset + this->chrlength) {
      /* return ((this->querylength - this->trim_right) - this->queryend) + (this->chroffset + this->chrlength) - this->alignend; */
      return (this->querylength - this->queryend) + (this->chroffset + this->chrlength) - this->alignend_trim;
    } else {
      return -1;
    }
  }
}


T
Substring_copy (T old) {
  T new;

  if (old == NULL) {
    return NULL;
  } else {
    new = (T) MALLOC_OUT(sizeof(*new));

    new->exactp = old->exactp;
    new->nmismatches_whole = old->nmismatches_whole;
    new->nmismatches_bothdiff = old->nmismatches_bothdiff;
    new->nmismatches_refdiff = old->nmismatches_refdiff;
    new->nmatches = old->nmatches;

    new->trim_left = old->trim_left;
    new->trim_right = old->trim_right;
    new->trim_left_splicep = old->trim_left_splicep;
    new->trim_right_splicep = old->trim_right_splicep;

    new->chrnum = old->chrnum;
    new->chroffset = old->chroffset;
    new->chrhigh = old->chrhigh;
    new->chrlength = old->chrlength;

    new->left_genomicseg = old->left_genomicseg;
    new->left = old->left;
    new->genomicstart = old->genomicstart;
    new->genomicend = old->genomicend;

    new->start_endtype = old->start_endtype;
    new->end_endtype = old->end_endtype;

    new->querystart_orig = old->querystart_orig;
    new->queryend_orig = old->queryend_orig;
    new->querystart = old->querystart;
    new->queryend = old->queryend;
    new->querylength = old->querylength;

    new->alignstart = old->alignstart;
    new->alignend = old->alignend;

    new->alignstart_trim = old->alignstart_trim;
    new->alignend_trim = old->alignend_trim;

    new->alignoffset = old->alignoffset;
    new->extraleft = old->extraleft;
    new->extraright = old->extraright;

    new->genomiclength = old->genomiclength;
    new->plusp = old->plusp;
    new->genestrand = old->genestrand;
    new->first_read_p = old->first_read_p;

    if (old->genomic_bothdiff == NULL) {
      new->genomic_bothdiff = (char *) NULL;
      new->genomic_refdiff = (char *) NULL;
    } else if (old->genomic_refdiff == old->genomic_bothdiff) {
      new->genomic_bothdiff = (char *) CALLOC_OUT(strlen(old->genomic_bothdiff)+1,sizeof(char));
      strcpy(new->genomic_bothdiff,old->genomic_bothdiff);
      new->genomic_refdiff = new->genomic_bothdiff;
    } else {
      new->genomic_bothdiff = (char *) CALLOC_OUT(strlen(old->genomic_bothdiff)+1,sizeof(char));
      strcpy(new->genomic_bothdiff,old->genomic_bothdiff);
      new->genomic_refdiff = (char *) CALLOC_OUT(strlen(old->genomic_refdiff)+1,sizeof(char));
      strcpy(new->genomic_refdiff,old->genomic_refdiff);
    }

    new->mapq_loglik = old->mapq_loglik;

    new->chimera_sensep = old->chimera_sensep;

    new->splicecoord = old->splicecoord;
    new->splicesites_knowni = old->splicesites_knowni;
    new->chimera_knownp = old->chimera_knownp;
    new->chimera_novelp = old->chimera_novelp;
    new->chimera_modelpos = old->chimera_modelpos;
    new->chimera_pos = old->chimera_pos;
    new->chimera_prob = old->chimera_prob;

    new->splicecoord_2 = old->splicecoord_2;
    new->splicesites_knowni_2 = old->splicesites_knowni_2;
    new->chimera_knownp_2 = old->chimera_knownp_2;
    new->chimera_novelp_2 = old->chimera_novelp_2;
    new->chimera_modelpos_2 = old->chimera_modelpos_2;
    new->chimera_pos_2 = old->chimera_pos_2;
    new->chimera_prob_2 = old->chimera_prob_2;

    return new;
  }
}



T
Substring_new_donor (Univcoord_T donor_coord, int donor_knowni, int donor_pos, int donor_nmismatches,
		     double donor_prob, Univcoord_T left, Compress_T query_compress,
		     int querylength, bool plusp, int genestrand, bool first_read_p, bool sensep,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength) {
  T new;
  int querystart, queryend, extraleft, extraright;
  Univcoord_T genomicstart, genomicend, alignstart, alignend;
  Endtype_T start_endtype, end_endtype;
  bool trim_left_p, trim_right_p;

  /* Previously checked if left >= chroffset + chrlength to exclude
     the duplicate length, but now excluding all translocations to
     circular chromosomes */

  if (chroffset + chrlength < chrhigh) {
    /* Don't splice to circular chromosomes */
    return (T) NULL;

  } else if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
    if (sensep == true) {
      start_endtype = END;
      end_endtype = DON;

      querystart = 0;
      queryend = donor_pos;
      extraleft = 0;
      extraright = 2;
      alignstart = genomicstart;
      alignend = genomicstart + donor_pos;
      trim_left_p = true;	/* querystart == 0 */
      trim_right_p = false;

    } else {
      extraleft = 2;
      extraright = 0;
      start_endtype = DON;
      end_endtype = END;

      querystart = donor_pos;
      queryend = querylength;
      alignstart = genomicstart + donor_pos;
      alignend = genomicend;
      trim_left_p = false;
      trim_right_p = true;	/* queryend == querylength */
    }

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    if (sensep == true) {
      start_endtype = END;
      end_endtype = DON;

      querystart = 0;
      queryend = querylength - donor_pos;
      extraleft = 0;
      extraright = 2;
      alignstart = genomicstart;
      alignend = genomicstart - (querylength - donor_pos);
      trim_left_p = true;	/* querystart == 0 */
      trim_right_p = false;

    } else {
      extraleft = 2;
      extraright = 0;
      start_endtype = DON;
      end_endtype = END;

      querystart = querylength - donor_pos;
      queryend = querylength;
      alignstart = genomicstart - (querylength - donor_pos);
      alignend = genomicend;
      trim_left_p = false;
      trim_right_p = true;	/* queryend == querylength */
    }
  }

  if ((new = Substring_new(donor_nmismatches,chrnum,chroffset,chrhigh,chrlength,
			   left,genomicstart,genomicend,query_compress,
			   start_endtype,end_endtype,querystart,queryend,querylength,
			   alignstart,alignend,/*genomiclength*/querylength,
			   extraleft,extraright,/*exactp*/false,plusp,genestrand,first_read_p,
			   trim_left_p,trim_right_p,/*minlength*/0)) == NULL) {
    return (T) NULL;
  }

  debug2(printf("Making new donor with splicesites_i %d, coord %u and left %u, plusp %d, sensep %d, query %d..%d, genome %u..%u\n",
		donor_knowni,donor_coord,left,plusp,sensep,querystart,queryend,alignstart - chroffset,alignend - chroffset));
  new->splicecoord = donor_coord;
  new->splicesites_knowni = donor_knowni;

  new->chimera_modelpos = left + donor_pos;
  assert(new->splicecoord == new->chimera_modelpos);
  new->chimera_sensep = sensep;
  if (donor_knowni >= 0) {
    new->chimera_knownp = true;
    /* new->chimera_novelp = false */
  } else {
    /* new->chimera_knownp = false; */
    new->chimera_novelp = true;
  }

  if (plusp == true) {
    new->chimera_pos = donor_pos;
  } else {
    new->chimera_pos = querylength - donor_pos;
  }
  new->chimera_prob = donor_prob;

  return new;
}


T
Substring_new_acceptor (Univcoord_T acceptor_coord, int acceptor_knowni, int acceptor_pos, int acceptor_nmismatches,
			double acceptor_prob, Univcoord_T left, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, bool first_read_p, bool sensep,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength) {
  T new;
  int querystart, queryend, extraleft, extraright;
  Univcoord_T genomicstart, genomicend, alignstart, alignend;
  Endtype_T start_endtype, end_endtype;
  bool trim_left_p, trim_right_p;

  /* Previously checked if left >= chroffset + chrlength to exclude
     the duplicate length, but now excluding all translocations to
     circular chromosomes */

  if (chroffset + chrlength < chrhigh) {
    /* Don't splice to circular chromosomes */
    return (T) NULL;

  } else if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
    if (sensep == true) {
      extraleft = 2;
      extraright = 0;
      start_endtype = ACC;
      end_endtype = END;

      querystart = acceptor_pos;
      queryend = querylength;
      alignstart = genomicstart + acceptor_pos;
      alignend = genomicend;
      trim_left_p = false;
      trim_right_p = true;	/* queryend == querylength */

    } else {
      extraleft = 0;
      extraright = 2;
      start_endtype = END;
      end_endtype = ACC;

      querystart = 0;
      queryend = acceptor_pos;
      alignstart = genomicstart;
      alignend = genomicstart + acceptor_pos;
      trim_left_p = true;	/* querystart == 0 */
      trim_right_p = false;
    }

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    if (sensep == true) {
      extraleft = 2;
      extraright = 0;
      start_endtype = ACC;
      end_endtype = END;

      querystart = querylength - acceptor_pos;
      queryend = querylength;
      alignstart = genomicstart - (querylength - acceptor_pos);
      alignend = genomicend;
      trim_left_p = false;
      trim_right_p = true;	/* queryend == querylength */

    } else {
      extraleft = 0;
      extraright = 2;
      start_endtype = END;
      end_endtype = ACC;

      querystart = 0;
      queryend = querylength - acceptor_pos;
      alignstart = genomicstart;
      alignend = genomicstart - (querylength - acceptor_pos);
      trim_left_p = true;	/* querystart == 0 */
      trim_right_p = false;
    }
  }

  if ((new = Substring_new(acceptor_nmismatches,chrnum,chroffset,chrhigh,chrlength,
			   left,genomicstart,genomicend,query_compress,
			   start_endtype,end_endtype,querystart,queryend,querylength,
			   alignstart,alignend,/*genomiclength*/querylength,
			   extraleft,extraright,/*exactp*/false,plusp,genestrand,first_read_p,
			   trim_left_p,trim_right_p,/*minlength*/0)) == NULL) {
    return (T) NULL;
  }

  debug2(printf("Making new acceptor with splicesites_i %d, coord %u and left %u, plusp %d, sensep %d, query %d..%d, genome %u..%u\n",
		acceptor_knowni,acceptor_coord,left,plusp,sensep,querystart,queryend,alignstart - chroffset,alignend - chroffset));
  new->splicecoord = acceptor_coord;
  new->splicesites_knowni = acceptor_knowni;

  new->chimera_modelpos = left + acceptor_pos;
  assert(new->splicecoord == new->chimera_modelpos);
  new->chimera_sensep = sensep;
  if (acceptor_knowni >= 0) {
    new->chimera_knownp = true;
    /* new->chimera_novelp = false */
  } else {
    /* new->chimera_knownp = false; */
    new->chimera_novelp = true;
  }

  if (plusp == true) {
    new->chimera_pos = acceptor_pos;
  } else {
    new->chimera_pos = querylength - acceptor_pos;
  }
  new->chimera_prob = acceptor_prob;

  return new;
}



T
Substring_new_shortexon (Univcoord_T acceptor_coord, int acceptor_knowni, Univcoord_T donor_coord, int donor_knowni,
			 int acceptor_pos, int donor_pos, int nmismatches,
			 double acceptor_prob, double donor_prob, Univcoord_T left,
			 Compress_T query_compress, int querylength,
			 bool plusp, int genestrand, bool first_read_p, bool sensep,
			 bool acceptor_ambp, bool donor_ambp,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength) {
  T new;
  int querystart, queryend;
  Univcoord_T genomicstart, genomicend, alignstart, alignend;
  Endtype_T start_endtype, end_endtype;

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
    if (sensep == true) {
      start_endtype = (acceptor_ambp == true) ? AMB_ACC : ACC;
      end_endtype = (donor_ambp == true) ? AMB_DON : DON;
      querystart = acceptor_pos;
      queryend = donor_pos;
      alignstart = genomicstart + acceptor_pos;
      alignend = genomicstart + donor_pos;
    } else {
      start_endtype = (donor_ambp == true) ? AMB_DON : DON;
      end_endtype = (acceptor_ambp == true) ? AMB_ACC : ACC;
      querystart = donor_pos;
      queryend = acceptor_pos;
      alignstart = genomicstart + donor_pos;
      alignend = genomicstart + acceptor_pos;
    }

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    if (sensep == true) {
      start_endtype = (acceptor_ambp == true) ? AMB_ACC : ACC;
      end_endtype = (donor_ambp == true) ? AMB_DON : DON;
      querystart = querylength - acceptor_pos;
      queryend = querylength - donor_pos;
      alignstart = genomicstart - (querylength - acceptor_pos);
      alignend = genomicstart - (querylength - donor_pos);
    } else {
      start_endtype = (donor_ambp == true) ? AMB_DON : DON;
      end_endtype = (acceptor_ambp == true) ? AMB_ACC : ACC;
      querystart = querylength - donor_pos;
      queryend = querylength - acceptor_pos;
      alignstart = genomicstart - (querylength - donor_pos);
      alignend = genomicstart - (querylength - acceptor_pos);
    }
  }

  if ((new = Substring_new(nmismatches,chrnum,chroffset,chrhigh,chrlength,
			   left,genomicstart,genomicend,query_compress,
			   start_endtype,end_endtype,querystart,queryend,querylength,
			   alignstart,alignend,/*genomiclength*/querylength,
			   /*extraleft*/2,/*extraright*/2,/*exactp*/false,plusp,genestrand,first_read_p,
			   /*trim_left_p*/false,/*trim_right_p*/false,
			   /*minlength*/0)) == NULL) {
    return (T) NULL;
  }

  debug2(printf("Making new middle with left %u, plusp %d\n",left,plusp));
  new->splicecoord = acceptor_coord;
  new->splicesites_knowni = acceptor_knowni;
  new->splicecoord_2 = donor_coord;
  new->splicesites_knowni_2 = donor_knowni;

  new->chimera_modelpos = left + acceptor_pos;
  new->chimera_modelpos_2 = left + donor_pos;
  new->chimera_sensep = sensep;

  if (acceptor_knowni >= 0) {
    new->chimera_knownp = true;
    /* new->chimera_novelp = false; */
  } else {
    /* new->chimera_knownp = false; */
    new->chimera_novelp = true;
  }

  if (donor_knowni >= 0) {
    new->chimera_knownp_2 = true;
    /* new->chimera_novelp_2 = false; */
  } else {
    /* new->chimera_knownp_2 = false; */
    new->chimera_novelp_2 = true;
  }

  if (plusp == true) {
    new->chimera_pos = acceptor_pos;
    new->chimera_pos_2 = donor_pos;
  } else {
    new->chimera_pos = querylength - acceptor_pos;
    new->chimera_pos_2 = querylength - donor_pos;
  }

  new->chimera_prob = acceptor_prob;
  new->chimera_prob_2 = donor_prob;

  return new;
}



void
Substring_assign_donor_prob (T donor) {

  if (donor == NULL) {
    return;
  }

  if (donor->chimera_knownp == false) {
    /* Prob already assigned */
  } else if (donor->plusp == donor->chimera_sensep) {
    donor->chimera_prob = Maxent_hr_donor_prob(donor->chimera_modelpos,donor->chroffset);
  } else {
    donor->chimera_prob = Maxent_hr_antidonor_prob(donor->chimera_modelpos,donor->chroffset);
  }

  return;
}

void
Substring_assign_acceptor_prob (T acceptor) {

  if (acceptor == NULL) {
    return;
  }

  if (acceptor->chimera_knownp == false) {
    /* Prob already assigned */
  } else if (acceptor->plusp == acceptor->chimera_sensep) {
    acceptor->chimera_prob = Maxent_hr_acceptor_prob(acceptor->chimera_modelpos,acceptor->chroffset);
  } else {
    acceptor->chimera_prob = Maxent_hr_antiacceptor_prob(acceptor->chimera_modelpos,acceptor->chroffset);
  }

  return;
}


void
Substring_assign_shortexon_prob (T shortexon) {

  if (shortexon->chimera_knownp == false) {
    /* Prob1 already assigned */
  } else if (shortexon->plusp == shortexon->chimera_sensep) {
    shortexon->chimera_prob = Maxent_hr_acceptor_prob(shortexon->chimera_modelpos,shortexon->chroffset);
  } else {
    shortexon->chimera_prob = Maxent_hr_antiacceptor_prob(shortexon->chimera_modelpos,shortexon->chroffset);
  }


  if (shortexon->chimera_knownp_2 == false) {
    /* Prob2 already assigned */
  } else if (shortexon->plusp == shortexon->chimera_sensep) {
    shortexon->chimera_prob_2 = Maxent_hr_donor_prob(shortexon->chimera_modelpos_2,shortexon->chroffset);
  } else {
    shortexon->chimera_prob_2 = Maxent_hr_antidonor_prob(shortexon->chimera_modelpos_2,shortexon->chroffset);
  }

  return;
}



static int
ascending_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->chimera_pos < y->chimera_pos) {
    return -1;
  } else if (x->chimera_pos > y->chimera_pos) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else if (x->chimera_knownp == true && y->chimera_knownp == false) {
    return -1;
  } else if (y->chimera_knownp == true && x->chimera_knownp == false) {
    return +1;
  } else {
    return 0;
  }
}

static int
descending_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->chimera_pos < y->chimera_pos) {
    return -1;
  } else if (x->chimera_pos > y->chimera_pos) {
    return +1;
  } else if (x->genomicstart > y->genomicstart) {
    return -1;
  } else if (x->genomicstart < y->genomicstart) {
    return +1;
  } else if (x->chimera_knownp == true && y->chimera_knownp == false) {
    return -1;
  } else if (y->chimera_knownp == true && x->chimera_knownp == false) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Substring_sort_chimera_halves (List_T hitlist, bool ascendingp) {
  List_T sorted = NULL, p;
  T x, *hits;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitlist);
  debug(printf("Checking %d spliceends for duplicates...",n));
  if (n == 0) {
    debug(printf("\n"));
    return NULL;
  }

  hits = (T *) MALLOCA(n * sizeof(T));
  List_fill_array_and_free((void **) hits,&hitlist);

  if (ascendingp == true) {
    qsort(hits,n,sizeof(T),ascending_pos_cmp);
  } else {
    qsort(hits,n,sizeof(T),descending_pos_cmp);
  }

  /* Check for duplicates */
  eliminate = (bool *) CALLOCA(n,sizeof(bool));
  for (i = 0; i < n; i++) {
    x = hits[i];
    j = i+1;
    while (j < n && hits[j]->chimera_pos == x->chimera_pos && hits[j]->genomicstart == x->genomicstart) {
      eliminate[j] = true;
      j++;
    }
  }

  debug(j = 0);
  for (i = n-1; i >= 0; i--) {
    x = hits[i];
    if (eliminate[i] == false) {
      sorted = List_push(sorted,x);
    } else {
      Substring_free(&x);
      debug(j++);
    }
  }
  debug(printf("%d eliminated\n",j));

  FREEA(hits);
  FREEA(eliminate);

  return sorted;
}


static void
print_snp_labels (FILE *fp, T this, Shortread_T queryseq) {
  int *snps, nsnps, querypos, i;
  char *label, *seq1, *seq2;
  bool allocp, printp = false;
  Interval_T interval;
  Chrpos_T position;

  if (this->plusp == true) {
    snps = IIT_get_with_divno(&nsnps,snps_iit,
			      snps_divint_crosstable[this->chrnum],
			      this->alignstart - this->chroffset + 1,this->alignend - this->chroffset,
			      /*sortp*/false);
  } else {
    snps = IIT_get_with_divno(&nsnps,snps_iit,
			      snps_divint_crosstable[this->chrnum],
			      this->alignend - this->chroffset + 1,this->alignstart - this->chroffset,
			      /*sortp*/false);
  }

  fprintf(fp,",snps:");

  seq1 = Shortread_fullpointer_uc(queryseq);
  if (this->genomic_bothdiff == NULL) {
    seq2 = Shortread_fullpointer(queryseq);
  } else {
    seq2 = this->genomic_bothdiff;
  }

  if (this->plusp) {

#if 0
    for (i = 0; i < nsnps; i++) {
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = position - (this->genomicstart-this->chroffset) - 1;
      printf("%d ",querypos);
    }
    printf("\n");
#endif

    for (i = 0; i < nsnps; i++) {
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = position - (this->genomicstart-this->chroffset) - 1;
      assert(querypos >= 0 && querypos < this->genomiclength);

#if 0
      alleles = IIT_typestring(snps_iit,Interval_type(interval));
      if (c == alleles[0] || c == alleles[1]) ;
#endif

      if (isupper(seq2[querypos]) && seq1[querypos] != seq2[querypos]) {
	label = IIT_label(snps_iit,snps[i],&allocp);
	if (printp) {
	  fprintf(fp,"|");
	}
	fprintf(fp,"%d@",querypos+1);
	fprintf(fp,"%s",label);
	printp = true;
	if (allocp) FREE(label);
      }
    }

  } else {

    for (i = nsnps-1; i >= 0; i--) {
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = (this->genomicstart-this->chroffset) - position;
      assert(querypos >= 0 && querypos < this->genomiclength);

#if 0
      /* printf("\n%d%c\n",querypos,c); */
      alleles = IIT_typestring(snps_iit,Interval_type(interval));
      if (c == alleles[0] || c == alleles[1]) ;
#endif

      if (isupper(seq2[querypos]) && seq1[querypos] != seq2[querypos]) {
	label = IIT_label(snps_iit,snps[i],&allocp);
	if (printp) {
	  fprintf(fp,"|");
	}
	fprintf(fp,"%d@",querypos+1);
	fprintf(fp,"%s",label);
	printp = true;
	if (allocp) FREE(label);
      }
    }

  }

  FREE(snps);

  return;
}



bool
Substring_contains_known_splicesite (T this) {
  Univcoord_T left;
  bool method2;
  int length;
#ifdef DEBUG7
  Chrpos_T start, end;
  bool method1;
#endif

  if (splicesites_iit == NULL) {
    return false;
  } else if (this->plusp == true) {
#ifdef DEBUG7
    start = this->alignstart+KNOWN_SPLICESITE_EDGE - this->chroffset + 1;
    end = this->alignend-KNOWN_SPLICESITE_EDGE - this->chroffset;
    method1 = IIT_exists_with_divno(splicesites_iit,splicesites_divint_crosstable[this->chrnum],start,end);
    debug7(printf("Checking if #%d:%u..%u has a known splice site: %d\n",this->chrnum,start,end,method1));
#endif

    left = this->alignstart;
    length = this->alignend - left;
    method2 = Splicetrie_splicesite_p(left,1,length-1);
    debug7(printf("Checking if #%d:%u (%d..%d) has a known splice site: %d\n",
		  this->chrnum,left - this->chroffset,1,length-1,method2));

    /* assert(method2 == method1); */
    /* if (method2 == true) exit(0); */
    return method2;

  } else {
#ifdef DEBUG7
    start = this->alignend+KNOWN_SPLICESITE_EDGE - this->chroffset + 1;
    end = this->alignstart-KNOWN_SPLICESITE_EDGE - this->chroffset;
    method1 = IIT_exists_with_divno(splicesites_iit,splicesites_divint_crosstable[this->chrnum],start,end);
    debug7(printf("Checking if #%d:%u..%u has a known splice site: %d\n",this->chrnum,start,end,method1));
#endif

    left = this->alignend;
    length = this->alignstart - left;
    method2 = Splicetrie_splicesite_p(left,1,length-1);
    debug7(printf("Checking if #%d:%u (%d..%d) has a known splice site: %d\n",
		  this->chrnum,left - this->chroffset,1,length-1,method2));

    /* assert(method2 == method1); */
    /* if (method2 == true) exit(0); */
    return method2;
  }
}



/* Has a copy in pair.c */
static void
print_splicesite_labels (FILE *fp, T this, int typeint, int chimera_pos, char *tag) {
  Chrpos_T splicesitepos;
  int *splicesites, nsplicesites, i;
  char *label;
  bool allocp;

  if (this->plusp == true) {
    splicesitepos = this->genomicstart - this->chroffset + chimera_pos;
  } else {
    splicesitepos = this->genomicstart - this->chroffset - chimera_pos;
  }

  /* Note: this->chimera_knownp might not be set for GMAP alignments */
  if (this->chimera_knownp == true) {
    /* Note: IIT_get_typed_signed_with_divno does not work here */
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[this->chrnum],
						    splicesitepos,splicesitepos+1U,typeint);
    if (nsplicesites == 0) {
#if 0
      /* Can happen if splicesites_iit are for introns, and not splicesites */
      fprintf(stderr,"Supposed to have a splicesite at chrnum %d, %u..%u, type %d\n",
	      this->chrnum,splicesitepos,splicesitepos+1U,typeint);
#endif
    } else {
      fprintf(fp,",%s:",tag);
      label = IIT_label(splicesites_iit,splicesites[0],&allocp);
      fprintf(fp,"%s",label);
      if (allocp) FREE(label);

      for (i = 1; i < nsplicesites; i++) {
	label = IIT_label(splicesites_iit,splicesites[i],&allocp);
	fprintf(fp,"|%s",label);
	if (allocp) FREE(label);
      }
      FREE(splicesites);
    }
  }

  return;
}


#if 0
/* replaced by intragenic_splice_p in stage1hr.c */
/* donor_genomicpos and acceptor_genomicpos are Univcoord_T */
bool
Substring_intragenic_splice_p (Genomicpos_T splicedistance, Chrnum_T chrnum,
			       Genomicpos_T donor_genomicpos, Genomicpos_T acceptor_genomicpos,
			       Genomicpos_T chroffset) {
  Genomicpos_T obsdistance;
  Genomicpos_T chrpos;
  int *splicesites, nsplicesites, i;
  char *annot, *restofheader;
  bool alloc_header_p;

  debug5(printf("Entered Substring_intragenic_splice_p with donor %u and acceptor %u, chroffset %u => chrpos %u and %u\n",
		donor_genomicpos,acceptor_genomicpos,chroffset,donor_genomicpos-chroffset,acceptor_genomicpos-chroffset));

  chrpos = donor_genomicpos - chroffset;
  splicesites = IIT_get_typed_with_divno(&nsplicesites,splicesites_iit,
					 splicesites_divint_crosstable[chrnum],
					 chrpos,chrpos,donor_typeint,/*sortp*/false);
  for (i = 0; i < nsplicesites; i++) {
    annot = IIT_annotation(&restofheader,splicesites_iit,splicesites[i],&alloc_header_p);
    debug5(printf("Comparing obsdistance %s at donor %u with splicedistance %u\n",restofheader,chrpos,splicedistance));
    if (sscanf(restofheader,"%u",&obsdistance) == 1 &&
	obsdistance >= splicedistance) {
      if (alloc_header_p == true) {
	FREE(restofheader);
      }
      FREE(splicesites);
      return true;
    }
    if (alloc_header_p == true) {
      FREE(restofheader);
    }
  }
  FREE(splicesites);

  chrpos = acceptor_genomicpos - chroffset;
  splicesites = IIT_get_typed_with_divno(&nsplicesites,splicesites_iit,
					 splicesites_divint_crosstable[chrnum],
					 acceptor_genomicpos,acceptor_genomicpos,donor_typeint,
					 /*sortp*/false);
  for (i = 0; i < nsplicesites; i++) {
    annot = IIT_annotation(&restofheader,splicesites_iit,splicesites[i],&alloc_header_p);
    debug5(printf("Comparing obsdistance %s at acceptor %u with splicedistance %u\n",restofheader,chrpos,splicedistance));
    if (sscanf(restofheader,"%u",&obsdistance) == 1 &&
	obsdistance >= splicedistance) {
      if (alloc_header_p == true) {
	FREE(restofheader);
      }
      FREE(splicesites);
      return true;
    }
    if (alloc_header_p == true) {
      FREE(restofheader);
    }
  }
  FREE(splicesites);

  return false;
}
#endif


/* Taken from NCBI Blast 2.2.29, algo/blast/core/blast_stat.c */
/* Karlin-Altschul formula: m n exp(-lambda * S + log k) = k m n exp(-lambda * S) */
/* Also in pair.c */

static double
blast_evalue (int alignlength, int nmismatches) {
  double k = 0.1;
  double lambda = 1.58;		/* For a +1, -1 scoring scheme */
  double score;
  
  score = (double) ((alignlength - nmismatches) /* scored as +1 */ - nmismatches /* scored as -1 */);

  return k * (double) alignlength * genomelength * exp(-lambda * score);
}

static double
blast_bitscore (int alignlength, int nmismatches) {
  double k = 0.1;
  double lambda = 1.58;		/* For a +1, -1 scoring scheme */
  double score;
  
  score = (double) ((alignlength - nmismatches) /* scored as +1 */ - nmismatches /* scored as -1 */);
  return (score * lambda - log(k)) / log(2.0);
}


void
Substring_print_m8 (FILE *fp, T substring, Shortread_T headerseq, char *acc_suffix,
		    char *chr, bool invertp) {
  double identity;
  int alignlength_trim;

  fprintf(fp,"%s%s",Shortread_accession(headerseq),acc_suffix); /* field 0: accession */

  fprintf(fp,"\t%s",chr);	/* field 1: chr */

  /* field 2: identity */
  if (substring->plusp == true) {
    alignlength_trim = (int) (substring->alignend_trim - substring->alignstart_trim);
  } else {
    alignlength_trim = (int) (substring->alignstart_trim - substring->alignend_trim);
  }

  identity = (double) (alignlength_trim - substring->nmismatches_bothdiff)/(double) alignlength_trim;
  fprintf(fp,"\t%.1f",100.0*identity);


  fprintf(fp,"\t%d",alignlength_trim); /* field 3: query length */

  fprintf(fp,"\t%d",substring->nmismatches_bothdiff); /* field 4: nmismatches */

  fprintf(fp,"\t0");		/* field 5: gap openings */

  fprintf(fp,"\t%d",substring->querystart + 1); /* field 6: query start */

  fprintf(fp,"\t%d",substring->queryend); /* field 7: query end */

  /* fields 8 and 9: chr start and end */
  if (substring->plusp == true) {
    if (invertp == false) {
      fprintf(fp,"\t%u\t%u",substring->alignstart_trim - substring->chroffset + 1U,
	      substring->alignend_trim - substring->chroffset);
    } else {
      fprintf(fp,"\t%u\t%u",substring->alignend_trim - substring->chroffset,
	      substring->alignstart_trim - substring->chroffset + 1U);
    }
  } else {
    if (invertp == false) {
      fprintf(fp,"\t%u\t%u",substring->alignstart_trim - substring->chroffset,
	      substring->alignend_trim - substring->chroffset + 1U);
    } else {
      fprintf(fp,"\t%u\t%u",substring->alignend_trim - substring->chroffset + 1U,
	      substring->alignstart_trim - substring->chroffset);
    }
  }

  /* field 10: E value */
  fprintf(fp,"\t%.2g",blast_evalue(alignlength_trim,substring->nmismatches_bothdiff));

 /* field 11: bit score */
  fprintf(fp,"\t%.1f",blast_bitscore(alignlength_trim,substring->nmismatches_bothdiff));
  
  fprintf(fp,"\n");

  return;
}



static void
print_forward (FILE *fp, char *string, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    fprintf(fp,"%c",string[i]);
  }
  return;
}




static void
print_lc (FILE *fp, char *string, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    fprintf(fp,"%c",(char) tolower(string[i]));
  }
  return;
}



static void
print_revcomp (FILE *fp, char *nt, int len) {
  int i;

  for (i = len-1; i >= 0; --i) {
    fprintf(fp,"%c",complCode[(int) nt[i]]);
  }
  return;
}

static void
print_revcomp_lc (FILE *fp, char *nt, int len) {
  int i;

  for (i = len-1; i >= 0; --i) {
    fprintf(fp,"%c",(char) tolower(complCode[(int) nt[i]]));
  }
  return;
}



static void
print_genomic (FILE *fp, T substring, char *deletion, int deletionlength, bool invertp,
	       Shortread_T queryseq) {
  int i;

  if (invertp == false) {
    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      Shortread_print_oneline_uc(fp,queryseq);

    } else if (show_refdiff_p == true) {
      print_forward(fp,substring->genomic_refdiff,substring->queryend);
      if (deletion != NULL) {
	print_lc(fp,deletion,deletionlength);
      }
      print_forward(fp,&(substring->genomic_refdiff[substring->queryend]),substring->querylength - substring->queryend);
    } else {
      print_forward(fp,substring->genomic_bothdiff,substring->queryend);
      if (deletion != NULL) {
	print_lc(fp,deletion,deletionlength);
      }
      print_forward(fp,&(substring->genomic_bothdiff[substring->queryend]),substring->querylength - substring->queryend);
    }

    for (i = 0; i < Shortread_choplength(queryseq); i++) {
      fprintf(fp,"*");
    }
    fprintf(fp,"\t");
    fprintf(fp,"%d..%d",1 + substring->querystart,substring->queryend);

  } else {
    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      Shortread_print_oneline_revcomp_uc(fp,queryseq);

    } else if (show_refdiff_p == true) {
      print_revcomp(fp,&(substring->genomic_refdiff[substring->querystart]),substring->querylength - substring->querystart);
      if (deletion != NULL) {
	print_revcomp_lc(fp,deletion,deletionlength);
      }
      print_revcomp(fp,substring->genomic_refdiff,substring->querystart);

    } else {
      print_revcomp(fp,&(substring->genomic_bothdiff[substring->querystart]),substring->querylength - substring->querystart);
      if (deletion != NULL) {
	print_revcomp_lc(fp,deletion,deletionlength);
      }
      print_revcomp(fp,substring->genomic_bothdiff,substring->querystart);
    }
    for (i = 0; i < Shortread_choplength(queryseq); i++) {
      fprintf(fp,"*");
    }
    fprintf(fp,"\t");
    fprintf(fp,"%d..%d",1 + substring->querylength - substring->queryend,
	   substring->querylength - substring->querystart);
  }
  return;
}


static void
print_coordinates (FILE *fp, T substring, char *chr, bool invertp) {

  if (substring->plusp == true) {
    if (invertp == false) {
      fprintf(fp,"+%s:%u..%u",chr,substring->alignstart_trim - substring->chroffset + 1U,
	     substring->alignend_trim - substring->chroffset);
    } else {
      fprintf(fp,"-%s:%u..%u",chr,substring->alignend_trim - substring->chroffset,
	     substring->alignstart_trim - substring->chroffset + 1U);
    }
  } else {
    if (invertp == false) {
      fprintf(fp,"-%s:%u..%u",chr,substring->alignstart_trim - substring->chroffset,
	     substring->alignend_trim - substring->chroffset + 1U);
    } else {
      fprintf(fp,"+%s:%u..%u",chr,substring->alignend_trim - substring->chroffset + 1U,
	     substring->alignstart_trim - substring->chroffset);
    }
  }

  return;
}



void
Substring_print_single (FILE *fp, T substring, Shortread_T queryseq,
			char *chr, bool invertp) {

  print_genomic(fp,substring,/*deletion*/(char *) NULL,/*deletionlength*/0,invertp,queryseq);
  fprintf(fp,"\t");
  print_coordinates(fp,substring,chr,invertp);

  fprintf(fp,"\t");
  if (invertp == false) {
    switch (substring->start_endtype) {
    case END: fprintf(fp,"start:%d",substring->trim_left); break;
    case TERM: fprintf(fp,"term:%d",substring->trim_left); break;
    default: fprintf(stderr,"start_endtype is %d\n",substring->start_endtype); abort(); break;
    }
  } else {
    switch (substring->end_endtype) {
    case END: fprintf(fp,"start:%d",substring->trim_right); break;
    case TERM: fprintf(fp,"term:%d",substring->trim_right); break;
    default: fprintf(stderr,"end_endtype is %d\n",substring->end_endtype); abort(); break;
    }
  }

  fprintf(fp,"..");

  if (invertp == false) {
    switch (substring->end_endtype) {
    case END: fprintf(fp,"end:%d",substring->trim_right); break;
    case TERM: fprintf(fp,"term:%d",substring->trim_right); break;
    default: fprintf(stderr,"end_endtype is %d\n",substring->end_endtype); abort(); break;
    }
  } else {
    switch (substring->start_endtype) {
    case END: fprintf(fp,"end:%d",substring->trim_left); break;
    case TERM: fprintf(fp,"term:%d",substring->trim_left); break;
    default: fprintf(stderr,"start_endtype is %d\n",substring->start_endtype); abort(); break;
    }
  }

  fprintf(fp,",matches:%d,sub:%d",substring->nmatches,substring->nmismatches_bothdiff);
  if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }

  return;
}


void
Substring_print_insertion_1 (FILE *fp, T substring1, T substring2, int nindels, 
			     Shortread_T queryseq, char *chr, bool invertp) {
  T substring;

  if (invertp == false) {
    substring = substring1;
    print_genomic(fp,substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false,
		  queryseq);

  } else {
    substring = substring2;
    print_genomic(fp,substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true,
		  queryseq);
  }

  fprintf(fp,"\t");

  print_coordinates(fp,substring,chr,invertp);


  fprintf(fp,"\t");
  if (invertp == false) {
    fprintf(fp,"start:%d..ins:%d,matches:%d,sub:%d",
	    substring->trim_left,nindels,substring->nmatches,substring->nmismatches_bothdiff);
  } else {
    fprintf(fp,"start:%d..ins:%d,matches:%d,sub:%d",
	    substring->trim_right,nindels,substring->nmatches,substring->nmismatches_bothdiff);
  }
  if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }


  return;
}

void
Substring_print_insertion_2 (FILE *fp, T substring1, T substring2, int nindels,
			     Shortread_T queryseq, char *chr, bool invertp) {
  T substring;

  if (invertp == false) {
    substring = substring2;
    print_genomic(fp,substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false,
		  queryseq);

  } else {
    substring = substring1;
    print_genomic(fp,substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true,
		  queryseq);
  }

  fprintf(fp,"\t");

  print_coordinates(fp,substring,chr,invertp);


  fprintf(fp,"\t");
  if (invertp == false) {
    fprintf(fp,"ins:%d..end:%d,matches:%d,sub:%d",
	    nindels,substring->trim_right,substring->nmatches,substring->nmismatches_bothdiff);
  } else {
    fprintf(fp,"ins:%d..end:%d,matches:%d,sub:%d",
	    nindels,substring->trim_left,substring->nmatches,substring->nmismatches_bothdiff);
  }
  if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }

  return;
}


void
Substring_print_deletion_1 (FILE *fp, T substring1, T substring2, int nindels, 
			    char *deletion, Shortread_T queryseq, char *chr, 
			    bool invertp) {
  T substring;

  if (invertp == false) {
    substring = substring1;
    print_genomic(fp,substring1,deletion,nindels,/*invertp*/false,queryseq);
  } else {
    substring = substring2;
    print_genomic(fp,substring2,deletion,nindels,/*invertp*/true,queryseq);
  }

  fprintf(fp,"\t");

  print_coordinates(fp,substring,chr,invertp);


  fprintf(fp,"\t");
  if (invertp == false) {
    fprintf(fp,"start:%d..del:%d,matches:%d,sub:%d",
	    substring->trim_left,nindels,substring->nmatches,substring->nmismatches_bothdiff);
  } else {
    fprintf(fp,"start:%d..del:%d,matches:%d,sub:%d",
	    substring->trim_right,nindels,substring->nmatches,substring->nmismatches_bothdiff);
  }
  if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }

  return;
}


void
Substring_print_deletion_2 (FILE *fp, T substring1, T substring2, int nindels, 
			    Shortread_T queryseq, char *chr, bool invertp) {
  T substring;

  if (invertp == false) {
    substring = substring2;
    print_genomic(fp,substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false,
		  queryseq);

  } else {
    substring = substring1;
    print_genomic(fp,substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true,
		  queryseq);
  }

  fprintf(fp,"\t");

  print_coordinates(fp,substring,chr,invertp);

  fprintf(fp,"\t");
  if (invertp == false) {
    fprintf(fp,"del:%d..end:%d,matches:%d,sub:%d",
	    nindels,substring->trim_right,substring->nmatches,substring->nmismatches_bothdiff);
  } else {
    fprintf(fp,"del:%d..end:%d,matches:%d,sub:%d",
	    nindels,substring->trim_left,substring->nmatches,substring->nmismatches_bothdiff);
  }
  if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }
  
  return;
}


/* This logic used in splice part of SAM_print */
static void
print_splice_distance (FILE *fp, T donor, T acceptor, Chrpos_T distance, bool sensep, char *tag) {
  bool normalp = true;

  if (donor == NULL || acceptor == NULL) {
    /* Don't print anything */
  } else if (distance == 0U) {
    fprintf(fp,",splice_type:%s",TRANSLOCATION_TEXT);
  } else {
    if (donor->plusp != acceptor->plusp) {
      fprintf(fp,",splice_type:%s",INVERSION_TEXT);
      normalp = false;
    } else if (donor->plusp == true) {
      if (sensep == true) {
	if (acceptor->genomicstart < donor->genomicstart) {
	  fprintf(fp,",splice_type:%s",SCRAMBLE_TEXT);
	  normalp = false;
	}
      } else {
	if (donor->genomicstart < acceptor->genomicstart) {
	  fprintf(fp,",splice_type:%s",SCRAMBLE_TEXT);
	  normalp = false;
	}
      }
    } else {
      if (sensep == true) {
	if (donor->genomicstart < acceptor->genomicstart) {
	  fprintf(fp,",splice_type:%s",SCRAMBLE_TEXT);
	  normalp = false;
	}
      } else {
	if (acceptor->genomicstart < donor->genomicstart) {
	  fprintf(fp,",splice_type:%s",SCRAMBLE_TEXT);
	  normalp = false;
	}
      }
    }
    if (normalp == true) {
      fprintf(fp,",splice_type:%s",CONSISTENT_TEXT);
    }
    fprintf(fp,",%s:%u",tag,distance);
  }

  return;
}

static void
print_shortexon_splice_distances (FILE *fp, Chrpos_T distance1, Chrpos_T distance2) {
  if (distance1 == 0U || distance2 == 0U) {
    /* Skip */
  } else {
    fprintf(fp,",splice_type:%s",CONSISTENT_TEXT);
    fprintf(fp,",splice_dist_1:%u",distance1);
    fprintf(fp,",splice_dist_2:%u",distance2);
  }

  return;
}



void
Substring_print_donor (FILE *fp, T donor, bool sensep, bool invertp, Shortread_T queryseq,
		       Univ_IIT_T chromosome_iit, T acceptor, Chrpos_T chimera_distance) {
  char *chr, *label_tag, *splice_dist_tag;
  bool allocp;

#ifdef CHECK_KNOWNI
  Chrpos_T splicesitepos;
  int *splicesites;
  int nsplicesites;
#endif

  print_genomic(fp,donor,/*deletion*/NULL,/*deletionlength*/0,invertp,queryseq);

  fprintf(fp,"\t");
  chr = Univ_IIT_label(chromosome_iit,donor->chrnum,&allocp);
  print_coordinates(fp,donor,chr,invertp);

  /* printf("donor chimera_pos is %d\n",donor->chimera_pos); */
  fprintf(fp,"\t");
  if (sensep == true && invertp == false) {
    fprintf(fp,"start:%d..donor:%.2f",donor->trim_left,donor->chimera_prob);
    label_tag = "label_2";
    splice_dist_tag = "splice_dist_2";
  } else if (sensep == true && invertp == true) {
    fprintf(fp,"donor:%.2f..end:%d",donor->chimera_prob,donor->trim_left);
    label_tag = "label_1";
    splice_dist_tag = "splice_dist_1";
  } else if (sensep == false && invertp == false) {
    fprintf(fp,"donor:%.2f..end:%d",donor->chimera_prob,donor->trim_right);
    label_tag = "label_1";
    splice_dist_tag = "splice_dist_1";
  } else if (sensep == false && invertp == true) {
    fprintf(fp,"start:%d..donor:%.2f",donor->trim_right,donor->chimera_prob);
    label_tag = "label_2";
    splice_dist_tag = "splice_dist_2";
  }

  fprintf(fp,",matches:%d,sub:%d",donor->nmatches,donor->nmismatches_bothdiff);
  if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",donor->nmismatches_refdiff - donor->nmismatches_bothdiff,donor->nmismatches_refdiff);
    if (print_snplabels_p && donor->nmismatches_refdiff > donor->nmismatches_bothdiff) {
      print_snp_labels(fp,donor,queryseq);
    }
  }

  if (sensep == true && invertp == false) {
    fprintf(fp,",dir:sense");
  } else if (sensep == true && invertp == true) {
    fprintf(fp,",dir:antisense");
  } else if (sensep == false && invertp == false) {
    fprintf(fp,",dir:antisense");
  } else if (sensep == false && invertp == true) {
    fprintf(fp,",dir:sense");
  }

  if (acceptor != NULL) {
    print_splice_distance(fp,donor,acceptor,chimera_distance,sensep,splice_dist_tag);
  }

#ifdef CHECK_KNOWNI
  if (donor->chimera_knownp == false && splicesites_iit) {
    if (donor->plusp == true) {
      splicesitepos = donor->genomicstart - donor->chroffset + donor->chimera_pos;
    } else {
      splicesitepos = donor->genomicstart - donor->chroffset - donor->chimera_pos;
    }
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[donor->chrnum],
						    splicesitepos,splicesitepos+1U,donor_typeint);
    assert(nsplicesites == 0);
  }
#endif

  if (donor->chimera_knownp && splicesites_iit) {
    print_splicesite_labels(fp,donor,donor_typeint,donor->chimera_pos,label_tag);
  }

  if (allocp == true) {
    FREE(chr);
  }

  return;
}

void 
Substring_print_acceptor (FILE *fp, T acceptor, bool sensep, bool invertp, Shortread_T queryseq, 
			  Univ_IIT_T chromosome_iit, T donor, Chrpos_T chimera_distance) {
  char *chr, *label_tag, *splice_dist_tag;
  bool allocp;

#ifdef CHECK_KNOWNI
  Chrpos_T splicesitepos;
  int *splicesites;
  int nsplicesites;
#endif

  print_genomic(fp,acceptor,/*deletion*/NULL,/*deletionlength*/0,invertp,queryseq);

  fprintf(fp,"\t");
  chr = Univ_IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
  print_coordinates(fp,acceptor,chr,invertp);

  /* printf("acceptor chimera_pos is %d\n",acceptor->chimera_pos); */
  fprintf(fp,"\t");
  if (sensep == true && invertp == false) {
    fprintf(fp,"acceptor:%.2f..end:%d",acceptor->chimera_prob,acceptor->trim_right);
    label_tag = "label_1";
    splice_dist_tag = "splice_dist_1";
  } else if (sensep == true && invertp == true) {
    fprintf(fp,"start:%d..acceptor:%.2f",acceptor->trim_right,acceptor->chimera_prob);
    label_tag = "label_2";
    splice_dist_tag = "splice_dist_2";
  } else if (sensep == false && invertp == false) {
    fprintf(fp,"start:%d..acceptor:%.2f",acceptor->trim_left,acceptor->chimera_prob);
    label_tag = "label_2";
    splice_dist_tag = "splice_dist_2";
  } else if (sensep == false && invertp == true) {
    fprintf(fp,"acceptor:%.2f..end:%d",acceptor->chimera_prob,acceptor->trim_left);
    label_tag = "label_1";
    splice_dist_tag = "splice_dist_1";
  }

  fprintf(fp,",matches:%d,sub:%d",acceptor->nmatches,acceptor->nmismatches_bothdiff);
  if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",acceptor->nmismatches_refdiff - acceptor->nmismatches_bothdiff,acceptor->nmismatches_refdiff);
    if (print_snplabels_p && acceptor->nmismatches_refdiff > acceptor->nmismatches_bothdiff) {
      print_snp_labels(fp,acceptor,queryseq);
    }
  }

  if (sensep == true && invertp == false) {
    fprintf(fp,",dir:sense");
  } else if (sensep == true && invertp == true) {
    fprintf(fp,",dir:antisense");
  } else if (sensep == false && invertp == false) {
    fprintf(fp,",dir:antisense");
  } else if (sensep == false && invertp == true) {
    fprintf(fp,",dir:sense");
  }

  if (donor != NULL) {
    print_splice_distance(fp,donor,acceptor,chimera_distance,sensep,splice_dist_tag);
  }

#ifdef CHECK_KNOWNI
  if (acceptor->chimera_knownp == false && splicesites_iit) {
    if (acceptor->plusp == true) {
      splicesitepos = acceptor->genomicstart - acceptor->chroffset + acceptor->chimera_pos;
    } else {
      splicesitepos = acceptor->genomicstart - acceptor->chroffset - acceptor->chimera_pos;
    }
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[acceptor->chrnum],
						    splicesitepos,splicesitepos+1U,acceptor_typeint);
    assert(nsplicesites == 0);
  }
#endif


  if (acceptor->chimera_knownp && splicesites_iit) {
    print_splicesite_labels(fp,acceptor,acceptor_typeint,acceptor->chimera_pos,label_tag);
  }


  if (allocp == true) {
    FREE(chr);
  }

  return;
}


void
Substring_print_shortexon (FILE *fp, T shortexon, bool sensep, bool invertp, Shortread_T queryseq,
			   Univ_IIT_T chromosome_iit, Chrpos_T distance1, Chrpos_T distance2) {
  char *chr;
  bool allocp;

#ifdef CHECK_KNOWNI
  Chrpos_T splicesitepos;
  int *splicesites;
  int nsplicesites;
#endif

#ifdef CHECK_KNOWNI
    if (shortexon->chimera_knownp == false && splicesites_iit) {
      if (shortexon->plusp == true) {
	splicesitepos = shortexon->genomicstart - shortexon->chroffset + shortexon->chimera_pos;
      } else {
	splicesitepos = shortexon->genomicstart - shortexon->chroffset - shortexon->chimera_pos;
      }
      splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						      splicesites_divint_crosstable[shortexon->chrnum],
						      splicesitepos,splicesitepos+1U,acceptor_typeint);
      assert(nsplicesites == 0);
    }

    if (shortexon->chimera_knownp_2 == false && splicesites_iit) {
      if (shortexon->plusp == true) {
	splicesitepos = shortexon->genomicstart - shortexon->chroffset + shortexon->chimera_pos_2;
      } else {
	splicesitepos = shortexon->genomicstart - shortexon->chroffset - shortexon->chimera_pos_2;
      }
      splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						      splicesites_divint_crosstable[shortexon->chrnum],
						      splicesitepos,splicesitepos+1U,donor_typeint);
      assert(nsplicesites == 0);
    }
#endif

  print_genomic(fp,shortexon,/*deletion*/NULL,/*deletionlength*/0,invertp,queryseq);

  fprintf(fp,"\t");

  chr = Univ_IIT_label(chromosome_iit,shortexon->chrnum,&allocp);
  print_coordinates(fp,shortexon,chr,invertp);

  fprintf(fp,"\t");
  if (sensep == true && invertp == false) {
    fprintf(fp,"acceptor:%.2f..donor:%.2f",shortexon->chimera_prob,shortexon->chimera_prob_2);
  } else if (sensep == true && invertp == true) {
    fprintf(fp,"donor:%.2f..acceptor:%.2f",shortexon->chimera_prob_2,shortexon->chimera_prob);
  } else if (sensep == false && invertp == false) {
    fprintf(fp,"donor:%.2f..acceptor:%.2f",shortexon->chimera_prob_2,shortexon->chimera_prob);
  } else if (sensep == false && invertp == true) {
    fprintf(fp,"acceptor:%.2f..donor:%.2f",shortexon->chimera_prob,shortexon->chimera_prob_2);
  }

  fprintf(fp,",matches:%d,sub:%d",shortexon->nmatches,shortexon->nmismatches_bothdiff);
  if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",shortexon->nmismatches_refdiff - shortexon->nmismatches_bothdiff,shortexon->nmismatches_refdiff);
    if (print_snplabels_p && shortexon->nmismatches_refdiff > shortexon->nmismatches_bothdiff) {
      print_snp_labels(fp,shortexon,queryseq);
    }
  }

  if (sensep == true && invertp == false) {
    fprintf(fp,",dir:sense");
    print_shortexon_splice_distances(fp,distance1,distance2);

    if (shortexon->chimera_knownp && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,acceptor_typeint,
			      shortexon->chimera_pos,/*tag*/"label_1");
    }
    if (shortexon->chimera_knownp_2 && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,donor_typeint,
			      shortexon->chimera_pos_2,/*tag*/"label_2");
    }

  } else if (sensep == true && invertp == true) {
    fprintf(fp,",dir:antisense");
    print_shortexon_splice_distances(fp,distance1,distance2);

    if (shortexon->chimera_knownp_2 && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,donor_typeint,
			      shortexon->chimera_pos_2,/*tag*/"label_1");
    }
    if (shortexon->chimera_knownp && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,acceptor_typeint,
			      shortexon->chimera_pos,/*tag*/"label_2");
    }

  } else if (sensep == false && invertp == false) {
    fprintf(fp,",dir:antisense");
    print_shortexon_splice_distances(fp,distance1,distance2);



    if (shortexon->chimera_knownp_2 && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,donor_typeint,
			      shortexon->chimera_pos_2,/*tag*/"label_1");
    }

    if (shortexon->chimera_knownp && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,acceptor_typeint,
			      shortexon->chimera_pos,/*tag*/"label_2");
    }

  } else if (sensep == false && invertp == true) {
    fprintf(fp,",dir:sense");
    print_shortexon_splice_distances(fp,distance1,distance2);
    if (shortexon->chimera_knownp && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,acceptor_typeint,
			      shortexon->chimera_pos,/*tag*/"label_1");
    }
    if (shortexon->chimera_knownp_2 && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,donor_typeint,
			      shortexon->chimera_pos_2,/*tag*/"label_2");
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  return;
}


/* Needs to be here to access splicesites_iit */
void
Substring_print_gmap (FILE *fp, struct Pair_T *pairs, int npairs, int nsegments, bool invertedp,
		      Endtype_T start_endtype, Endtype_T end_endtype,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		      int querylength, bool watsonp, int cdna_direction, int score,
		      int insertlength, int pairscore, int mapq_score, Univ_IIT_T chromosome_iit) {
  Pair_print_gsnap(fp,pairs,npairs,nsegments,invertedp,
		   start_endtype,end_endtype,chrnum,chroffset,chrhigh,querylength,
		   watsonp,cdna_direction,score,insertlength,pairscore,mapq_score,
		   chromosome_iit,splicesites_iit,splicesites_divint_crosstable,
		   donor_typeint,acceptor_typeint);
  return;
}


/************************************************************************
 *   Resolve against known genes
 ************************************************************************/

Overlap_T
Substring_gene_overlap (T this, bool favor_multiexon_p) {
  Chrpos_T start, end;

  if (genes_iit == NULL) {
    return false;
  } else if (this->plusp == true) {
    start = this->alignstart - this->chroffset + 1;
    end = this->alignend - this->chroffset;
    debug7(printf("Checking if #%d:%u..%u overlaps an exon in a known gene\n",this->chrnum,start,end));
    return IIT_gene_overlap(genes_iit,genes_divint_crosstable[this->chrnum],start,end,favor_multiexon_p);

  } else {
    start = this->alignend - this->chroffset + 1;
    end = this->alignstart - this->chroffset;
    debug7(printf("Checking if #%d:%u..%u overlaps an exon in a known gene\n",this->chrnum,start,end));
    return IIT_gene_overlap(genes_iit,genes_divint_crosstable[this->chrnum],start,end,favor_multiexon_p);
  }
}



/************************************************************************
 *   Multiclean
 ************************************************************************/

static char *
get_total_tally (long int *tally, char *ptr) {
  int n;
  char *end;

  if ((end = index(ptr,'\n')) == NULL) {
    fprintf(stderr,"Premature end of line %s\n",ptr);
    return 0;
  }
  /* fprintf(stderr,"Getting tally for %.*s\n",end-ptr,ptr); */

  while (ptr < end) {
    while (ptr < end && !isdigit((int) *ptr)) {
      ptr++;
    }
    if (ptr < end) {
      sscanf(ptr,"%d",&n);
#if 0
      debug(if (n > 0) printf(" %d",n));
#endif
      (*tally) += n;
      while (ptr < end && !isspace(*ptr)) {
	ptr++;
      }
      while (ptr < end && isspace(*ptr)) {
	ptr++;
      }
    }
  }

  return ptr;
}

long int
Substring_tally (T this, IIT_T tally_iit, int *tally_divint_crosstable) {
  long int total = 0U;
  Interval_T interval;
  char *annotation, *restofheader, *ptr;
#if 0
  bool alloc_chr_p;
#endif
  bool allocp;
  Chrpos_T chrpos, intervalend;

  char *chr;
  Chrpos_T coordstart, coordend, pos5, pos3;
  int *matches;
  int nmatches, i;
  
  pos5 = this->alignstart_trim - this->chroffset;
  pos3 = this->alignend_trim - this->chroffset;

  if (pos5 < pos3) {
    coordstart = pos5;
    coordend = pos3;
  } else {
    coordstart = pos3;
    coordend = pos5;
  }
  coordstart += 1U;		/* Because tally IIT is 1-based */
  debug(printf("coordstart = %u, coordend = %u\n",coordstart,coordend));

#if 0
  chr = Univ_IIT_label(chromosome_iit,this->chrnum,&alloc_chr_p);
#endif
  matches = IIT_get_with_divno(&nmatches,tally_iit,tally_divint_crosstable[this->chrnum],
			       coordstart,coordend,/*sortp*/false);

  for (i = 0; i < nmatches; i++) {
    annotation = IIT_annotation(&restofheader,tally_iit,matches[i],&allocp);

    interval = IIT_interval(tally_iit,matches[i]);
    chrpos = Interval_low(interval);
    intervalend = Interval_high(interval);

    ptr = annotation;

    while (chrpos < coordstart) {
      if ((ptr = index(ptr,'\n')) == NULL) {
	fprintf(stderr,"Premature end of tally from %u to %u\n",
		Interval_low(interval),Interval_high(interval));
	return total;
      } else {
	ptr++;
      }
      chrpos++;
    }

    while (chrpos <= intervalend && chrpos <= coordend) {
      ptr = get_total_tally(&total,ptr);
      ptr++;
      chrpos++;
    }

    if (allocp == true) {
      FREE(restofheader);
    }
  }

  FREE(matches);

#if 0
  if (alloc_chr_p) {
    FREE(chr);
  }
#endif

  debug(printf("Subtotal = %ld\n",total));
  return total;
}


bool
Substring_runlength_p (T this, IIT_T runlength_iit, int *runlength_divint_crosstable) {
  Chrpos_T coordstart, coordend, pos5, pos3;
  
  pos5 = this->alignstart_trim - this->chroffset;
  pos3 = this->alignend_trim - this->chroffset;

  if (pos5 < pos3) {
    coordstart = pos5;
    coordend = pos3;
  } else {
    coordstart = pos3;
    coordend = pos5;
  }
  coordstart += 1U;		/* Because runlength IIT is 1-based */
  debug(printf("coordstart = %u, coordend = %u\n",coordstart,coordend));

  /* chr = Univ_IIT_label(chromosome_iit,this->chrnum,&alloc_chr_p); */
  return IIT_exists_with_divno(runlength_iit,runlength_divint_crosstable[this->chrnum],
				coordstart,coordend);
}


int
Substring_count_mismatches_region (T this, int trim_left, int trim_right,
				   Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  int left_bound, right_bound;

  if (this == NULL) {
    return 0;
  }

  left_bound = trim_left;
  right_bound = this->querylength - trim_right;

  if (this->queryend_orig < left_bound) {
    return 0;
  } else if (this->querystart_orig > right_bound) {
    return 0;
  } else {
    if (this->querystart_orig > left_bound) {
      left_bound = this->querystart_orig;
    }
    if (this->queryend_orig < right_bound) {
      right_bound = this->queryend_orig;
    }

    if (this->plusp) {
      return Genome_count_mismatches_substring(query_compress_fwd,this->left,/*pos5*/left_bound,
					       /*pos3*/right_bound,/*plusp*/true,this->genestrand,this->first_read_p);
    } else {
      return Genome_count_mismatches_substring(query_compress_rev,this->left,
					       /*pos5*/this->querylength - right_bound,
					       /*pos3*/this->querylength - left_bound,
					       /*plusp*/false,this->genestrand,this->first_read_p);
    }
  }
}


/************************************************************************
 *   Conversion to Pair_T format
 ************************************************************************/

List_T
Substring_convert_to_pairs (List_T pairs, T substring, Shortread_T queryseq,
			    int clipdir, int hardclip_low, int hardclip_high,
			    bool first_read_p, int queryseq_offset) {
  int querystart, queryend, querypos, i;
  Chrpos_T chrpos;
  char *seq1;
  char genome;

  if (substring == NULL) {
    return pairs;
  }


  seq1 = Shortread_fullpointer_uc(queryseq);
  if (substring->plusp == true) {
    if (hardclip_low > substring->querystart) {
      querystart = hardclip_low;
    } else {
      querystart = substring->querystart;
    }

    if (substring->querylength - hardclip_high < substring->queryend) {
      queryend = substring->querylength - hardclip_high;
    } else {
      queryend = substring->queryend;
    }
    debug6(printf("querystart %d, queryend %d, plusp %d\n",querystart,queryend,substring->plusp));
    debug6(printf("alignstart %u, alignend %u\n",substring->alignstart_trim - substring->chroffset,
		  substring->alignend_trim - substring->chroffset));

    chrpos = substring->genomicstart + querystart - substring->chroffset + 1U;
    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++) {
	pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos++,
							  seq1[i],/*comp*/MATCH_COMP,seq1[i]));
      }
    } else if (show_refdiff_p == true) {
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++) {
	if (isupper(genome = substring->genomic_refdiff[i])) {
	  pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos++,
							    seq1[i],/*comp*/MATCH_COMP,genome));
	} else {
	  pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos++,
							    seq1[i],/*comp*/MISMATCH_COMP,toupper(genome)));
	}
      }
    } else {
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++) {
	if (isupper(genome = substring->genomic_bothdiff[i])) {
	  pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos++,
							    seq1[i],/*comp*/MATCH_COMP,genome));
	} else {
	  pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos++,
							    seq1[i],/*comp*/MISMATCH_COMP,toupper(genome)));
	}
      }
    }

  } else {
    if (hardclip_high > substring->querystart) {
      querystart = hardclip_high;
    } else {
      querystart = substring->querystart;
    }

    if (substring->querylength - hardclip_low < substring->queryend) {
      queryend = substring->querylength - hardclip_low;
    } else {
      queryend = substring->queryend;
    }
    debug6(printf("querystart %d, queryend %d, plusp %d\n",querystart,queryend,substring->plusp));
    debug6(printf("alignstart %u, alignend %u\n",substring->alignstart_trim - substring->chroffset,
		  substring->alignend_trim - substring->chroffset));

    chrpos = substring->genomicstart - querystart - substring->chroffset;
    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++) {
	pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos--,
							  seq1[i],/*comp*/MATCH_COMP,seq1[i]));
      }
    } else if (show_refdiff_p == true) {
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++) {
	if (isupper(genome = substring->genomic_refdiff[i])) {
	  pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos--,
							    seq1[i],/*comp*/MATCH_COMP,genome));
	} else {
	  pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos--,
							    seq1[i],/*comp*/MISMATCH_COMP,toupper(genome)));
	}
      }
    } else {
      for (i = querystart, querypos = queryseq_offset + querystart; i < queryend; i++) {
	if (isupper(genome = substring->genomic_bothdiff[i])) {
	  pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos--,
							    seq1[i],/*comp*/MATCH_COMP,genome));
	} else {
	  pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrpos--,
							    seq1[i],/*comp*/MISMATCH_COMP,toupper(genome)));
	}
      }
    }
  }

  return pairs;
}


List_T
Substring_add_insertion (List_T pairs, T substringA, T substringB, int insertionlength, Shortread_T queryseq,
			 int clipdir, int hardclip_low, int hardclip_high, bool first_read_p, int queryseq_offset) {
  int querystartA, queryendA, querystartB, queryendB, querypos, i;
  Chrpos_T chrendA;
  char *seq1;


  if (substringA->plusp == true) {
    if (hardclip_low > substringA->querystart) {
      querystartA = hardclip_low;
    } else {
      querystartA = substringA->querystart;
    }

    if (substringA->querylength - hardclip_high < substringA->queryend) {
      queryendA = substringA->querylength - hardclip_high;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_low > substringB->querystart) {
      querystartB = hardclip_low;
    } else {
      querystartB = substringB->querystart;
    }

    if (substringB->querylength - hardclip_high < substringB->queryend) {
      queryendB = substringB->querylength - hardclip_high;
    } else {
      queryendB = substringB->queryend;
    }

    chrendA = substringA->genomicstart + queryendA - substringA->chroffset + 1U;

  } else {
    if (hardclip_high > substringA->querystart) {
      querystartA = hardclip_high;
    } else {
      querystartA = substringA->querystart;
    }

    if (substringA->querylength - hardclip_low < substringA->queryend) {
      queryendA = substringA->querylength - hardclip_low;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_high > substringB->querystart) {
      querystartB = hardclip_high;
    } else {
      querystartB = substringB->querystart;
    }

    if (substringB->querylength - hardclip_low < substringB->queryend) {
      queryendB = substringB->querylength - hardclip_low;
    } else {
      queryendB = substringB->queryend;
    }

    chrendA = substringA->genomicstart - queryendA - substringA->chroffset;
  }

  if (querystartA <= queryendA && querystartB <= queryendB) {
    seq1 = Shortread_fullpointer_uc(queryseq);
    querypos = queryseq_offset + queryendA;
    i = queryendA;
    while (--insertionlength >= 0) {
      pairs = List_push_out(pairs,(void *) Pair_new_out(querypos++,/*genomepos*/chrendA,
							seq1[i++],/*comp*/INDEL_COMP,' '));
    }
  }

  return pairs;
}


List_T
Substring_add_deletion (List_T pairs, T substringA, T substringB, char *deletion, int deletionlength,
			int clipdir, int hardclip_low, int hardclip_high, bool first_read_p, int queryseq_offset) {
  int querystartA, queryendA, querystartB, queryendB, k;
  Chrpos_T chrendA;

  if (substringA->plusp == true) {
    if (hardclip_low > substringA->querystart) {
      querystartA = hardclip_low;
    } else {
      querystartA = substringA->querystart;
    }

    if (substringA->querylength - hardclip_high < substringA->queryend) {
      queryendA = substringA->querylength - hardclip_high;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_low > substringB->querystart) {
      querystartB = hardclip_low;
    } else {
      querystartB = substringB->querystart;
    }

    if (substringB->querylength - hardclip_high < substringB->queryend) {
      queryendB = substringB->querylength - hardclip_high;
    } else {
      queryendB = substringB->queryend;
    }

    chrendA = substringA->genomicstart + queryendA - substringA->chroffset + 1U;

    if (querystartA < queryendA && querystartB < queryendB) {
      queryendA += queryseq_offset;
      for (k = 0; k < deletionlength; k++) {
	pairs = List_push_out(pairs,(void *) Pair_new_out(queryendA,/*genomepos*/chrendA++,
							  ' ',/*comp*/INDEL_COMP,deletion[k]));
      }
    }

  } else {
    if (hardclip_high > substringA->querystart) {
      querystartA = hardclip_high;
    } else {
      querystartA = substringA->querystart;
    }

    if (substringA->querylength - hardclip_low < substringA->queryend) {
      queryendA = substringA->querylength - hardclip_low;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_high > substringB->querystart) {
      querystartB = hardclip_high;
    } else {
      querystartB = substringB->querystart;
    }

    if (substringB->querylength - hardclip_low < substringB->queryend) {
      queryendB = substringB->querylength - hardclip_low;
    } else {
      queryendB = substringB->queryend;
    }

    chrendA = substringA->genomicstart - queryendA - substringA->chroffset;

    if (querystartA <= queryendA && querystartB <= queryendB) {
      queryendA += queryseq_offset;
      for (k = 0; k < deletionlength; k++) {
	pairs = List_push_out(pairs,(void *) Pair_new_out(queryendA,/*genomepos*/chrendA--,
							  ' ',/*comp*/INDEL_COMP,deletion[k]));
      }
    }
  }

  return pairs;
}



List_T
Substring_add_intron (List_T pairs, T substringA, T substringB,
		      int clipdir, int hardclip_low, int hardclip_high,
		      bool first_read_p, int queryseq_offset) {
  int querystartA, queryendA, querystartB, queryendB;
  Chrpos_T chrendA;

  if (substringA->plusp == true) {
    if (hardclip_low > substringA->querystart) {
      querystartA = hardclip_low;
    } else {
      querystartA = substringA->querystart;
    }

    if (substringA->querylength - hardclip_high < substringA->queryend) {
      queryendA = substringA->querylength - hardclip_high;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_low > substringB->querystart) {
      querystartB = hardclip_low;
    } else {
      querystartB = substringB->querystart;
    }

    if (substringB->querylength - hardclip_high < substringB->queryend) {
      queryendB = substringB->querylength - hardclip_high;
    } else {
      queryendB = substringB->queryend;
    }

    chrendA = substringA->genomicstart + queryendA - substringA->chroffset + 1U;

  } else {
    if (hardclip_high > substringA->querystart) {
      querystartA = hardclip_high;
    } else {
      querystartA = substringA->querystart;
    }

    if (substringA->querylength - hardclip_low < substringA->queryend) {
      queryendA = substringA->querylength - hardclip_low;
    } else {
      queryendA = substringA->queryend;
    }

    if (hardclip_high > substringB->querystart) {
      querystartB = hardclip_high;
    } else {
      querystartB = substringB->querystart;
    }

    if (substringB->querylength - hardclip_low < substringB->queryend) {
      queryendB = substringB->querylength - hardclip_low;
    } else {
      queryendB = substringB->queryend;
    }

    chrendA = substringA->genomicstart - queryendA - substringA->chroffset;
  }

  if (querystartA <= queryendA && querystartB <= queryendB) {
    /* Add gapholder */
    /* All we really need for Pair_print_sam is to set gapp to be true */
    pairs = List_push_out(pairs,(void *) Pair_new_out(queryseq_offset + queryendA,/*genomepos*/chrendA,
						      ' ',/*comp*/FWD_CANONICAL_INTRON_COMP,' '));
  }

  return pairs;
}

