static char rcsid[] = "$Id: pair.c 158535 2015-02-12 21:33:37Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "pair.h"
#include "pairdef.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include <math.h>		/* For rint(), abs() */
#include <ctype.h>		/* For toupper */

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "comp.h"
#include "complement.h"
#include "intron.h"
#include "listdef.h"
#include "intlist.h"
#include "separator.h"
#include "scores.h"
#include "segmentpos.h"
#include "translation.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "mapq.h"

#ifndef PMAP
#include "substring.h"		/* For Endtype_T */
#include "stage3hr.h"
#include "samflags.h"
#endif

#ifdef GSNAP
#include "samprint.h"
#endif


#define ONEBASEDP 1		/* 1-based coordinates.  Also defined in segmentpos.c */

#define MIN_INTRONLEN 20	/* For deciding between N and D in cigar string */



/* Check for ANSI mode, which does not include rint */
#ifdef __STRICT_ANSI__
#define rint(x) floor(0.5+(x))
#endif

#define DEFAULT_MARGIN 14

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Print pointer information in Pair_dump_one */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* PSL indels */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Pair_fracidentity_max */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* compute_md_string */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Phase information */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* trimming */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* end_bound and start_bound */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* binary search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif



#define TRIM_MATCH_SCORE 1
static int trim_mismatch_score;
static int trim_indel_score;
static bool gff3_separators_p;
static bool sam_insert_0M_p = false;
static bool force_xs_direction_p;
static bool md_lowercase_variant_p;
static bool snps_p;
static double genomelength;	/* For BLAST E-value */
static Cigar_action_T cigar_action;


void
Pair_setup (int trim_mismatch_score_in, int trim_indel_score_in,
	    bool gff3_separators_p_in, bool sam_insert_0M_p_in, bool force_xs_direction_p_in,
	    bool md_lowercase_variant_p_in, bool snps_p_in, Univcoord_T genomelength_in,
	    Cigar_action_T cigar_action_in) {
  trim_mismatch_score = trim_mismatch_score_in;
  trim_indel_score = trim_indel_score_in;
  gff3_separators_p = gff3_separators_p_in;
  sam_insert_0M_p = sam_insert_0M_p_in;
  force_xs_direction_p = force_xs_direction_p_in;
  md_lowercase_variant_p = md_lowercase_variant_p_in;
  snps_p = snps_p_in;
  genomelength = (double) genomelength_in;
  cigar_action = cigar_action_in;

  return;
}



#define T Pair_T

int
Pair_querypos (T this) {
  return this->querypos;
}

Chrpos_T
Pair_genomepos (T this) {
  return this->genomepos;
}

char
Pair_cdna (T this) {
  return this->cdna;
}

char
Pair_comp (T this) {
  return this->comp;
}

char
Pair_genome (T this) {
  return this->genome;
}

char
Pair_genomealt (T this) {
  return this->genomealt;
}

bool
Pair_gapp (T this) {
  return this->gapp;
}

bool
Pair_shortexonp (T this) {
  return this->shortexonp;
}


void
Pair_print_ends (List_T pairs) {
  List_T p;
  T start, end;

  if (pairs == NULL) {
    printf("0..0, 0..0\n");
  } else {
    start = (T) pairs->first;
    for (p = pairs; p != NULL; p = p->rest) {
      end = (T) p->first;
    }
    printf("%d..%d %u..%u",start->querypos,end->querypos,start->genomepos,end->genomepos);
  }
  return;
}


void
Pair_set_genomepos (struct Pair_T *pairarray, int npairs,
		    Univcoord_T chroffset, Univcoord_T chrhigh, bool watsonp) {
  int i;
  Chrpos_T chraliaslength;

  if (watsonp == true) {
    /* No need to adjust, since we are using chromosomal coordinates already */
  } else {
    chraliaslength = chrhigh - chroffset;
    for (i = 0; i < npairs; i++) {
      pairarray[i].genomepos = chraliaslength - pairarray[i].genomepos;
    }
  }
  return;
}


#if 0
/* Don't change list, just pairarray */
void
Pair_set_genomepos_list (List_T pairs, Univcoord_T chroffset,
			 Univcoord_T chrhigh, bool watsonp) {
  List_T p;
  T pair;
  Chrpos_T chraliaslength;

  if (watsonp == true) {
    /* No need to adjust, since we are using chromosomal coordinates already */
  } else {
    chraliaslength = chrhigh - chroffset;
    for (p = pairs; p != NULL; p = p->rest) {
      pair = (T) p->first;
      pair->genomepos = chraliaslength - pair->genomepos;
    }
  }

  return;
}
#endif


/* For outbuffer usage (e.g., truncation), use Pair_clip_bounded_array instead */
/* Note: This code is designed to handle source, which may still have
   gaps with querypos undefined */
List_T
Pair_clip_bounded_list (List_T source, int minpos, int maxpos) {
  List_T dest, *prev, p;
  T pair;
  int starti = -1, endi = -1, i;

  if (source == NULL) {
    return (List_T) NULL;
  } else {
    for (p = source, i = 0; p != NULL; p = p->rest, i++) {
      pair = (Pair_T) List_head(p);
      if (pair->querypos == minpos) {
	starti = i;		/* Advances in case of ties */
      } else if (pair->querypos > minpos && starti < 0) {
	starti = i;		/* Handles case where minpos was skipped */
      }

      if (pair->querypos == maxpos && endi < 0) {
	endi = i + 1;		/* Does not advance in case of tie */
      } else if (pair->querypos > maxpos && endi < 0) {
	endi = i;	   /* Handles case where maxpos was skipped */
      }
    }

    if (starti < 0 && endi < 0) {
      /* None of the pairs fall within bounds */
      return (List_T) NULL;
    } else {
      if (starti < 0) {
	starti = 0;
      }
      if (endi < 0) {
	endi = i;
      }
    }

    p = source;
    i = 0;
    while (i < starti) {
      p = p->rest;
      i++;
    }

    dest = p;
    prev = &p->rest;
    while (i < endi) {
      prev = &p->rest;
      p = p->rest;
      i++;
    }

    *prev = NULL;		/* Clip rest of list */
    return dest;
  }
}


int
Pair_clip_bounded_array (struct T *source, int npairs, int minpos, int maxpos) {
  T pair;
  int starti = -1, endi = -1, i, k;

#if 0
  printf("Pair_clip_bounded_array called with %d pairs, minpos %d, maxpos %d\n",npairs,minpos,maxpos);
  Pair_dump_array(source,npairs,true);
#endif

  for (i = 0; i < npairs; i++) {
    pair = &(source[i]);
    if (pair->querypos == minpos) {
      starti = i;		/* Advances in case of ties */
    } else if (pair->querypos > minpos && starti < 0) {
      starti = i;		/* Handles case where minpos was skipped */
    }

    if (pair->querypos == maxpos && endi < 0) {
      endi = i + 1;		/* Does not advance in case of tie */
    } else if (pair->querypos > maxpos && endi < 0) {
      endi = i;	   /* Handles case where maxpos was skipped */
    }
  }

  if (starti < 0 && endi < 0) {
    /* None of the pairs fall within bounds.  Don't do anything. */
    return npairs;
  } else {
    if (starti < 0) {
      starti = 0;
    }
    if (endi < 0) {
      endi = i;
    }
  }

  k = 0;
  for (i = starti; i < endi; i++) {
    memcpy((void *) &(source[k++]),(void *) &(source[i]),sizeof(struct T));
  }

  return endi - starti;
}




/* Head of list is the medial part of the read */
List_T
Pair_protect_end5 (List_T pairs, Pairpool_T pairpool) {
  List_T p;
  T pair;

  p = pairs;

  /* Go until known splice is seen */
  while (p != NULL && ((T) p->first)->gapp == false) {
    pair = (T) p->first;
    pair->protectedp = true;
    p = p->rest;
  }

  /* Handle known splice */
  if (p != NULL) {
    pair = (T) p->first;
    pair->protectedp = true;
    p = p->rest;
  }

  /* Continue until distal indel is seen */
  while (p != NULL && ((T) p->first)->cdna != ' ' && ((T) p->first)->genome != ' ') {
    pair = (T) p->first;
    pair->protectedp = true;
    p = p->rest;
  }

  /* Do not protect the sequence after the distal indel */
  while (p != NULL) {
    pair = (T) p->first;
    pair->protectedp = false;
    p = p->rest;
  }

  return pairs;
}


/* Head of list is the 3' distal end of the read */
List_T
Pair_protect_end3 (List_T pairs, Pairpool_T pairpool) {
  List_T result = NULL, p;
  T pair;

  p = pairs = List_reverse(pairs); /* Now head is medial end */

  /* Go until known splice is seen */
  while (p != NULL && ((T) p->first)->gapp == false) {
    pair = (T) p->first;
    pair->protectedp = true;
    result = Pairpool_push_existing(result,pairpool,pair);
    p = p->rest;
  }

  /* Handle known splice */
  if (p != NULL) {
    pair = (T) p->first;
    pair->protectedp = true;
    result = Pairpool_push_existing(result,pairpool,pair);
    p = p->rest;
  }
    
  /* Continue until distal indel is seen */
  while (p != NULL && ((T) p->first)->cdna != ' ' && ((T) p->first)->genome != ' ') {
    pair = (T) p->first;
    pair->protectedp = true;
    result = Pairpool_push_existing(result,pairpool,pair);
    p = p->rest;
  }
    
  /* Do not protect the sequence after the distal indel */
  while (p != NULL) {
    pair = (T) p->first;
    pair->protectedp = false;
    result = Pairpool_push_existing(result,pairpool,pair);
    p = p->rest;
  }

  return result;
}


void
Pair_protect_list (List_T pairs) {
  List_T p;
  T pair;

  for (p = pairs; p != NULL; p = p->rest) {
    pair = (T) p->first;
    pair->protectedp = true;
  }

  return;
}



/* For output thread only.  Pairs needed by worker threads are made in pairpool.c */
T
Pair_new_out (int querypos, Chrpos_T genomepos, char cdna, char comp, char genome) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  
  new->querypos = querypos;
  new->genomepos = genomepos;
  new->aapos = 0;
  new->aaphase_g = -1;
  new->aaphase_e = -1;
  new->cdna = cdna;
  new->comp = comp;
  new->genome = genome;
  new->genomealt = genome;
  new->dynprogindex = 0;

  new->aa_g = ' ';
  new->aa_e = ' ';
  new->shortexonp = false;
  new->introntype = NONINTRON;

  switch (comp) {
  case FWD_CANONICAL_INTRON_COMP /*'>'*/:
  case REV_CANONICAL_INTRON_COMP /*'<'*/:
  case NONINTRON_COMP /*'='*/:
  case SHORTGAP_COMP /*'~'*/:
  case INTRONGAP_COMP /*'.'*/:
  case FWD_GCAG_INTRON_COMP /*')'*/:
  case REV_GCAG_INTRON_COMP /*'('*/:
  case FWD_ATAC_INTRON_COMP /*']'*/:
  case REV_ATAC_INTRON_COMP /*'['*/:
  case DUALBREAK_COMP /*'#'*/:
    new->gapp = true; break;
  default: new->gapp = false;
  }

  new->extraexonp = false;

  new->protectedp = false;
  new->disallowedp = false;

  return new;
}

void
Pair_free_out (T *old) {
  if (*old) {
    FREE_OUT(*old);
  }
  return;
}


/* Print routines */

static char *RULER = "    .    :    .    :    .    :    .    :    .    :";
static void
print_top_ruler (FILE *fp, int n, int npairs, int margin, int wraplength) {
  fprintf(fp,"%*d ",margin,n);
  if (n + wraplength < npairs) {
    fprintf(fp,"%s\n",RULER);
  } else {
    fprintf(fp,"%.*s\n",npairs-n,RULER);
  }
  return;
}

/*
static void
print_bottom_ruler (int n, int npairs, int margin, int wraplength) {
  printf("%*s ",margin,"");
  if (n + wraplength < npairs) {
    printf("%s\n",RULER);
  } else {
    printf("%.*s\n",npairs-n,RULER);
  }
  return;
}
*/


static void
print_cdna_sequence (FILE *fp, struct T *ptr, int n, int npairs, int margin, int wraplength) {
  struct T *this;
  int i;

  this = ptr;
  fprintf(fp,"%*u ",margin,this->querypos + ONEBASEDP);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    putc(this->cdna,fp);
  }
  putc('\n',fp);
  return;
}

static int
find_aapos_in_line (struct T *ptr, int n, int npairs, int wraplength, 
		    bool genomep) {
  struct T *this, *last;

  if (npairs - n < wraplength) {
    last = &ptr[npairs - n - 1];
  } else {
    last = &ptr[wraplength - 1];
  }
  this = ptr;
  while (this <= last && (genomep ? this->aa_g : this->aa_e) == ' ') {
    this++;
  }

  if (this > last) {
    /* No aa found */
    return -1;
  } else {
    return this->aapos;
  }
}


static void
print_peptide (FILE *fp, struct T *ptr, int n, int npairs, int margin,
	       int wraplength, bool genomep) {
  struct T *this;
  int aapos, i;

  if ((aapos = find_aapos_in_line(ptr,n,npairs,wraplength,genomep)) < 0) {
    fprintf(fp,"%*s ",margin,"");
  } else {
    /* 4 is length of "aa.c" and "aa.g" */
    if (genomep == true) {
      fprintf(fp,"aa.g%*d ",margin-4,aapos);
    } else {
      fprintf(fp,"aa.c%*d ",margin-4,aapos);
    }
  }

  if (genomep == true) {
    for (i = 0; n < npairs && i < wraplength; n++, i++) {
      this = ptr++;
      putc(this->aa_g,fp);
    }
  } else {
    for (i = 0; n < npairs && i < wraplength; n++, i++) {
      this = ptr++;
      putc(this->aa_e,fp);
    }
  }

  putc('\n',fp);
  return;
}

static void
print_alignment (FILE *fp, struct T *ptr, int n, int npairs, bool diagnosticp, 
		 int margin, int wraplength) {
  struct T *this;
  int i;

  fprintf(fp,"%*s ",margin,"");
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    
    if (diagnosticp == true) {
      /* Subtract 1 because dynprogindices start at +1 and -1 */
      if (this->comp == DYNPROG_MATCH_COMP) {
	if (this->dynprogindex > 0) {
	  fprintf(fp,"%c",(this->dynprogindex-1)%26+'a');
	} else if (this->dynprogindex < 0) {
	  fprintf(fp,"%c",(-this->dynprogindex-1)%26+'A');
	} else {
	  putc(DYNPROG_MATCH_COMP,fp);
	}
      } else if (this->shortexonp == true) {
	putc(DIAGNOSTIC_SHORTEXON_COMP,fp);
      } else {
	putc(this->comp,fp);
      }

    } else if (this->comp == DYNPROG_MATCH_COMP) {
      putc(MATCH_COMP,fp);
    } else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
      putc(AMBIGUOUS_COMP,fp);
#else
      putc(MISMATCH_COMP,fp);
#endif
    } else if (this->comp == SHORTGAP_COMP) {
      putc(INDEL_COMP,fp);
    } else if (this->comp == EXTRAEXON_COMP) {
      putc(INTRONGAP_COMP,fp);
    } else {
      putc(this->comp,fp);
    }
  }

  putc('\n',fp);
  return;
}


static void
print_genomic_sequence (FILE *fp, struct T *ptr, int n, int npairs,
			char *chrstring, Univcoord_T chroffset,
			int margin, int wraplength) {
  struct T *this;
  int i;
  char Buffer[100];

  this = ptr;
  if (chrstring == NULL) {
    sprintf(Buffer,"%u",chroffset+this->genomepos + ONEBASEDP);
  } else {
    sprintf(Buffer,"%s:%u",chrstring,this->genomepos + ONEBASEDP);
  }
  fprintf(fp,"%*s ",margin,Buffer);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    if (this->comp == EXTRAEXON_COMP) {
      putc(INTRONGAP_CHAR,fp);
    } else {
      putc(this->genome,fp);
    }
  }
  putc('\n',fp);
  return;
}

static void
print_genomicalt_sequence (FILE *fp, struct T *ptr, int n, int npairs,
			   char *chrstring, Univcoord_T chroffset,
			   int margin, int wraplength) {
  struct T *this;
  int i;
  char Buffer[100];

  this = ptr;
  if (chrstring == NULL) {
    sprintf(Buffer,"%u",chroffset+this->genomepos + ONEBASEDP);
  } else {
    sprintf(Buffer,"%s:%u",chrstring,this->genomepos + ONEBASEDP);
  }
  fprintf(fp,"%*s ",margin,Buffer);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    if (this->comp == EXTRAEXON_COMP) {
      putc(INTRONGAP_CHAR,fp);
    } else if (this->genomealt == this->genome) {
      putc(' ',fp);
    } else {
      putc(this->genomealt,fp);
    }
  }
  putc('\n',fp);
  return;
}


static int
compute_margin (struct T *start, struct T *end, char *chrstring,
		Univcoord_T chroffset) {
  int margin;
  char Buffer[100];

  if (chrstring == NULL) {
    sprintf(Buffer,"%u",chroffset + start->genomepos + ONEBASEDP);
  } else {
    sprintf(Buffer,"%s:%u",chrstring,start->genomepos + ONEBASEDP);
  }
  margin = (int) strlen(Buffer) + 1;

  if (chrstring == NULL) {
    sprintf(Buffer,"%u",chroffset + end->genomepos + ONEBASEDP);
  } else {
    sprintf(Buffer,"%s:%u",chrstring,end->genomepos + ONEBASEDP);
  }
  if ((int) strlen(Buffer) + 1 > margin) {
    margin = (int) strlen(Buffer) + 1;
  }

  if (margin < DEFAULT_MARGIN) {
    margin = DEFAULT_MARGIN;
  }

  return margin;
}


/*
static char
intron_symbol_rev (char c) {
  switch (c) {
  case '>': return '<';
  case ')': return '(';
  case ']': return '[';
  case '<': return '>';
  case '(': return ')';
  case '[': return ']';
  default: return c;
  }
}
*/

static char complCode[128] = COMPLEMENT_LC;

static struct T *
invert_path (struct T *old, int npairs) {
  struct T *new;
  int i, j;

  new = (struct T *) MALLOC(npairs*sizeof(struct T));
  for (i = 0, j = npairs-1; i < npairs; i++, j--) {
    memcpy(&(new[j]),&(old[i]),sizeof(struct T));
    new[j].comp = complCode[(int) old[i].comp];
  }
  return new;
}

static struct T *
invert_and_revcomp_path (struct T *old, int npairs) {
  struct T *new;
  int i, j;

  new = (struct T *) MALLOC(npairs*sizeof(struct T));
  for (i = 0, j = npairs-1; i < npairs; i++, j--) {
    memcpy(&(new[j]),&(old[i]),sizeof(struct T));
    new[j].cdna = complCode[(int) old[i].cdna];
    new[j].genome = complCode[(int) old[i].genome];
    new[j].genomealt = complCode[(int) old[i].genomealt];
    new[j].comp = complCode[(int) old[i].comp];
  }
  return new;
}


static struct T *
invert_and_revcomp_path_and_coords (struct T *old, int npairs, int querylength) {
  struct T *new;
  int i, j;

  new = (struct T *) MALLOC(npairs*sizeof(struct T));
  for (i = 0, j = npairs-1; i < npairs; i++, j--) {
    memcpy(&(new[j]),&(old[i]),sizeof(struct T));
    new[j].querypos = (querylength - 1) - old[i].querypos;
    new[j].cdna = complCode[(int) old[i].cdna];
    new[j].genome = complCode[(int) old[i].genome];
    new[j].genomealt = complCode[(int) old[i].genomealt];
    new[j].comp = complCode[(int) old[i].comp];
  }
  return new;
}


static void
add_intronlengths (struct T *pairs, int npairs) {
  struct T *this = NULL, *ptr;
  int space, margin, i, j, k, gapstart;
  char intronstring[20], cdnabreak[20], genomicbreak[20], comp;
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;

  i = 0;
  while (i < npairs) {
    /* prev = this; */
    this = &(pairs[i++]);

    if (this->extraexonp == true) {
      /* Don't add any lengths */
    } else if (this->gapp) {
      comp = this->comp;
      gapstart = i-1;
      space = 0;
      while (this->gapp) {
	this = &(pairs[i++]);
	space++;
      }

      if (comp == DUALBREAK_COMP || comp == EXTRAEXON_COMP) {
	sprintf(cdnabreak,"%d",abs(this->querypos - last_querypos)-1);
	if (this->genomepos < last_genomepos) {
	  sprintf(genomicbreak,"%d",last_genomepos - this->genomepos - 1);
	} else {
	  sprintf(genomicbreak,"%d",this->genomepos - last_genomepos - 1);
	}

	margin = (space - strlen(cdnabreak))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < (int) strlen(cdnabreak); k++) {
	  ptr = &(pairs[j++]);
	  ptr->cdna = cdnabreak[k];
	}

	margin = (space - strlen(genomicbreak))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < (int) strlen(genomicbreak); k++) {
	  ptr = &(pairs[j++]);
	  ptr->genome = genomicbreak[k];
	  /* ptr->genomealt = ' '; */
	}

      } else {			/* Intron */
	if (this->genomepos < last_genomepos) {
	  sprintf(intronstring,"%d",last_genomepos - this->genomepos - 1);
	} else {
	  sprintf(intronstring,"%d",this->genomepos - last_genomepos - 1);
	}
	margin = (space - strlen(intronstring))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < (int) strlen(intronstring); k++) {
	  ptr = &(pairs[j++]);
	  ptr->cdna = intronstring[k];
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }	
  return;
}


/* Needed to recompute translation_length in parts of chimeras */
int
Pair_translation_length (struct T *pairs, int npairs) {
  int translation_length = 0;
  int i;

  for (i = 0; i < npairs; i++) {
    if (pairs[i].aa_e == ' ') {
    } else if (pairs[i].aa_e == '*') {
    } else {
      translation_length++;
    }
  }
  return translation_length;
}


void
Pair_print_continuous (FILE *fp, struct T *pairs, int npairs, bool watsonp,
		       bool diagnosticp, bool genomefirstp, int invertmode,
		       bool nointronlenp) {
  T this;
  struct T *save = NULL, *ptr;
  int n = 0;

  if (watsonp == true) {
    ptr = pairs;
  } else if (invertmode == 0) {
    ptr = pairs;
  } else if (invertmode == 1) {
    save = ptr = invert_path(pairs,npairs);
  } else if (invertmode == 2) {
    save = ptr = invert_and_revcomp_path(pairs,npairs);
  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }
  if (nointronlenp == false) {
    add_intronlengths(ptr,npairs);
  }

  if (genomefirstp == true) {
    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putc(this->genome,fp);
    }
    putc('\n',fp);

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      if (this->comp == MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	putc(AMBIGUOUS_COMP,fp);
#else
	putc(MISMATCH_COMP,fp);
#endif
      } else {
	putc(this->comp,fp);
      }
    }
    putc('\n',fp);

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putc(this->cdna,fp);
    }
    putc('\n',fp);

  } else {
    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putc(this->cdna,fp);
    }
    putc('\n',fp);

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      if (this->comp == MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	putc(AMBIGUOUS_COMP,fp);
#else
	putc(MISMATCH_COMP,fp);
#endif
      } else {
	putc(this->comp,fp);
      }
    }
    putc('\n',fp);

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putc(this->genome,fp);
    }
    putc('\n',fp);
  }

  if (save != NULL) {
    FREE(save);
  }
  return;
}  



void
Pair_print_continuous_byexon (FILE *fp, struct T *pairs, int npairs, bool watsonp, bool diagnosticp, int invertmode) {
  T this;
  struct T *save = NULL, *ptr;
  int i = 0, j;

  if (watsonp == true) {
    ptr = pairs;
  } else if (invertmode == 0) {
    ptr = pairs;
  } else if (invertmode == 1) {
    save = ptr = invert_path(pairs,npairs);
  } else if (invertmode == 2) {
    save = ptr = invert_and_revcomp_path(pairs,npairs);
  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }

  ptr = pairs;
  while (i < npairs) {
    j = i;
    this = ptr;

    while (j < npairs && this->gapp == false) {
      putc(this->genome,fp);
      this++;
      j++;
    }
    putc('\n',fp);

    j = i;
    this = ptr;
    while (j < npairs && this->gapp == false) {
      if (this->comp == MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	putc(AMBIGUOUS_COMP,fp);
#else
	putc(MISMATCH_COMP,fp);
#endif
      } else {
	putc(this->comp,fp);
      }
      this++;
      j++;
    }
    putc('\n',fp);

    j = i;
    this = ptr;
    while (j < npairs && this->gapp == false) {
      putc(this->cdna,fp);
      this++;
      j++;
    }
    fprintf(fp,"\n\n");

    i = j;
    while (i < npairs && this->gapp == true) {
      this++;
      i++;
    }
    ptr = this;
  }

  if (save != NULL) {
    FREE(save);
  }
  return;
}  


void
Pair_print_alignment (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
		      Univcoord_T chroffset, Univ_IIT_T chromosome_iit, bool watsonp,
		      bool diagnosticp, int invertmode, bool nointronlenp,
		      int wraplength) {
  struct T *save = NULL, *ptr;
  int n = 0, i;
  char *chrstring = NULL;
  int margin;

  if (watsonp == true) {
    ptr = pairs;

  } else if (invertmode == 0) {
    /* Given cDNA sequence, use minus genome strand */
    ptr = pairs;

  } else if (invertmode == 1) {
    /* Invert cDNA sequence, use minus genome strand */
    save = ptr = invert_path(pairs,npairs);

  } else if (invertmode == 2) {
    /* Invert cDNA sequence, use plus genome strand */
    save = ptr = invert_and_revcomp_path(pairs,npairs);

  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }
  
  if (nointronlenp == false) {
    add_intronlengths(ptr,npairs);
  }
  if (chrnum != 0) {
    if (invertmode == 2) {
      chrstring = Chrnum_to_string(chrnum,chromosome_iit);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  margin = compute_margin(&(pairs[0]),&(pairs[npairs-1]),chrstring,chroffset);

  while (n < npairs) {
    print_top_ruler(fp,n,npairs,margin,wraplength);
    print_peptide(fp,ptr,n,npairs,margin,wraplength,/*genomep*/true);
    if (snps_p) {
      print_genomicalt_sequence(fp,ptr,n,npairs,chrstring,
				chroffset,margin,wraplength);
    }
    print_genomic_sequence(fp,ptr,n,npairs,chrstring,
			   chroffset,margin,wraplength);
    print_alignment(fp,ptr,n,npairs,diagnosticp,margin,wraplength);
    print_cdna_sequence(fp,ptr,n,npairs,margin,wraplength);
    print_peptide(fp,ptr,n,npairs,margin,wraplength,/*genomep*/false);
    putc('\n',fp);
    for (i = 0; n < npairs && i < wraplength; n++, i++) {
      ptr++;
    }
  }
  if (chrstring != NULL) {
    FREE(chrstring);
  }
  if (save != NULL) {
    FREE(save);
  }
  return;
}  

void
Pair_print_pathsummary (FILE *fp, int pathnum, T start, T end, Chrnum_T chrnum,
			Univcoord_T chroffset, Univ_IIT_T chromosome_iit, bool referencealignp, 
			IIT_T altstrain_iit, char *strain, Univ_IIT_T contig_iit, char *dbversion,
			int querylength_given, int skiplength, int trim_start, int trim_end,
			int nexons, int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels, int goodness,
			bool watsonp, int cdna_direction,
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend, bool maponlyp,
			bool diagnosticp, int stage2_source, int stage2_indexsize) {
  int querypos1, querypos2, den;
  double fracidentity, coverage, trimmed_coverage;
  Univcoord_T position1, position2;
  Chrpos_T chrpos1, chrpos2;
  char *refstrain, *comma1, *comma2, *chr;

  querypos1 = start->querypos;
  querypos2 = end->querypos;

  fprintf(fp,"  Path %d: ",pathnum);
  fprintf(fp,"query %d%s%d (%d bp) => ",
	 querypos1 + ONEBASEDP,SEPARATOR,querypos2 + ONEBASEDP,querypos2-querypos1+1);

  chrpos1 = start->genomepos;
  chrpos2 = end->genomepos;

  comma1 = Genomicpos_commafmt(chrpos1 + ONEBASEDP);
  comma2 = Genomicpos_commafmt(chrpos2 + ONEBASEDP);
  if (chrnum == 0) {
    if (watsonp) {
      fprintf(fp,"genome %s%s%s (%d bp)\n",
	     comma1,SEPARATOR,comma2,chrpos2-chrpos1+1);
    } else {
      fprintf(fp,"genome %s%s%s (%d bp)\n",
	     comma1,SEPARATOR,comma2,chrpos2-chrpos1-1);
    }
  } else {
    chr = Chrnum_to_string(chrnum,chromosome_iit);
    if (watsonp) {
      fprintf(fp,"genome %s:%s%s%s (%d bp)\n",chr,comma1,SEPARATOR,comma2,chrpos2-chrpos1+1);
    } else {
      fprintf(fp,"genome %s:%s%s%s (%d bp)\n",chr,comma1,SEPARATOR,comma2,chrpos2-chrpos1-1);
    }
    FREE(chr);
  }
  FREE(comma2);
  FREE(comma1);

  if (maponlyp == false) {

    if (diagnosticp == true) {
      /* fprintf(fp,"    Stage 2 diag runtime: %.3f sec\n",stage2_diag_runtime); */
      /* fprintf(fp,"    Stage 2 align runtime: %.3f sec\n",stage2_align_runtime); */
      fprintf(fp,"    Stage 2 source: %d\n",stage2_source);
      fprintf(fp,"    Stage 2 indexsize: %d\n",stage2_indexsize);
      /* fprintf(fp,"    Stage 3 runtime: %.3f sec\n",stage3_runtime); */
      /* fprintf(fp,"    Stage 3 defectrate: %f\n",stage3_defectrate); */
      fprintf(fp,"    Goodness: %d\n",goodness);
    }

    fprintf(fp,"    cDNA direction: ");
    if (cdna_direction > 0) {
      fprintf(fp,"sense\n");
    } else if (cdna_direction < 0) {
      fprintf(fp,"antisense\n");
    } else {
      fprintf(fp,"indeterminate\n");
    }
  }

  if (altstrain_iit != NULL) {
    if (strain == NULL) {
      refstrain = IIT_typestring(altstrain_iit,/*straintype*/0);
      if (refstrain[0] == '\0') {
	/* Backward compatibility with old altstrain_iit */
	fprintf(fp,"    Strain: reference\n");
      } else {
	fprintf(fp,"    Strain: %s (reference)\n",refstrain);
      }
    } else {
      fprintf(fp,"    Strain: %s\n",strain);
    }
  }

  position1 = chroffset + chrpos1;
  position2 = chroffset + chrpos2;
  comma1 = Genomicpos_commafmt(position1 + ONEBASEDP);
  comma2 = Genomicpos_commafmt(position2 + ONEBASEDP);
  if (dbversion == NULL) {
    fprintf(fp,"    Genomic pos: %s%s%s",comma1,SEPARATOR,comma2);
  } else {
    fprintf(fp,"    Genomic pos: %s:%s%s%s",dbversion,comma1,SEPARATOR,comma2);
  }
  if (chrpos1 <= chrpos2) {
    fprintf(fp," (+ strand)\n");
  } else {
    fprintf(fp," (- strand)\n");
  }
  FREE(comma2);
  FREE(comma1);

  if (contig_iit != NULL) {
    if (position1 <= position2) {
      Segmentpos_print_accessions(fp,contig_iit,position1,position2,referencealignp,strain);
    } else {
      Segmentpos_print_accessions(fp,contig_iit,position2,position1,referencealignp,strain);
    }
  }
    
  if (maponlyp == false) {
    fprintf(fp,"    Number of exons: %d\n",nexons);

#ifdef PMAP
    coverage = (double) (querypos2 - querypos1 + 1)/(double) (3*(querylength_given + skiplength));
    /* coverage = (double) (matches + mismatches + qindels)/(double) (3*(querylength_given + skiplength)); */

    /* Can have coverage greater than given querylength because of added '*' at end */
    if (coverage > 1.0) {
      coverage = 1.0;
    }
#else
    /* coverage = (double) (matches + mismatches + qindels)/(double) (querylength_given + skiplength); */
    coverage = (double) (querypos2 - querypos1 + 1)/(double) (querylength_given + skiplength);
#endif
    fprintf(fp,"    Coverage: %.1f",((double) rint(1000.0*coverage))/10.0);
#ifdef PMAP
    fprintf(fp," (query length: %d aa)\n",querylength_given);
#else
    fprintf(fp," (query length: %d bp)\n",querylength_given);
    if (querypos2 + 1 > trim_end) {
      trim_end = querypos2 + 1;
    }
    if (querypos1 < trim_start) {
      trim_start = querypos1;
    }

    trimmed_coverage = (double) (querypos2 - querypos1 + 1)/(double) (trim_end - trim_start + skiplength);
    fprintf(fp,"    Trimmed coverage: %.1f",((double) rint(1000.0*trimmed_coverage))/10.0);
    fprintf(fp," (trimmed length: %d bp, trimmed region: %d..%d)",
	   trim_end-trim_start,trim_start+ONEBASEDP,trim_end-1+ONEBASEDP);
    putc('\n',fp);
#endif

    if ((den = matches + mismatches + qindels + tindels) == 0) {
      fracidentity = 1.0;
    } else {
      fracidentity = (double) matches/(double) den;
    }

    /* The definition of indels here should be consistent with Stage3_indels */
    fprintf(fp,"    Percent identity: %.1f (%d matches, %d mismatches, %d indels, %d unknowns)\n",
	   ((double) rint(1000.0*fracidentity))/10.0,matches,mismatches,qindels+tindels,unknowns);
    if (qindels + tindels > 0) {
      fprintf(fp,"    Non-intron gaps: %d openings, %d bases in cdna; %d openings, %d bases in genome\n",
	     qopens,qindels,topens,tindels);
    } 

#ifndef PMAP
    if (translation_length > 0) {
      if (cdna_direction >= 0) {
	fprintf(fp,"    Translation: %d..%d (%d aa)\n",
	       translation_start+ONEBASEDP,translation_end+ONEBASEDP,translation_length);
      } else {
	fprintf(fp,"    Translation: %d..%d (%d aa)\n",
	       translation_end+ONEBASEDP,translation_start+ONEBASEDP,translation_length);
      }
    } else if (relaastart > 0) {
      if (relaastart < relaaend) {
	fprintf(fp,"    Protein coords: %d..%d\n",relaastart,relaaend);
      } else {
	fprintf(fp,"    Protein coords: %d..%d\n",relaaend,relaastart);
      }
    }
#endif

    /* fprintf(fp,"    Defect rate (percent): %.1f\n",defect_rate*100.0); */

    /* putc('\n',fp); -- Done by caller */
  }

  return;
}


void
Pair_print_coordinates (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
			Univcoord_T chroffset, Univ_IIT_T chromosome_iit,
			bool watsonp, int invertmode) {
  T this;
  struct T *save = NULL;
  int i;
  char *chrstring = NULL;

  if (watsonp == true) {
    /* ptr = pairs; */

  } else if (invertmode == 0) {
    /* Given cDNA sequence, use minus genome strand */
    /* ptr = pairs; */

  } else if (invertmode == 1) {
    /* Invert cDNA sequence, use minus genome strand */
    save = invert_path(pairs,npairs);

  } else if (invertmode == 2) {
    /* Invert cDNA sequence, use plus genome strand */
    save = invert_and_revcomp_path(pairs,npairs);

  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }
  
  if (chrnum != 0) {
    if (invertmode == 2) {
      chrstring = Chrnum_to_string(chrnum,chromosome_iit);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  for (i = 0; i < npairs; i++) {
    this = pairs++;
    if (this->gapp == false) {
#ifdef DEBUG5
      fprintf(fp,"%d %d %c\t",this->aapos,this->aaphase_e,this->aa_e);
#else
      if (this->aaphase_e != 0) {
	fprintf(fp,"%d\t",this->aapos);
      } else {
	fprintf(fp,"%d %c\t",this->aapos,this->aa_e);
      }
#endif
      fprintf(fp,"%d %c\t",this->querypos + ONEBASEDP,this->cdna);
      if (chrstring == NULL) {
	fprintf(fp,"%u %u %c",this->genomepos + ONEBASEDP,
		chroffset + this->genomepos + ONEBASEDP,
		this->genome);
      } else {
	fprintf(fp,"%s:%u %u %c",chrstring,
		this->genomepos + ONEBASEDP,
		chroffset + this->genomepos + ONEBASEDP,
		this->genome);
      }
      if (this->genomealt != this->genome) {
	fprintf(fp," %c",this->genomealt);
      }

#ifdef DEBUG5
      fprintf(fp,"\t%d %c",this->aaphase_g,this->aa_g);
#else
      if (this->aaphase_g != 0) {
	fprintf(fp,"\t");
      } else {
	fprintf(fp,"\t%c",this->aa_g);
      }
#endif
      putc('\n',fp);
    }
  }

  if (chrstring != NULL) {
    FREE(chrstring);
  }
  if (save != NULL) {
    FREE(save);
  }
  return;
}  


void
Pair_dump_one (T this, bool zerobasedp) {

  debug1(printf("%p ",this));

  if (this->gapp == true && this->extraexonp == false) {
    printf("*** Gap: queryjump = %d, genomejump = %d, type: ",this->queryjump,this->genomejump);
    switch (this->comp) {
    case FWD_CANONICAL_INTRON_COMP: printf("> GT-AG"); break;
    case FWD_GCAG_INTRON_COMP: printf(") GC-AG"); break;
    case FWD_ATAC_INTRON_COMP: printf("] AT-AC"); break;
    case REV_ATAC_INTRON_COMP: printf("[ AT-AC"); break;
    case REV_GCAG_INTRON_COMP: printf("( GC-AG"); break;
    case REV_CANONICAL_INTRON_COMP: printf("< GT-AG"); break;
    case SHORTGAP_COMP: printf("~ shortgap"); break;
    case NONINTRON_COMP: printf("= nonintron"); break;
    default: printf("? unknown"); break;
    }

    if (this->knowngapp == true) {
      printf(" known");
    }

    printf(" donor:%f acceptor:%f",this->donor_prob,this->acceptor_prob);
    printf(" ***");

  } else {
    printf("%d %d %c ",
	   this->querypos + !zerobasedp,this->genomepos + !zerobasedp,this->cdna);

    /* Subtract 1 because dynprogindices start at +1 and -1 */
    if (this->dynprogindex > 0) {
      printf("%c%c",this->comp,(this->dynprogindex-1)%26+'a');
    } else if (this->dynprogindex < 0) {
      printf("%c%c",this->comp,(-this->dynprogindex-1)%26+'A');
    } else {
      putchar(this->comp);
    }
    printf(" %c",this->genome);
    if (this->genomealt != this->genome) {
      printf(" alt:%c",this->genomealt);
    }
  }

  if (this->protectedp == true) {
    printf(" protected");
  }

  if (this->disallowedp == true) {
    printf(" disallowed");
  }

  if (this->shortexonp == true) {
    printf(" shortexon");
  }

  if (this->state == BAD) {
    printf(" bad");
  }

  return;
}


/* Useful for debugging */
void
Pair_dump_list (List_T pairs, bool zerobasedp) {
  T this;
  List_T p;

  printf("***Start of list***\n");
  for (p = pairs; p != NULL; p = List_next(p)) {
    this = List_head(p);
    Pair_dump_one(this,zerobasedp);
    printf("\n");
  }
  printf("***End of list***\n");
  return;
}  

void
Pair_dump_array (struct T *pairs, int npairs, bool zerobasedp) {
  struct T *this;
  int i;

  for (i = 0; i < npairs; i++) {
    this = pairs++;
    printf("%d: %d %d %d %c ",
	   i,this->querypos + !zerobasedp,this->genomepos + !zerobasedp,this->aapos,
	   this->cdna);

    /* Subtract 1 because dynprogindices start at +1 and -1 */
    if (this->dynprogindex > 0) {
      printf("%c%c",this->comp,(this->dynprogindex-1)%26+'a');
    } else if (this->dynprogindex < 0) {
      printf("%c%c",this->comp,(-this->dynprogindex-1)%26+'A');
    } else {
      putchar(this->comp);
    }
    printf(" %c",this->genome);
    if (this->genomealt != this->genome) {
      printf(" alt:%c",this->genomealt);
    }

    if (this->aaphase_g == 0 || this->aaphase_e == 0) {
      printf(" => %c %c",this->aa_g,this->aa_e);
    }
    printf("\n");
  }
  return;
}  


Chrpos_T
Pair_genomicpos (struct T *pairs, int npairs, int querypos, bool headp) {
  struct T *this;
  int i;

  if (headp == true) {
    for (i = 0; i < npairs; i++) {
      this = pairs++;
      if (this->querypos == querypos) {
	return this->genomepos;
      } else if (this->querypos > querypos) {
	return 0;
      }
    }
  } else {
    pairs += npairs;
    for (i = npairs-1; i >= 0; --i) {
      this = --pairs;
      if (this->querypos == querypos) {
	return this->genomepos;
      } else if (this->querypos < querypos) {
	return 0;
      }
    }
  }
  return 0;
}

int
Pair_codon_changepos (struct T *pairs, int npairs, int aapos, int cdna_direction) {
  struct T *this, *start, *end;
  int changepos = 0, i, ngenome = 0, ncdna = 0;

  i = 0;
  this = pairs;
  while (i < npairs && this->aapos != aapos) {
    this++;
    i++;
  }
  start = this;

  while (i < npairs && (ngenome < 3 || ncdna < 3)) {
    if (this->gapp == false) {
      if (this->genome != ' ') {
	ngenome++;
      }
      if (this->cdna != ' ') {
	ncdna++;
      }
    }
    this++;
    i++;
  }
  end = --this;

  if (cdna_direction < 0) {
    for (this = end; this >= start; --this) {
      if (this->gapp == true) {
      } else if (this->genome == ' ') {
      } else if (this->cdna == ' ') {
      } else if (this->genome != this->cdna) {
	return changepos;
      } else {
	changepos++;
      }
    }
  } else {
    for (this = start; this <= end; this++) {
      if (this->gapp == true) {
      } else if (this->genome == ' ') {
      } else if (this->cdna == ' ') {
      } else if (this->genome != this->cdna) {
	return changepos;
      } else {
	changepos++;
      }
    }
  }

  return changepos;
}  


void
Pair_check_list (List_T pairs) {
  T this;
  List_T p;
  int prev_querypos;

  if (pairs == NULL) {
    return;
  } else {
    this = List_head(pairs);
    prev_querypos = this->querypos;
    /* prev_genomepos = this->genomepos; */

    for (p = List_next(pairs); p != NULL; p = List_next(p)) {
      this = List_head(p);
      if (this->gapp == false) {
	if (this->querypos > prev_querypos) {
	  printf("Problem at querypos %d\n",this->querypos);
	}
#if 0
	/* No longer a valid check after genomepos converted to chrpos */
	if (this->genomepos > prev_genomepos) {
	  printf("Problem at genomepos %d\n",this->genomepos);
	}
#endif
	prev_querypos = this->querypos;
	/* prev_genomepos = this->genomepos; */
      }
    }
  }
  return;
}  


bool
Pair_check_array (struct T *pairs, int npairs) {
  bool result = false;
  struct T *this;
  int prev_querypos;
  int i;

  if (npairs == 0) {
    return false;
  } else {
    this = pairs++;
    prev_querypos = this->querypos;
    /* prev_genomepos = this->genomepos; */

    for (i = 1; i < npairs; i++) {
      this = pairs++;
      if (this->querypos < prev_querypos) {
	fprintf(stderr,"Problem at querypos %d\n",this->querypos);
	result = true;
      } else if (this->querypos - prev_querypos > 1) {
	/* Could be the result of a dual break */
	fprintf(stderr,"Jump at querypos %d\n",this->querypos);
	result = false;
      }
#if 0
      /* No longer a valid check after genomepos converted to chrpos */
      if (this->genomepos < prev_genomepos) {
	fprintf(stderr,"Problem at genomepos %d\n",this->genomepos);
	result = true;
      }
#endif
      prev_querypos = this->querypos;
      /* prev_genomepos = this->genomepos; */
    }
  }
  return result;
}  


/* Called by output thread for --merge-overlap feature.  Modeled after Substring_convert_to_pairs. */
List_T
Pair_convert_array_to_pairs (List_T pairs, struct T *pairarray, int npairs, bool plusp, int querylength,
			     int clipdir, int hardclip_low, int hardclip_high, bool first_read_p, int queryseq_offset) {
  T pair;
  int querystart, queryend, i;

  if (plusp == true) {
    querystart = hardclip_low;
    queryend = querylength - hardclip_high;

  } else {
    querystart = hardclip_high;
    queryend = querylength - hardclip_low;
  }

  for (i = 0; i < npairs; i++) {
    pair = &(pairarray[i]);
    if (pair->querypos >= querystart && pair->querypos < queryend) {
      pairs = List_push_out(pairs,(void *) Pair_new_out(pair->querypos + queryseq_offset,/*genomepos*/pair->genomepos,
							pair->cdna,pair->comp,pair->genome));
    }
  }
      
  return pairs;
}



#if 0
static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[(int) j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}
#endif

static void
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

  return;
}


static double
donor_score (Univcoord_T genomicpos, Univcoord_T chroffset, bool revcomp, Genome_T genome,
	     Univ_IIT_T chromosome_iit) {
  Univcoord_T left;
  Chrnum_T chrnum;
  int nunknowns;
  char gbuffer[MAXENT_MAXLENGTH];
  Genomecomp_T *genome_blocks;

  if (revcomp == false) {
    if ((genome_blocks = Genome_blocks(genome)) != NULL) {
      /* Add 1 to get from exon end to intron start */
      return Maxent_hr_donor_prob(genomicpos + 1U,chroffset);
    } else {
      left = genomicpos + 1 - DONOR_MODEL_LEFT_MARGIN; /* Add 1 to get from exon end to intron start */
      Genome_fill_buffer(&chrnum,&nunknowns,genome,left,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer,chromosome_iit);
#if 0
      printf("\n");
      printf("%s donor truestrand:+ left:%u\n",gbuffer,left);
      printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,"");
#endif
      return Maxent_donor_prob(gbuffer);
    }

  } else {
    if ((genome_blocks = Genome_blocks(genome)) != NULL) {
      return Maxent_hr_antidonor_prob(genomicpos,chroffset);
    } else {
      left = genomicpos - DONOR_MODEL_RIGHT_MARGIN - 1;
      Genome_fill_buffer(&chrnum,&nunknowns,genome,left,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer,chromosome_iit);
      make_complement_inplace(gbuffer,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1);
#if 0
      printf("\n");
      printf("%s donor truestrand:- left:%u\n",gbuffer,left);
      printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,"");
#endif
      return Maxent_donor_prob(gbuffer);
    }
  }
}


static double
acceptor_score (Univcoord_T genomicpos, Univcoord_T chroffset, bool revcomp, Genome_T genome,
		Univ_IIT_T chromosome_iit) {
  Univcoord_T left;
  Chrnum_T chrnum;
  int nunknowns;
  char gbuffer[MAXENT_MAXLENGTH];
  Genomecomp_T *genome_blocks;

  if (revcomp == false) {
    /* sense on plus strand, or antisense on minus strand */
    if ((genome_blocks = Genome_blocks(genome)) != NULL) {
      return Maxent_hr_acceptor_prob(genomicpos,chroffset);
    } else {
      left = genomicpos - ACCEPTOR_MODEL_LEFT_MARGIN;
      Genome_fill_buffer(&chrnum,&nunknowns,genome,left,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer,chromosome_iit);
#if 0
      printf("\n");
      printf("%s acceptor truestrand:+ left:%u\n",gbuffer,left);
      printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,"");
#endif
      return Maxent_acceptor_prob(gbuffer);
    }

  } else {
    if ((genome_blocks = Genome_blocks(genome)) != NULL) {
      /* Add 1 to get from exon end to intron start */
      return Maxent_hr_antiacceptor_prob(genomicpos + 1U,chroffset);
    } else {
      left = genomicpos - ACCEPTOR_MODEL_RIGHT_MARGIN;
      Genome_fill_buffer(&chrnum,&nunknowns,genome,left,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer,chromosome_iit);
      make_complement_inplace(gbuffer,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1);
#if 0
      printf("\n");
      printf("%s acceptor truestrand:- left:%u\n",gbuffer,left);
      printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,"");
#endif
      return Maxent_acceptor_prob(gbuffer);
    }
  }
}



static bool
unknown_base (char c) {
  switch (c) {
  case 'A': case 'C': case 'G': case 'T': case 'U':
  case 'a': case 'c': case 'g': case 't': case 'u': return false;
  default: return true;
  }
}
    
void
Pair_print_exonsummary (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
			Univcoord_T chroffset, Genome_T genome, Univ_IIT_T chromosome_iit,
			bool watsonp, int cdna_direction, bool genomefirstp, int invertmode) {
  bool in_exon = false;
  struct T *save = NULL, *ptr, *this = NULL;
  int exon_querystart = -1, exon_queryend;
  Chrpos_T exon_genomestart = -1U, exon_genomeend, intron_start, intron_end;
  int num = 0, den = 0, i;
  char *chrstring = NULL;
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;


  if (watsonp == true) {
    ptr = pairs;
  } else if (invertmode == 0) {
    ptr = pairs;
  } else if (invertmode == 1) {
    save = ptr = invert_path(pairs,npairs);
  } else if (invertmode == 2) {
    save = ptr = invert_and_revcomp_path(pairs,npairs);
  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }
  
  if (chrnum != 0) {
    if (invertmode == 2) {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,/*watsonp*/true);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  debug(Pair_dump_array(pairs,npairs,/*zerobasedp*/true));

  for (i = 0; i < npairs; i++) {
    /* prev = this; */
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + ONEBASEDP;
	exon_genomeend = last_genomepos + ONEBASEDP;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
	if (genomefirstp == true) {
	  fprintf(fp,"    ");
	  if (chrnum == 0) {
	    fprintf(fp,"%u-%u",chroffset+exon_genomestart,chroffset+exon_genomeend);
	  } else {
	    fprintf(fp,"%s:%d-%d",chrstring,exon_genomestart,exon_genomeend);
	  }
	  fprintf(fp,"  (%d-%d)",exon_querystart,exon_queryend);
	} else {
	  fprintf(fp,"    %d-%d",exon_querystart,exon_queryend);
	  fprintf(fp,"  ");
	  if (chrnum == 0) {
	    fprintf(fp,"(%u-%u)",chroffset+exon_genomestart,chroffset+exon_genomeend);
	  } else {
	    fprintf(fp,"(%s:%d-%d)",chrstring,exon_genomestart,exon_genomeend);
	  }
	}
	if (den == 0) {
	  fprintf(fp,"   %d%%",100);
	} else {
	  fprintf(fp,"   %d%%",(int) floor(100.0*(double) num/(double) den));
	}
	if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	  fprintf(fp," ->");
	  /* sensep = true; */
	} else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	  fprintf(fp," <-");
	  /* sensep = false; */
	} else if (this->comp == FWD_GCAG_INTRON_COMP) {
	  fprintf(fp," -)");
	  /* sensep = true; */
	} else if (this->comp == REV_GCAG_INTRON_COMP) {
	  fprintf(fp," (-");
	  /* sensep = false; */
	} else if (this->comp == FWD_ATAC_INTRON_COMP) {
	  fprintf(fp," -]");
	  /* sensep = true; */
	} else if (this->comp == REV_ATAC_INTRON_COMP) {
	  fprintf(fp," [-");
	  /* sensep = false; */
	} else if (this->comp == NONINTRON_COMP) {
	  fprintf(fp," ==");
	  /* sensep = true; */
	} else {
	  fprintf(fp," ##");
	  /* sensep = true; */
	}
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + ONEBASEDP;
	exon_genomestart = this->genomepos + ONEBASEDP;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}
	if (i > 0) {
	  if (intron_end > intron_start) {
	    fprintf(fp,"   ...%d...",intron_end - intron_start + 1);
	  } else {
	    fprintf(fp,"   ...%d...",intron_start - intron_end + 1);
	  }

	  if (exon_querystart > exon_queryend + 1) {
	    fprintf(fp,"   ***query_skip:%d***",exon_querystart-(exon_queryend+1));
	  }

	  if (genome != NULL) {
	    if (cdna_direction >= 0) {
	      fprintf(fp,"  %.3f, %.3f",
		      donor_score(chroffset+exon_genomeend-1,chroffset,!watsonp,genome,chromosome_iit),
		      acceptor_score(chroffset+exon_genomestart-1,chroffset,!watsonp,genome,chromosome_iit));
	    } else {
	      fprintf(fp,"  %.3f, %.3f",
		      acceptor_score(chroffset+exon_genomeend-1,chroffset,watsonp,genome,chromosome_iit),
		      donor_score(chroffset+exon_genomestart-1,chroffset,watsonp,genome,chromosome_iit));
	    }
	  }

	  putc('\n',fp);
	}
	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Previously not counted in numerator or denominator */
	den++;
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	/* Comp must be a space */
	/* Don't count in numerator or denominator */
#endif
      } else {
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	  num++;
	} else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	  num++;
#else
	  den--;
#endif
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  /* prev = this; */
  exon_queryend = last_querypos + ONEBASEDP;
  exon_genomeend = last_genomepos + ONEBASEDP;
  if (genomefirstp == true) {
    fprintf(fp,"    ");
    if (chrnum == 0) {
      fprintf(fp,"%u-%u",chroffset+exon_genomestart,chroffset+exon_genomeend);
    } else {
      fprintf(fp,"%s:%d-%d",chrstring,exon_genomestart,exon_genomeend);
    }
    fprintf(fp,"  (%d-%d)",exon_querystart,exon_queryend);
  } else {
    fprintf(fp,"    %d-%d",exon_querystart,exon_queryend);
    fprintf(fp,"  ");
    if (chrnum == 0) {
      fprintf(fp,"(%u-%u)",chroffset+exon_genomestart,chroffset+exon_genomeend);
    } else {
      fprintf(fp,"(%s:%d-%d)",chrstring,exon_genomestart,exon_genomeend);
    }
  }
  if (den == 0) {
    fprintf(fp,"   %d%%",100);
  } else {
    fprintf(fp,"   %d%%",(int) floor(100.0*(double) num/(double) den));
  }
  fprintf(fp,"\n\n");

  if (chrstring != NULL) {
    FREE(chrstring);
  }
  if (save != NULL) {
    FREE(save);
  }

  return;
}

static void
tokens_free (List_T *tokens) {
  List_T p;
  char *token;

  for (p = *tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FREE(token);
  }
  List_free(&(*tokens));

  return;
}


/* Tokens used by compressed and gff3 formats */

static void
print_tokens_compressed (FILE *fp, List_T tokens) {
  List_T p;
  int tokencount = 1;
  char *token, *lasttoken = NULL;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    if (lasttoken == NULL) {
      fprintf(fp,"\t%s",token);
      lasttoken = token;
    } else if (!strcmp(token,lasttoken)) {
      tokencount++;
    } else {
      if (tokencount > 1) {
	fprintf(fp,"!%d",tokencount);
      }
      fprintf(fp," %s",token);
      lasttoken = token;
      tokencount = 1;
    }
  }
  if (tokencount > 1) {
    fprintf(fp,"!%d",tokencount);
  }

  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FREE(token);
  }

  return;
}

static void
print_tokens_gff3 (FILE *fp, List_T tokens) {
  List_T p;
  char *token;
  
  if (tokens != NULL) {
    p = tokens;
    token = (char *) List_head(p);
    fprintf(fp,"%s",token);

    for (p = List_next(p); p != NULL; p = List_next(p)) {
      token = (char *) List_head(p);
      fprintf(fp," %s",token);
    }
  }

  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FREE(token);
  }

  return;
}

static List_T
push_token (List_T tokens, char *token) {
  char *copy;

  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);
  return List_push(tokens,(void *) copy);
}


/* Definition of GFF3 format is at http://song.sourceforge.net/gff3.shtml */

static void
print_gff3_gene (FILE *fp, int pathnum, char *sourcename, char *accession, char *chrstring, Chrpos_T start_genomepos, 
		 Chrpos_T end_genomepos, bool watsonp, int cdna_direction) {

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"gene\t");		/* 3: type */

  if (start_genomepos < end_genomepos) {
    fprintf(fp,"%u\t%u\t",start_genomepos,end_genomepos); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",end_genomepos,start_genomepos); /* 4,5: start, end */
  }

  fprintf(fp,".\t");		/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      fprintf(fp,"+\t");
    } else {
      fprintf(fp,"-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      fprintf(fp,"-\t");		/* 7: strand */
    } else {
      fprintf(fp,"+\t");
    }
  }

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.path%d;Name=%s\n",accession,pathnum,accession);

  return;
}

static void
print_gff3_mrna (FILE *fp, int pathnum, T start, T end,
		 char *sourcename, char *accession, char *chrstring, Chrpos_T start_genomepos, 
		 Chrpos_T end_genomepos, int querylength_given, int skiplength,
		 int matches, int mismatches, int qindels, int tindels, int unknowns,
		 bool watsonp, int cdna_direction) {
  int den;
  int querypos1, querypos2;
  double coverage, fracidentity;

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"mRNA\t");		/* 3: type */
  if (start_genomepos < end_genomepos) {
    fprintf(fp,"%u\t%u\t",start_genomepos,end_genomepos); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",end_genomepos,start_genomepos); /* 4,5: start, end */
  }

  fprintf(fp,".\t");		/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      fprintf(fp,"+\t");
    } else {
      fprintf(fp,"-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      fprintf(fp,"-\t");		/* 7: strand */
    } else {
      fprintf(fp,"+\t");
    }
  }

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.mrna%d;Name=%s;Parent=%s.path%d;",
	 accession,pathnum,accession,accession,pathnum);

  querypos1 = start->querypos;
  querypos2 = end->querypos;

#ifdef PMAP
  coverage = (double) (querypos2 - querypos1 + 1)/(double) (3*(querylength_given + skiplength));
  /* Can have coverage greater than given querylength because of added '*' at end */
  if (coverage > 1.0) {
    coverage = 1.0;
  }
#else
  coverage = (double) (querypos2 - querypos1 + 1)/(double) (querylength_given + skiplength);
#endif
  fprintf(fp,"coverage=%.1f;",((double) rint(1000.0*coverage))/10.0);

  if ((den = matches + mismatches + qindels + tindels) == 0) {
    fracidentity = 1.0;
  } else {
    fracidentity = (double) matches/(double) den;
  }
  fprintf(fp,"identity=%.1f;",((double) rint(1000.0*fracidentity))/10.0);
  fprintf(fp,"matches=%d;mismatches=%d;indels=%d;unknowns=%d",
	  matches,mismatches,qindels+tindels,unknowns);

  putc('\n',fp);

  return;
}


static void
print_gff3_exon (FILE *fp, int exonno, int pathnum, char *sourcename, char *accession, char *chrstring,
		 int exon_genomestart, int exon_genomeend,
		 int exon_querystart, int exon_queryend, bool watsonp, int cdna_direction,
		 int pctidentity) {

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"exon\t");		/* 3: type */
  if (exon_genomestart < exon_genomeend) {
    fprintf(fp,"%u\t%u\t",exon_genomestart,exon_genomeend); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",exon_genomeend,exon_genomestart); /* 4,5: start, end */
  }
  fprintf(fp,"%d\t",pctidentity);	/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      fprintf(fp,"+\t");
    } else {
      fprintf(fp,"-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      fprintf(fp,"-\t");		/* 7: strand */
    } else {
      fprintf(fp,"+\t");
    }
  }

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.mrna%d.exon%d;",accession,pathnum,exonno);
  fprintf(fp,"Name=%s;",accession);
  fprintf(fp,"Parent=%s.mrna%d;",accession,pathnum);
  if (cdna_direction >= 0) {
    fprintf(fp,"Target=%s %d %d +\n",accession,exon_querystart,exon_queryend);
  } else {
    fprintf(fp,"Target=%s %d %d -\n",accession,exon_queryend,exon_querystart);
  }

  return;
}

static void
print_gff3_cds (FILE *fp, int cdsno, int pathnum, char *sourcename, char *accession, char *chrstring,
		int cds_genomestart, int cds_genomeend,
		int cds_querystart, int cds_queryend, bool watsonp, int cdna_direction,
		int pctidentity, int cds_phase) {

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"CDS\t");		/* 3: type */
  if (cds_genomestart < cds_genomeend) {
    fprintf(fp,"%u\t%u\t",cds_genomestart,cds_genomeend); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",cds_genomeend,cds_genomestart); /* 4,5: start, end */
  }
  fprintf(fp,"%d\t",pctidentity);	/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      fprintf(fp,"+\t");
    } else {
      fprintf(fp,"-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      fprintf(fp,"-\t");		/* 7: strand */
    } else {
      fprintf(fp,"+\t");
    }
  }

  fprintf(fp,"%d\t",cds_phase);	/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.mrna%d.cds%d;",accession,pathnum,cdsno);
  fprintf(fp,"Name=%s;",accession);
  fprintf(fp,"Parent=%s.mrna%d;",accession,pathnum);
  if (cdna_direction >= 0) {
    fprintf(fp,"Target=%s %d %d +\n",accession,cds_querystart,cds_queryend);
  } else {
    fprintf(fp,"Target=%s %d %d -\n",accession,cds_queryend,cds_querystart);
  }

  return;
}


static void
print_gff3_cdna_match (FILE *fp, int pathnum, char *sourcename, char *accession, char *chrstring,
		       int exon_genomestart, int exon_genomeend,
		       int exon_querystart, int exon_queryend, bool watsonp,
		       int pctidentity, List_T tokens) {
  
  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"cDNA_match\t");		/* 3: type */
  if (exon_genomestart < exon_genomeend) {
    fprintf(fp,"%u\t%u\t",exon_genomestart,exon_genomeend); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",exon_genomeend,exon_genomestart); /* 4,5: start, end */
  }
  fprintf(fp,"%d\t",pctidentity);	/* 6: score */

  /* 7: strand */
  if (watsonp == true) {
    fprintf(fp,"+\t");
  } else {
    fprintf(fp,"-\t");
  }

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.path%d;",accession,pathnum);
  fprintf(fp,"Name=%s;",accession);
  fprintf(fp,"Target=%s %d %d;Gap=",accession,exon_querystart,exon_queryend);
  print_tokens_gff3(fp,tokens);
  putc('\n',fp);

  return;
}


static char
strand_char (int strand) {
  switch (strand) {
    case  1: return '+';
    case -1: return '-';
      /* case  0: return '?'; -- Now returning '.' for unknown strand */
    default: return '.';
  }
}


static void
print_gff3_est_match (FILE *fp, int pathnum, T start, T end,
		      char *sourcename, char *accession, char *chrstring,
		      int exon_genomestart, int exon_genomeend,
		      int exon_querystart, int exon_queryend,
		      int querylength_given, int skiplength, int matches, int mismatches, int qindels, int tindels,
		      int unknowns, bool watsonp, int cdna_direction, int pctidentity, List_T tokens) {
  int feature_strand, target_strand;
  double coverage, fracidentity;
  int den;
  int querypos1, querypos2;

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"EST_match\t");	/* 3: type */
  if (exon_genomestart < exon_genomeend) {
    fprintf(fp,"%u\t%u\t",exon_genomestart,exon_genomeend); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",exon_genomeend,exon_genomestart); /* 4,5: start, end */
  }
  fprintf(fp,"%d\t",pctidentity);	/* 6: score */

  /* 7: strand */
  feature_strand = watsonp ? cdna_direction : -cdna_direction;
  fprintf(fp,"%c\t",strand_char(feature_strand));

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.path%d;",accession,pathnum);
  fprintf(fp,"Name=%s;",accession);
  target_strand = cdna_direction != 0 ? cdna_direction : (watsonp ? 1 : -1);
  fprintf(fp,"Target=%s %d %d %c;Gap=",accession,exon_querystart,exon_queryend,
      strand_char(target_strand));
  print_tokens_gff3(fp,tokens);

  querypos1 = start->querypos;
  querypos2 = end->querypos;

#ifdef PMAP
  coverage = (double) (querypos2 - querypos1 + 1)/(double) (3*(querylength_given + skiplength));
  /* Can have coverage greater than given querylength because of added '*' at end */
  if (coverage > 1.0) {
    coverage = 1.0;
  }
#else
  coverage = (double) (querypos2 - querypos1 + 1)/(double) (querylength_given + skiplength);
#endif
  fprintf(fp,";coverage=%.1f",((double) rint(1000.0*coverage))/10.0);

  if ((den = matches + mismatches + qindels + tindels) == 0) {
    fracidentity = 1.0;
  } else {
    fracidentity = (double) matches/(double) den;
  }
  fprintf(fp,";identity=%.1f",((double) rint(1000.0*fracidentity))/10.0);
  fprintf(fp,";matches=%d;mismatches=%d;indels=%d;unknowns=%d",
	  matches,mismatches,qindels+tindels,unknowns);

  putc('\n',fp);
}


static void
print_gff3_exons_forward (FILE *fp, struct T *pairs, int npairs, int pathnum, T start, T end,
			  char *sourcename, char *accession, char *chrstring,
			  int querylength_given, int skiplength, int matches, int mismatches,
			  int qindels, int tindels, int unknowns, bool watsonp, int cdna_direction,
			  bool gff_introns_p, bool gff_gene_format_p, bool gff_estmatch_format_p) {
  bool in_exon = false;
  struct T *ptr, *this = NULL;
  int exon_querystart = -1, exon_queryend;
  Chrpos_T exon_genomestart = -1, exon_genomeend, intron_start, intron_end;
  int pctidentity, num = 0, den = 0, exonno = 0, i;
  int Mlength = 0, Ilength = 0, Dlength = 0;
  List_T tokens = NULL;
  char token[10];
#if 0
  int intronno = 0;
#endif
  int estmatch_querystart, estmatch_queryend, estmatch_genomestart, estmatch_genomeend;
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    /* prev = this; */
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
	
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	if (gff_gene_format_p == true) {
	  print_gff3_exon(fp,++exonno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
			  exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
	} else {
	  if (Mlength > 0) {
	    sprintf(token,"M%d",Mlength);
	    tokens = push_token(tokens,token);
	  } else if (Ilength > 0) {
	    sprintf(token,"I%d",Ilength);
	    tokens = push_token(tokens,token);
	  } else if (Dlength > 0) {
	    sprintf(token,"D%d",Dlength);
	    tokens = push_token(tokens,token);
	  }
	  if (gff_estmatch_format_p == false) {
	    tokens = List_reverse(tokens);
	    /* ++exonno; */
	    print_gff3_cdna_match(fp,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
				  exon_querystart,exon_queryend,watsonp,pctidentity,tokens);
	    List_free(&tokens);
	  }
	}

	Mlength = Ilength = Dlength = 0;
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}

	if (gff_estmatch_format_p == true && i > 0) {
	  sprintf(token,"N%u",abs(intron_end - intron_start) + 1);
	  tokens = push_token(tokens,token);
	} else if (gff_introns_p == true) {
	  if (i > 0) {
#if 0
	    printf_gff3_intron(++intronno,pathnum,sourcename,accession,chrstring,?,?,intron_start,intron_end,watsonp);
#endif
	  }
	  putc('\n',fp);
	}

	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (gff_gene_format_p == true) {
	  /* Don't deal with tokens */
	} else if (this->genome == ' ') {
	  if (Mlength > 0) {
	    sprintf(token,"M%d",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Dlength > 0) {
	    /* unlikely */
	    sprintf(token,"D%d",Dlength);
	    tokens = push_token(tokens,token);
	    Dlength = 0;
	  }
	  Ilength++;
	} else if (this->cdna == ' ') {
	  if (Mlength > 0) {
	    sprintf(token,"M%d",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Ilength > 0) {
	    sprintf(token,"I%d",Ilength);
	    tokens = push_token(tokens,token);
	    Ilength = 0;
	  }
	  Dlength++;
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

	/* Previously not counted in numerator or denominator */
	den++;

      } else {
	/* Count in token even if unknown base */

	if (gff_gene_format_p == true) {
	  /* Don't deal with tokens */
	} else if (Ilength > 0) {
	  sprintf(token,"I%d",Ilength);
	  tokens = push_token(tokens,token);
	  Ilength = 0;
	} else if (Dlength > 0) {
	  sprintf(token,"D%d",Dlength);
	  tokens = push_token(tokens,token);
	  Dlength = 0;
	}
	Mlength++;

#ifdef PMAP	
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	  num++;
	} else if (this->comp == AMBIGUOUS_COMP) {
	  num++;
	}
#else
	if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	  /* Comp must be a space */
	  /* Don't count in numerator or denominator */
	} else {
	  den++;
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	    num++;
	  } else if (this->comp == AMBIGUOUS_COMP) {
	    den--;
	  }
	}
#endif

      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  /* prev = this; */
  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;

  if (den == 0) {
    pctidentity = 100;
  } else {
    pctidentity = (int) floor(100.0*(double) num/(double) den);
  }
	
  if (gff_gene_format_p == true) {
    print_gff3_exon(fp,++exonno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		    exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
  } else {
    if (Mlength > 0) {
      sprintf(token,"M%d",Mlength);
      tokens = push_token(tokens,token);
    } else if (Ilength > 0) {
      sprintf(token,"I%d",Ilength);
      tokens = push_token(tokens,token);
    } else if (Dlength > 0) {
      sprintf(token,"D%d",Dlength);
      tokens = push_token(tokens,token);
    }
    if (gff_estmatch_format_p == true) {
      estmatch_querystart = pairs->querypos + 1;
      estmatch_queryend = exon_queryend;
      estmatch_genomestart = pairs->genomepos + 1;
      estmatch_genomeend = exon_genomeend;
      if (watsonp) {
	tokens = List_reverse(tokens);
      }
      print_gff3_est_match(fp,pathnum,start,end,sourcename,accession,chrstring,
			   estmatch_genomestart,estmatch_genomeend,
			   estmatch_querystart,estmatch_queryend,
			   querylength_given,skiplength,matches,mismatches,qindels,tindels,unknowns,
			   watsonp,cdna_direction,pctidentity,tokens);
    } else {
      tokens = List_reverse(tokens);
      /* ++exonno; */
      print_gff3_cdna_match(fp,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
			    exon_querystart,exon_queryend,watsonp,pctidentity,tokens);
    }
    List_free(&tokens);
  }

  return;
}

static void
print_gff3_exons_backward (FILE *fp, struct T *pairs, int npairs, int pathnum, char *sourcename, char *accession, char *chrstring,
			   bool watsonp, int cdna_direction, bool gff_introns_p) {
  bool in_exon = false;
  struct T *ptr, *this = NULL;
  int exon_querystart = -1, exon_queryend;
  Chrpos_T exon_genomestart = -1, exon_genomeend;
  int pctidentity, num = 0, den = 0, exonno = 0, i;
#if 0
  int intronno = 0;
  Chrpos_T intron_start, intron_end;
#endif
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;

  ptr = &(pairs[npairs-1]);
  for (i = npairs-1; i >= 0; i--) {
    /* prev = this; */
    this = ptr--;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_exon(fp,++exonno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
			exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);

	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;

	if (gff_introns_p == true) {
	  if (i > 0) {
#if 0
	    printf_gff3_intron(++intronno,pathnum,sourcename,accession,chrstring,?,?,intron_start,intron_end,watsonp);
#endif
	  }
	  putc('\n',fp);
	}

	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Previously not counted in numerator or denominator */
	den++;

#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	/* Comp must be a space */
	/* Don't count in numerator or denominator */
#endif
      } else {
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	  num++;
	} else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	  num++;
#else
	  den--;
#endif
	}
      }
    }
    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  /* prev = this; */
  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;

  if (den == 0) {
    pctidentity = 100;
  } else {
    pctidentity = (int) floor(100.0*(double) num/(double) den);
  }
	
  print_gff3_exon(fp,++exonno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		  exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
  return;
}


static void
print_gff3_cdss_forward (FILE *fp, struct T *pairs, int npairs, int pathnum, char *sourcename, char *accession, char *chrstring,
			 bool watsonp, int cdna_direction) {
  bool in_cds = false;
  struct T *ptr, *this = NULL;
  int exon_querystart = -1, exon_queryend, exon_phase;
  Chrpos_T exon_genomestart = -1, exon_genomeend;
  int pctidentity, num = 0, den = 0, cdsno = 0;
#if 0
  Chrpos_T intron_start, intron_end;
#endif
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;

  ptr = pairs;
  while (ptr < &(pairs[npairs])) {
    /* prev = this; */
    this = ptr++;

    if (in_cds == true) {
      if (this->aaphase_e == -1) { /* was aaphase_g */
	/* End of cds */
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	  
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_cds(fp,++cdsno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		       exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);

	in_cds = false;

      } else {
	/* Continuation of cds */
	if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	  /* Previously not counted in numerator or denominator */
	  den++;
	  
#ifndef PMAP
	} else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	  /* Comp must be a space */
	  /* Don't count in numerator or denominator */
#endif
	} else {
	  den++;
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	    num++;
	  } else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	    num++;
#else
	    den--;
#endif
	  }
	}
      }

    } else {
      if (this->aaphase_e == -1) {
	/* Continuation of non-cds */
      } else {
	/* Start of cds */
	exon_querystart = this->querypos + 1;
	exon_phase = this->aaphase_e; /* ? was aaphase_g */
	exon_genomestart = this->genomepos + 1;

	num = den = 0;
	in_cds = true;
      }
    }
    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  if (in_cds == true) {
    exon_queryend = last_querypos + 1;
    exon_genomeend = last_genomepos + 1;

    if (den == 0) {
      pctidentity = 100;
    } else {
      pctidentity = (int) floor(100.0*(double) num/(double) den);
    }
	
    print_gff3_cds(fp,++cdsno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		   exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);
  }

  return;
}

static void
print_gff3_cdss_backward (FILE *fp, struct T *pairs, int npairs, int pathnum, char *sourcename, char *accession, char *chrstring,
			  bool watsonp, int cdna_direction) {
  bool in_cds = false;
  struct T *ptr, *this = NULL;
  int exon_querystart = -1, exon_queryend, exon_phase;
  Chrpos_T exon_genomestart = -1, exon_genomeend;
  int pctidentity, num = 0, den = 0, cdsno = 0;
#if 0
  Chrpos_T intron_start, intron_end;
#endif
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;


  ptr = &(pairs[npairs-1]);
  while (ptr >= &(pairs[0])) {
    /* prev = this; */
    this = ptr--;

    if (in_cds == true) {
      if (this->aaphase_e == -1) { /* was aaphase_g */
	/* End of cds */
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	  
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_cds(fp,++cdsno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		       exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);

	in_cds = false;

      } else {
	/* Continuation of cds */
	if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	  /* Previously not counted in numerator or denominator */
	  den++;

#ifndef PMAP
	} else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	  /* Comp must be a space */
	  /* Don't count in numerator or denominator */
#endif
	} else {
	  den++;
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	    num++;
	  } else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	    num++;
#else
	    den--;
#endif
	  }
	}
      }

    } else {
      if (this->aaphase_e == -1) { /* was aaphase_g */
	/* Continuation of non-cds */
      } else {
	/* Start of cds */
	exon_querystart = this->querypos + 1;
	exon_phase = this->aaphase_e; /* ? was aaphase_g */
	exon_genomestart = this->genomepos + 1;

	num = den = 0;
	in_cds = true;
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  if (in_cds == true) {
    exon_queryend = last_querypos + 1;
    exon_genomeend = last_genomepos + 1;

    if (den == 0) {
      pctidentity = 100;
    } else {
      pctidentity = (int) floor(100.0*(double) num/(double) den);
    }

    print_gff3_cds(fp,++cdsno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		   exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);
  }

  return;
}


void
Pair_print_gff3 (FILE *fp, struct T *pairs, int npairs, int pathnum, char *accession, 
		 T start, T end, Chrnum_T chrnum, Univ_IIT_T chromosome_iit, Sequence_T usersegment,
		 int translation_end,
		 int querylength_given, int skiplength, int matches, int mismatches, 
		 int qindels, int tindels, int unknowns, bool watsonp, int cdna_direction,
		 bool gff_gene_format_p, bool gff_estmatch_format_p, char *sourcename) {
  char *chrstring = NULL;
  Chrpos_T chrpos1, chrpos2;

  if (chrnum == 0) {
    chrstring = Sequence_accession(usersegment);
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  if (sourcename == NULL) {
    sourcename = "NA";
  }

  if (gff_gene_format_p == true) {
    chrpos1 = start->genomepos;
    chrpos2 = end->genomepos;

    print_gff3_gene(fp,pathnum,sourcename,accession,chrstring,chrpos1+1,chrpos2+1,watsonp,cdna_direction);
    print_gff3_mrna(fp,pathnum,start,end,sourcename,accession,chrstring,chrpos1+1,chrpos2+1,
		    querylength_given,skiplength,matches,mismatches,qindels,tindels,unknowns,
		    watsonp,cdna_direction);

    if (cdna_direction >= 0) {
      print_gff3_exons_forward(fp,pairs,npairs,pathnum,start,end,sourcename,accession,chrstring,
			       querylength_given,skiplength,matches,mismatches,qindels,tindels,unknowns,
			       watsonp,cdna_direction,/*gff_introns_p*/false,/*gff_gene_format_p*/true,
			       /*gff_estmatch_format_p*/false);
      if (translation_end > 0) {
	print_gff3_cdss_forward(fp,pairs,npairs,pathnum,sourcename,accession,chrstring,watsonp,
				cdna_direction);
      }
    } else {
      print_gff3_exons_backward(fp,pairs,npairs,pathnum,sourcename,accession,chrstring,watsonp,
				cdna_direction,/*gff_introns_p*/false);
      if (translation_end > 0) {
	print_gff3_cdss_backward(fp,pairs,npairs,pathnum,sourcename,accession,chrstring,watsonp,
				 cdna_direction);
      }
    }

  } else {
    print_gff3_exons_forward(fp,pairs,npairs,pathnum,start,end,sourcename,accession,chrstring,
			     querylength_given,skiplength,matches,mismatches,qindels,tindels,unknowns,
			     watsonp,cdna_direction,/*gff_introns_p*/false,/*gff_gene_format_p*/false,
			     gff_estmatch_format_p);
  }

  if (gff3_separators_p == true) {
    fprintf(fp,"###\n");		/* Terminates alignment */
  }

  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


int
Pair_circularpos (int *alias, struct T *pairs, int npairs, Chrpos_T chrlength, bool plusp, int querylength) {
  int i;
  struct T *ptr;

  if (plusp == true) {
    i = 0;
    ptr = pairs;
    if (ptr->genomepos >= chrlength) {
      /* All of read after trimming is in circular alias */
      *alias = +1;
      return -1;
    } else {
      while (i < npairs && ptr->genomepos < chrlength) {
	i++;
	ptr++;
      }
      if (i >= npairs) {
	/* All of read after trimming is in circular proper */
	*alias = -1;
	return -1;
      } else {
	/* Some of read is in circular proper and some is in circular alias */
	*alias = 0;
	return ptr->querypos;
      }
    }

  } else {
    i = npairs - 1;
    ptr = &(pairs[i]);
    if (ptr->genomepos >= chrlength) {
      /* All of read after trimming is in circular alias */
      *alias = +1;
      return -1;
    } else {
      while (i >= 0 && ptr->genomepos < chrlength) {
	i--;
	ptr--;
      }
      if (i < 0) {
	/* All of read after trimming is in circular proper */
	*alias = -1;
	return -1;
      } else {
	/* Some of read is in circular proper and some is in circular alias */
	*alias = 0;
	return (querylength - ptr->querypos - 1);
      }
    }
  }
}


#ifndef PMAP
/************************************************************************
 *   GSNAP
 ************************************************************************/

/* Based on procedure in substring.c */
static void
print_splicesite_labels (FILE *fp, Chrnum_T chrnum, Chrpos_T splicesitepos,
			 char *tag, IIT_T splicesites_iit, int *splicesites_divint_crosstable,
			 int typeint) {
  int *splicesites, nsplicesites, i;
  char *label;
  bool allocp;

  /* splicesitepos of 0U indicates it was not set */
  if (splicesites_iit != NULL && splicesitepos > 0U) {
    /* Note: IIT_get_typed_signed_with_divno does not work here */
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[chrnum],
						    splicesitepos,splicesitepos+1U,typeint);
    if (nsplicesites == 0) {
#if 0
      /* GMAP can have novel splice sites in the middle */
      fprintf(stderr,"Supposed to have a splicesite at chrnum %d, %u..%u, type %d\n",
	      chrnum,splicesitepos,splicesitepos+1U,typeint);
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

static void
print_endtypes (FILE *fp,
		Endtype_T endtype1, int ntrim1, int nindels1, Chrpos_T prev_splice_dist,
		Endtype_T endtype2, int ntrim2, int nindels2, Chrpos_T splice_dist,
		int nmatches, int nmismatches_refdiff, int nmismatches_bothdiff,
		Chrnum_T chrnum, Univcoord_T chroffset,
		Chrpos_T exon_genomestart, Chrpos_T exon_genomeend,
		bool watsonp, int cdna_direction,
		IIT_T splicesites_iit, int *splicesites_divint_crosstable,
		int donor_typeint, int acceptor_typeint) {
  double prob;
  int prev_splicesitepos = 0, splicesitepos = 0;
  int typeint1, typeint2;

  if (endtype1 == END) {
    fprintf(fp,"start:%d",ntrim1);
  } else if (endtype1 == INS) {
    fprintf(fp,"ins:%d",nindels1);
  } else if (endtype1 == DEL) {
    fprintf(fp,"del:%d",nindels1);
  } else if (endtype1 == DON || endtype1 == AMB_DON) {
    typeint1 = donor_typeint;
    if (watsonp == true) {
      prev_splicesitepos = exon_genomestart-1U;
      if (cdna_direction > 0) {
	prob = Maxent_hr_donor_prob(chroffset+exon_genomestart,chroffset);
      } else if (cdna_direction < 0) {
	prob = Maxent_hr_antidonor_prob(chroffset+exon_genomestart-1U,chroffset);
      } else {
	abort();
      }
    } else {
      prev_splicesitepos = exon_genomestart;
      if (cdna_direction > 0) {
	prob = Maxent_hr_antidonor_prob(chroffset+exon_genomestart-1U,chroffset);
      } else if (cdna_direction < 0) {
	prob = Maxent_hr_donor_prob(chroffset+exon_genomestart,chroffset);
      } else {
	abort();
      }
    }
    fprintf(fp,"donor:%.2f",prob);
  } else if (endtype1 == ACC || endtype1 == AMB_ACC) {
    typeint1 = acceptor_typeint;
    if (watsonp == true) {
      prev_splicesitepos = exon_genomestart-1U;
      if (cdna_direction > 0) {
	prob = Maxent_hr_acceptor_prob(chroffset+exon_genomestart-1U,chroffset);
      } else if (cdna_direction < 0) {
	prob = Maxent_hr_antiacceptor_prob(chroffset+exon_genomestart,chroffset);
      } else {
	abort();
      }
    } else {
      prev_splicesitepos = exon_genomestart;
      if (cdna_direction > 0) {
	prob = Maxent_hr_antiacceptor_prob(chroffset+exon_genomestart,chroffset);
      } else if (cdna_direction < 0) {
	prob = Maxent_hr_acceptor_prob(chroffset+exon_genomestart-1U,chroffset);
      } else {
	abort();
      }
    }
    fprintf(fp,"acceptor:%.2f",prob);
  } else {
    fprintf(fp,"unknown");
  }

  fprintf(fp,"..");

  if (endtype2 == END) {
    fprintf(fp,"end:%d",ntrim2);
  } else if (endtype2 == INS) {
    fprintf(fp,"ins:%d",nindels2);
  } else if (endtype2 == DEL) {
    fprintf(fp,"del:%d",nindels2);
  } else if (endtype2 == DON || endtype2 == AMB_DON) {
    typeint2 = donor_typeint;
    if (watsonp == true) {
      splicesitepos = exon_genomeend;
      if (cdna_direction > 0) {
	prob = Maxent_hr_donor_prob(chroffset+exon_genomeend,chroffset);
      } else if (cdna_direction < 0) {
	prob = Maxent_hr_antidonor_prob(chroffset+exon_genomeend-1U,chroffset);
      } else {
	abort();
      }
    } else {
      splicesitepos = exon_genomeend-1U;
      if (cdna_direction > 0) {
	prob = Maxent_hr_antidonor_prob(chroffset+exon_genomeend-1U,chroffset);
      } else if (cdna_direction < 0) {
	prob = Maxent_hr_donor_prob(chroffset+exon_genomeend,chroffset);
      } else {
	abort();
      }
    }
    fprintf(fp,"donor:%.2f",prob);
  } else if (endtype2 == ACC || endtype2 == AMB_ACC) {
    typeint2 = acceptor_typeint;
    if (watsonp == true) {
      splicesitepos = exon_genomeend;
      if (cdna_direction > 0) {
	prob = Maxent_hr_acceptor_prob(chroffset+exon_genomeend-1U,chroffset);
      } else if (cdna_direction < 0) {
	prob = Maxent_hr_antiacceptor_prob(chroffset+exon_genomeend,chroffset);
      } else {
	abort();
      }
    } else {
      splicesitepos = exon_genomeend-1U;
      if (cdna_direction > 0) {
	prob = Maxent_hr_antiacceptor_prob(chroffset+exon_genomeend,chroffset);
      } else if (cdna_direction < 0) {
	prob = Maxent_hr_acceptor_prob(chroffset+exon_genomeend-1U,chroffset);
      } else {
	abort();
      }
    }
    fprintf(fp,"acceptor:%.2f",prob);
  } else {
    fprintf(fp,"unknown");
  }

  fprintf(fp,",matches:%d,sub:%d",nmatches,nmismatches_bothdiff);
  fprintf(fp,"+%d=%d",nmismatches_refdiff - nmismatches_bothdiff,nmismatches_refdiff);

  if (prev_splice_dist != 0 && splice_dist != 0) {
    /* Double introns */
    if (cdna_direction > 0) {
      fprintf(fp,",dir:sense,splice_type:consistent");
    } else {
      fprintf(fp,",dir:antisense,splice_type:consistent");
    }
    fprintf(fp,",splice_dist_1:%u,splice_dist_2:%u",prev_splice_dist,splice_dist);
    print_splicesite_labels(fp,chrnum,prev_splicesitepos,"label_1",splicesites_iit,
			    splicesites_divint_crosstable,typeint1);
    print_splicesite_labels(fp,chrnum,splicesitepos,"label_2",splicesites_iit,
			    splicesites_divint_crosstable,typeint2);

  } else if (prev_splice_dist != 0) {
    /* Prev intron */
    if (cdna_direction > 0) {
      fprintf(fp,",dir:sense,splice_type:consistent");
    } else {
      fprintf(fp,",dir:antisense,splice_type:consistent");
    }
    fprintf(fp,",splice_dist_1:%u",prev_splice_dist);
    print_splicesite_labels(fp,chrnum,prev_splicesitepos,"label_1",splicesites_iit,
			    splicesites_divint_crosstable,typeint1);
    if (endtype2 == AMB_DON || endtype2 == AMB_ACC) {
      print_splicesite_labels(fp,chrnum,splicesitepos,"label_2",splicesites_iit,
			      splicesites_divint_crosstable,typeint2);
    }

  } else if (splice_dist != 0) {
    /* Next intron */
    if (cdna_direction > 0) {
      fprintf(fp,",dir:sense,splice_type:consistent");
    } else {
      fprintf(fp,",dir:antisense,splice_type:consistent");
    }
    fprintf(fp,",splice_dist_2:%u",splice_dist);
    if (endtype1 == AMB_DON || endtype1 == AMB_ACC) {
      print_splicesite_labels(fp,chrnum,splicesitepos,"label_1",splicesites_iit,
			      splicesites_divint_crosstable,typeint1);
    }
    print_splicesite_labels(fp,chrnum,splicesitepos,"label_2",splicesites_iit,
			    splicesites_divint_crosstable,typeint2);
  } else {
    /* No introns */
    if (endtype1 == AMB_DON || endtype1 == AMB_ACC) {
      print_splicesite_labels(fp,chrnum,splicesitepos,"label_1",splicesites_iit,
			      splicesites_divint_crosstable,typeint1);
    }
    if (endtype2 == AMB_DON || endtype2 == AMB_ACC) {
      print_splicesite_labels(fp,chrnum,splicesitepos,"label_2",splicesites_iit,
			      splicesites_divint_crosstable,typeint2);
    }
  }

  return;
}




/* Based on print_pair_info in stage3hr.c */
static void
print_pair_info (FILE *fp, int insertlength, int pairscore, Pairtype_T pairtype) {
  fprintf(fp,"pair_score:%d",pairscore);
  fprintf(fp,",insert_length:%d",insertlength);

  switch (pairtype) {
  case CONCORDANT: break;
  case PAIRED_SCRAMBLE: fprintf(fp,",pairtype:scramble"); break;
  case PAIRED_INVERSION: fprintf(fp,",pairtype:inversion"); break;
  case PAIRED_TOOLONG: fprintf(fp,",pairtype:toolong"); break;
  case CONCORDANT_TRANSLOCATIONS: break;
  case CONCORDANT_TERMINAL: break;
  case PAIRED_UNSPECIFIED: abort();
  case UNSPECIFIED: break;
  case UNPAIRED: abort();
  }

  return;
}



void
Pair_print_gsnap (FILE *fp, struct T *pairs_querydir, int npairs, int nsegments, bool invertedp,
		  Endtype_T start_endtype, Endtype_T end_endtype,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		  int querylength, bool watsonp, int cdna_direction, int score,
		  int insertlength, int pairscore, int mapq_score,
		  Univ_IIT_T chromosome_iit, IIT_T splicesites_iit,
		  int *splicesites_divint_crosstable, int donor_typeint, int acceptor_typeint) {
  bool in_exon = true;
  struct T *pairs, *ptr, *ptr0, *this = NULL;
  int exon_querystart = -1, exon_queryend;
  Chrpos_T exon_genomestart = -1U, exon_genomeend;
  int querypos, nmismatches_refdiff, nmismatches_bothdiff, nmatches,
    ntrim_start, ntrim_end, nindels, prev_nindels, i;
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;
  Univcoord_T pos;
  Chrpos_T splice_dist = 0U, prev_splice_dist;
  char *chr, strand, c, c_alt;
  Endtype_T endtype, prev_endtype;
  bool allocp, firstp = true;


  if (invertedp == true) {
    pairs = invert_and_revcomp_path_and_coords(pairs_querydir,npairs,querylength);
    watsonp = !watsonp;
    cdna_direction = -cdna_direction;
  } else {
    pairs = pairs_querydir;
  }


  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
  if (watsonp == true) {
    strand = '+';
  } else {
    strand = '-';
  }

  fprintf(fp," ");		/* Beginning of GSNAP line */

  ptr = pairs;
  exon_querystart = ptr->querypos + 1;
  exon_genomestart = ptr->genomepos + 1;
  endtype = start_endtype;
  nmismatches_refdiff = nmismatches_bothdiff = nmatches = 0;

  ntrim_start = ptr->querypos;

  /* Print leading trimmed characters */
  if (watsonp == true) {
    if (ntrim_start >= (exon_genomestart - 1)) {
      for (querypos = 0; querypos < ntrim_start - exon_genomestart + 1; querypos++) {
	fprintf(fp,"*");
      }
      pos = chroffset;
      for ( ; querypos < ntrim_start; querypos++) {
	c = Genome_get_char_blocks(&c_alt,pos++);
	fprintf(fp,"%c",tolower(c));
      }

    } else {
      pos = chroffset + (exon_genomestart - 1) - ntrim_start;
      for (querypos = 0; querypos < ntrim_start; querypos++) {
	c = Genome_get_char_blocks(&c_alt,pos++);
	fprintf(fp,"%c",tolower(c));
      }
    }

  } else {
    if ((pos = chroffset + (exon_genomestart - 1) + ntrim_start) >= chrhigh) {
      assert(ntrim_start - (int) (chrhigh - chroffset - exon_genomestart + 1) < querylength);
      for (querypos = 0; querypos <= ntrim_start - (int) (chrhigh - chroffset - exon_genomestart + 1); querypos++) {
	fprintf(fp,"*");
      }
      pos = chrhigh - 1;
      for ( ; querypos < ntrim_start; querypos++) {
	c = Genome_get_char_blocks(&c_alt,pos--);
	fprintf(fp,"%c",tolower(complCode[(int) c]));
      }

    } else {
      for (querypos = 0; querypos < ntrim_start; querypos++) {
	c = Genome_get_char_blocks(&c_alt,pos--);
	fprintf(fp,"%c",tolower(complCode[(int) c]));
      }
    }
  }

  i = 0;
  while (i < npairs) {
    /* prev = this; */
    this = ptr++;
    i++;

    if (this->gapp) {
      if (in_exon == true) {
	/* SPLICE START */
	ptr0 = ptr;
	while (ptr0->gapp) {
	  ptr0++;
	}
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	prev_endtype = endtype;

	if (prev_endtype == INS || prev_endtype == DEL) {
	  prev_nindels = nindels;
	  prev_splice_dist = 0U;
	} else if (prev_endtype == DON) {
	  prev_splice_dist = splice_dist;
	  prev_endtype = ACC;
	} else if (prev_endtype == ACC) {
	  prev_splice_dist = splice_dist;
	  prev_endtype = DON;
	} else {
	  prev_splice_dist = 0U;
	}

	if (cdna_direction > 0) {
	  endtype = DON;
	} else {
	  endtype = ACC;
	}

	if (watsonp == true) {
	  splice_dist = ptr0->genomepos - last_genomepos - 1;
	  /* prev_splicesitepos = exon_genomestart-1U; */
	  /* splicesitepos = exon_genomeend; */
	} else {
	  splice_dist = last_genomepos - ptr0->genomepos - 1;
	  /* prev_splicesitepos = exon_genomestart; */
	  /* splicesitepos = exon_genomeend-1U; */
	}


	fprintf(fp,"%c",tolower(ptr[-1].genome)); /* dinucleotide */
	fprintf(fp,"%c",tolower(ptr[0].genome));
	for (querypos = exon_queryend+2; querypos < querylength; querypos++) {
	  fprintf(fp,"-");
	}

	fprintf(fp,"\t%d..%d",exon_querystart,exon_queryend);
	fprintf(fp,"\t%c%s:%u..%u",strand,chr,exon_genomestart,exon_genomeend);
	fprintf(fp,"\t");
	print_endtypes(fp,prev_endtype,ntrim_start,prev_nindels,prev_splice_dist,
		       endtype,ntrim_end,/*nindels*/0,splice_dist,
		       nmatches,nmismatches_refdiff,nmismatches_bothdiff,chrnum,chroffset,
		       exon_genomestart,exon_genomeend,watsonp,cdna_direction,
		       splicesites_iit,splicesites_divint_crosstable,
		       donor_typeint,acceptor_typeint);

	if (firstp == true) {
	  fprintf(fp,"\tsegs:%d,align_score:%d,mapq:%d",nsegments,score,mapq_score);
	  fprintf(fp,",method:gmap");
	  fprintf(fp,"\t");
	  print_pair_info(fp,insertlength,pairscore,/*pairtype*/CONCORDANT);
	  firstp = false;
	}

	nmismatches_refdiff = nmismatches_bothdiff = nmatches = 0;
	fprintf(fp,"\n");

	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* May want to print dinucleotides */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* SPLICE CONTINUATION */
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;

	fprintf(fp,",");
	for (querypos = 0; querypos < this->querypos - 2; querypos++) {
	  fprintf(fp,"-");
	}
	fprintf(fp,"%c",tolower(ptr[-3].genome)); /* dinucleotide */
	fprintf(fp,"%c",tolower(ptr[-2].genome));

	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
	  /* INSERTION */
	  exon_queryend = last_querypos + 1;
	  exon_genomeend = last_genomepos + 1;
	  prev_endtype = endtype;
	  endtype = INS;

	  if (prev_endtype == INS || prev_endtype == DEL) {
	    prev_nindels = nindels;
	    prev_splice_dist = 0U;
	  } else if (prev_endtype == DON) {
	    prev_splice_dist = splice_dist;
	    prev_endtype = ACC;
	  } else if (prev_endtype == ACC) {
	    prev_splice_dist = splice_dist;
	    prev_endtype = DON;
	  } else {
	    prev_splice_dist = 0U;
	  }

	  /* indel_pos = this->querypos; */
	  nindels = 0;
	  while (i < npairs && this->gapp == false && this->genome == ' ') {
	    nindels++;
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  /* Finish rest of this line */
	  for (querypos = exon_queryend; querypos < querylength; querypos++) {
	    fprintf(fp,"-");
	  }
	  fprintf(fp,"\t%d..%d",exon_querystart,exon_queryend);
	  fprintf(fp,"\t%c%s:%u..%u",strand,chr,exon_genomestart,exon_genomeend);
	  fprintf(fp,"\t");
	  print_endtypes(fp,prev_endtype,ntrim_start,prev_nindels,prev_splice_dist,
			 endtype,ntrim_end,nindels,/*splice_dist*/0U,
			 nmatches,nmismatches_refdiff,nmismatches_bothdiff,chrnum,chroffset,
			 exon_genomestart,exon_genomeend,watsonp,cdna_direction,
			 splicesites_iit,splicesites_divint_crosstable,
			 donor_typeint,acceptor_typeint);

	  if (firstp == true) {
	    fprintf(fp,"\tsegs:%d,align_score:%d,mapq:%d",nsegments,score,mapq_score);
	    fprintf(fp,",method:gmap");
	    fprintf(fp,"\t");
	    print_pair_info(fp,insertlength,pairscore,/*pairtype*/CONCORDANT);
	    firstp = false;
	  }

	  fprintf(fp,"\n,");

	  this = ptr;
	  exon_querystart = this->querypos + 1;
	  exon_genomestart = this->genomepos + 1;
	  nmismatches_refdiff = nmismatches_bothdiff = nmatches = 0;

	  /* Start of next line */
	  for (querypos = 1; querypos < exon_querystart; querypos++) {
	    fprintf(fp,"-");
	  }

	} else if (this->cdna == ' ') {
	  /* DELETION */
	  exon_queryend = last_querypos + 1;
	  exon_genomeend = last_genomepos + 1;
	  prev_endtype = endtype;
	  endtype = DEL;

	  if (prev_endtype == INS || prev_endtype == DEL) {
	    prev_nindels = nindels;
	    prev_splice_dist = 0U;
	  } else if (prev_endtype == DON) {
	    prev_endtype = ACC;
	    prev_splice_dist = splice_dist;
	  } else if (prev_endtype == ACC) {
	    prev_endtype = DON;
	    prev_splice_dist = splice_dist;
	  } else {
	    prev_splice_dist = 0U;
	  }

	  /* indel_pos = this->querypos; */
	  nindels = 0;
	  while (i < npairs && this->gapp == false && this->cdna == ' ') {
	    fprintf(fp,"%c",tolower(this->genome));
	    nindels++;
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  /* Finish rest of this line */
	  for (querypos = exon_queryend + nindels; querypos < querylength; querypos++) {
	    fprintf(fp,"-");
	  }
	  fprintf(fp,"\t%d..%d",exon_querystart,exon_queryend);
	  fprintf(fp,"\t%c%s:%u..%u",strand,chr,exon_genomestart,exon_genomeend);
	  fprintf(fp,"\t");
	  print_endtypes(fp,prev_endtype,ntrim_start,prev_nindels,prev_splice_dist,
			 endtype,ntrim_end,nindels,/*splice_dist*/0U,
			 nmatches,nmismatches_refdiff,nmismatches_bothdiff,chrnum,chroffset,
			 exon_genomestart,exon_genomeend,watsonp,cdna_direction,
			 splicesites_iit,splicesites_divint_crosstable,
			 donor_typeint,acceptor_typeint);

	  if (firstp == true) {
	    fprintf(fp,"\tsegs:%d,align_score:%d,mapq:%d",nsegments,score,mapq_score);
	    fprintf(fp,",method:gmap");
	    fprintf(fp,"\t");
	    print_pair_info(fp,insertlength,pairscore,/*pairtype*/CONCORDANT);
	    firstp = false;
	  }

	  fprintf(fp,"\n,");

	  this = ptr;
	  exon_querystart = this->querypos + 1;
	  exon_genomestart = this->genomepos + 1;
	  nmismatches_refdiff = nmismatches_bothdiff = nmatches = 0;

	  /* Start of next line */
	  for (querypos = 1; querypos < exon_querystart; querypos++) {
	    fprintf(fp,"-");
	  }

	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	c = this->genome;
	if (this->genome == this->cdna) {
	  fprintf(fp,"%c",c);
	  nmatches++;
	} else if (this->genomealt == this->cdna) {
	  fprintf(fp,"%c",c);
	  nmismatches_refdiff++;
	} else {
	  fprintf(fp,"%c",tolower(c));
	  nmismatches_bothdiff++;
	  nmismatches_refdiff++;
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  /* prev = this; */
  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;
  prev_endtype = endtype;
  endtype = end_endtype;

  if (prev_endtype == INS || prev_endtype == DEL) {
    prev_nindels = nindels;
    prev_splice_dist = 0U;
  } else if (prev_endtype == DON) {
    prev_endtype = ACC;
    prev_splice_dist = splice_dist;
  } else if (prev_endtype == ACC) {
    prev_endtype = DON;
    prev_splice_dist = splice_dist;
  } else {
    prev_splice_dist = 0U;
  }

  ntrim_end = querylength - exon_queryend;

  /* Print trailing trimmed characters */
  if (watsonp == true) {
    pos = chroffset + (exon_genomeend - 1);
    if (pos + ntrim_end >= chrhigh) {
      assert((int) (chrhigh - chroffset - exon_genomeend) < querylength);
      for (i = 0; i < (int) (chrhigh - chroffset - exon_genomeend); i++) {
	c = Genome_get_char_blocks(&c_alt,++pos);
	fprintf(fp,"%c",tolower(c));
      }
      for ( ; i < ntrim_end; i++) {
	fprintf(fp,"*");
      }

    } else {
      for (i = 0; i < ntrim_end; i++) {
	c = Genome_get_char_blocks(&c_alt,++pos);
	fprintf(fp,"%c",tolower(c));
      }
    }

  } else {
    pos = chroffset + (exon_genomeend - 1);
    if (ntrim_end >= (exon_genomeend - 1)) {
      for (i = 0; i < exon_genomeend - 1; i++) {
	c = Genome_get_char_blocks(&c_alt,--pos);
	fprintf(fp,"%c",tolower(complCode[(int) c]));
      }
      for ( ; i < ntrim_end; i++) {
	fprintf(fp,"*");
      }

    } else {
      for (i = 0; i < ntrim_end; i++) {
	c = Genome_get_char_blocks(&c_alt,--pos);
	fprintf(fp,"%c",tolower(complCode[(int) c]));
      }
    }
  }

  fprintf(fp,"\t%d..%d",exon_querystart,exon_queryend);
  fprintf(fp,"\t%c%s:%u..%u",strand,chr,exon_genomestart,exon_genomeend);
  fprintf(fp,"\t");
  print_endtypes(fp,prev_endtype,ntrim_start,prev_nindels,prev_splice_dist,
		 /*endtype2*/end_endtype,ntrim_end,/*nindels*/0,/*splice_dist*/0U,
		 nmatches,nmismatches_refdiff,nmismatches_bothdiff,chrnum,chroffset,
		 exon_genomestart,exon_genomeend,watsonp,cdna_direction,
		 splicesites_iit,splicesites_divint_crosstable,
		 donor_typeint,acceptor_typeint);

  if (firstp == true) {
    fprintf(fp,"\tsegs:%d,align_score:%d,mapq:%d",nsegments,score,mapq_score);
    fprintf(fp,",method:gmap");
    fprintf(fp,"\t");
    print_pair_info(fp,insertlength,pairscore,/*pairtype*/CONCORDANT);
    firstp = false;
  }

  fprintf(fp,"\n");

  if (allocp) {
    FREE(chr);
  }

  if (invertedp == true) {
    FREE(pairs);
  }

  return;
}


#ifdef GSNAP
/* Taken from NCBI Blast 2.2.29, algo/blast/core/blast_stat.c */
/* Karlin-Altschul formula: m n exp(-lambda * S + log k) = k m n exp(-lambda * S) */
/* Also in substring.c */

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


static void
print_m8_line (FILE *fp, int exon_querystart, int exon_queryend,
	       char *chr, Chrpos_T exon_genomestart, Chrpos_T exon_genomeend,
	       int nmismatches_bothdiff, Shortread_T headerseq, char *acc_suffix) {
  double identity;
  int alignlength_trim;

  fprintf(fp,"%s%s",Shortread_accession(headerseq),acc_suffix); /* field 0: accession */

  fprintf(fp,"\t%s",chr);	/* field 1: chr */

  /* field 2: identity */
  alignlength_trim = exon_queryend - exon_querystart;
  identity = (double) (alignlength_trim - nmismatches_bothdiff)/(double) alignlength_trim;
  fprintf(fp,"\t%.1f",100.0*identity);


  fprintf(fp,"\t%d",alignlength_trim); /* field 3: query length */

  fprintf(fp,"\t%d",nmismatches_bothdiff); /* field 4: nmismatches */

  fprintf(fp,"\t0");		/* field 5: gap openings */

  /* fields 6 and 7: query start and end */
  fprintf(fp,"\t%d\t%d",exon_querystart,exon_queryend);

  /* fields 8 and 9: chr start and end */
  fprintf(fp,"\t%u\t%u",exon_genomestart,exon_genomeend);

  /* field 10: E value */
  fprintf(fp,"\t%.2g",blast_evalue(alignlength_trim,nmismatches_bothdiff));

 /* field 11: bit score */
  fprintf(fp,"\t%.1f",blast_bitscore(alignlength_trim,nmismatches_bothdiff));
  
  fprintf(fp,"\n");

  return;
}


void
Pair_print_m8 (FILE *fp, struct T *pairs_querydir, int npairs, bool invertedp,
	       Chrnum_T chrnum, Shortread_T queryseq, Shortread_T headerseq,
	       char *acc_suffix, Univ_IIT_T chromosome_iit) {
  bool in_exon = true;
  struct T *pairs, *ptr, *ptr0, *this = NULL;
  int exon_querystart = -1, exon_queryend;
  Chrpos_T exon_genomestart = -1U, exon_genomeend;
  int nmismatches_refdiff, nmismatches_bothdiff, nmatches, i;
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;
  char *chr;
  int querylength;
  bool allocp;

  querylength = Shortread_fulllength(queryseq);

  if (invertedp == true) {
    pairs = invert_and_revcomp_path_and_coords(pairs_querydir,npairs,querylength);
  } else {
    pairs = pairs_querydir;
  }


  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);

  ptr = pairs;
  exon_querystart = ptr->querypos + 1;
  exon_genomestart = ptr->genomepos + 1;
  nmismatches_refdiff = nmismatches_bothdiff = nmatches = 0;

  i = 0;
  while (i < npairs) {
    this = ptr++;
    i++;

    if (this->gapp) {
      if (in_exon == true) {
	/* SPLICE START */
	ptr0 = ptr;
	while (ptr0->gapp) {
	  ptr0++;
	}
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;

	print_m8_line(fp,exon_querystart,exon_queryend,chr,exon_genomestart,exon_genomeend,
		      nmismatches_bothdiff,headerseq,acc_suffix);

	nmismatches_refdiff = nmismatches_bothdiff = nmatches = 0;

	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* May want to print dinucleotides */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* SPLICE CONTINUATION */
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;

	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
	  /* INSERTION */
	  exon_queryend = last_querypos + 1;
	  exon_genomeend = last_genomepos + 1;

	  /* indel_pos = this->querypos; */
	  while (i < npairs && this->gapp == false && this->genome == ' ') {
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  this = ptr;
	  exon_querystart = this->querypos + 1;
	  exon_genomestart = this->genomepos + 1;
	  nmismatches_refdiff = nmismatches_bothdiff = nmatches = 0;

	} else if (this->cdna == ' ') {
	  /* DELETION */
	  exon_queryend = last_querypos + 1;
	  exon_genomeend = last_genomepos + 1;

	  /* indel_pos = this->querypos; */
	  while (i < npairs && this->gapp == false && this->cdna == ' ') {
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  /* Finish rest of this line */
	  print_m8_line(fp,exon_querystart,exon_queryend,chr,exon_genomestart,exon_genomeend,
			nmismatches_bothdiff,headerseq,acc_suffix);

	  this = ptr;
	  exon_querystart = this->querypos + 1;
	  exon_genomestart = this->genomepos + 1;
	  nmismatches_refdiff = nmismatches_bothdiff = nmatches = 0;

	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	/* c = this->genome; */
	if (this->genome == this->cdna) {
	  nmatches++;
	} else if (this->genomealt == this->cdna) {
	  nmismatches_refdiff++;
	} else {
	  nmismatches_bothdiff++;
	  nmismatches_refdiff++;
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;

  print_m8_line(fp,exon_querystart,exon_queryend,chr,exon_genomestart,exon_genomeend,
		nmismatches_bothdiff,headerseq,acc_suffix);

  if (allocp) {
    FREE(chr);
  }

  if (invertedp == true) {
    FREE(pairs);
  }

  return;
}
#endif


/* Modified from print_endtypes */
static void
splice_site_probs (double *sense_prob, double *antisense_prob,
		   bool prev_splicesitep, bool splicesitep, Univcoord_T chroffset,
		   int exon_genomestart, int exon_genomeend, bool watsonp) {

  if (prev_splicesitep == true) {
    if (watsonp == true) {
      /* printf("watsonp is true, so looking up acceptor/antidonor at %u+%u-1\n",chroffset,exon_genomestart); */
      *sense_prob += Maxent_hr_acceptor_prob(chroffset+exon_genomestart-1U,chroffset);
      *antisense_prob += Maxent_hr_antidonor_prob(chroffset+exon_genomestart-1U,chroffset);
    } else {
      /* printf("watsonp is false, so looking up antiacceptor/donor at %u+%u\n",chroffset,exon_genomestart); */
      *sense_prob += Maxent_hr_antiacceptor_prob(chroffset+exon_genomestart,chroffset);
      *antisense_prob += Maxent_hr_donor_prob(chroffset+exon_genomestart,chroffset);
    }
  }

  if (splicesitep == true) {
    if (watsonp == true) {
      /* printf("watsonp is true, so looking up donor/antiacceptor at %u+%u\n",chroffset,exon_genomeend); */
      *sense_prob += Maxent_hr_donor_prob(chroffset+exon_genomeend,chroffset);
      *antisense_prob += Maxent_hr_antiacceptor_prob(chroffset+exon_genomeend,chroffset);
    } else {
      /* printf("watsonp is false, so looking up antiacceptor/donor at %u+%u-1\n",chroffset,exon_genomeend); */
      *sense_prob += Maxent_hr_antidonor_prob(chroffset+exon_genomeend-1U,chroffset);
      *antisense_prob += Maxent_hr_acceptor_prob(chroffset+exon_genomeend-1U,chroffset);
    }
  }
  /* printf("sense %g, antisense %g\n",*sense_prob,*antisense_prob); */

  return;
}


/* Modified from Pair_print_gsnap */
int
Pair_guess_cdna_direction_array (int *sensedir, struct T *pairs_querydir, int npairs, bool invertedp,
				 Univcoord_T chroffset, bool watsonp) {
  double sense_prob = 0.0, antisense_prob = 0.0;
  bool in_exon = true;
  struct T *pairs, *ptr, *this = NULL;
  int i;
  Chrpos_T exon_genomestart = -1U, exon_genomeend;
  Chrpos_T last_genomepos = -1U;
  bool splicesitep, prev_splicesitep;


  if (invertedp == true) {
    fprintf(stderr,"Pair_guess_cdna_direction cannot handle invertedp\n");
    /* pairs = invert_and_revcomp_path_and_coords(pairs_querydir,npairs,querylength); */
    /* watsonp = !watsonp; */
    abort();
  } else {
    pairs = pairs_querydir;
  }

  if (pairs == NULL) {
    *sensedir = SENSE_NULL;
    return 0;
  } else {
    ptr = pairs;
    exon_genomestart = ptr->genomepos + 1;
    splicesitep = false;
  }

  i = 0;
  while (i < npairs) {
    this = ptr++;
    i++;

    if (this->gapp) {
      if (in_exon == true) {
	/* SPLICE START */
#if 0
	ptr0 = ptr;
	while (ptr0->gapp) {
	  ptr0++;
	}
#endif
	exon_genomeend = last_genomepos + 1;

	prev_splicesitep = splicesitep;
	splicesitep = true;

	splice_site_probs(&sense_prob,&antisense_prob,
			  prev_splicesitep,splicesitep,chroffset,
			  exon_genomestart,exon_genomeend,watsonp);

	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* May want to print dinucleotides */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* SPLICE CONTINUATION */
	exon_genomestart = this->genomepos + 1;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
	  /* INSERTION */
	  exon_genomeend = last_genomepos + 1;
	  prev_splicesitep = splicesitep;
	  splicesitep = false;

	  while (i < npairs && this->gapp == false && this->genome == ' ') {
	    this = ptr++;
	    i++;
	  }
	  /* ptr--; */
	  /* i--; */

	  splice_site_probs(&sense_prob,&antisense_prob,
			    prev_splicesitep,splicesitep,chroffset,
			    exon_genomestart,exon_genomeend,watsonp);

	  this = ptr;
	  exon_genomestart = this->genomepos + 1;

	} else if (this->cdna == ' ') {
	  /* DELETION */
	  exon_genomeend = last_genomepos + 1;
	  prev_splicesitep = splicesitep;
	  splicesitep = false;

	  while (i < npairs && this->gapp == false && this->cdna == ' ') {
	    this = ptr++;
	    i++;
	  }
	  /* ptr--; */
	  /* i--; */

	  splice_site_probs(&sense_prob,&antisense_prob,
			    prev_splicesitep,splicesitep,chroffset,
			    exon_genomestart,exon_genomeend,watsonp);

	  this = ptr;
	  exon_genomestart = this->genomepos + 1;

	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      }
    }

    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  exon_genomeend = last_genomepos + 1;
  prev_splicesitep = splicesitep;
  splicesitep = false;

  splice_site_probs(&sense_prob,&antisense_prob,
		    prev_splicesitep,splicesitep,chroffset,
		    exon_genomestart,exon_genomeend,watsonp);

  if (invertedp == true) {
    FREE(pairs);
  }

  if (sense_prob == 0.0 && antisense_prob == 0.0) {
    *sensedir = SENSE_NULL;
    return 0;
  } else if (sense_prob >= antisense_prob) {
    *sensedir = SENSE_FORWARD;
    return +1;
  } else {
    *sensedir = SENSE_ANTI;
    return -1;
  }
}


#if 0
static char
get_genomic_nt_array (char *g_alt, int genomicpos, Univcoord_T chroffset, Univcoord_T chrhigh,
		      bool watsonp) {
  char c2, c2_alt;
  Univcoord_T pos;

  if (watsonp) {
    if ((pos = chroffset + genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      *g_alt = '*';
      return '*';
      
    } else if (pos >= chrhigh) {
      *g_alt = '*';
      return '*';
      
    } else {
      return Genome_get_char_blocks(&(*g_alt),pos);
    }

  } else {
    /* coordinates already processed by Pair_set_genomepos */
    if ((pos = chroffset + genomicpos) < chroffset) { /* Must be <, and not <=, or dynamic programming will fail */
      return '*';
      
    } else if (pos >= chrhigh) {
      return '*';
      
    } else {
      c2 = Genome_get_char_blocks(&c2_alt,pos);
    }
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}
#endif


void
Pair_fix_cdna_direction_array (struct T *pairs_querydir, int npairs, int cdna_direction) {
  struct T *ptr, *this = NULL;
  int i;

  ptr = pairs_querydir;
  i = 0;

  while (i < npairs) {
    this = ptr++;
    i++;

    if (this->gapp && this->comp == NONINTRON_COMP) {
      if (cdna_direction > 0) {
	switch (this->introntype) {
	case GTAG_FWD: this->comp = FWD_CANONICAL_INTRON_COMP; break;
	case GCAG_FWD: this->comp = FWD_GCAG_INTRON_COMP; break;
	case ATAC_FWD: this->comp = FWD_ATAC_INTRON_COMP; break;
	default: this->comp = NONINTRON_COMP;
	}
      } else if (cdna_direction < 0) {
	switch (this->introntype) {
	case ATAC_REV: this->comp = REV_ATAC_INTRON_COMP; break;
	case GCAG_REV: this->comp = REV_GCAG_INTRON_COMP; break;
	case GTAG_REV: this->comp = REV_CANONICAL_INTRON_COMP; break;
	default: this->comp = NONINTRON_COMP; break;
	}
      }
    }
  }

  return;
}




int
Pair_gsnap_nsegments (int *total_nmismatches, int *total_nindels, int *nintrons,
		      int *nindelbreaks, struct T *pairs, int npairs) {
  int nsegments = 0;
  bool in_exon = true;
  struct T *ptr, *ptr0, *this = NULL;
  int i;

  ptr = pairs;
  *total_nmismatches = 0;
  *total_nindels = 0;
  *nintrons = 0;
  *nindelbreaks = 0;

  i = 0;
  while (i < npairs) {
    this = ptr++;
    i++;

    if (this->gapp) {
      if (in_exon == true) {
	/* SPLICE START */
	ptr0 = ptr;
	while (ptr0->gapp) {
	  ptr0++;
	}

	(*nintrons) += 1;
	nsegments++;

	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* May want to print dinucleotides */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* SPLICE CONTINUATION */
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
	  /* INSERTION */
	  while (i < npairs && this->genome == ' ') {
	    (*total_nindels) += 1;
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  (*nindelbreaks) += 1;
	  nsegments++;

	} else if (this->cdna == ' ') {
	  /* DELETION */
	  while (i < npairs && this->cdna == ' ') {
	    (*total_nindels) -= 1;
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  (*nindelbreaks) += 1;
	  nsegments++;

	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	if (this->genome != this->cdna) {
	  (*total_nmismatches) += 1;
	}
      }
    }
  }

  nsegments++;

  return nsegments;
}



/************************************************************************
 *   SAM
 ************************************************************************/

/* Derived from print_tokens_gff3 */
static void
print_tokens_sam (FILE *fp, List_T tokens) {
  List_T p;
  char *token;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    fprintf(fp,"%s",token);
    FREE(token);
  }

  return;
}

/* Derived from print_tokens_gff3 */
static int
tokens_cigarlength (List_T tokens) {
  int length = 0, tokenlength;
  List_T p;
  char *token;
  char type;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    type = token[strlen(token)-1];
    /* Should include 'H', but that gets added according to hardclip_low and hardclip_high */
    if (type == 'S' || type == 'I' || type == 'M') {
      sscanf(token,"%d",&tokenlength);
      length += tokenlength;
    }
  }

  return length;
}



/* Only for GMAP program */
static unsigned int
compute_sam_flag_nomate (int pathnum, int npaths, bool first_read_p, bool watsonp, bool sam_paired_p) {
  unsigned int flag = 0U;

  if (sam_paired_p == true) {
    flag |= PAIRED_READ;
    if (first_read_p == true) {
      flag |= FIRST_READ_P;
    } else {
      flag |= SECOND_READ_P;
    }
  }

  if (npaths == 0) {
    flag |= QUERY_UNMAPPED;
  } else if (watsonp == false) {
    flag |= QUERY_MINUSP;
  }

#if 0
  /* Will let external program decide what is primary */
  if (pathnum > 1) {
    flag |= NOT_PRIMARY;
  }
#endif

  return flag;
}


/* Modeled after Shortread_print_chopped */
static void
print_chopped (FILE *fp, char *contents, int querylength,
	       int hardclip_start, int hardclip_end) {
  int i;

  for (i = hardclip_start; i < querylength - hardclip_end; i++) {
    putc(contents[i],fp);
  }
  return;
}

/* Differs from Shortread version, in that hardclip_high and hardclip_low are not reversed */
static void
print_chopped_revcomp (FILE *fp, char *contents, int querylength,
		       int hardclip_start, int hardclip_end) {
  int i;

  for (i = querylength - 1 - hardclip_end; i >= hardclip_start; --i) {
    putc(complCode[(int) contents[i]],fp);
  }
  return;
}


static void
print_chopped_end (FILE *fp, char *contents, int querylength,
		   int hardclip_start, int hardclip_end) {
  int i;

  if (hardclip_start > 0) {
    for (i = 0; i < hardclip_start; i++) {
      putc(contents[i],fp);
    }
    return;

  } else {
    for (i = querylength - hardclip_end; i < querylength; i++) {
      putc(contents[i],fp);
    }
    return;
  }
}

/* Differs from Shortread version, in that hardclip_high and hardclip_low are not reversed */
static void
print_chopped_end_revcomp (FILE *fp, char *contents, int querylength,
			   int hardclip_start, int hardclip_end) {
  int i;

  if (hardclip_start > 0) {
    for (i = hardclip_start - 1; i >= 0; --i) {
      putc(complCode[(int) contents[i]],fp);
    }
    return;

  } else {
    for (i = querylength - 1; i >= querylength - hardclip_end; --i) {
      putc(complCode[(int) contents[i]],fp);
    }
    return;
  }
}



/* Modeled after Shortread_print_quality */
static void
print_quality (FILE *fp, char *quality, int querylength,
	       int hardclip_start, int hardclip_end, int shift) {
  int i;
  int c;

  if (quality == NULL) {
    putc('*',fp);
  } else {
    for (i = hardclip_start; i < querylength - hardclip_end; i++) {
      if ((c = quality[i] + shift) <= 32) {
	fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		shift,quality[i]);
	abort();
      } else {
	putc(c,fp);
      }
    }
  }
  return;
}


static void
print_quality_revcomp (FILE *fp, char *quality, int querylength,
		       int hardclip_start, int hardclip_end, int shift) {
  int i;
  int c;

  if (quality == NULL) {
    putc('*',fp);
  } else {
    for (i = querylength - 1 - hardclip_end; i >= hardclip_start; --i) {
      if ((c = quality[i] + shift) <= 32) {
	fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		shift,quality[i]);
	abort();
      } else {
	putc(c,fp);
      }
    }
  }

  return;
}



static int
sensedir_from_cdna_direction (int cdna_direction) {
  if (cdna_direction > 0) {
    return SENSE_FORWARD;
  } else if (cdna_direction < 0) {
    return SENSE_ANTI;
  } else {
    return SENSE_NULL;
  }
}


/* Derived from print_gff3_cdna_match */
/* Assumes pairarray has been hard clipped already */
static void
print_sam_line (FILE *fp, char *abbrev, bool first_read_p, char *acc1, char *acc2, char *chrstring,
		bool watsonp, int cdna_direction, List_T cigar_tokens, List_T md_tokens,
		int nmismatches_refdiff, int nmismatches_bothdiff, int nindels,
		bool intronp, char *queryseq_ptr, char *quality_string,
		int hardclip_start, int hardclip_end, int querylength, Chimera_T chimera, int quality_shift,
		int pathnum, int npaths, int absmq_score, int first_absmq, int second_absmq, unsigned int flag,
		Chrnum_T chrnum, Univ_IIT_T chromosome_iit, Chrpos_T chrpos, Chrpos_T chrlength,
#ifdef GSNAP
		Shortread_T queryseq, Resulttype_T resulttype, int pair_mapq_score, int end_mapq_score,
		char *mate_chrstring, Chrnum_T mate_chrnum, Chrnum_T mate_effective_chrnum,
		Chrpos_T mate_chrpos, Chrpos_T mate_chrlength, int mate_cdna_direction, int pairedlength,
#else
		int mapq_score, struct T *pairarray, int npairs,
#endif
		char *sam_read_group_id, bool invertp, bool merged_overlap_p) {
  int sensedir;

  if (cigar_action == CIGAR_ACTION_IGNORE) {
    /* Don't check */
  } else if (tokens_cigarlength(cigar_tokens) + hardclip_start + hardclip_end == querylength) {
    /* Okay */
  } else if (cigar_action == CIGAR_ACTION_WARNING) {
    fprintf(stderr,"Warning: for %s, CIGAR length %d plus hardclips %d and %d do not match sequence length %d\n",
	    acc1,tokens_cigarlength(cigar_tokens),hardclip_start,hardclip_end,querylength);
  } else {
    /* CIGAR_ACTION_ABORT */
    fprintf(stderr,"Error: for %s, CIGAR length %d plus hardclips %d and %d do not match sequence length %d\n",
	    acc1,tokens_cigarlength(cigar_tokens),hardclip_start,hardclip_end,querylength);
    abort();
  }

  /* 1. QNAME or Accession */
  if (acc2 == NULL) {
    fprintf(fp,"%s\t",acc1);
  } else {
    fprintf(fp,"%s,%s\t",acc1,acc2);
  }

  /* 2. Flags */
  fprintf(fp,"%u\t",flag);

  /* 3. RNAME or Chrstring */
  /* 4. POS or Chrlow */
  /* Taken from GMAP part of SAM_chromosomal_pos */
  if (chrpos > chrlength) {
    fprintf(fp,"%s\t%u\t",chrstring,chrpos - chrlength /*+ 1U*/);
  } else {
    fprintf(fp,"%s\t%u\t",chrstring,chrpos /*+ 1U*/);
  }

  /* 5. MAPQ or Mapping quality */
#ifdef GSNAP
  fprintf(fp,"%d\t",pair_mapq_score);
#else
  fprintf(fp,"%d\t",mapq_score);
#endif

  /* 6. CIGAR */
  print_tokens_sam(fp,cigar_tokens);

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
#ifdef GSNAP
  if (mate_chrpos == 0U) {
    fprintf(fp,"\t*\t0");
  } else if (mate_chrpos > mate_chrlength) {
    fprintf(fp,"\t%s\t%u",mate_chrstring,mate_chrpos - mate_chrlength /* +1U*/);
  } else {
    fprintf(fp,"\t%s\t%u",mate_chrstring,mate_chrpos /* +1U*/);
  }
#else
  fprintf(fp,"\t*\t0");
#endif

  /* 9. ISIZE: Insert size */
#ifdef GSNAP
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (watsonp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#else
  fprintf(fp,"\t0");
#endif

  /* 10. SEQ: queryseq and 11. QUAL: quality_scores */
  fprintf(fp,"\t");
  if (watsonp == true) {
    print_chopped(fp,queryseq_ptr,querylength,hardclip_start,hardclip_end);
    fprintf(fp,"\t");
    print_quality(fp,quality_string,querylength,hardclip_start,hardclip_end,
		  quality_shift);
  } else {
    print_chopped_revcomp(fp,queryseq_ptr,querylength,hardclip_start,hardclip_end);
    fprintf(fp,"\t");
    print_quality_revcomp(fp,quality_string,querylength,hardclip_start,hardclip_end,
			  quality_shift);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH */
  if (hardclip_start > 0 || hardclip_end > 0) {
    fprintf(fp,"\tXH:Z:");
    if (watsonp == true) {
      print_chopped_end(fp,queryseq_ptr,querylength,hardclip_start,hardclip_end);
    } else {
      print_chopped_end_revcomp(fp,queryseq_ptr,querylength,hardclip_start,hardclip_end);
    }
  }

#ifdef GSNAP
  if (queryseq != NULL) {
    /* 12. TAGS: XB */
    Shortread_print_barcode(fp,queryseq);

    /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
    Shortread_print_chop(fp,queryseq,invertp);
  }
#endif

  /* 12. TAGS: MD string */
  fprintf(fp,"\tMD:Z:");
  print_tokens_sam(fp,md_tokens);

  /* 12. TAGS: NH */
  fprintf(fp,"\tNH:i:%d",npaths);
  
  /* 12. TAGS: HI */
  fprintf(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\tNM:i:%d",nmismatches_refdiff + nindels);

  if (snps_p) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\tXW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }


  /* 12. TAGS: SM */
#ifdef GSNAP
  fprintf(fp,"\tSM:i:%d",end_mapq_score);
#else
  fprintf(fp,"\tSM:i:%d",40);
#endif

  /* 12. TAGS: XQ */
  fprintf(fp,"\tXQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  if (intronp == true) {
#ifdef GSNAP
    if ((sensedir = sensedir_from_cdna_direction(cdna_direction)) == SENSE_NULL) {
      sensedir = sensedir_from_cdna_direction(mate_cdna_direction);
    }
#else
    sensedir = sensedir_from_cdna_direction(cdna_direction);
#endif

    if (sensedir == SENSE_FORWARD) {
      if (watsonp == true) {
	fprintf(fp,"\tXS:A:+");
      } else {
	fprintf(fp,"\tXS:A:-");
      }

    } else if (sensedir == SENSE_ANTI) {
      if (watsonp == true) {
	fprintf(fp,"\tXS:A:-");
      } else {
	fprintf(fp,"\tXS:A:+");
      }

    } else if (force_xs_direction_p == true) {
      /* Could not determine sense, so just report arbitrarily as + */
      /* This option provided for users of Cufflinks, which cannot handle XS:A:? */
      fprintf(fp,"\tXS:A:+");

    } else {
      /* Non-canonical, so report as such */
      fprintf(fp,"\tXS:A:?");
    }
  }

  /* 12. TAGS: XT */
  if (chimera != NULL) {
    fprintf(fp,"\tXT:Z:");
    Chimera_print_sam_tag(fp,chimera,chromosome_iit);
  }

  /* 12. TAGS: XG */
  if (merged_overlap_p) {
    fprintf(fp,"\tXG:Z:O");
  } else {
    fprintf(fp,"\tXG:Z:M");
  }

  fprintf(fp,"\n");

  return;
}


void
Pair_alias_circular (struct T *pairs, int npairs, Chrpos_T chrlength) {
  int i;
  struct T *ptr;

  i = 0;
  ptr = pairs;
  while (i < npairs) {
    assert(ptr->genomepos < chrlength);
    ptr->genomepos += chrlength;
    i++;
    ptr++;
  }

  return;
}

void
Pair_unalias_circular (struct T *pairs, int npairs, Chrpos_T chrlength) {
  int i;
  struct T *ptr;

  i = 0;
  ptr = pairs;
  while (i < npairs) {
    assert(ptr->genomepos >= chrlength);
    ptr->genomepos -= chrlength;
    i++;
    ptr++;
  }

  return;
}


static struct T *
hardclip_pairs (int *clipped_npairs, int hardclip_start, int hardclip_end,
		struct T *pairs, int npairs, int querylength) {
  struct T *clipped_pairs, *ptr;
  int i, starti;

  debug10(printf("Entered hardclip_pairs with hardclip_start %d, hardclip_end %d, querylength %d\n",
		 hardclip_start,hardclip_end,querylength));
  debug10(Pair_dump_array(pairs,npairs,true));
  debug10(printf("Starting with %d pairs\n",npairs));

  i = 0;
  ptr = pairs;
  while (i < npairs && ptr->querypos < hardclip_start) {
    i++;
    ptr++;
  }
  while (i < npairs && (ptr->gapp == true || ptr->cdna == ' ' || ptr->genome == ' ')) {
    i++;
    ptr++;
  }

  if (i >= npairs) {
    /* hardclip_start passes right end of read, so invalid */
    hardclip_start = 0;
  } else if (hardclip_start > 0) {
    hardclip_start = ptr->querypos;
  }

  starti = i;
  debug10(printf("starti is %d\n",starti));

  clipped_pairs = ptr;

  while (i < npairs && ptr->querypos < querylength - hardclip_end) {
    i++;
    ptr++;
  }

  i--;
  ptr--;
  while (i >= starti && (ptr->gapp == true || ptr->cdna == ' ' || ptr->genome == ' ')) {
    i--;
    ptr--;
  }
  
  if (i < 0) {
    /* hardclip_end passes left end of read, so invalid */
    hardclip_end = 0;
  } else if (hardclip_end > 0) {
    hardclip_end = querylength - 1 - ptr->querypos;
  }

  if (hardclip_start == 0 && hardclip_end == 0) {
    debug10(printf("Unable to hard clip\n"));
    *clipped_npairs = npairs;
    clipped_pairs = pairs;
  } else {
    *clipped_npairs = i - starti + 1;
  }

  debug10(printf("Ending with %d pairs\n",*clipped_npairs));
  debug10(printf("Exiting hardclip_pairs with hardclip_start %d, hardclip_end %d\n",
		 hardclip_start,hardclip_end));

  return clipped_pairs;
}


List_T
Pair_clean_cigar (List_T tokens, bool watsonp) {
  List_T clean, unique = NULL, p;
  char token[10], *curr_token, *last_token;
  int length = 0;
  char type, last_type = ' ';
  bool duplicatep = false;

  for (p = tokens; p != NULL; p = List_next(p)) {
    curr_token = (char *) List_head(p);
    type = curr_token[strlen(curr_token)-1];
    if (type == last_type) {
      length += atoi(last_token);
      FREE(last_token);
      duplicatep = true;
    } else {
      if (last_type == ' ') {
	/* Skip */
      } else if (duplicatep == false) {
	unique = List_push(unique,(void *) last_token);
      } else {
	length += atoi(last_token);
	FREE(last_token);
	sprintf(token,"%d%c",length,last_type);
	unique = push_token(unique,token);
      }
      last_type = type;
      duplicatep = false;
      length = 0;
    }
    last_token = curr_token;
  }
  if (last_type == ' ') {
    /* Skip */
  } else if (duplicatep == false) {
    unique = List_push(unique,(void *) last_token);
  } else {
    length += atoi(last_token);
    FREE(last_token);
    sprintf(token,"%d%c",length,last_type);
    unique = push_token(unique,token);
  }
  List_free(&tokens);


  if (sam_insert_0M_p == false) {
    /* Return result */
    if (watsonp) {
      /* Put tokens in forward order */
      return unique;
    } else {
      /* Keep tokens in reverse order */
      return List_reverse(unique);
    }

  } else {
    /* Insert "0M" between adjacent I and D operations */
    last_type = ' ';
    clean = (List_T) NULL;
    for (p = unique; p != NULL; p = List_next(p)) {
      curr_token = (char *) List_head(p);
      type = curr_token[strlen(curr_token)-1];
      if (last_type == 'I' && type == 'D') {
	clean = push_token(clean,"0M");
      } else if (last_type == 'D' && type == 'I') {
	clean = push_token(clean,"0M");
      }
      clean = List_push(clean,(void *) curr_token);
      last_type = type;
    }
    List_free(&unique);

    /* Return result */
    if (watsonp) {
      /* Put tokens in forward order */
      return List_reverse(clean);
    } else {
      /* Keep tokens in reverse order */
      return clean;
    }
  }
}


static List_T
compute_cigar (bool *intronp, int *hardclip_start, int *hardclip_end, struct T *pairs, int npairs, int querylength_given,
	       bool watsonp, int cdna_direction, int chimera_part) {
  List_T tokens = NULL;
  char token[10];
  int Mlength = 0, Ilength = 0, Dlength = 0;
  bool in_exon = false, deletionp;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1;
  int exon_queryend = -1;
  Chrpos_T exon_genomestart = -1;
  Chrpos_T exon_genomeend, genome_gap;
  Chrpos_T intron_start, intron_end;
  int query_gap;
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;
  int i;

  /* *chimera_hardclip_start = *chimera_hardclip_high = 0; */
  *intronp = false;

  ptr = pairs;

  if (chimera_part == +1) {
    if (ptr->querypos > *hardclip_start) {
      if (ptr->querypos > 0) {
	/* Clip to beginning */
	*hardclip_start = ptr->querypos;
	sprintf(token,"%dH",*hardclip_start);
	tokens = push_token(tokens,token);
      }
    } else {
      if (*hardclip_start > 0) {
	/* Clip to hard clip boundary */
	sprintf(token,"%dH",*hardclip_start);
	tokens = push_token(tokens,token);
      }
    }
  } else {
    if (*hardclip_start > 0) {
      sprintf(token,"%dH",*hardclip_start);
      tokens = push_token(tokens,token);
    }
    if (ptr->querypos > (*hardclip_start)) {
      sprintf(token,"%dS",ptr->querypos - (*hardclip_start));
      tokens = push_token(tokens,token);
    }
  }

  this = (T) NULL;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

#if 0
    print_tokens_sam(stdout,tokens);
    printf("querypos %d, %c %c, exon %u..%u, intron %u..%u\n",
	   this->querypos,this->cdna,this->genome,exon_genomestart,exon_genomeend,
	   intron_start,intron_end);
#endif

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
	
	if (Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(tokens,token);
	} else if (Ilength > 0) {
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(tokens,token);
	} else if (Dlength > 0) {
	  sprintf(token,"%dD",Dlength);
	  tokens = push_token(tokens,token);
	}

	Mlength = Ilength = Dlength = 0;

	in_exon = false;
      }

    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}

	if (prev != NULL) {
	  /* Gap */
	  genome_gap = abs(intron_end - intron_start) + 1;

	  deletionp = false;
#ifdef CONVERT_INTRONS_TO_DELETIONS
	  if (cdna_direction > 0) {
	    if (prev->comp == FWD_CANONICAL_INTRON_COMP ||
		prev->comp == FWD_GCAG_INTRON_COMP ||
		prev->comp == FWD_ATAC_INTRON_COMP) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else {
	      sprintf(token,"%uD",genome_gap);
	      deletionp = true;
	    }
	  } else if (cdna_direction < 0) {
	    if (prev->comp == REV_CANONICAL_INTRON_COMP ||
		prev->comp == REV_GCAG_INTRON_COMP ||
		prev->comp == REV_ATAC_INTRON_COMP) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else {
	      sprintf(token,"%uD",genome_gap);
	      deletionp = true;
	    }
	  } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN){
	    sprintf(token,"%uN",genome_gap);
	    *intronp = true;
	  } else {
	    sprintf(token,"%uD",genome_gap);
	    deletionp = true;
	  }
#else
	  sprintf(token,"%uN",genome_gap);
	  *intronp = true;
#endif
	  tokens = push_token(tokens,token);

	  /* Check for dual gap.  Doesn't work for hard clipping. */
	  assert(exon_queryend >= 0);

	  query_gap = this->querypos - exon_queryend;
	  assert(query_gap >= 0);
	  if (query_gap > 0) {
	    if (deletionp == true && sam_insert_0M_p == true) {
	      /* Put zero matches between deletion and insertion, since some programs will complain */
	      sprintf(token,"0M");
	      tokens = push_token(tokens,token);
	    }

	    sprintf(token,"%uI",query_gap);
	    tokens = push_token(tokens,token);
	  }
	}

	in_exon = true;
      }

      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (this->genome == ' ') {
	  /* Insertion relative to genome */
	  if (Mlength > 0) {
	    sprintf(token,"%dM",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Dlength > 0) {
	    /* unlikely */
	    sprintf(token,"%dD",Dlength);
	    tokens = push_token(tokens,token);
	    Dlength = 0;
	  }
	  Ilength++;
	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (Mlength > 0) {
	    sprintf(token,"%dM",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Ilength > 0) {
	    sprintf(token,"%dI",Ilength);
	    tokens = push_token(tokens,token);
	    Ilength = 0;
	  }
	  Dlength++;
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	/* Count even if unknown base */

	if (Ilength > 0) {
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(tokens,token);
	  Ilength = 0;
	} else if (Dlength > 0) {
	  sprintf(token,"%dD",Dlength);
	  tokens = push_token(tokens,token);
	  Dlength = 0;
	}
	Mlength++;
      }
    }

    if (this != NULL) {
      if (this->cdna != ' ') {
	last_querypos = this->querypos;
      }
      if (this->genome != ' ') {
	last_genomepos = this->genomepos;
      }
    }
  }

  /* prev = this; */
  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;

  if (Mlength > 0) {
    sprintf(token,"%dM",Mlength);
    tokens = push_token(tokens,token);
  } else if (Ilength > 0) {
    sprintf(token,"%dI",Ilength);
    tokens = push_token(tokens,token);
  } else if (Dlength > 0) {
    sprintf(token,"%dD",Dlength);
    tokens = push_token(tokens,token);
  }


  /* Terminal clipping */
  if (chimera_part == -1) {
    if (last_querypos < querylength_given - 1 - (*hardclip_end)) {
      if (last_querypos < querylength_given - 1) {
	/* Clip to end */
	*hardclip_end = querylength_given - 1 - last_querypos;
	sprintf(token,"%dH",*hardclip_end);
	tokens = push_token(tokens,token);
      }
    } else {
      if (*hardclip_end > 0) {
	/* Clip to hard clip boundary */
	sprintf(token,"%dH",*hardclip_end);
	tokens = push_token(tokens,token);
      }
    }
  } else {
    if (last_querypos < querylength_given - 1 - (*hardclip_end)) {
      sprintf(token,"%dS",querylength_given - 1 - (*hardclip_end) - last_querypos);
      tokens = push_token(tokens,token);
    }
    if (*hardclip_end > 0) {
      sprintf(token,"%dH",*hardclip_end);
      tokens = push_token(tokens,token);
    }
  }

  return Pair_clean_cigar(tokens,watsonp);
}


/* Copied from samprint.c */
static bool
check_cigar_types (Intlist_T cigar_types) {
  Intlist_T p;
  int type, last_type = 'M';
  bool M_present_p = false;

  for (p = cigar_types; p != NULL; p = Intlist_next(p)) {
    type = Intlist_head(p);
    if (type == 'M') {
      M_present_p = true;
#if 0
    } else if (type == 'H' && last_type == 'S') {
      debug1(printf("check_cigar_types detects adjacent S and H, so returning false\n"));
      return false;
    } else if (type == 'S' && last_type == 'H') {
      debug1(printf("check_cigar_types detects adjacent S and H, so returning false\n"));
      return false;
#endif
    }
  }

  return M_present_p;
}


bool
Pair_check_cigar (struct T *pairs, int npairs, int querylength_given,
		  int clipdir, int hardclip5, int hardclip3,
		  bool watsonp, int cdna_direction, bool first_read_p, bool circularp) {
  bool result;
  Intlist_T cigar_types = NULL;
  int hardclip_low, hardclip_high;
  int Mlength = 0, Ilength = 0, Dlength = 0;
  bool in_exon = false, deletionp;
  struct T *ptr, *prev, *this = NULL;
  int exon_queryend;
  int query_gap;
  int last_querypos = -1;
  int i;

  if (circularp == true) {
    if (watsonp == true) {
      hardclip_low = hardclip5;
      hardclip_high = hardclip3;
    } else {
      hardclip_low = hardclip3;
      hardclip_high = hardclip5;
    }
  } else {
    /* Incoming hardclip5 and hardclip3 are due to overlaps, not chimera */
    if (clipdir >= 0) {
      if (watsonp == true) {
	if (first_read_p == true) {
	  hardclip_high = hardclip5;
	  hardclip_low = 0;
	} else {
	  hardclip_high = 0;
	  hardclip_low = hardclip3;
	}
      } else {
	if (first_read_p == true) {
	  hardclip_low = hardclip5;
	  hardclip_high = 0;
	} else {
	  hardclip_low = 0;
	  hardclip_high = hardclip3;
	}
      }
    } else {
      if (watsonp == true) {
	if (first_read_p == true) {
	  hardclip_low = hardclip5;
	  hardclip_high = 0;
	} else {
	  hardclip_low = 0;
	  hardclip_high = hardclip3;
	}
      } else {
	if (first_read_p == true) {
	  hardclip_high = hardclip5;
	  hardclip_low = 0;
	} else {
	  hardclip_high = 0;
	  hardclip_low = hardclip3;
	}
      }
    }
  }


  ptr = pairs;

#if 0
  /* This procedure is used to check circular alignments */
  if (chimera_part == +1) {
    if (ptr->querypos > hardclip_low) {
      if (ptr->querypos > 0) {
	/* Clip to beginning */
	hardclip_low = ptr->querypos;
	cigar_types = Intlist_push(cigar_types,'H');
      }
    } else {
      if (hardclip_low > 0) {
	/* Clip to hard clip boundary */
	cigar_types = Intlist_push(cigar_types,'H');
      }
    }
  } else {
#endif
    if (hardclip_low > 0) {
      cigar_types = Intlist_push(cigar_types,'H');
    }
    if (ptr->querypos > hardclip_low) {
      cigar_types = Intlist_push(cigar_types,'S');
    }
#if 0
  }
#endif

  this = (T) NULL;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
#if 0
	exon_genomeend = last_genomepos + 1;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
#endif
	
	if (Mlength > 0) {
	  cigar_types = Intlist_push(cigar_types,'M');
	} else if (Ilength > 0) {
	  cigar_types = Intlist_push(cigar_types,'I');
	} else if (Dlength > 0) {
	  cigar_types = Intlist_push(cigar_types,'D');
	}

	Mlength = Ilength = Dlength = 0;

	in_exon = false;
      }

    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
#if 0
	/* Needed only for full token */
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}
#endif

	if (prev != NULL) {
	  /* Gap */
	  /* genome_gap = abs(intron_end - intron_start) + 1; */

	  deletionp = false;
#ifdef CONVERT_INTRONS_TO_DELETIONS
	  if (cdna_direction > 0) {
	    if (prev->comp == FWD_CANONICAL_INTRON_COMP ||
		prev->comp == FWD_GCAG_INTRON_COMP ||
		prev->comp == FWD_ATAC_INTRON_COMP) {
	      cigar_types = Intlist_push(cigar_types,'N');
	      /* *intronp = true; */
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      cigar_types = Intlist_push(cigar_types,'N');
	      /* *intronp = true; */
	    } else {
	      cigar_types = Intlist_push(cigar_types,'D');
	      deletionp = true;
	    }
	  } else if (cdna_direction < 0) {
	    if (prev->comp == REV_CANONICAL_INTRON_COMP ||
		prev->comp == REV_GCAG_INTRON_COMP ||
		prev->comp == REV_ATAC_INTRON_COMP) {
	      cigar_types = Intlist_push(cigar_types,'N');
	      /* *intronp = true; */
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      cigar_types = Intlist_push(cigar_types,'N');
	      /* *intronp = true; */
	    } else {
	      cigar_types = Intlist_push(cigar_types,'D');
	      deletionp = true;
	    }
	  } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN){
	    cigar_types = Intlist_push(cigar_types,'N');
	    /* *intronp = true; */
	  } else {
	    cigar_types = Intlist_push(cigar_types,'D');
	    deletionp = true;
	  }
#else
	  cigar_types = Intlist_push(cigar_types,'N');
	  /* *intronp = true; */
#endif

	  /* Check for dual gap.  Doesn't work for hard clipping. */
	  assert(exon_queryend >= 0);

	  query_gap = this->querypos - exon_queryend;
	  assert(query_gap >= 0);
	  if (query_gap > 0) {
	    if (deletionp == true && sam_insert_0M_p == true) {
	      /* Put zero matches between deletion and insertion, since some programs will complain */
	      cigar_types = Intlist_push(cigar_types,'M');
	    }

	    cigar_types = Intlist_push(cigar_types,'I');
	  }
	}

	in_exon = true;
      }

      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (this->genome == ' ') {
	  /* Insertion relative to genome */
	  if (Mlength > 0) {
	    cigar_types = Intlist_push(cigar_types,'M');
	    Mlength = 0;
	  } else if (Dlength > 0) {
	    /* unlikely */
	    cigar_types = Intlist_push(cigar_types,'D');
	    Dlength = 0;
	  }
	  Ilength++;
	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (Mlength > 0) {
	    cigar_types = Intlist_push(cigar_types,'M');
	    Mlength = 0;
	  } else if (Ilength > 0) {
	    cigar_types = Intlist_push(cigar_types,'I');
	    Ilength = 0;
	  }
	  Dlength++;
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	/* Count even if unknown base */

	if (Ilength > 0) {
	  cigar_types = Intlist_push(cigar_types,'I');
	  Ilength = 0;
	} else if (Dlength > 0) {
	  cigar_types = Intlist_push(cigar_types,'D');
	  Dlength = 0;
	}
	Mlength++;
      }
    }

    if (this != NULL) {
      if (this->cdna != ' ') {
	last_querypos = this->querypos;
      }
#if 0
      if (this->genome != ' ') {
	last_genomepos = this->genomepos;
      }
#endif
    }
  }

  /* prev = this; */
  exon_queryend = last_querypos + 1;
  /* exon_genomeend = last_genomepos + 1; */

  if (Mlength > 0) {
    cigar_types = Intlist_push(cigar_types,'M');
  } else if (Ilength > 0) {
    cigar_types = Intlist_push(cigar_types,'I');
  } else if (Dlength > 0) {
    cigar_types = Intlist_push(cigar_types,'D');
  }


  /* Terminal clipping */
#if 0
  /* This procedure is used to check circular alignments */
  if (chimera_part == -1) {
    if (last_querypos < querylength_given - 1 - hardclip_high) {
      if (last_querypos < querylength_given - 1) {
	/* Clip to end */
	hardclip_high = querylength_given - 1 - last_querypos;
	cigar_types = Intlist_push(cigar_types,'H');
      }
    } else {
      if (hardclip_high > 0) {
	/* Clip to hard clip boundary */
	cigar_types = Intlist_push(cigar_types,'H');
      }
    }
  } else {
#endif
    if (last_querypos < querylength_given - 1 - hardclip_high) {
      cigar_types = Intlist_push(cigar_types,'S');
    }
    if (hardclip_high > 0) {
      cigar_types = Intlist_push(cigar_types,'H');
    }
#if 0
  }
#endif

  result = check_cigar_types(cigar_types);

  Intlist_free(&cigar_types);
  return result;
}



typedef enum {IN_MATCHES, IN_MISMATCHES, IN_DELETION} MD_state_T;

#if 0
static void
state_print (MD_state_T state) {
  switch (state) {
  case IN_MATCHES: printf("IN_MATCHES"); break;
  case IN_MISMATCHES: printf("IN_MISMATCHES"); break;
  case IN_DELETION: printf("IN_DELETION"); break;
  default: abort();
  }
  return;
}
#endif


#if 0
static List_T
compute_md_string_old (int *nmismatches, struct T *pairs, int npairs, bool watsonp) {
  List_T tokens = NULL;
  char token[10], *first_token;
  int nmatches = 0;
  struct T *ptr, *prev, *this = NULL;
  MD_state_T state = IN_MISMATCHES;
  int i;

  ptr = pairs;
  *nmismatches = 0;

  /* Ignore initial soft clipping */

  if (watsonp == true) {
    for (i = 0; i < npairs; i++) {
      prev = this;
      this = ptr++;

      if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	nmatches++;
	state = IN_MATCHES;

      } else if (this->comp == MISMATCH_COMP) {
	*nmismatches += 1;
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    tokens = push_token(tokens,token);
	    nmatches = 0;
	  }

	} else if (state == IN_DELETION) {
	  tokens = push_token(tokens,"0");
	}
	state = IN_MISMATCHES;

	sprintf(token,"%c",watsonp ? this->genome : complCode[(int) this->genome]);
	tokens = push_token(tokens,token);

      } else if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
#if 0
	  /* Insertion relative to genome.  Ignored in MD string (but not in cigar). */
	  nmatches++;
	  state = IN_MATCHES;
#endif

	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (state == IN_MATCHES) {
	    if (nmatches > 0) {
	      sprintf(token,"%d",nmatches);
	      tokens = push_token(tokens,token);
	      nmatches = 0;
	    }
	    tokens = push_token(tokens,"^");

	  } else if (state == IN_MISMATCHES) {
	    tokens = push_token(tokens,"^");

	  }
	  state = IN_DELETION;

	  sprintf(token,"%c",watsonp ? this->genome : complCode[(int) this->genome]);
	  tokens = push_token(tokens,token);
	}

      } else {
	/* Ignore */
      }
    }

    /* Ignore terminal soft clipping */

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      tokens = push_token(tokens,token);
    }

    /* Put tokens in forward order */
    tokens = List_reverse(tokens);

  } else {

    for (i = 0; i < npairs; i++) {
      prev = this;
      this = ptr++;

      if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	if (state == IN_DELETION) {
	  tokens = push_token(tokens,"^");
	}
	nmatches++;
	state = IN_MATCHES;

      } else if (this->comp == MISMATCH_COMP) {
	*nmismatches += 1;
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    tokens = push_token(tokens,token);
	    nmatches = 0;
	  }

	} else if (state == IN_DELETION) {
	  tokens = push_token(tokens,"^");
	}
	state = IN_MISMATCHES;

	sprintf(token,"%c",watsonp ? this->genome : complCode[(int) this->genome]);
	tokens = push_token(tokens,token);

      } else if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
#if 0
	  /* Insertion relative to genome.  Ignored in MD string, but not in cigar string. */
	  if (state == IN_DELETION) {
	    tokens = push_token(tokens,"^");
	  }
	  nmatches++;
	  state = IN_MATCHES;
#endif

	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (state == IN_MATCHES) {
	    if (nmatches > 0) {
	      sprintf(token,"%d",nmatches);
	      tokens = push_token(tokens,token);
	      nmatches = 0;
	    }

	  } else if (state == IN_MISMATCHES) {
	    tokens = push_token(tokens,"0");

	  }
	  state = IN_DELETION;

	  sprintf(token,"%c",watsonp ? this->genome : complCode[(int) this->genome]);
	  tokens = push_token(tokens,token);
	}

      } else {
	/* Ignore */
      }
    }

    /* Ignore terminal soft clipping */

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      tokens = push_token(tokens,token);
    }
    
    /* Keep tokens in reverse order */
  }


  /* Insert initial 0 token if necessary */
  if (tokens != NULL) {
    first_token = (char *) List_head(tokens);
    if (!isdigit(first_token[0])) {
      tokens = push_token(tokens,"0");
    }
  }
  
  return tokens;
}
#endif


static List_T
compute_md_string (int *nmismatches_refdiff, int *nmismatches_bothdiff, int *nindels,
		   struct T *pairs, int npairs, bool watsonp, List_T cigar_tokens) {
  List_T md_tokens = NULL, p;
  char *cigar_token, token[10], *first_token, type;
  Pair_T this;
  int nmatches = 0, length;
  MD_state_T state = IN_MISMATCHES;
  int i, k = 0;

  *nmismatches_refdiff = *nmismatches_bothdiff = *nindels = 0;

  debug4(Pair_dump_array(pairs,npairs,true));
  debug4(printf("watsonp %d\n",watsonp));

  if (watsonp == true) {
    for (p = cigar_tokens; p != NULL; p = List_next(p)) {
      cigar_token = (char *) List_head(p);
      debug4(printf("token is %s\n",cigar_token));
      type = cigar_token[strlen(cigar_token)-1];
      length = atoi(cigar_token);
    
      if (type == 'H') {
	/* k += length; */

      } else if (type == 'S') {
	/* k += length; */

      } else if (type == 'M') {
	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  debug4(printf("M %d/%d comp %c\n",i,length,this->comp));
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	    nmatches++;
	    state = IN_MATCHES;

	  } else if (this->comp == MISMATCH_COMP) {
	    if (state == IN_MATCHES) {
	      sprintf(token,"%d",nmatches);
	      md_tokens = push_token(md_tokens,token);
	      nmatches = 0;
	    } else if (state == IN_DELETION) {
	      md_tokens = push_token(md_tokens,"0");
	    }
	    state = IN_MISMATCHES;

	    *nmismatches_refdiff += 1;
	    if (md_lowercase_variant_p && this->cdna == this->genomealt) {
	      /* A mismatch against the reference only => alternate variant */
	      sprintf(token,"%c",tolower(this->genome));
	    } else {
	      /* A true mismatch against both variants */
	      *nmismatches_bothdiff += 1;
	      sprintf(token,"%c",this->genome);
	    }
	    md_tokens = push_token(md_tokens,token);

	  } else {
	    fprintf(stderr,"Unexpected comp '%c'\n",this->comp);
	    exit(9);
	  }
	}

      } else if (type == 'I') {
	while (k < npairs && pairs[k].comp == INDEL_COMP && pairs[k].genome == ' ') {
	  *nindels += 1;
	  k++;
	}
	state = IN_MATCHES;

      } else if (type == 'N') {
	while (k < npairs && pairs[k].gapp == true) {
	  k++;
	}

      } else if (type == 'D') {
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    md_tokens = push_token(md_tokens,token);
	    nmatches = 0;
	  }
	}

	if (state != IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}
	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  sprintf(token,"%c",this->genome);
	  md_tokens = push_token(md_tokens,token);
	  *nindels += 1;
	}

	state = IN_DELETION;

      } else {
	fprintf(stderr,"Don't recognize type %c\n",type);
	abort();
      }
    }

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      md_tokens = push_token(md_tokens,token);
    }

    md_tokens = List_reverse(md_tokens);

  } else {
    cigar_tokens = List_reverse(cigar_tokens);
    for (p = cigar_tokens; p != NULL; p = List_next(p)) {
      cigar_token = (char *) List_head(p);
      debug4(printf("token is %s\n",cigar_token));
      type = cigar_token[strlen(cigar_token)-1];
      length = atoi(cigar_token);
    
      if (type == 'H') {
	/* k += length; */

      } else if (type == 'S') {
	/* k += length; */

      } else if (type == 'M') {
	if (state == IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}

	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  debug4(printf("M %d/%d comp %c\n",i,length,this->comp));
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	    nmatches++;
	    state = IN_MATCHES;

	  } else if (this->comp == MISMATCH_COMP) {
	    if (state == IN_MATCHES) {
	      sprintf(token,"%d",nmatches);
	      md_tokens = push_token(md_tokens,token);
	      nmatches = 0;
	    }
	    state = IN_MISMATCHES;

	    *nmismatches_refdiff += 1;

	    if (md_lowercase_variant_p && this->cdna == this->genomealt) {
	      /* A mismatch against the reference only => alternate variant */
	      sprintf(token,"%c",tolower(complCode[(int) this->genome]));
	    } else {
	      *nmismatches_bothdiff += 1;
	      sprintf(token,"%c",complCode[(int) this->genome]);
	    }
	    md_tokens = push_token(md_tokens,token);


	  } else {
	    fprintf(stderr,"Unexpected comp '%c'\n",this->comp);
	    abort();
	  }
	}

      } else if (type == 'I') {
	if (state == IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}

	while (k < npairs && pairs[k].comp == INDEL_COMP && pairs[k].genome == ' ') {
	  *nindels += 1;
	  k++;
	}
	state = IN_MATCHES;

      } else if (type == 'N') {
#if 0
	/* Ignore deletion adjacent to intron, to avoid double ^^ */
	if (state == IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}
#endif

	while (k < npairs && pairs[k].gapp == true) {
	  k++;
	}

      } else if (type == 'D') {
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    md_tokens = push_token(md_tokens,token);
	    nmatches = 0;
	  }
	} else if (state == IN_MISMATCHES) {
	  md_tokens = push_token(md_tokens,"0");
	}

	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  sprintf(token,"%c",complCode[(int) this->genome]);
	  md_tokens = push_token(md_tokens,token);
	  *nindels += 1;
	}
	state = IN_DELETION;

      } else {
	fprintf(stderr,"Don't recognize type %c\n",type);
	abort();
      }
    }

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      md_tokens = push_token(md_tokens,token);
    }

    /* Restore cigar_tokens */
    cigar_tokens = List_reverse(cigar_tokens);
  }

  assert(k == npairs);

  /* Insert initial 0 token if necessary */
  if (md_tokens != NULL) {
    first_token = (char *) List_head(md_tokens);
    if (!isdigit(first_token[0])) {
      md_tokens = push_token(md_tokens,"0");
    }
  }

  return md_tokens;
}




void
Pair_print_sam (FILE *fp, char *abbrev, struct T *pairs, int npairs,
		char *acc1, char *acc2, Chrnum_T chrnum, Univ_IIT_T chromosome_iit, Sequence_T usersegment,
		char *queryseq_ptr, char *quality_string,
		int clipdir, int hardclip_low, int hardclip_high, int querylength_given,
		bool watsonp, int cdna_direction, int chimera_part, Chimera_T chimera,
		int quality_shift, bool first_read_p, int pathnum, int npaths,
		int absmq_score, int first_absmq, int second_absmq, Chrpos_T chrpos, Chrpos_T chrlength,
#ifdef GSNAP
		Shortread_T queryseq, Resulttype_T resulttype, unsigned int flag,
		int pair_mapq_score, int end_mapq_score,
		Chrnum_T mate_chrnum, Chrnum_T mate_effective_chrnum,
		Chrpos_T mate_chrpos, Chrpos_T mate_chrlength,
		int mate_cdna_direction, int pairedlength,
#else
		int mapq_score, bool sam_paired_p,
#endif
		char *sam_read_group_id, bool invertp, bool circularp, bool merged_overlap_p) {
  char *chrstring = NULL;
#ifdef GSNAP
  char *mate_chrstring, *mate_chrstring_alloc = NULL;
#else
  unsigned int flag;
#endif

  List_T cigar_tokens = NULL, md_tokens = NULL;
  int nmismatches_refdiff, nmismatches_bothdiff, nindels;
  bool intronp, ignore_intronp;
  int hardclip_start, hardclip_end;
  int hardclip_start_zero = 0, hardclip_end_zero = 0;
  struct T *clipped_pairs;
  int clipped_npairs;


  if (chrnum == 0) {
    chrstring = Sequence_accession(usersegment);
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

#ifdef GSNAP
  if (mate_chrpos == 0U) {
    mate_chrstring = "*";
  } else if (mate_chrnum == 0) {
    if (/* chrpos > 0U && chrnum > 0 && */ mate_effective_chrnum == chrnum) {
      mate_chrstring = "=";
    } else {
      mate_chrstring = mate_chrstring_alloc = Chrnum_to_string(mate_effective_chrnum,chromosome_iit);
    }
  } else {
    if (/* chrpos > 0U && chrnum > 0 && */ mate_chrnum == chrnum) {
      mate_chrstring = "=";
    } else {
      mate_chrstring = mate_chrstring_alloc = Chrnum_to_string(mate_chrnum,chromosome_iit);
    }
  }
#else
  flag = compute_sam_flag_nomate(pathnum,npaths,first_read_p,watsonp,sam_paired_p);
#endif

  debug4(printf("Entered Pair_print_sam with clipdir %d, watsonp %d, first_read_p %d, hardclip5 %d, and hardclip3 %d\n",
		clipdir,watsonp,first_read_p,hardclip5,hardclip3));

  if (watsonp == true) {
    hardclip_start = hardclip_low;
    hardclip_end = hardclip_high;
  } else {
    hardclip_start = hardclip_high;
    hardclip_end = hardclip_low;
  }
  debug4(printf("hardclip_start %d, hardclip_end %d\n",hardclip_start,hardclip_end));


  /* Get CIGAR and intronp for entire read */
  cigar_tokens = compute_cigar(&intronp,&hardclip_start_zero,&hardclip_end_zero,pairs,npairs,querylength_given,
			       watsonp,cdna_direction,chimera_part);
  if (hardclip_start == 0 && hardclip_end == 0) {
    clipped_pairs = pairs;
    clipped_npairs = npairs;
  } else {
    clipped_pairs = hardclip_pairs(&clipped_npairs,hardclip_start,hardclip_end,
				   pairs,npairs,querylength_given);
  }
  tokens_free(&cigar_tokens);

  /* Cigar updates hardclip5 and hardclip3 for chimeras */
  cigar_tokens = compute_cigar(&ignore_intronp,&hardclip_start,&hardclip_end,clipped_pairs,clipped_npairs,querylength_given,
			       watsonp,cdna_direction,chimera_part);

  md_tokens = compute_md_string(&nmismatches_refdiff,&nmismatches_bothdiff,&nindels,
				clipped_pairs,clipped_npairs,watsonp,cigar_tokens);

  print_sam_line(fp,abbrev,first_read_p,acc1,acc2,chrstring,
		 watsonp,cdna_direction,cigar_tokens,md_tokens,
		 nmismatches_refdiff,nmismatches_bothdiff,nindels,
		 intronp,queryseq_ptr,quality_string,hardclip_start,hardclip_end,
		 querylength_given,chimera,quality_shift,pathnum,npaths,
		 absmq_score,first_absmq,second_absmq,flag,
		 chrnum,chromosome_iit,chrpos,chrlength,
#ifdef GSNAP
		 queryseq,resulttype,pair_mapq_score,end_mapq_score,mate_chrstring,
		 mate_chrnum,mate_effective_chrnum,mate_chrpos,mate_chrlength,
		 mate_cdna_direction,pairedlength,
#else
		 mapq_score,clipped_pairs,clipped_npairs,
#endif
		 sam_read_group_id,invertp,merged_overlap_p);

  /* Print procedures free the character strings */
  List_free(&md_tokens);
  List_free(&cigar_tokens);

#ifdef GSNAP
  if (mate_chrstring_alloc != NULL) {
    FREE(mate_chrstring_alloc);
  }
#endif
  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Pair_print_sam_nomapping (FILE *fp, char *abbrev, char *acc1, char *acc2, char *queryseq_ptr,
			  char *quality_string, int querylength, int quality_shift,
			  bool first_read_p, bool sam_paired_p, char *sam_read_group_id) {
  unsigned int flag;

#ifdef GSNAP
  fprintf(stderr,"Unexpected call to Pair_print_sam_nomapping in GSNAP\n");
  abort();
#endif

  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }
  
  /* 2. FLAG */
  flag = compute_sam_flag_nomate(/*pathnum*/0,/*npaths*/0,first_read_p,/*watsonp*/true,sam_paired_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  fprintf(fp,"\t*");

  /* 4. POS: chrpos */
  fprintf(fp,"\t0");

  /* 5. MAPQ: Mapping quality */
  /* Picard says MAPQ should be 0 for an unmapped read */
  fprintf(fp,"\t0");

  /* 6. CIGAR */
  fprintf(fp,"\t*");

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  /* 9. ISIZE: Insert size */
  fprintf(fp,"\t*\t0\t0\t");

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  print_chopped(fp,queryseq_ptr,querylength,/*hardclip_start*/0,/*hardclip_end*/0);
  fprintf(fp,"\t");
  print_quality(fp,quality_string,querylength,/*hardclip_start*/0,/*hardclip_end*/0,
		quality_shift);

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }
  
  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  fprintf(fp,"\n");

  return;
}


#endif



Uintlist_T
Pair_exonbounds (struct T *pairs, int npairs, Univcoord_T chroffset) {
  Uintlist_T exonbounds = NULL;
  struct T *ptr, *this = NULL;
  bool in_exon = false;
  int i;
  Chrpos_T last_genomepos = -1U;

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    /* prev = this; */
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* exon genomeend */
	exonbounds = Uintlist_push(exonbounds,chroffset + last_genomepos);
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* exon genomestart */
	exonbounds = Uintlist_push(exonbounds,chroffset + this->genomepos);
	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  /* prev = this; */
  exonbounds = Uintlist_push(exonbounds,chroffset + last_genomepos);

  return Uintlist_reverse(exonbounds);
}


static int
count_psl_blocks_nt (Intlist_T *blockSizes, Intlist_T *qStarts, Uintlist_T *tStarts, struct T *pairs_directional,
		     int npairs, int querylength, bool watsonp) {
  int nblocks = 0, i;
  int block_querystart, block_queryend;
  struct T *ptr = pairs_directional, *this = NULL;
  bool in_block = false;
  int last_querypos = -1;
  /* Chrpos_T last_genomepos = -1U; */

  for (i = 0; i < npairs; i++) {
    /* prev = this; */
    this = ptr++;

    if (this->gapp) {
      if (in_block == true) {
	nblocks++;
	block_queryend = last_querypos;
	debug2(fprintf(fp,"Block size: %d\n",abs(block_queryend-block_querystart)+1));
	*blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);
	in_block = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else if (this->cdna == ' ' || this->genome == ' ') {
      if (in_block == true) {
	nblocks++;
	block_queryend = last_querypos;
	debug2(fprintf(fp,"Block size: %d\n",abs(block_queryend-block_querystart)+1));
	*blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);
	in_block = false;
      }

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
         or SHORTGAP_COMP */
      if (in_block == false) {
	block_querystart = this->querypos;
	if (watsonp == true) {
	  debug2(fprintf(fp,"Pushing qstart: %d\n",block_querystart));
	  *qStarts = Intlist_push(*qStarts,block_querystart);
	} else {
	  debug2(fprintf(fp,"Pushing qstart: %d\n",querylength-block_querystart-1));
	  *qStarts = Intlist_push(*qStarts,querylength-block_querystart-1);
	}
	*tStarts = Uintlist_push(*tStarts,this->genomepos);
	in_block = true;
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
#if 0
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
#endif
  }

  if (in_block == true) {
    /* prev = this; */
    nblocks++;
    block_queryend = last_querypos;
    debug2(fprintf(fp,"Block size: %d\n",abs(block_queryend-block_querystart)+1));
    *blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);
  }

  *blockSizes = Intlist_reverse(*blockSizes);
  *qStarts = Intlist_reverse(*qStarts);
  *tStarts = Uintlist_reverse(*tStarts);

  return nblocks;
}


static int
count_psl_blocks_pro (Intlist_T *blockSizes, Intlist_T *qStarts, Uintlist_T *tStarts, struct T *pairs_directional,
		      int npairs, bool watsonp, Chrpos_T chrlength) {
  int nblocks = 0, i;
  int naminoacids = 0;
  int block_querystart;
  struct T *ptr = pairs_directional, *this = NULL;
  bool in_block = false;
#ifdef NOGAPSINBLOCK
  struct T *prev;
#endif

  for (i = 0; i < npairs; i++) {
#ifdef NOGAPSINBLOCK
    prev = this;
#endif
    this = ptr++;

    if (this->gapp) {
      if (in_block == true) {
	nblocks++;
	*blockSizes = Intlist_push(*blockSizes,naminoacids);
	in_block = false;
	naminoacids = 0;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

#ifdef NOGAPSINBLOCK
    } else if (this->cdna == ' ' || this->genome == ' ') {
      if (in_block == true) {
	nblocks++;
	block_queryend = last_querypos;
	*blockSizes = Intlist_push(*blockSizes,block_queryend/3-(block_querystart+2)/3+1);
	in_block = false;
      }
#endif

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
         or SHORTGAP_COMP */
      if (this->aa_e != ' ') {
	naminoacids++;
      }
      if (in_block == false) {
	block_querystart = this->querypos;
	*qStarts = Intlist_push(*qStarts,(block_querystart+2)/3);
	if (watsonp == true) {
	  *tStarts = Uintlist_push(*tStarts,this->genomepos);
	} else {
#if 0
	  /* Should be this */
	  *tStarts = Uintlist_push(*tStarts,this->genomepos);
#else
	  /* But is actually this */
	  *tStarts = Uintlist_push(*tStarts,chrlength - this->genomepos - 1);
#endif
	}
	in_block = true;
      }
    }
  }

  if (in_block == true) {
#ifdef NOGAPSINBLOCK
    prev = this;
#endif
    nblocks++;
    *blockSizes = Intlist_push(*blockSizes,naminoacids);
  }

  *blockSizes = Intlist_reverse(*blockSizes);
  *qStarts = Intlist_reverse(*qStarts);
  *tStarts = Uintlist_reverse(*tStarts);

  return nblocks;
}


static void
compute_gap_lengths_int (int *nbreaks, int *length, Intlist_T blockSizes, Intlist_T Starts, int nblocks) {
  int i;
  int start, end;
  /* Intlist_T p = blockSizes, q = Starts; */

  debug2(fprintf(fp,"Entered compute_gap_lengths_int with nblocks = %d, and Starts having length %d\n",
		nblocks,Intlist_length(Starts)));
  *nbreaks = *length = 0;
  for (i = 0; i < nblocks - 1; i++) {
    if (i > 0) {
      start = Intlist_head(Starts);
      if (start - end > 0) {
	*nbreaks += 1;
	*length += (start - end);
      }
      debug2(fprintf(fp,"%d - %d = %d, gap = %d\n",start,end,start-end,*length));
    }
    end = Intlist_head(Starts) + Intlist_head(blockSizes);
    blockSizes = Intlist_next(blockSizes);
    Starts = Intlist_next(Starts);
  }

  if (i > 0) {
    start = Intlist_head(Starts);
    if (start - end > 0) {
      *nbreaks += 1;
      *length += (start - end);
    }
    debug2(fprintf(fp,"%d - %d = %d, gap = %d\n",start,end,start-end,*length));
  }

  return;
}

static void
compute_gap_lengths_uint (int *nbreaks, int *length, Intlist_T blockSizes, Uintlist_T Starts, int nblocks) {
  int i;
  int start, end;
  /*
  Intlist_T p = blockSizes;
  Uintlist_T q = Starts;
  */

  *nbreaks = *length = 0;
  for (i = 0; i < nblocks - 1; i++) {
    if (i > 0) {
      start = Uintlist_head(Starts);
      if (start - end > 0) {
	*nbreaks += 1;
	*length += (start - end);
      }
      debug2(fprintf(fp,"%d - %d = %d, gap = %d\n",start,end,start-end,*length));
    }
    end = Uintlist_head(Starts) + Intlist_head(blockSizes);
    blockSizes = Intlist_next(blockSizes);
    Starts = Uintlist_next(Starts);
  }

  if (i > 0) {
    start = Uintlist_head(Starts);
    if (start - end > 0) {
      *nbreaks += 1;
      *length += (start - end);
    }
    debug2(fprintf(fp,"%d - %d = %d, gap = %d\n",start,end,start-end,*length));
  }

  return;
}



static void
count_matches_pro (int *matches, int *mismatches, int *unknowns, 
		   struct T *pairs, int npairs) {
  struct T *this = NULL;
  int i;

  i = 0;
  while (i < npairs) {
    /* prev = this; */
    this = &(pairs[i++]);

    if (this->gapp == false) {
      if (this->aa_g != ' ' && this->aa_e != ' ') {
	if (this->aa_g == this->aa_e) {
	  *matches += 1;
	} else if (this->aa_e == 'X') {
	  *unknowns += 1;
	} else {
	  *mismatches += 1;
	}
      }
    }
  }	
  return;
}



void
Pair_print_pslformat_nt (FILE *fp, struct T *pairs, int npairs, T start, T end,
			 Sequence_T queryseq, Chrnum_T chrnum,
			 Univ_IIT_T chromosome_iit, Sequence_T usersegment,
			 int matches, int unknowns, int mismatches, 
			 bool watsonp) {
  Chrpos_T chrpos1, chrpos2;
  struct T *pairs_directional = NULL;
  Intlist_T blockSizes = NULL, qStarts = NULL, p;
  Uintlist_T tStarts = NULL, q;
  int nblocks;
  int qnbreaks, qlength, tnbreaks, tlength, querylength;
  char *chr;

#ifdef PMAP
    querylength = 3*Sequence_fulllength(queryseq);
#else
    querylength = Sequence_fulllength(queryseq);
#endif

  if (watsonp == true) {
    pairs_directional = pairs;
  } else {
    pairs_directional = invert_and_revcomp_path(pairs,npairs);
  }

  nblocks = count_psl_blocks_nt(&blockSizes,&qStarts,&tStarts,pairs_directional,npairs,
				querylength,watsonp);
  compute_gap_lengths_int(&qnbreaks,&qlength,blockSizes,qStarts,nblocks);
  compute_gap_lengths_uint(&tnbreaks,&tlength,blockSizes,tStarts,nblocks);

  fprintf(fp,"%d\t%d\t%d\t%d\t",matches,mismatches,/*repeatmatches*/0,unknowns);
  fprintf(fp,"%d\t%d\t%d\t%d\t",qnbreaks,qlength,tnbreaks,tlength);

  if (watsonp == true) {
    fprintf(fp,"+");
  } else {
    fprintf(fp,"-");
  }
  fprintf(fp,"\t%s\t%d",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));

  fprintf(fp,"\t%d\t%d",start->querypos,end->querypos+1);

  /* T name and T size */
  if (chrnum == 0) {
    fprintf(fp,"\t%s\t%u",Sequence_accession(usersegment),Sequence_fulllength(usersegment));
  } else {
    chr = Chrnum_to_string(chrnum,chromosome_iit);
    fprintf(fp,"\t%s\t%u",chr,Chrnum_length(chrnum,chromosome_iit));
    FREE(chr);
  }

  /* T start and T end */
  chrpos1 = start->genomepos;
  chrpos2 = end->genomepos;
  if (watsonp) {
    fprintf(fp,"\t%u\t%u",chrpos1,chrpos2+1U);
  } else {
    fprintf(fp,"\t%u\t%u",chrpos2,chrpos1+1U);
  }

  fprintf(fp,"\t%d",nblocks);

  fprintf(fp,"\t");
  for (p = blockSizes; p != NULL; p = Intlist_next(p)) {
    fprintf(fp,"%d,",Intlist_head(p));
  }

  fprintf(fp,"\t");
  for (p = qStarts; p != NULL; p = Intlist_next(p)) {
    fprintf(fp,"%d,",Intlist_head(p));
  }

  fprintf(fp,"\t");
  for (q = tStarts; q != NULL; q = Uintlist_next(q)) {
    fprintf(fp,"%u,",Uintlist_head(q));
  }

  Intlist_free(&blockSizes);
  Intlist_free(&qStarts);
  Uintlist_free(&tStarts);

  if (watsonp == false) {
    FREE(pairs_directional);
  }

  putc('\n',fp);
  return;
}

void
Pair_print_pslformat_pro (FILE *fp, struct T *pairs, int npairs, T start, T end,
			  Sequence_T queryseq, Chrnum_T chrnum,
			  Univ_IIT_T chromosome_iit, Sequence_T usersegment,
			  bool watsonp, int cdna_direction) {
  Chrpos_T chrpos1, chrpos2;
  Chrpos_T chrlength;
  Intlist_T blockSizes = NULL, qStarts = NULL, p;
  Uintlist_T tStarts = NULL, q;
  int nblocks, matches = 0, mismatches = 0, unknowns = 0;
  int qnbreaks, qlength, tnbreaks, tlength;
  char *chr;

  chrlength = Chrnum_length(chrnum,chromosome_iit);
  nblocks = count_psl_blocks_pro(&blockSizes,&qStarts,&tStarts,pairs,npairs,
				 watsonp,chrlength);
  compute_gap_lengths_int(&qnbreaks,&qlength,blockSizes,qStarts,nblocks);
  compute_gap_lengths_uint(&tnbreaks,&tlength,blockSizes,tStarts,nblocks);

  count_matches_pro(&matches,&mismatches,&unknowns,pairs,npairs);

  fprintf(fp,"%d\t%d\t%d\t%d\t",matches,mismatches,/*repeatmatches*/0,unknowns);
  fprintf(fp,"%d\t%d\t%d\t%d\t",qnbreaks,qlength,tnbreaks,tlength);

  if (cdna_direction >= 0) {
    fprintf(fp,"+");
  } else {
    fprintf(fp,"-");
  }

  if (watsonp == true) {
    fprintf(fp,"+");
  } else {
    fprintf(fp,"-");
  }
  fprintf(fp,"\t%s\t%d",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));

  fprintf(fp,"\t%d\t%d",(start->querypos+2)/3,end->querypos/3+1);

  /* T name and T size */
  if (chrnum == 0) {
    fprintf(fp,"\t%s\t%u",Sequence_accession(usersegment),Sequence_fulllength(usersegment));
  } else {
    chr = Chrnum_to_string(chrnum,chromosome_iit);
    fprintf(fp,"\tchr%s\t%u",chr,Chrnum_length(chrnum,chromosome_iit));
    FREE(chr);
  }

  /* T start and T end */
  chrpos1 = start->genomepos;
  chrpos2 = end->genomepos;
  if (watsonp) {
    fprintf(fp,"\t%u\t%u",chrpos1,chrpos2+1U);
  } else {
    fprintf(fp,"\t%u\t%u",chrpos2,chrpos1+1U);
  }

  nblocks = count_psl_blocks_pro(&blockSizes,&qStarts,&tStarts,pairs,npairs,
				 watsonp,chrlength);
  fprintf(fp,"\t%d",nblocks);
  fprintf(fp,"\t");

  for (p = blockSizes; p != NULL; p = Intlist_next(p)) {
    fprintf(fp,"%d,",Intlist_head(p));
  }

  fprintf(fp,"\t");
  for (p = qStarts; p != NULL; p = Intlist_next(p)) {
    fprintf(fp,"%d,",Intlist_head(p));
  }

  fprintf(fp,"\t");

  for (q = tStarts; q != NULL; q = Uintlist_next(q)) {
    fprintf(fp,"%u,",Uintlist_head(q));
  }

  Intlist_free(&blockSizes);
  Intlist_free(&qStarts);
  Uintlist_free(&tStarts);

  putc('\n',fp);
  return;
}

void
Pair_print_exons (FILE *fp, struct T *pairs, int npairs, int wraplength, int ngap, bool cdnap) {
  bool in_exon = false;
  struct T *ptr, *this = NULL;
  int i, exonno = 0, column = 0;

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	if (column != 0) {
	  putc('\n',fp);
	  column = 0;
	}
	fprintf(fp,"</exon>\n");
	in_exon = false;
	if (ngap > 0) {
	  fprintf(fp,"<intron %d>\n",exonno);
	  putc(this->genome,fp);
	  column = 1;
	}
      } else {
	if (ngap > 0) {
	  putc(this->genome,fp);
	  if (++column % wraplength == 0) {
	    putc('\n',fp);
	    column = 0;
	  }
	}
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP,
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	if (ngap > 0) {
	  if (exonno > 0) {
	    if (column != 0) {
	      putc('\n',fp);
	      column = 0;
	    }
	    fprintf(fp,"</intron>\n");
	  }
	}
	fprintf(fp,"<exon %d",++exonno);
	if (cdnap == true) {
	  if (this->aaphase_e >= 0) {
	    fprintf(fp,", phase %d",this->aaphase_e);
	  }
	} else {
	  if (this->aaphase_g >= 0) {
	    fprintf(fp,", phase %d",this->aaphase_g);
	  }
	}
	fprintf(fp,">\n");
	in_exon = true;
      }
      if (cdnap == true) {
	if (this->cdna != ' ') {
	  putc(this->cdna,fp);
	  if (++column % wraplength == 0) {
	    putc('\n',fp);
	    column = 0;
	  }
	}
      } else {
	if (this->genome != ' ') {
	  putc(this->genome,fp);
	  if (++column % wraplength == 0) {
	    putc('\n',fp);
	    column = 0;
	  }
	}
      }
    }
  }
  if (column != 0) {
    putc('\n',fp);
  }
  fprintf(fp,"</exon>\n");

  return;
}


int
Pair_nmatches_posttrim (int *max_match_length, List_T pairs, int pos5, int pos3) {
  int nmatches = 0, match_length;
  bool in_intron = false, indelp = false;
  List_T p;
  T this;

  *max_match_length = match_length = 0;
  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      if (!in_intron) {
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	indelp = true;
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->querypos < pos5) {
	/* Don't count match or mismatch */
      } else if (this->querypos >= pos3) {
	/* Don't count match or mismatch */
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	nmatches++;
	match_length++;
      } else if (this->comp == MISMATCH_COMP) {
	/* (*mismatches)++; */
	if (match_length > *max_match_length) {
	  *max_match_length = match_length;
	}
	match_length = 0;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
  }

  if (match_length > *max_match_length) {
    *max_match_length = match_length;
  }

  return nmatches;
}


int
Pair_array_nmatches_posttrim (struct T *pairarray, int npairs, int pos5, int pos3) {
  int nmatches = 0;
  bool in_intron = false, indelp = false;
  int i;
  T this;

  for (i = 0; i < npairs; i++) {
    this = &(pairarray[i]);
    if (this->gapp) {
      if (!in_intron) {
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	indelp = true;
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->querypos < pos5) {
	/* Don't count match or mismatch */
      } else if (this->querypos >= pos3) {
	/* Don't count match or mismatch */
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	nmatches++;
      } else if (this->comp == MISMATCH_COMP) {
	/* (*mismatches)++; */
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
  }

  return nmatches;
}


int
Pair_nmismatches_region (int *nindelbreaks, struct T *pairs, int npairs,
			 int trim_left, int trim_right, int start_amb_nmatches, int end_amb_nmatches,
			 int querylength) {
  int nmismatches = 0;
  bool in_intron = false, indelp = false;
  int i = 0;
  T this;

  *nindelbreaks = 0;

  /* Handle GMAP alignments that are not extended to the end */
  this = &(pairs[0]);
  if (this->querypos - start_amb_nmatches < trim_left) {
    /* Skip */
  } else {
    nmismatches += (this->querypos - start_amb_nmatches) - trim_left;
  }

  while (i < npairs) {
    this = &(pairs[i]);
    if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
      /* Count indelbreaks, even if outside of trimmed region */
      if (this->genome == ' ') {
	/* INSERTION */
	while (i < npairs && this->genome == ' ') {
	  /* (*total_nindels) += 1; */
	  this = &(pairs[i++]);
	}
	i--;
	(*nindelbreaks) += 1;
	
      } else if (this->cdna == ' ') {
	/* DELETION */
	while (i < npairs && this->cdna == ' ') {
	  /* (*total_nindels) -= 1; */
	  this = &(pairs[i++]);
	}
	i--;
	(*nindelbreaks) += 1;
      }

    } else if (this->querypos < trim_left) {
      /* Skip for counting mismatches */
    } else if (this->querypos >= querylength - trim_right) {
      /* Skip for counting mismatches */
    } else if (this->comp == MISMATCH_COMP) {
      nmismatches++;
    }
    i++;
  }

  /* Handle GMAP alignments that are not extended to the end */
  this = &(pairs[npairs-1]);
  if (this->querypos + end_amb_nmatches >= (querylength - 1) - trim_right) {
    /* Skip */
  } else {
    nmismatches += (querylength - 1 - trim_right) - (this->querypos + end_amb_nmatches);
  }

  return nmismatches;
}
	


void
Pair_fracidentity_simple (int *matches, int *unknowns, int *mismatches, List_T pairs) {
  bool in_intron = false;
  List_T p;
  T this;

  *matches = *unknowns = *mismatches = 0;
  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      if (!in_intron) {
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	(*unknowns)++;
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	(*matches)++;
      } else if (this->comp == MISMATCH_COMP) {
	(*mismatches)++;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
  }

  return;
}


void
Pair_fracidentity (int *matches, int *unknowns, int *mismatches, int *qopens, int *qindels, 
		   int *topens, int *tindels, int *ncanonical, int *nsemicanonical, int *nnoncanonical,
		   double *min_splice_prob, List_T pairs, int cdna_direction) {
  bool in_intron = false;
  List_T p;
  T this, prev = NULL;

  *matches = *unknowns = *mismatches = *qopens = *qindels = *topens = *tindels = 
    *ncanonical = *nsemicanonical = *nnoncanonical = 0;
  *min_splice_prob = 1.0;

  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      if (this->donor_prob < *min_splice_prob) {
	*min_splice_prob = this->donor_prob;
      }
      if (this->acceptor_prob < *min_splice_prob) {
	*min_splice_prob = this->acceptor_prob;
      }
      if (!in_intron) {
	if (cdna_direction > 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	    in_intron = true;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	    in_intron = true;
	  } else if (this->genomejump - this->queryjump < 50) {
	    (*topens)++;
	    (*tindels) += this->genomejump - this->queryjump;
	    /* in_intron = false */
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	    in_intron = true;
	  }

	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	    in_intron = true;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	    in_intron = true;
	  } else if (this->genomejump - this->queryjump < 50) {
	    (*topens)++;
	    (*tindels) += this->genomejump - this->queryjump;
	    /* in_intron = false */
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	    in_intron = true;
	  }
	}
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  (*tindels)++;		/* If genome has extra char, count it as a genome skip */
	  if (prev && prev->cdna != ' ') {
	    (*topens)++;
	  }
	} else if (this->genome == ' ') {
	  (*qindels)++;
	  if (prev && prev->genome != ' ') {
	    (*qopens)++;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	(*unknowns)++;
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	(*matches)++;
      } else if (this->comp == MISMATCH_COMP) {
	(*mismatches)++;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    prev = this;
  }

  return;
}


int
Pair_fracidentity_array (int *matches, int *unknowns, int *mismatches, int *qopens, int *qindels, 
			 int *topens, int *tindels, int *ncanonical, int *nsemicanonical, int *nnoncanonical,
			 double *min_splice_prob, struct T *ptr, int npairs, int cdna_direction) {
  bool in_intron = false;
  int i;
  T this, prev = NULL;

  *matches = *unknowns = *mismatches = *qopens = *qindels = *topens = *tindels = 
    *ncanonical = *nsemicanonical = *nnoncanonical = 0;
  *min_splice_prob = 1.0;

  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->gapp) {
      if (this->donor_prob < *min_splice_prob) {
	*min_splice_prob = this->donor_prob;
      }
      if (this->acceptor_prob < *min_splice_prob) {
	*min_splice_prob = this->acceptor_prob;
      }
      if (!in_intron) {
	if (cdna_direction > 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	    in_intron = true;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	    in_intron = true;
	  } else if (this->genomejump - this->queryjump < 50) {
	    (*topens)++;
	    (*tindels) += this->genomejump - this->queryjump;
	    /* in_intron = false */
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	    in_intron = true;
	  }

	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	    in_intron = true;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	    in_intron = true;
	  } else if (this->genomejump - this->queryjump < 50) {
	    (*topens)++;
	    (*tindels) += this->genomejump - this->queryjump;
	    /* in_intron = false */
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	    in_intron = true;
	  }
	}
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  (*tindels)++;		/* If genome has extra char, count it as a genome skip */
	  if (prev && prev->cdna != ' ') {
	    (*topens)++;
	  }
	} else if (this->genome == ' ') {
	  (*qindels)++;
	  if (prev && prev->genome != ' ') {
	    (*qopens)++;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	(*unknowns)++;
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	(*matches)++;
      } else if (this->comp == MISMATCH_COMP) {
	(*mismatches)++;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    prev = this;
  }

  return (*matches) + MISMATCH*(*mismatches)
    + QOPEN*(*qopens) + QINDEL*(*qindels) + TOPEN*(*topens) + TINDEL*(*tindels)
    - CANONICAL_POINTS*(*nnoncanonical);
}


#if 0
/* Called on first and last exons during distal/medial calculation */
/* Procedure seems to give random results */
int
Pair_fracidentity_changepoint (List_T pairs, int cdna_direction) {
  int changepoint = 0, maxscore = 0, score = 0;
  int i = 0;

  bool in_intron = false;
  List_T p;
  T this, prev = NULL;

  for (p = pairs; p != NULL; p = p->rest) {
    i++;
    this = p->first;
    debug3(fprintf(fp,"%d: ",i));
    debug3(Pair_dump_one(this,/*zerobasedp*/false));
    if (this->gapp) {
      if (!in_intron) {
#if 0
	/* Don't expect an intron */
	if (cdna_direction > 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }

	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }
	}
#endif
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  score += TINDEL;
	  if (prev && prev->cdna != ' ') {
	    score += TOPEN;
	  }
	} else if (this->genome == ' ') {
	  score += QINDEL;
	  if (prev && prev->genome != ' ') {
	    score += QOPEN;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
#if 0
	score += (MATCH + MATCH); /* Give more weight to matches to allow for poor quality at ends */
#else
	score += MATCH;
#endif
	if (score > maxscore) {
	  maxscore = score;
	  changepoint = i;
	  debug3(fprintf(fp," => maxscore %d",maxscore));
	}
      } else if (this->comp == MISMATCH_COMP) {
	score += MISMATCH;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    debug3(fprintf(fp,"\n"));
    prev = this;
  }

  return changepoint;
}
#endif


int
Pair_fracidentity_score (List_T pairs, int cdna_direction) {
  int score = 0;
  int i = 0;

  bool in_intron = false;
  List_T p;
  T this, prev = NULL;

  for (p = pairs; p != NULL; p = p->rest) {
    i++;
    this = p->first;
    debug3(fprintf(fp,"%d: ",i));
    debug3(Pair_dump_one(this,/*zerobasedp*/false));
    if (this->gapp) {
      if (!in_intron) {
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  score += TINDEL;
	  if (prev && prev->cdna != ' ') {
	    score += TOPEN;
	  }
	} else if (this->genome == ' ') {
	  score += QINDEL;
	  if (prev && prev->genome != ' ') {
	    score += QOPEN;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	score += MATCH;
      } else if (this->comp == MISMATCH_COMP) {
	score += MISMATCH;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    debug3(fprintf(fp,"\n"));
    prev = this;
  }

  return score;
}


double
Pair_frac_error (List_T pairs, int cdna_direction) {
  int matches, unknowns, mismatches, qopens, qindels,
    topens, tindels, ncanonical, nsemicanonical, nnoncanonical;
  int den;
  double min_splice_prob;

  Pair_fracidentity(&matches,&unknowns,&mismatches,&qopens,&qindels, 
		    &topens,&tindels,&ncanonical,&nsemicanonical,&nnoncanonical,
		    &min_splice_prob,pairs,cdna_direction);

  if ((den = matches + mismatches + qindels + tindels) == 0) {
    return 1.0;
  } else {
    return (double) (mismatches + qindels + tindels)/(double) den;
  }
}

void
Pair_fracidentity_bounded (int *matches, int *unknowns, int *mismatches, 
			   int *qopens, int *qindels, int *topens, int *tindels,
			   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
			   struct T *ptr, int npairs, 
			   int cdna_direction, int minpos, int maxpos) {
  bool in_intron = false;
  T this, prev = NULL;
  int i;

  *matches = *unknowns = *mismatches = *qopens = *qindels = *topens = *tindels = 
    *ncanonical = *nsemicanonical = *nnoncanonical = 0;
  
  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->gapp) {
      if (!in_intron) {
	if (this->querypos >= minpos && this->querypos <= maxpos) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->querypos >= minpos && this->querypos <= maxpos) {
	if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	  if (this->cdna == ' ') {
	    (*tindels)++;		/* If genome has extra char, count it as a genome skip */
	    if (prev && prev->cdna != ' ') {
	      (*topens)++;
	    }
	  } else if (this->genome == ' ') {
	    (*qindels)++;
	    if (prev && prev->genome != ' ') {
	      (*qopens)++;
	    }
	  } else {
	    fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		    this->comp,this->cdna,this->genome);
	    abort();
	  }
#ifndef PMAP
	} else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	  (*unknowns)++;
#endif
	} else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	  (*matches)++;
	} else if (this->comp == MISMATCH_COMP) {
	  (*mismatches)++;
	} else {
	  fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	  abort();
	}
      }
    }
    prev = this;
  }
  return;
}

static const Except_T Array_bounds_error = { "Exceeded array bounds" };


void
Pair_matchscores (int *matchscores, struct T *ptr, int npairs, int querylength) {
  T this;
  int querypos;
  int i;

  for (i = 0; i < npairs; i++) {
    this = ptr++;
    querypos = this->querypos;

    if (this->gapp) {
      matchscores[querypos] = 0;	/* Count as mismatch; make evidence support the gap */
    } else if (this->comp == MISMATCH_COMP) {
      matchscores[querypos] = 0; /* For mismatch */
    } else if (this->comp == INDEL_COMP) {
      matchscores[querypos] = -1;	/* Ignore indels */
    } else {
      matchscores[querypos] = 1; /* For match */
    }
  }

  return;
}


#if 0
/* Called only by chop_ends_by_changepoint in stage3.c, which is not used anymore */
int *
Pair_matchscores_list (int *nmatches, int *ntotal, int *length, List_T pairs) {
  int *matchscores;
  T this;
  List_T p;
  int i = 0;

  matchscores = (int *) CALLOC(List_length(pairs),sizeof(int));
  *nmatches = *ntotal = *length = 0;

  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      matchscores[i++] = 0;	/* Count as mismatch; make evidence support the gap */
      (*ntotal) += 1;
    } else if (this->comp == MISMATCH_COMP) {
      matchscores[i++] = 0; /* For mismatch */
      (*ntotal) += 1;
#ifndef PMAP
    } else if (this->comp == AMBIGUOUS_COMP) {
      matchscores[i++] = 0; /* For cases involving 'N' */
      (*ntotal) += 1;
#endif
    } else if (this->comp == INDEL_COMP) {
      matchscores[i++] = -1;	/* Ignore indels */
    } else {
      matchscores[i++] = 1; /* For match */
      (*nmatches) += 1;
      (*ntotal) += 1;
    }
    (*length) += 1;
  }

  return matchscores;
}
#endif


void
Pair_pathscores (bool *gapp, int *pathscores, struct T *ptr, int npairs, 
		 int cdna_direction, int querylength, cDNAEnd_T cdnaend, int pre_extension_slop) {
  int querypos, querystart, queryend;
  int basescore;
  bool in_intron = false;
  T this, prev = NULL;
  int i;

  /* Determine these before ptr changes */
  this = &(ptr[0]);
  querystart = this->querypos;
  this = &(ptr[npairs-1]);
  queryend = this->querypos;
  /* printf("Entered Pair_pathscores with querystart %d and queryend %d\n",querystart,queryend); */

  /* Allow transitions slightly outside of the ends
     (pre_extension_slop) when finding non-extended paths to pair, but
     not when finding the breakpoint for the final pair, which has
     been extended */
  if (cdnaend == FIVE) {
    /* left part of chimera */
    for (querypos = 0; querypos < querystart; querypos++) {
      gapp[querypos] = true;
    }
    for (querypos = queryend + pre_extension_slop; querypos < querylength; querypos++) {
      gapp[querypos] = true;
    }
  } else {
    /* right part of chimera */
    for (querypos = 0; querypos < querystart - pre_extension_slop; querypos++) {
      gapp[querypos] = true;
    }
    for (querypos = queryend; querypos < querylength; querypos++) {
      gapp[querypos] = true;
    }
  }

  /* Initialize to cover the ends that aren't aligned */
  for (querypos = 0; querypos < querylength; querypos++) {
    pathscores[querypos] = QINDEL;
  }

  for (i = 0; i < npairs; i++) {
    this = ptr++;

    querypos = this->querypos;
    if (querypos >= querylength) {
      fprintf(stderr,"Pair_pathscores: querypos %d >= querylength %d\n",querypos,querylength);
      Pair_dump_array(ptr,npairs,/*zerobasedp*/true);
      fflush(stdout);
      abort();
      RAISE(Array_bounds_error);
    }

    if (this->gapp) {
      gapp[querypos] = true;
      if (in_intron == false) {
	/* Adds only a single reward/penalty per intron */
	if (cdna_direction > 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    pathscores[querypos] = CANONICAL_POINTS;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    pathscores[querypos] = SEMICANONICAL_POINTS;
	  } else {
	    pathscores[querypos] = NONCANONICAL_POINTS; /* noncanonical */
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    pathscores[querypos] = CANONICAL_POINTS;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    pathscores[querypos] = SEMICANONICAL_POINTS;
	  } else {
	    pathscores[querypos] = NONCANONICAL_POINTS; /* noncanonical */
	  }
	} else {
	  pathscores[querypos] = NONCANONICAL_POINTS; /* indeterminate */
	}
	in_intron = true;
      }

    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  pathscores[querypos] = TINDEL;
	  if (prev && prev->cdna != ' ') {
	    pathscores[querypos] = TOPEN;
	  }
	} else if (this->genome == ' ') {
	  pathscores[querypos] = QINDEL;
	  if (prev && prev->genome != ' ') {
	    pathscores[querypos] = QOPEN;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	pathscores[querypos] = +1; /* For match */
      } else if (this->comp == MISMATCH_COMP) {
	pathscores[querypos] = MISMATCH;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    prev = this;
  }

#if 0
  /* Gets querystart to queryend inclusive */
  if (0 && querystart == 0) {
    for (i = 1; i <= queryend; i++) {
      pathscores[i] += pathscores[i-1];
    }
  } else {
    for (i = querystart; i <= queryend; i++) {
      pathscores[i] += pathscores[i-1];
    }
  }
#endif

#if 0
  if (cdnaend == FIVE) {
    for (i = queryend + 1; i < querylength; i++) {
      pathscores[i] = pathscores[i-1] + QINDEL;
    }
  } else if (cdnaend == THREE) {
    for (i = querystart - 1; i >= 0; --i) {
      pathscores[i] = pathscores[i+1] - QINDEL;
    }
    for (i = queryend + 1; i < querylength; i++) {
      pathscores[i] = pathscores[i-1];
    }
  }
#endif

  if (cdnaend == FIVE) {
    for (i = 1; i < querylength; i++) {
      pathscores[i] += pathscores[i-1];
    }
    basescore = pathscores[querystart];
  } else if (cdnaend == THREE) {
    for (i = querylength-2; i >= 0; --i) {
      pathscores[i] += pathscores[i+1];
    }
    basescore = pathscores[queryend];
  }

  for (i = 0; i < querylength; i++) {
    pathscores[i] -= basescore;
  }

  return;
}


int
Pair_nexons_approx (List_T pairs) {
  int nexons = 0;
  bool in_exon = false;
  T this;
  List_T p;
  
  for (p = pairs; p != NULL; p = List_next(p)) {
    this = List_head(p);
    if (this->gapp) {
      if (in_exon) {
	in_exon = false;
      }
    } else {
      if (!in_exon) {
	nexons++;
	in_exon = true;
      }
    }
  }

  return nexons;
}


int
Pair_nexons (struct T *pairs, int npairs) {
  int nexons = 0;
  struct T *ptr, *this = NULL;
  bool in_exon = false;
  int i;
  
  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->gapp) {
      if (in_exon) {
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      if (!in_exon) {
	nexons++;
	in_exon = true;
      }
    }
  }

  return nexons;
}


bool
Pair_consistentp (int *ncanonical, struct T *pairs, int npairs, int cdna_direction) {
  bool in_intron = false;
  struct T *this;
  int i;

  *ncanonical = 0;
  for (i = 0; i < npairs; i++) {
    this = pairs++;
    if (this->gapp) {
      if (!in_intron) {
	if (cdna_direction > 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP || 
	      this->comp == REV_GCAG_INTRON_COMP || 
	      this->comp == REV_ATAC_INTRON_COMP) {
	    return false;
	  } else if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP || 
	      this->comp == FWD_GCAG_INTRON_COMP || 
	      this->comp == FWD_ATAC_INTRON_COMP) {
	    return false;
	  } else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  }
	} else if (cdna_direction == 0) {
	  /* Set cdna_direction for next time */
	  if (this->comp == FWD_CANONICAL_INTRON_COMP || 
	      this->comp == FWD_GCAG_INTRON_COMP || 
	      this->comp == FWD_ATAC_INTRON_COMP) {
	    cdna_direction = +1;
	  } else if (this->comp == REV_CANONICAL_INTRON_COMP || 
		     this->comp == REV_GCAG_INTRON_COMP || 
		     this->comp == REV_ATAC_INTRON_COMP) {
	    cdna_direction = -1;
	  } 
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
    }
  }

  return true;
}


static void
invert_intron (char *donor, char *acceptor) {
  char temp;

  temp = donor[0];
  donor[0] = complCode[(int) acceptor[1]];
  acceptor[1] = complCode[(int) temp];

  temp = donor[1];
  donor[1] = complCode[(int) acceptor[0]];
  acceptor[0] = complCode[(int) temp];
  
  return;
}


void
Pair_print_protein_genomic (FILE *fp, struct T *ptr, int npairs, int wraplength, bool forwardp) {
  struct T *this;
  int xpos = 0, i;

  if (forwardp == true) {
    for (i = 0; i < npairs; i++) {
      this = ptr++;
      if (this->aa_g != ' ') {
	if (xpos == wraplength) {
	  putc('\n',fp);
	  xpos = 0;
	}
#ifdef PMAP
	putc(this->aa_g,fp);
	xpos++;
#else
	if (this->aa_g != '*') {
	  putc(this->aa_g,fp);
	  xpos++;
	}
#endif
      }
    }
    putc('\n',fp);

  } else {
    for (i = npairs-1; i >= 0; i--) {
      this = ptr--;
      if (this->aa_g != ' ') {
	if (xpos == wraplength) {
	  putc('\n',fp);
	  xpos = 0;
	}
#ifdef PMAP
	abort();
	putc(this->aa_g,fp);
	xpos++;
#else
	if (this->aa_g != '*') {
	  putc(this->aa_g,fp);
	  xpos++;
	}
#endif
      }
    }
    putc('\n',fp);

  }

  return;
}

#ifdef PMAP
void
Pair_print_nucleotide_cdna (FILE *fp, struct T *ptr, int npairs, int wraplength) {
  struct T *this;
  int xpos = 0, i;

  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->cdna != ' ') {
      if (xpos == wraplength) {
	putc('\n',fp);
	xpos = 0;
      }
      putc(this->cdna,fp);
      xpos++;
    }
  }
  putc('\n',fp);
  return;
}
#else
void
Pair_print_protein_cdna (FILE *fp, struct T *ptr, int npairs, int wraplength, bool forwardp) {
  struct T *this;
  int xpos = 0, i;

  if (forwardp == true) {
    for (i = 0; i < npairs; i++) {
      this = ptr++;
      if (this->aa_e != ' ') {
	if (xpos == wraplength) {
	  putc('\n',fp);
	  xpos = 0;
	}
	if (this->aa_e != '*') {
	  putc(this->aa_e,fp);
	  xpos++;
	}
      }
    }
    putc('\n',fp);

  } else {
    for (i = npairs-1; i >= 0; i--) {
      this = ptr--;
      if (this->aa_e != ' ') {
	if (xpos == wraplength) {
	  putc('\n',fp);
	  xpos = 0;
	}
	if (this->aa_e != '*') {
	  putc(this->aa_e,fp);
	  xpos++;
	}
      }
    }
    putc('\n',fp);
  }

  return;
}
#endif


void
Pair_print_compressed (FILE *fp, int pathnum, int npaths, T start, T end, Sequence_T queryseq, char *dbversion,
		       Sequence_T usersegment, int nexons, double fracidentity,
		       struct T *pairs, int npairs, Chrnum_T chrnum,
		       Univcoord_T chroffset, Univ_IIT_T chromosome_iit, int querylength_given,
		       int skiplength, int trim_start, int trim_end, bool checksump,
		       int chimerapos, int chimeraequivpos, double donor_prob, double acceptor_prob,
		       int chimera_cdna_direction, char *strain, bool watsonp, int cdna_direction) {
  Chrpos_T chrpos1, chrpos2;
  Univcoord_T position1, position2;

  bool in_exon = false;
  List_T tokens = NULL;
  struct T *ptr = pairs, *this = NULL;
  int querypos1, querypos2;
  int exon_querystart = -1, exon_queryend;
  Chrpos_T exon_genomestart = -1, exon_genomeend, intron_start, intron_end;
  int num = 0, den = 0, runlength = 0, i;
  int print_dinucleotide_p;
  char token[10], donor[3], acceptor[3], *chr;
  double coverage;
  /* double trimmed_coverage; */
  int last_querypos = -1;
  Chrpos_T last_genomepos = -1U;

  donor[0] = donor[1] = donor[2] = '\0';
  acceptor[0] = acceptor[1] = acceptor[2] = '\0';

  querypos1 = start->querypos;
  querypos2 = end->querypos;

  fprintf(fp,">%s ",Sequence_accession(queryseq));
  if (dbversion != NULL) {
    fprintf(fp,"%s ",dbversion);
  } else if (usersegment != NULL && Sequence_accession(usersegment) != NULL) {
    fprintf(fp,"%s ",Sequence_accession(usersegment));
  } else {
    fprintf(fp,"user-provided ");
  }
#ifdef PMAP
  fprintf(fp,"%d/%d %d %d",pathnum,npaths,(querylength_given+skiplength)*3,nexons);
  coverage = (double) (querypos2 - querypos1 + 1)/(double) ((querylength_given+skiplength)*3);
  fprintf(fp," %.1f",((double) rint(1000.0*coverage)));
#else
  coverage = (double) (querypos2 - querypos1 + 1)/(double) (querylength_given+skiplength);
  if (end->querypos + 1 > trim_end) {
    trim_end = end->querypos + 1;
  }
  if (start->querypos < trim_start) {
    trim_start = start->querypos;
  }
  /*
  trimmed_coverage = (double) (end->querypos - start->querypos + 1)/(double) (trim_end - trim_start + skiplength);
  fprintf(fp,">%s %s %d/%d %d(%d) %d",
	 Sequence_accession(queryseq),dbversion,pathnum,npaths,
	 querylength_given+skiplength,trim_end-trim_start,nexons);
  fprintf(fp," %.1f(%.1f)",((double) rint(1000.0*coverage))/10.0,((double) rint(1000.0*trimmed_coverage))/10.0);
  */
  fprintf(fp,"%d/%d %d %d",pathnum,npaths,querylength_given+skiplength,nexons);
  fprintf(fp," %.1f",((double) rint(1000.0*coverage))/10.0);
#endif
  fprintf(fp," %.1f",((double) rint(1000.0*fracidentity))/10.0);

  start = &(pairs[0]);
  end = &(pairs[npairs-1]);
  fprintf(fp," %d%s%d",start->querypos + ONEBASEDP,"..",end->querypos + ONEBASEDP);

  chrpos1 = start->genomepos;
  chrpos2 = end->genomepos;
  position1 = chroffset + chrpos1;
  position2 = chroffset + chrpos2;
  fprintf(fp," %u%s%u",position1 + ONEBASEDP,"..",position2 + ONEBASEDP);

  if (chrnum == 0) {
    fprintf(fp," %u%s%u",chrpos1 + ONEBASEDP,"..",chrpos2 + ONEBASEDP);
  } else {
    chr = Chrnum_to_string(chrnum,chromosome_iit);
    fprintf(fp," %s:%u%s%u",chr,chrpos1 + ONEBASEDP,"..",chrpos2 + ONEBASEDP);
    FREE(chr);
  }

  if (chrpos1 <= chrpos2) {
    fprintf(fp," +");
  } else {
    fprintf(fp," -");
  }

  if (cdna_direction > 0) {
    fprintf(fp," dir:sense");
  } else if (cdna_direction < 0) {
    fprintf(fp," dir:antisense");
  } else {
    fprintf(fp," dir:indet");
  }

  if (checksump == true) {
    fprintf(fp," md5:");
    Sequence_print_digest(fp,queryseq);
  }

  if (chimerapos >= 0) {
    if (chimeraequivpos == chimerapos) {
      if (donor_prob > 0.0 && acceptor_prob > 0.0) {
	if (chimera_cdna_direction >= 0) {
	  fprintf(fp," chimera:%d(>)/%.3f/%.3f",chimerapos + ONEBASEDP,donor_prob,acceptor_prob);
	} else {
	  fprintf(fp," chimera:%d(<)/%.3f/%.3f",chimerapos + ONEBASEDP,donor_prob,acceptor_prob);
	}
      } else {
	fprintf(fp," chimera:%d",chimerapos + ONEBASEDP);
      }
    } else {
      fprintf(fp," chimera:%d..%d",chimerapos + ONEBASEDP,chimeraequivpos + ONEBASEDP);
    }
  }

  if (strain != NULL) {
    fprintf(fp," strain:%s",strain);
  }

  putc('\n',fp);

  for (i = 0; i < npairs; i++) {
    /* prev = this; */
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_queryend = last_querypos + ONEBASEDP;
	exon_genomeend = last_genomepos + ONEBASEDP;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}

	fprintf(fp,"\t%u %u",exon_genomestart,exon_genomeend);
	fprintf(fp," %d %d",exon_querystart,exon_queryend);
	if (den == 0) {
	  fprintf(fp," 100");
	} else {
	  fprintf(fp," %d",(int) floor(100.0*(double) num/(double) den));
	}
	print_dinucleotide_p = 1;
	if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	  sprintf(token,"%d>",runlength);
	} else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	  sprintf(token,"%d<",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == NONINTRON_COMP) {
	  sprintf(token,"%d=",runlength);
	} else if (this->comp == FWD_GCAG_INTRON_COMP) {
	  sprintf(token,"%d)",runlength);
	} else if (this->comp == REV_GCAG_INTRON_COMP) {
	  sprintf(token,"%d(",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == FWD_ATAC_INTRON_COMP) {
	  sprintf(token,"%d]",runlength);
	} else if (this->comp == REV_ATAC_INTRON_COMP) {
	  sprintf(token,"%d[",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == DUALBREAK_COMP) {
	  sprintf(token,"%d#",runlength);
	  print_dinucleotide_p = 0;
	} else if (this->comp == EXTRAEXON_COMP) {
	  sprintf(token,"%d#",runlength);
	  print_dinucleotide_p = 0;
	} else {
	  fprintf(stderr,"Can't parse comp '%c' in compression for %s\n",
		  this->comp,Sequence_accession(queryseq));
	  abort();
	}
	tokens = push_token(tokens,token);
	tokens = List_reverse(tokens);
	print_tokens_compressed(fp,tokens);
	List_free(&tokens);
	fprintf(fp,"\t%d",exon_queryend - exon_querystart + 1);

	runlength = 0;
	donor[0] = this->genome;
	donor[1] = '\0';
	in_exon = false;
      } else if (donor[1] == '\0') {
	donor[1] = this->genome;
      } else {
	acceptor[0] = acceptor[1];
	acceptor[1] = this->genome;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + ONEBASEDP;
	exon_genomestart = this->genomepos + ONEBASEDP;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}
	if (i > 0) {
	  if (intron_end > intron_start) {
	    fprintf(fp,"\t%d",intron_end - intron_start + 1);
	  } else {
	    fprintf(fp,"\t%d",intron_start - intron_end + 1);
	  }
	  if (print_dinucleotide_p == -1) {
	    invert_intron(donor,acceptor);
	  }
	  if (print_dinucleotide_p != 0) {
	    if ((donor[0] == 'G' || donor[0] == 'g') &&
		(donor[1] == 'T' || donor[1] == 't') &&
		(acceptor[0] == 'A' || acceptor[0] == 'a') &&
		(acceptor[1] == 'G' || acceptor[1] == 'g')) {
	      /* Do nothing */
	    } else {
	      fprintf(fp,"\t%c%c-%c%c",toupper(donor[0]),toupper(donor[1]),toupper(acceptor[0]),toupper(acceptor[1]));
	    }
	  }
#if 0
	  if (exon_querystart > exon_queryend + 1) {
	    fprintf(fp,"***");
	  }
#endif
	  putc('\n',fp);
	}

	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (this->genome == ' ') {
	  sprintf(token,"%d^%c",runlength,this->cdna);
	} else if (this->cdna == ' ') {
	  sprintf(token,"%dv",runlength);
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}
	tokens = push_token(tokens,token);
	runlength = 0;
	/* Don't increment den */

      } else if (this->comp == MISMATCH_COMP) {
	sprintf(token,"%dx%c",runlength,this->cdna);
	tokens = push_token(tokens,token);
	runlength = 0;
	den++;

#ifndef PMAP
      } else if (this->comp == AMBIGUOUS_COMP) {
	sprintf(token,"%d:%c",runlength,this->cdna);
	tokens = push_token(tokens,token);
	runlength = 0;
	den++;
	num++;
#endif

      } else {
	runlength++;
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	  /* AMBIGUOUS_COMP handled above */
	  num++;
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }
  
  /* prev = this; */
  exon_queryend = last_querypos + ONEBASEDP;
  exon_genomeend = last_genomepos + ONEBASEDP;
  
  fprintf(fp,"\t%d %d",exon_genomestart,exon_genomeend);
  fprintf(fp," %d %d",exon_querystart,exon_queryend);
  if (den == 0) {
    fprintf(fp," 100");
  } else {
    fprintf(fp," %d",(int) floor(100.0*(double) num/(double) den));
  }

  sprintf(token,"%d*",runlength);
  tokens = push_token(tokens,token);
  tokens = List_reverse(tokens);
  print_tokens_compressed(fp,tokens);
  List_free(&tokens);

  fprintf(fp,"\t%d",exon_queryend - exon_querystart + 1);
  putc('\n',fp);

  return;
}


void
Pair_print_iit_map (FILE *fp, Sequence_T queryseq, char *accession,
		    T start, T end, Chrnum_T chrnum, Univ_IIT_T chromosome_iit) {
  char *chrstring = NULL;
  Chrpos_T chrpos1, chrpos2;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  /* Made identical to code for Pair_print_iit_exon_map */
  chrpos1 = start->genomepos + ONEBASEDP;
  chrpos2 = end->genomepos + ONEBASEDP;
  fprintf(fp,">%s %s:%u..%u\n",accession,chrstring,chrpos1,chrpos2);
  Sequence_print_header(fp,queryseq,/*checksump*/false);

  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Pair_print_iit_exon_map (FILE *fp, struct T *pairs, int npairs, Sequence_T queryseq, char *accession,
			 T start, T end, Chrnum_T chrnum, Univ_IIT_T chromosome_iit) {
  int i;
  bool in_exon = false;
  struct T *ptr = pairs, *this = NULL;
  Chrpos_T exon_genomestart = -1, exon_genomeend;
  char *chrstring = NULL;
  Chrpos_T chrpos1, chrpos2;
  Chrpos_T last_genomepos = -1U;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  chrpos1 = start->genomepos + ONEBASEDP;
  chrpos2 = end->genomepos + ONEBASEDP;
  fprintf(fp,">%s %s:%u..%u\n",accession,chrstring,chrpos1,chrpos2);
  Sequence_print_header(fp,queryseq,/*checksump*/false);

  for (i = 0; i < npairs; i++) {
    /* prev = this; */
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_genomeend = last_genomepos + ONEBASEDP;
	fprintf(fp,"%u %u\n",exon_genomestart,exon_genomeend);
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_genomestart = this->genomepos + ONEBASEDP;
	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }
  
  /* prev = this; */
  exon_genomeend = last_genomepos + ONEBASEDP;
  
  fprintf(fp,"%u %u\n",exon_genomestart,exon_genomeend);

  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Pair_print_splicesites (FILE *fp, struct T *pairs, int npairs, char *accession,
			int nexons, Chrnum_T chrnum, Univ_IIT_T chromosome_iit, bool watsonp) {
  int exoni = 0, i;
  bool in_exon = false;
  struct T *ptr = pairs, *this = NULL;
  Chrpos_T exon_genomestart = -1U, exon_genomeend;
  char *chrstring = NULL;
  Chrpos_T last_genomepos = -1U, intron_length;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  for (i = 0; i < npairs; i++) {
    /* prev = this; */
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_genomeend = last_genomepos + ONEBASEDP;
	if (watsonp) {
	  fprintf(fp,">%s.exon%d/%d %s:%u..%u donor",accession,exoni,nexons,chrstring,exon_genomeend,exon_genomeend+1U);
	} else {
	  fprintf(fp,">%s.exon%d/%d %s:%u..%u donor",accession,exoni,nexons,chrstring,exon_genomeend,exon_genomeend-1U);
	}
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exoni++;
	if (exoni > 1) {
	  exon_genomestart = this->genomepos + ONEBASEDP;
	  if (watsonp) {
	    intron_length = exon_genomestart - exon_genomeend - 1U;
	    fprintf(fp," %u\n",intron_length); /* For previous donor */
	    fprintf(fp,">%s.exon%d/%d %s:%u..%u acceptor",accession,exoni,nexons,chrstring,exon_genomestart-1U,exon_genomestart);
	    fprintf(fp," %u\n",intron_length);
	  } else {
	    intron_length = exon_genomeend - exon_genomestart - 1U;
	    fprintf(fp," %u\n",intron_length); /* For previous donor */
	    fprintf(fp,">%s.exon%d/%d %s:%u..%u acceptor",accession,exoni,nexons,chrstring,exon_genomestart+1U,exon_genomestart);
	    fprintf(fp," %u\n",intron_length);
	  }
	}

	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }
  
  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Pair_print_introns (FILE *fp, struct T *pairs, int npairs, char *accession,
		    int nexons, Chrnum_T chrnum, Univ_IIT_T chromosome_iit) {
  int exoni = 0, i;
  bool in_exon = false;
  struct T *ptr = pairs, *this = NULL;
  Chrpos_T exon_genomestart = -1, exon_genomeend;
  char *chrstring = NULL;
  Chrpos_T last_genomepos = -1U;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  for (i = 0; i < npairs; i++) {
    /* prev = this; */
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_genomeend = last_genomepos + ONEBASEDP;
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exoni++;
	if (exoni > 1) {
	  exon_genomestart = this->genomepos + ONEBASEDP;
	  fprintf(fp,">%s.intron%d/%d %s:%u..%u\n",accession,exoni-1,nexons-1,chrstring,exon_genomeend,exon_genomestart);
	}

	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }
  
  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}



/* goal_start < goal_end */
Chrpos_T
Pair_binary_search_ascending (int *querypos, int lowi, int highi, struct T *pairarray,
			      Chrpos_T goal_start, Chrpos_T goal_end) {
  int middlei;

  debug10(printf("entered binary search_ascending with lowi=%d, highi=%d, goal=%u..%u\n",
		 lowi,highi,goal_start,goal_end));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    while (middlei < highi && pairarray[middlei].cdna == ' ') {
      /* Go forward past pairs corresponding to gaps */
      middlei++;
    }
    if (middlei >= highi) {
      middlei = lowi + ((highi - lowi) / 2);
      while (middlei >= lowi && pairarray[middlei].cdna == ' ') {
	/* Go backward past pairs corresponding to gaps */
	middlei--;
      }
      if (middlei < lowi) {
	debug10(printf("all intermediate pairs are gaps\n"));
#if 0
	*querypos = pairarray[lowi].querypos;
	return pairarray[lowi].genomepos;
#else
	return 0U;
#endif
      }
    }

    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u..%u\n",
		   lowi,pairarray[lowi].genomepos,middlei,pairarray[middlei].genomepos,
		   highi,pairarray[highi].genomepos,goal_start,goal_end));
    if (goal_end < pairarray[middlei].genomepos) {
      highi = middlei;
    } else if (goal_start > pairarray[middlei].genomepos) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      *querypos = pairarray[middlei].querypos;
      return pairarray[middlei].genomepos;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return 0U;
}

/* goal_start > goal_end */
Chrpos_T
Pair_binary_search_descending (int *querypos, int lowi, int highi, struct T *pairarray,
			       Chrpos_T goal_start, Chrpos_T goal_end) {
  int middlei;

  debug10(printf("entered binary search_descending with lowi=%d, highi=%d, goal=%u..%u\n",
		 lowi,highi,goal_start,goal_end));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    while (middlei < highi && pairarray[middlei].cdna == ' ') {
      /* Go forward past pairs corresponding to gaps */
      middlei++;
    }
    if (middlei >= highi) {
      middlei = lowi + ((highi - lowi) / 2);
      while (middlei >= lowi && pairarray[middlei].cdna == ' ') {
	/* Go backward past pairs corresponding to gaps */
	middlei--;
      }
      if (middlei < lowi) {
	debug10(printf("all intermediate pairs are gaps\n"));
#if 0
	*querypos = pairarray[lowi].querypos;
	return pairarray[lowi].genomepos;
#else
	return 0U;
#endif
      }
    }

    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u..%u\n",
		   lowi,pairarray[lowi].genomepos,middlei,pairarray[middlei].genomepos,
		   highi,pairarray[highi].genomepos,goal_start,goal_end));
    if (goal_end > pairarray[middlei].genomepos) {
      highi = middlei;
    } else if (goal_start < pairarray[middlei].genomepos) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      *querypos = pairarray[middlei].querypos;
      return pairarray[middlei].genomepos;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return 0U;
}


#ifndef PMAP

bool
Pairarray_contains_p (struct T *pairarray, int npairs, int querypos) {
  int i;
  T pair;

  for (i = 0; i < npairs; i++) {
    pair = &(pairarray[i]);
    if (pair->querypos > querypos) {
      return false;
    } else if (pair->querypos < querypos) {
      /* continue */
    } else if (pair->gapp == true) {
      return false;
    } else if (pair->cdna == ' ') {
      return false;
    } else if (pair->genome == ' ') {
      return false;
    } else {
      return true;
    }
  }

  return false;
}


Chrpos_T
Pair_genomicpos_low (int hardclip_low, int hardclip_high, struct T *pairarray, int npairs, int querylength,
		     bool watsonp, bool hide_soft_clips_p) {
  struct T *clipped_pairs;
  int clipped_npairs;
  T pair;

#if 0
  if (clipdir >= 0) {
    if (watsonp == true) {
      if (first_read_p == true) {
	hardclip_high = hardclip5;
	hardclip_low = 0;
      } else {
	hardclip_high = 0;
	hardclip_low = hardclip3;
      }
    } else {
      if (first_read_p == true) {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
    }
  } else {
    if (watsonp == true) {
      if (first_read_p == true) {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
    } else {
      if (first_read_p == true) {
	hardclip_high = hardclip5;
	hardclip_low = 0;
      } else {
	hardclip_high = 0;
	hardclip_low = hardclip3;
      }
    }
  }
#endif

  if (watsonp == true) {
    clipped_pairs = hardclip_pairs(&clipped_npairs,hardclip_low,hardclip_high,
				   pairarray,npairs,querylength);
    pair = &(clipped_pairs[0]);
    if (hide_soft_clips_p == true) {
      assert(pair->querypos == 0);
      return pair->genomepos + 1U - pair->querypos;
    } else {
      return pair->genomepos + 1U;
    }
  } else {
    /* Swap hardclip_low and hardclip_high */
    clipped_pairs = hardclip_pairs(&clipped_npairs,hardclip_high,hardclip_low,
				   pairarray,npairs,querylength);
    pair = &(clipped_pairs[clipped_npairs-1]);
    if (hide_soft_clips_p == true) {
      assert(pair->querypos == querylength - 1);
      return pair->genomepos + 1U + (querylength - 1 - pair->querypos);
    } else {
      return pair->genomepos + 1U;
    }
  }
}

#endif


Chrpos_T
Pairarray_genomicbound_from_start (struct T *pairarray, int npairs, int overlap) {
  int i;
  struct T pair;

  i = 0;
  pair = pairarray[i];
  while (i < npairs && overlap > 0) {
    pair = pairarray[i];
    if (pair.cdna != ' ') {
      overlap--;
    }
    i++;
  }

  return pair.genomepos;
}

Chrpos_T
Pairarray_genomicbound_from_end (struct T *pairarray, int npairs, int overlap) {
  int i;
  struct T pair;

  i = npairs-1;
  pair = pairarray[i];
  while (i >= 0 && overlap > 0) {
    pair = pairarray[i];
    if (pair.cdna != ' ') {
      overlap--;
    }
    i--;
  }

  return pair.genomepos;
}


int
Pair_cdna_direction (List_T pairs) {
  int cdna_direction = 0;
  bool in_intron = false;
  T this;
  List_T p;
  
  for (p = pairs; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    if (this->gapp) {
      if (!in_intron) {
	if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	  cdna_direction += 1;
	} else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	  cdna_direction -= 1;
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
    }
  }

  return cdna_direction;
}


/* Returns first pair that exceeds breakpoint */
T
Pair_start_bound (int *cdna_direction, List_T pairs, int breakpoint) {
  T start = NULL, this;
  bool in_intron = false;
  List_T p;

  debug9(printf("Entering Pair_start_bound with breakpoint %d\n",breakpoint));

  *cdna_direction = 0;

  if ((p = pairs) != NULL) {
    start = this = (T) p->first;
  }

  while (p != NULL) {
    this = (T) p->first;
    debug9(Pair_dump_one(this,true));
    debug9(printf("\n"));


    if (this->gapp == true) {
      /* Skip */
    } else if (this->querypos > breakpoint) {
      while (p != NULL) {
	this = (T) List_head(p);

	if (this->gapp) {
	  debug9(printf("For start bound, saw gap with comp %c\n",this->comp));
	  if (!in_intron) {
	    if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	      *cdna_direction += 1;
	    } else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	      *cdna_direction -= 1;
	    }
	    in_intron = true;
	  }
	} else {
	  if (in_intron) {
	    in_intron = false;
	  }
	}

	p = p->rest;
      }

      if (*cdna_direction > 0) {
	*cdna_direction = +1;
      } else if (*cdna_direction < 0) {
	*cdna_direction = -1;
      }
      return start;

    } else {
      start = this;
    }

    p = p->rest;
  }

#if 0
  /* Found no gap beyond start */
  if (*cdna_direction > 0) {
    *cdna_direction = +1;
  } else if (*cdna_direction < 0) {
    *cdna_direction = -1;
  }
#endif

  return start;
}


/* Returns last pair that exceeds breakpoint */
T
Pair_end_bound (int *cdna_direction, List_T pairs, int breakpoint) {
  T end = NULL, this;
  bool in_intron = false;
  List_T p;

  debug9(printf("Entering Pair_end_bound with breakpoint %d\n",breakpoint));

  *cdna_direction = 0;

  if ((p = pairs) != NULL) {
    end = this = (T) p->first;
  }

  while (p != NULL) {
    this = (T) p->first;
    debug9(Pair_dump_one(this,true));
    debug9(printf("\n"));
    if (this->gapp) {
      debug9(printf("For end bound, saw gap with comp %c\n",this->comp));
      if (!in_intron) {
	if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	  *cdna_direction += 1;
	} else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	  *cdna_direction -= 1;
	}
	in_intron = true;
      }

    } else {
      if (in_intron) {
	in_intron = false;
      }

      if (this->querypos > breakpoint) {

	if (*cdna_direction > 0) {
	  *cdna_direction = +1;
	} else if (*cdna_direction < 0) {
	  *cdna_direction = -1;
	}
	return end;

      } else {
	end = this;
      }
    }

    p = p->rest;
  }

  if (*cdna_direction > 0) {
    *cdna_direction = +1;
  } else if (*cdna_direction < 0) {
    *cdna_direction = -1;
  }
  return end;
}



List_T
Pair_trim_ends (bool *trim5p, bool *trim3p, List_T pairs, int ambig_end_length_5, int ambig_end_length_3) {
  List_T trimmed = NULL;
  int trim_right = 0, trim_left = -1; /* Needs to be -1 to avoid trimming when pairs is NULL */
  int bestscore, score;
  int pairi;
  List_T p, pairptr;
  T this;
  int i;
  bool in_indelp;

  debug8(printf("Entered trim_ends\n"));
  if (pairs == NULL) {
    *trim5p = *trim3p = 0;
    return (List_T) NULL;
  }


  /* Find trim_right */
  bestscore = 0;
  score = 0;
  in_indelp = false;
  this = (T) NULL;
  for (p = pairs, pairi = 0; p != NULL; p = p->rest, pairi++) {
    this = p->first;

    if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
      if (in_indelp == false) {
	score += trim_indel_score;
	if (score < 0) {
	  score = 0;
	}
	in_indelp = true;
      }

    } else {
      in_indelp = false;
      if (this->gapp) {
	/* Don't count */
	
      } else if (this->comp == INTRONGAP_COMP) {
	/* Do nothing */

      } else if (this->cdna == 'N' || this->comp == MISMATCH_COMP) {
	score += trim_mismatch_score;
	if (score < 0) {
	  score = 0;
	} else if (score >= bestscore) { /* Want >= and not >, so extend to ends */
	  bestscore = score;
	  trim_right = pairi;
	}

      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	score += TRIM_MATCH_SCORE;
	if (score >= bestscore) { /* Want >= and not >, so extend to ends */
	  bestscore = score;
	  trim_right = pairi;
	}

      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }

    debug8(printf("pairi %d, querypos %d, genomepos %u, comp %c: Trim right score %d, trim_right %d, protectedp %d\n",
		  pairi,this->querypos,this->genomepos,this->comp,score,trim_right,this->protectedp));
  }

  if (this == NULL) {
    fprintf(stderr,"check for trim_right yields this == NULL\n");
    abort();
  } else if (ambig_end_length_3 > 0) {
    debug8(printf("Not disturbing ambiguous end on right\n"));
    trim_right = 0;
  } else if (this->protectedp == true) {
    debug8(printf("Protected against trim_right\n"));
    trim_right = 0;
  } else {
    trim_right = pairi - 1 - trim_right;
    debug8(printf("Final: Trim right pairi %d, score %d, trim_right %d\n",pairi,score,trim_right));
  }
  debug8(printf("\n"));


  /* Find trim_left */
  pairs = List_reverse(pairs);
  bestscore = 0;
  score = 0;
  in_indelp = false;
  this = (T) NULL;
  for (p = pairs, pairi = 0; p != NULL; p = p->rest, pairi++) {
    this = p->first;

    if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
      if (in_indelp == false) {
	score += trim_indel_score;
	if (score < 0) {
	  score = 0;
	}
	in_indelp = true;
      }

    } else {
      in_indelp = false;

      if (this->gapp) {
	/* Don't count */
	
      } else if (this->comp == INTRONGAP_COMP) {
	/* Do nothing */

      } else if (this->cdna == 'N' || this->comp == MISMATCH_COMP) {
	score += trim_mismatch_score;
	if (score < 0) {
	  score = 0;
	} else if (score >= bestscore) { /* Want >= and not >, so extend to ends */
	  bestscore = score;
	  trim_left = pairi;
	}

      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	score += TRIM_MATCH_SCORE;
	if (score >= bestscore) { /* Want >= and not >, so extend to ends */
	  bestscore = score;
	  trim_left = pairi;
	}
	
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }

    debug8(printf("pairi %d, querypos %d, genomepos %u, comp %c: Trim left score %d, trim_left %d, protectedp %d\n",
		  pairi,this->querypos,this->genomepos,this->comp,score,trim_left,this->protectedp));
  }

  if (this == NULL) {
    fprintf(stderr,"check for trim_left yields this == NULL\n");
    abort();
  } else if (ambig_end_length_5 > 0) {
    debug8(printf("Not disturbing ambiguous end on left\n"));
    trim_left = pairi - 1;
  } else if (this->protectedp == true) {
    debug8(printf("Protected against trim_left\n"));
    trim_left = pairi - 1;
  } else {
    debug8(printf("Final: Trim left pairi %d, score %d, trim_left %d\n",pairi,score,trim_left));
  }
  debug8(printf("\n"));


  /* trim */
  if (trim_right == 0) {
    *trim3p = false;
  } else {
    *trim3p = true;
  }

  if (trim_left == 0) {
    *trim5p = false;
  } else {
    *trim5p = true;
  }

  i = 0;
  while (i < trim_right) {
    pairs = Pairpool_pop(pairs,&this);
    i++;
  }

  while (i <= trim_left) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&this);
#ifdef WASTE
    path = Pairpool_push_existing(path,pairpool,pair);
#else
    trimmed = List_push_existing(trimmed,pairptr);
#endif
    i++;
  }

  debug8(Pair_dump_list(trimmed,/*zerobasedp*/true));

  return trimmed;
}


#ifdef GSNAP
static int quality_score_adj = 33; /* Default is Sanger */

/* Taken from mapq.c */

static float
mismatch_logprob[MAX_QUALITY_SCORE+1] =
  /* log(1/3*10^(-Q/10)) */
  {-1.098612,
   -1.328871, -1.559129, -1.789388, -2.019646, -2.249905,
   -2.480163, -2.710422, -2.940680, -3.170939, -3.401197,
   -3.631456, -3.861714, -4.091973, -4.322231, -4.552490,
   -4.782748, -5.013007, -5.243265, -5.473524, -5.703782,
   -5.934041, -6.164299, -6.394558, -6.624817, -6.855075,
   -7.085334, -7.315592, -7.545851, -7.776109, -8.006368,
   -8.236626, -8.466885, -8.697143, -8.927402, -9.157660,
   -9.387919, -9.618177, -9.848436, -10.078694, -10.308953};



/* Look also at Substring_compute_mapq */
float
Pair_compute_mapq (struct T *pairarray, int npairs, int trim_left, int trim_right, int querylength,
		   char *quality_string, bool trim_terminals_p) {
  float loglik = 0.0;
  int Q;
  T pair;
  int querypos, i;

  for (i = 0; i < npairs; i++) {
    pair = &(pairarray[i]);
    if (trim_terminals_p == true && pair->querypos < trim_left) {
      /* Skip */
    } else if (trim_terminals_p == true && pair->querypos >= querylength - trim_right) {
      /* Skip */
    } else if (pair->comp == MISMATCH_COMP) {
      /* printf("Got a mismatch at querypos %d, cdna %c\n",pair->querypos,pair->cdna); */
      querypos = pair->querypos;
      Q = (quality_string == NULL) ? MAX_QUALITY_SCORE : quality_string[querypos] - quality_score_adj;
      if (Q < 0) {
	fprintf(stderr,"Warning: quality score %c (ASCII %d) - %d (quality-zero-score) = %d, which is less than 0.  May need to specify --quality-protocol or --quality-zero-score\n",
		quality_string[querypos],(int) quality_string[querypos],quality_score_adj, Q);
	/* fprintf(stderr,"Position %d in %d-%d of %s\n",querypos,querystart,queryend,quality_string); */
	Q = 0;
      } else if (Q > MAX_QUALITY_SCORE_INPUT) {
	fprintf(stderr,"Warning: quality score %c (ASCII %d) - %d (quality-zero-score) = %d, which exceeds %d.  May need to specify --quality-protocol or --quality-zero-score\n",
		quality_string[querypos],(int) quality_string[querypos],quality_score_adj,Q,MAX_QUALITY_SCORE_INPUT);
	/* fprintf(stderr,"Position %d in %d-%d of %s\n",querypos,querystart,queryend,quality_string); */
	Q = MAX_QUALITY_SCORE;
      } else if (Q > MAX_QUALITY_SCORE) {
	Q = MAX_QUALITY_SCORE;
      }

      loglik += mismatch_logprob[Q];
    }
  }

  return loglik;
}



Overlap_T
Pair_gene_overlap (struct T *ptr /* = pairarray */, int npairs, IIT_T genes_iit, int divno,
		   bool favor_multiexon_p) {
  struct T *this = NULL;
  bool in_exon = false;
  int i;
  Chrpos_T low, high, exon_genomestart, exon_genomeend, last_genomepos = -1U;
  bool foundp = false;
  Overlap_T overlap;

  for (i = 0; i < npairs; i++) {
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_genomeend = last_genomepos;
	if (exon_genomestart <= exon_genomeend) {
	  low = exon_genomestart;
	  high = exon_genomeend;
	} else {
	  low = exon_genomeend;
	  high = exon_genomestart;
	}

	overlap = IIT_gene_overlap(genes_iit,divno,low,high,favor_multiexon_p);
	if (overlap == NO_KNOWN_GENE) {
	  /* Keep searching */
	} else if (overlap == KNOWN_GENE) {
	  if (favor_multiexon_p == true) {
	    /* Keep looking for a multiexon gene */
	    foundp = true;
	  } else {
	    return KNOWN_GENE;
	  }
	} else {
	  return KNOWN_GENE_MULTIEXON;
	}
	  
	in_exon = false;
      }

    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_genomestart = this->genomepos;
	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  exon_genomeend = last_genomepos;
  if (exon_genomestart <= exon_genomeend) {
    low = exon_genomestart;
    high = exon_genomeend;
  } else {
    low = exon_genomeend;
    high = exon_genomestart;
  }

  overlap = IIT_gene_overlap(genes_iit,divno,low,high,favor_multiexon_p);
  if (overlap == NO_KNOWN_GENE) {
    /* Keep searching */
  } else if (overlap == KNOWN_GENE) {
    if (favor_multiexon_p == true) {
      /* Keep looking for a multiexon gene */
      foundp = true;
    } else {
      return KNOWN_GENE;
    }
  } else {
    return KNOWN_GENE_MULTIEXON;
  }

  if (foundp == true) {
    return KNOWN_GENE;
  } else {
    return NO_KNOWN_GENE;
  }
}


void
Pair_init (int quality_score_adj_in) {
  quality_score_adj = quality_score_adj_in;
  return;
}

#endif


