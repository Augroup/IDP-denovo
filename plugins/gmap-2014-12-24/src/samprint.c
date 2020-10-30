static char rcsid[] = "$Id: samprint.c 160877 2015-03-13 00:31:23Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "samprint.h"
#include "samflags.h"
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "mem.h"
#include "complement.h"
#include "stage3hr.h"
#include "mapq.h"
#include "assert.h"


#define SANGER_ILLUMINA_DIFF 31
/* #define PRINT_AMBIG_COORDS 1 */

/* BAM appears to truncate the H information on the ends of a cigar */
/* Also, this provides the information needed for getting term information */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* compute_cigar */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* print_md_string */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


/* overlap */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* compute_chrpos */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


static bool quiet_if_excessive_p;
static int maxpaths_report;
static char *failedinput_root;
static bool fastq_format_p;
static bool hide_soft_clips_p;

static bool sam_multiple_primaries_p;
static bool force_xs_direction_p;
static bool md_lowercase_variant_p;
static IIT_T snps_iit;

void
SAM_setup (bool quiet_if_excessive_p_in, int maxpaths_report_in,
	   char *failedinput_root_in, bool fastq_format_p_in, bool hide_soft_clips_p_in,
	   bool sam_multiple_primaries_p_in,
	   bool force_xs_direction_p_in, bool md_lowercase_variant_p_in, IIT_T snps_iit_in) {
  quiet_if_excessive_p = quiet_if_excessive_p_in;
  failedinput_root = failedinput_root_in;
  fastq_format_p = fastq_format_p_in;
  hide_soft_clips_p = hide_soft_clips_p_in;
  maxpaths_report = maxpaths_report_in;
  sam_multiple_primaries_p = sam_multiple_primaries_p_in;
  force_xs_direction_p = force_xs_direction_p_in;
  md_lowercase_variant_p = md_lowercase_variant_p_in;
  snps_iit = snps_iit_in;
  return;
}


static FILE *fp_failedinput_1;
static FILE *fp_failedinput_2;

static FILE *fp_nomapping;
static FILE *fp_unpaired_uniq;
static FILE *fp_unpaired_circular;
static FILE *fp_unpaired_transloc;
static FILE *fp_unpaired_mult;
static FILE *fp_unpaired_mult_xs_1;
static FILE *fp_unpaired_mult_xs_2;
static FILE *fp_halfmapping_uniq;
static FILE *fp_halfmapping_circular;
static FILE *fp_halfmapping_transloc;
static FILE *fp_halfmapping_mult;
static FILE *fp_halfmapping_mult_xs_1;
static FILE *fp_halfmapping_mult_xs_2;
static FILE *fp_paired_uniq_circular;
static FILE *fp_paired_uniq_inv;
static FILE *fp_paired_uniq_scr;
static FILE *fp_paired_uniq_long;
static FILE *fp_paired_mult;
static FILE *fp_paired_mult_xs_1;
static FILE *fp_paired_mult_xs_2;
static FILE *fp_concordant_uniq;
static FILE *fp_concordant_circular;
static FILE *fp_concordant_transloc;
static FILE *fp_concordant_mult;
static FILE *fp_concordant_mult_xs_1;
static FILE *fp_concordant_mult_xs_2;


void
SAM_file_setup_single (FILE *fp_failedinput_in, FILE *fp_nomapping_in,
		       FILE *fp_unpaired_uniq_in, FILE *fp_unpaired_circular_in, FILE *fp_unpaired_transloc_in,
		       FILE *fp_unpaired_mult_in, FILE *fp_unpaired_mult_xs_1_in) {

  fp_failedinput_1 = fp_failedinput_in;

  fp_nomapping = fp_nomapping_in;
  fp_unpaired_uniq = fp_unpaired_uniq_in;
  fp_unpaired_circular = fp_unpaired_circular_in;
  fp_unpaired_transloc = fp_unpaired_transloc_in;
  fp_unpaired_mult = fp_unpaired_mult_in;
  fp_unpaired_mult_xs_1 = fp_unpaired_mult_xs_1_in;

  return;
}

void
SAM_file_setup_paired (FILE *fp_failedinput_1_in, FILE *fp_failedinput_2_in, FILE *fp_nomapping_in, 
		       FILE *fp_halfmapping_uniq_in, FILE *fp_halfmapping_circular_in, FILE *fp_halfmapping_transloc_in,
		       FILE *fp_halfmapping_mult_in, FILE *fp_halfmapping_mult_xs_1_in, FILE *fp_halfmapping_mult_xs_2_in,
		       FILE *fp_paired_uniq_circular_in, FILE *fp_paired_uniq_inv_in, FILE *fp_paired_uniq_scr_in,
		       FILE *fp_paired_uniq_long_in, FILE *fp_paired_mult_in, FILE *fp_paired_mult_xs_1_in, FILE *fp_paired_mult_xs_2_in,
		       FILE *fp_concordant_uniq_in, FILE *fp_concordant_circular_in, FILE *fp_concordant_transloc_in, 
		       FILE *fp_concordant_mult_in, FILE *fp_concordant_mult_xs_1_in, FILE *fp_concordant_mult_xs_2_in) {

  fp_failedinput_1 = fp_failedinput_1_in;
  fp_failedinput_2 = fp_failedinput_2_in;

  fp_nomapping = fp_nomapping_in;
  fp_halfmapping_uniq = fp_halfmapping_uniq_in;
  fp_halfmapping_circular = fp_halfmapping_circular_in;
  fp_halfmapping_transloc = fp_halfmapping_transloc_in;
  fp_halfmapping_mult = fp_halfmapping_mult_in;
  fp_halfmapping_mult_xs_1 = fp_halfmapping_mult_xs_1_in;
  fp_halfmapping_mult_xs_2 = fp_halfmapping_mult_xs_2_in;
  fp_paired_uniq_circular = fp_paired_uniq_circular_in;
  fp_paired_uniq_inv = fp_paired_uniq_inv_in;
  fp_paired_uniq_scr = fp_paired_uniq_scr_in;
  fp_paired_uniq_long = fp_paired_uniq_long_in;
  fp_paired_mult = fp_paired_mult_in;
  fp_paired_mult_xs_1 = fp_paired_mult_xs_1_in;
  fp_paired_mult_xs_2 = fp_paired_mult_xs_2_in;
  fp_concordant_uniq = fp_concordant_uniq_in;
  fp_concordant_circular = fp_concordant_circular_in;
  fp_concordant_transloc = fp_concordant_transloc_in;
  fp_concordant_mult = fp_concordant_mult_in;
  fp_concordant_mult_xs_1 = fp_concordant_mult_xs_1_in;
  fp_concordant_mult_xs_2 = fp_concordant_mult_xs_2_in;

  return;
}

void
SAM_file_setup_all (FILE *fp_failedinput_1_in, FILE *fp_failedinput_2_in, FILE *fp_nomapping_in,
		    FILE *fp_unpaired_uniq_in, FILE *fp_unpaired_circular_in, FILE *fp_unpaired_transloc_in,
		    FILE *fp_unpaired_mult_in, FILE *fp_unpaired_mult_xs_1_in, FILE *fp_unpaired_mult_xs_2_in,
		    FILE *fp_halfmapping_uniq_in, FILE *fp_halfmapping_circular_in, FILE *fp_halfmapping_transloc_in,
		    FILE *fp_halfmapping_mult_in, FILE *fp_halfmapping_mult_xs_1_in, FILE *fp_halfmapping_mult_xs_2_in,
		    FILE *fp_paired_uniq_circular_in, FILE *fp_paired_uniq_inv_in, FILE *fp_paired_uniq_scr_in,
		    FILE *fp_paired_uniq_long_in, FILE *fp_paired_mult_in, FILE *fp_paired_mult_xs_1_in, FILE *fp_paired_mult_xs_2_in,
		    FILE *fp_concordant_uniq_in, FILE *fp_concordant_circular_in, FILE *fp_concordant_transloc_in, 
		    FILE *fp_concordant_mult_in, FILE *fp_concordant_mult_xs_1_in, FILE *fp_concordant_mult_xs_2_in) {
  
  fp_failedinput_1 = fp_failedinput_1_in;
  fp_failedinput_2 = fp_failedinput_2_in;

  fp_nomapping = fp_nomapping_in;
  fp_unpaired_uniq = fp_unpaired_uniq_in;
  fp_unpaired_circular = fp_unpaired_circular_in;
  fp_unpaired_transloc = fp_unpaired_transloc_in;
  fp_unpaired_mult = fp_unpaired_mult_in;
  fp_unpaired_mult_xs_1 = fp_unpaired_mult_xs_1_in;
  fp_unpaired_mult_xs_2 = fp_unpaired_mult_xs_2_in;
  fp_halfmapping_uniq = fp_halfmapping_uniq_in;
  fp_halfmapping_circular = fp_halfmapping_circular_in;
  fp_halfmapping_transloc = fp_halfmapping_transloc_in;
  fp_halfmapping_mult = fp_halfmapping_mult_in;
  fp_halfmapping_mult_xs_1 = fp_halfmapping_mult_xs_1_in;
  fp_halfmapping_mult_xs_2 = fp_halfmapping_mult_xs_2_in;
  fp_paired_uniq_circular = fp_paired_uniq_circular_in;
  fp_paired_uniq_inv = fp_paired_uniq_inv_in;
  fp_paired_uniq_scr = fp_paired_uniq_scr_in;
  fp_paired_uniq_long = fp_paired_uniq_long_in;
  fp_paired_mult = fp_paired_mult_in;
  fp_paired_mult_xs_1 = fp_paired_mult_xs_1_in;
  fp_paired_mult_xs_2 = fp_paired_mult_xs_2_in;
  fp_concordant_uniq = fp_concordant_uniq_in;
  fp_concordant_circular = fp_concordant_circular_in;
  fp_concordant_transloc = fp_concordant_transloc_in;
  fp_concordant_mult = fp_concordant_mult_in;
  fp_concordant_mult_xs_1 = fp_concordant_mult_xs_1_in;
  fp_concordant_mult_xs_2 = fp_concordant_mult_xs_2_in;

  return;
}


unsigned int
SAM_compute_flag (bool plusp, Stage3end_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, int npaths_mate,
		  int absmq_score, int first_absmq, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;

  debug(printf("Resulttype: %s\n",Resulttype_string(resulttype)));

  if (npaths == 0) {
    debug(printf("npaths = 0, so QUERY_UNMAPPED %d\n",QUERY_UNMAPPED));
    flag |= QUERY_UNMAPPED;
  } else if (plusp == invertp) {
    debug(printf("plusp %d and invertp %d, so QUERY_MINUSP %d\n",
		 plusp,invertp,QUERY_MINUSP));
    flag |= QUERY_MINUSP;
  }

  if (resulttype == SINGLEEND_NOMAPPING || resulttype == SINGLEEND_UNIQ || resulttype == SINGLEEND_TRANSLOC || resulttype == SINGLEEND_MULT) {
    /* No first or second read or mate */
  } else {
    debug(printf("PAIRED_READ %d\n",PAIRED_READ));
    flag |= PAIRED_READ;
    if (first_read_p == true) {
      debug(printf("FIRST_READ %d\n",FIRST_READ_P));
      flag |= FIRST_READ_P;
    } else {
      debug(printf("SECOND_READ %d\n",SECOND_READ_P));
      flag |= SECOND_READ_P;
    }
    if (npaths_mate == 0) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (quiet_if_excessive_p && npaths_mate > maxpaths_report) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (mate == NULL) {
      /* Unpaired; no mate.  Not clear if should be MATE_UNMAPPED. */

    } else if (npaths == 0) {
      /* Need to check npaths == 0 in case clipping of overlaps results in a nomapping */
      if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
      /* Can distinguish concordant mappings by presence of insert length */
      debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
      flag |= PAIRED_MAPPING;
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else if (resulttype == PAIRED_UNIQ || resulttype == PAIRED_MULT) {
      /* Note: We are counting PAIRED_UNIQ and PAIRED_MULT as "paired" mappings.
	 However, we are no longer counting UNPAIRED_UNIQ as a "paired" mapping. */
      debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
      flag |= PAIRED_MAPPING;
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else {
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }
    }
  }

  if (pathnum > 1) {
    if (sam_multiple_primaries_p == false) {
      debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
      flag |= NOT_PRIMARY;
    } else if (absmq_score != first_absmq) {
      debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
      flag |= NOT_PRIMARY;
    } else {
      /* Just as good as first alignment, so don't mark as secondary */
    }
  }

  return flag;
}


#if 0
/* Replaced by new adjust_hardclips procedure in stage3hr.c */

/* Shifts low_querypos and high_querypos upward until a matching
   nucleotide is found from both hits.  If not found, the shifts
   low_querypos and high_querypos downward until a matching nucleotide
   is found. */

static void
adjust_hardclips (int *hardclip_low, Stage3end_T hit_low, int low_querylength,
		  int *hardclip_high, Stage3end_T hit_high, int high_querylength) {
  int orig_hardclip_low, orig_hardclip_high;
  Substring_T low_substring, high_substring;
  struct Pair_T *low_pairarray, *high_pairarray;
  int low_querystart, low_queryend, low_npairs, high_npairs;
  int low_querypos, high_querypos;
  bool plusp;

  debug3(printf("Entering adjust_hardclips with hardclip_low %d, hardclip_high %d\n",
		*hardclip_low,*hardclip_high));
  orig_hardclip_low = *hardclip_low;
  orig_hardclip_high = *hardclip_high;

  plusp = Stage3end_plusp(hit_low);

  if (Stage3end_hittype(hit_low) == GMAP && Stage3end_hittype(hit_high) == GMAP) {
    debug3(printf("Dual GMAP\n"));
    low_pairarray = Stage3end_pairarray(hit_low);
    low_npairs = Stage3end_npairs(hit_low);
    high_pairarray = Stage3end_pairarray(hit_high);
    high_npairs = Stage3end_npairs(hit_high);

    if (plusp == true) {
      if (hide_soft_clips_p == true) {
	low_querystart = 0;
      } else {
	low_querystart = Stage3end_gmap_querystart(hit_low);
      }
      if (*hardclip_low > low_querystart) {
	low_querypos = *hardclip_low;
	high_querypos = high_querylength - 1 - (*hardclip_high);
	while (low_querypos < low_querylength && high_querypos < high_querylength &&
	       (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false || 
		Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false)) {
	  (*hardclip_low)++;
	  (*hardclip_high)--;
	  low_querypos++;
	  high_querypos++;
	}
	if (low_querypos >= low_querylength || high_querypos >= high_querylength) {
	  debug3(printf("Querypos increase failed.  Trying querypos decrease.\n"));
	  (*hardclip_low)--;
	  (*hardclip_high)++;
	  low_querypos--;
	  high_querypos--;
	  while (low_querypos > 0 && high_querypos > 0 &&
		 (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false || 
		  Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false)) {
	    (*hardclip_low)--;
	    (*hardclip_high)++;
	    low_querypos--;
	    high_querypos--;
	  }
	  if (low_querypos <= 0 || high_querypos <= 0) {
	    *hardclip_low = orig_hardclip_low;
	    *hardclip_high = orig_hardclip_high;
	  }
	}
      }

    } else {
      if (hide_soft_clips_p == true) {
	low_queryend = low_querylength - 1;
      } else {
	low_queryend = Stage3end_gmap_queryend(hit_low);
      }
      if (low_querylength - *hardclip_low < low_queryend) {
	low_querypos = low_querylength - 1 - (*hardclip_low);
	high_querypos = *hardclip_high;
	while (low_querypos < low_querylength && high_querypos < high_querylength &&
	       (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false || 
		Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false)) {
	  (*hardclip_low)--;
	  (*hardclip_high)++;
	  low_querypos++;
	  high_querypos++;
	}
	if (low_querypos >= low_querylength || high_querypos >= high_querylength) {
	  debug3(printf("Querypos increase failed.  Trying querypos decrease.\n"));
	  (*hardclip_low)++;
	  (*hardclip_high)--;
	  low_querypos--;
	  high_querypos--;
	  while (low_querypos > 0 && high_querypos > 0 &&
		 (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false || 
		  Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false)) {
	    (*hardclip_low)++;
	    (*hardclip_high)--;
	    low_querypos--;
	    high_querypos--;
	  }
	  if (low_querypos <= 0 || high_querypos <= 0) {
	    *hardclip_low = orig_hardclip_low;
	    *hardclip_high = orig_hardclip_high;
	  }
	}
      }
    }

  } else if (Stage3end_hittype(hit_low) == GMAP) {
    debug3(printf("Low GMAP\n"));
    low_pairarray = Stage3end_pairarray(hit_low);
    low_npairs = Stage3end_npairs(hit_low);

    if (plusp == true) {
      if (hide_soft_clips_p == true) {
	low_querystart = 0;
      } else {
	low_querystart = Stage3end_gmap_querystart(hit_low);
      }
      if (*hardclip_low > low_querystart) {
	low_querypos = *hardclip_low;
	high_querypos = high_querylength - 1 - (*hardclip_high);
	high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	while (low_querypos < low_querylength && high_querypos < high_querylength &&
	       (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false || high_substring == NULL)) {
	  (*hardclip_low)++;
	  (*hardclip_high)--;
	  low_querypos++;
	  high_querypos++;
	  high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	}
	if (low_querypos >= low_querylength || high_querypos >= high_querylength) {
	  debug3(printf("Querypos increase failed.  Trying querypos decrease.\n"));
	  (*hardclip_low)--;
	  (*hardclip_high)++;
	  low_querypos--;
	  high_querypos--;
	  while (low_querypos > 0 && high_querypos > 0 &&
		 (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false || high_substring == NULL)) {
	    (*hardclip_low)--;
	    (*hardclip_high)++;
	    low_querypos--;
	    high_querypos--;
	    high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	  }
	  if (low_querypos <= 0 || high_querypos <= 0) {
	    *hardclip_low = orig_hardclip_low;
	    *hardclip_high = orig_hardclip_high;
	  }
	}
      }

    } else {
      if (hide_soft_clips_p == true) {
	low_queryend = low_querylength - 1;
      } else {
	low_queryend = Stage3end_gmap_queryend(hit_low);
      }
      if (low_querylength - *hardclip_low < low_queryend) {
	low_querypos = low_querylength - 1 - (*hardclip_low);
	high_querypos = *hardclip_high;
	high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	while (low_querypos < low_querylength && high_querypos < high_querylength &&
	       (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false || high_substring == NULL)) {
	  (*hardclip_low)--;
	  (*hardclip_high)++;
	  low_querypos++;
	  high_querypos++;
	  high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	}
	if (low_querypos >= low_querylength || high_querypos >= high_querylength) {
	  debug3(printf("Querypos increase failed.  Trying querypos decrease.\n"));
	  (*hardclip_low)++;
	  (*hardclip_high)--;
	  low_querypos--;
	  high_querypos--;
	  while (low_querypos > 0 && high_querypos > 0 &&
		 (Pairarray_contains_p(low_pairarray,low_npairs,low_querypos) == false || high_substring == NULL)) {
	    (*hardclip_low)++;
	    (*hardclip_high)--;
	    low_querypos--;
	    high_querypos--;
	    high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	  }
	  if (low_querypos <= 0 || high_querypos <= 0) {
	    *hardclip_low = orig_hardclip_low;
	    *hardclip_high = orig_hardclip_high;
	  }
	}
      }
    }

  } else if (Stage3end_hittype(hit_high) == GMAP) {
    debug3(printf("High GMAP\n"));
    high_pairarray = Stage3end_pairarray(hit_high);
    high_npairs = Stage3end_npairs(hit_high);

    if (plusp == true) {
      if (hide_soft_clips_p == true) {
	low_querystart = Substring_querystart_orig(Stage3end_substring_low(hit_low));
      } else {
	low_querystart = Substring_querystart(Stage3end_substring_low(hit_low));
      }
      if (*hardclip_low > low_querystart) {
	low_querypos = *hardclip_low;
	high_querypos = high_querylength - 1 - (*hardclip_high);
	low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	while (low_querypos < low_querylength && high_querypos < high_querylength &&
	       (low_substring == NULL || Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false)) {
	  (*hardclip_low)++;
	  (*hardclip_high)--;
	  low_querypos++;
	  high_querypos++;
	  low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	}
	if (low_querypos >= low_querylength || high_querypos >= high_querylength) {
	  debug3(printf("Querypos increase failed.  Trying querypos decrease.\n"));
	  (*hardclip_low)--;
	  (*hardclip_high)++;
	  low_querypos--;
	  high_querypos--;
	  while (low_querypos > 0 && high_querypos > 0 &&
	       (low_substring == NULL || Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false)) {
	    (*hardclip_low)--;
	    (*hardclip_high)++;
	    low_querypos--;
	    high_querypos--;
	    low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	  }
	  if (low_querypos <= 0 || high_querypos <= 0) {
	    *hardclip_low = orig_hardclip_low;
	    *hardclip_high = orig_hardclip_high;
	  }
	}
      }

    } else {
      if (hide_soft_clips_p == true) {
	low_queryend = Substring_queryend_orig(Stage3end_substring_low(hit_low));
      } else {
	low_queryend = Substring_queryend(Stage3end_substring_low(hit_low));
      }
      if (low_querylength - *hardclip_low < low_queryend) {
	low_querypos = low_querylength - 1 - (*hardclip_low);
	high_querypos = *hardclip_high;
	low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	while (low_querypos < low_querylength && high_querypos < high_querylength &&
	       (low_substring == NULL || Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false)) {
	  (*hardclip_low)--;
	  (*hardclip_high)++;
	  low_querypos++;
	  high_querypos++;
	  low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	}
	if (low_querypos >= low_querylength || high_querypos >= high_querylength) {
	  debug3(printf("Querypos increase failed.  Trying querypos decrease.\n"));
	  (*hardclip_low)++;
	  (*hardclip_high)--;
	  low_querypos--;
	  high_querypos--;
	  while (low_querypos > 0 && high_querypos > 0 &&
		 (low_substring == NULL || Pairarray_contains_p(high_pairarray,high_npairs,high_querypos) == false)) {
	    (*hardclip_low)++;
	    (*hardclip_high)--;
	    low_querypos--;
	    high_querypos--;
	    low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	  }
	  if (low_querypos <= 0 || high_querypos <= 0) {
	    *hardclip_low = orig_hardclip_low;
	    *hardclip_high = orig_hardclip_high;
	  }
	}
      }
    }

  } else {
    if (plusp == true) {
      debug3(printf("Both substrings, plus\n"));

      if (hide_soft_clips_p == true) {
	low_querystart = Substring_querystart_orig(Stage3end_substring_low(hit_low));
      } else {
	low_querystart = Substring_querystart(Stage3end_substring_low(hit_low));
      }

      if (*hardclip_low > low_querystart) {
	low_querypos = *hardclip_low;
	high_querypos = high_querylength - 1 - *hardclip_high;
	low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	while (low_querypos < low_querylength && high_querypos < high_querylength &&
	       (low_substring == NULL || high_substring == NULL)) {
	  (*hardclip_low)++;
	  (*hardclip_high)--;
	  low_querypos++;
	  high_querypos++;
	  low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	  high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	}
	if (low_querypos >= low_querylength || high_querypos >= high_querylength) {
	  debug3(printf("Querypos increase failed.  Tryiing querypos decrease.\n"));
	  (*hardclip_low)--;
	  (*hardclip_high)++;
	  low_querypos--;
	  high_querypos--;
	  while (low_querypos > 0 && high_querypos > 0 &&
		 (low_substring == NULL || high_substring == NULL)) {
	    (*hardclip_low)--;
	    (*hardclip_high)++;
	    low_querypos--;
	    high_querypos--;
	    low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	    high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	  }
	  if (low_querypos <= 0 || high_querypos <= 0) {
	    *hardclip_low = orig_hardclip_low;
	    *hardclip_high = orig_hardclip_high;
	  }
	}
      }

    } else {
      debug3(printf("Both substrings, minus\n"));

      if (hide_soft_clips_p == true) {
	low_queryend = Substring_queryend_orig(Stage3end_substring_low(hit_low));
      } else {
	low_queryend = Substring_queryend(Stage3end_substring_low(hit_low));
      }

      if (low_querylength - *hardclip_low < low_queryend) {
	low_querypos = low_querylength - 1 - (*hardclip_low);
	high_querypos = *hardclip_high;
	low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	while (low_querypos < low_querylength && high_querypos < high_querylength &&
	       (low_substring == NULL || high_substring == NULL)) {
	  (*hardclip_low)--;
	  (*hardclip_high)++;
	  low_querypos++;
	  high_querypos++;
	  low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	  high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	}
	if (low_querypos >= low_querylength || high_querypos >= high_querylength) {
	  debug3(printf("Querypos increase failed.  Trying querypos decrease.\n"));
	  (*hardclip_low)++;
	  (*hardclip_high)--;
	  low_querypos--;
	  high_querypos--;
	  while (low_querypos > 0 && high_querypos > 0 &&
		 (low_substring == NULL || high_substring == NULL)) {
	    (*hardclip_low)++;
	    (*hardclip_high)--;
	    low_querypos--;
	    high_querypos--;
	    low_substring = Stage3end_substring_containing(hit_low,low_querypos);
	    high_substring = Stage3end_substring_containing(hit_high,high_querypos);
	  }
	  if (low_querypos <= 0 || high_querypos <= 0) {
	    *hardclip_low = orig_hardclip_low;
	    *hardclip_high = orig_hardclip_high;
	  }
	}
      }
    }
  }
    
  debug3(printf("Exiting adjust_hardclips with hardclip_low %d, hardclip_high %d\n",
		*hardclip_low,*hardclip_high));

  return;
}
#endif



Chrpos_T
SAM_compute_chrpos (int hardclip_low, int hardclip_high, Stage3end_T this, int querylength) {
  Substring_T substring_hardclipped;
  Chrpos_T chrpos;
  int querystart, queryend;

  if (this == NULL) {
    return 0U;

  } else if (Stage3end_hittype(this) == GMAP) {
    chrpos = Pair_genomicpos_low(hardclip_low,hardclip_high,Stage3end_pairarray(this),Stage3end_npairs(this),
				 querylength,/*watsonp*/Stage3end_plusp(this),hide_soft_clips_p);

  } else if (Stage3end_plusp(this) == true) {
    substring_hardclipped = Stage3end_substring_containing(this,hardclip_low);
    debug4(printf("Plus: Substring containing hardclip_low %d is %d..%d => %u..%u\n",
		  hardclip_low,Substring_querystart(substring_hardclipped),Substring_queryend(substring_hardclipped),
		  Substring_chrstart(substring_hardclipped),Substring_chrend(substring_hardclipped)));

    /* Add 1 to report in 1-based coordinates */
    if (hide_soft_clips_p == true) {
      chrpos = Substring_alignstart(substring_hardclipped) - Substring_chroffset(substring_hardclipped) + 1U;
      querystart = Substring_querystart_orig(substring_hardclipped);
      /* queryend = Substring_queryend_orig(Stage3end_substring_high(this)); */
    } else {
      chrpos = Substring_alignstart_trim(substring_hardclipped) - Substring_chroffset(substring_hardclipped) + 1U;
      querystart = Substring_querystart(substring_hardclipped);
      /* queryend = Substring_queryend(Stage3end_substring_high(this)); */
    }

    debug4(printf("Incrementing chrpos by %d\n",hardclip_low - querystart));
    chrpos += hardclip_low - querystart;

  } else {
    substring_hardclipped = Stage3end_substring_containing(this,querylength - 1 - hardclip_low);
    debug4(printf("Minus: Substring containing querylength - 1 - hardclip_low %d is %d..%d => %u..%u\n",
		  querylength - 1 - hardclip_low,
		  Substring_querystart(substring_hardclipped),Substring_queryend(substring_hardclipped),
		  Substring_chrstart(substring_hardclipped),Substring_chrend(substring_hardclipped)));

    /* Add 1 to report in 1-based coordinates */
    if (hide_soft_clips_p == true) {
      chrpos = Substring_alignend(substring_hardclipped) - Substring_chroffset(substring_hardclipped) + 1U;
      /* querystart = Substring_querystart_orig(Stage3end_substring_high(this)); */
      queryend = Substring_queryend_orig(substring_hardclipped);
    } else {
      chrpos = Substring_alignend_trim(substring_hardclipped) - Substring_chroffset(substring_hardclipped) + 1U;
      /* querystart = Substring_querystart(Stage3end_substring_high(this)); */
      queryend = Substring_queryend(substring_hardclipped);
    }

    debug4(printf("Incrementing chrpos by %d\n",queryend - (querylength - hardclip_low)));
    chrpos += queryend - (querylength - hardclip_low);
  }
    
  return chrpos;
}

static void
print_chromosomal_pos (FILE *fp, Chrnum_T chrnum, Chrpos_T chrpos, Chrpos_T chrlength,
		       Univ_IIT_T chromosome_iit) {
  bool allocp;
  char *chr;

#if 0
  if (chrpos == 0U) {
    /* No mapping */
    fprintf(fp,"\t*\t0");
    return;
  }
#endif

  if (chrnum == 0) {
    /* Interchromosomal splice */
    fprintf(stderr,"Trying to print interchrosomal splice in one line\n");
    abort();

  } else {
    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
    fprintf(fp,"\t%s",chr);
    if (allocp == true) {
      FREE(chr);
    }

    /* chrpos already in 1-based coordinates */
    if (chrpos > chrlength) {
      fprintf(fp,"\t%u",chrpos - chrlength /*+1U*/);
    } else {
      fprintf(fp,"\t%u",chrpos /*+1U*/);
    }
    return;
  }
}

static void
print_mate_chromosomal_pos (FILE *fp, Chrnum_T mate_chrnum, Chrnum_T mate_effective_chrnum,
			    Chrpos_T mate_chrpos, Chrpos_T mate_chrlength, Chrnum_T anchor_chrnum, Chrpos_T anchor_chrpos,
			    Univ_IIT_T chromosome_iit) {
  bool allocp;
  char *chr;

  if (mate_chrpos == 0U) {
    /* No mapping */
    fprintf(fp,"\t*\t0");
    return;

  } else if (mate_chrnum == 0) {
    /* Interchromosomal splice.  Choose effective chrnum. */
    if (anchor_chrpos > 0U && anchor_chrnum > 0 && mate_effective_chrnum == anchor_chrnum) {
      fprintf(fp,"\t=");
    } else {
      chr = Univ_IIT_label(chromosome_iit,mate_effective_chrnum,&allocp);
      fprintf(fp,"\t%s",chr);
      if (allocp == true) {
	FREE(chr);
      }
    }
    
    /* chrpos already in 1-based coordinates */
    if (mate_chrpos > mate_chrlength) {
      fprintf(fp,"\t%u",mate_chrpos - mate_chrlength /*+1U*/);
    } else {
      fprintf(fp,"\t%u",mate_chrpos /*+1U*/);
    }
    return;

  } else {
    if (anchor_chrpos > 0U && anchor_chrnum > 0 && mate_chrnum == anchor_chrnum) {
      fprintf(fp,"\t=");
    } else {
      chr = Univ_IIT_label(chromosome_iit,mate_chrnum,&allocp);
      fprintf(fp,"\t%s",chr);
      if (allocp == true) {
	FREE(chr);
      }
    }
    
    /* chrpos already in 1-based coordinates */
    if (mate_chrpos > mate_chrlength) {
      fprintf(fp,"\t%u",mate_chrpos - mate_chrlength /*+1U*/);
    } else {
      fprintf(fp,"\t%u",mate_chrpos /*+1U*/);
    }
    return;
  }
}





static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}



/* npaths could be non-zero, if user selected --quiet-if-excessive */
void
SAM_print_nomapping (FILE *fp, char *abbrev, Shortread_T queryseq, Stage3end_T mate, char *acc1, char *acc2,
		     Univ_IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int npaths, int npaths_mate, Chrpos_T mate_chrpos, int quality_shift,
		     char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag;


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }
  
  /* 2. FLAG */
  flag = SAM_compute_flag(/*plusp (NA)*/true,mate,resulttype,first_read_p,
			  /*pathnum*/0,/*npaths*/0,npaths_mate,
			  /*absmq_score*/0,/*first_absmq*/0,invertp,invert_mate_p);
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
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     /*anchor_chrnum*/0,/*anchor_chrpos*/0U,chromosome_iit);


  /* 9. ISIZE: Insert size */
  fprintf(fp,"\t0");

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Since there is no mapping, we print the original query sequence. */
  fprintf(fp,"\t");
  if (invertp == false) {
    Shortread_print_chopped(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }
  
  /* 12. TAGS: NH */
  if (npaths > 0) {
    fprintf(fp,"\tNH:i:%d",npaths);
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  fprintf(fp,"\n");

  return;
}


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

static List_T
push_token (List_T tokens, char *token) {
  char *copy;

  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);
  return List_push(tokens,(void *) copy);
}


#if 0
/* Currently used for insertions and deletions */
static List_T
compute_cigar_old (List_T tokens, char type, int stringlength, int querypos, int querylength,
		   int hardclip_low, int hardclip_high, bool plusp, bool firstp, bool lastp) {
  char token[10];
  
  debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plusp %d\n",
		type,stringlength,querypos,querylength,hardclip_low,hardclip_high,plusp));

  if (firstp == true) {
    debug1(printf("firstp is true\n"));
    if (plusp == true) {
      if (hardclip_low > 0) {
	sprintf(token,"%dH",hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (querypos > hardclip_low) {
	sprintf(token,"%dS",querypos - hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    } else {
      if (hardclip_high > 0) {
	sprintf(token,"%dH",hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (querypos < querylength - hardclip_high) {
	sprintf(token,"%dS",querypos - hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  if (type == 'D' || type == 'N') {
    if (querypos < hardclip_low || querypos >= querylength - hardclip_high) {
      stringlength = 0;
    }

  } else if (plusp == true) {
    debug1(printf("Comparing querypos %d..%d against %d..%d\n",
		  querypos,querypos + stringlength,hardclip_low,querylength - hardclip_high));
    if (/* querypos < hardclip_low && */querypos + stringlength < hardclip_low) {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 1: stringlength 0\n"));
    } else if (querypos < hardclip_low) {
      if (querypos + stringlength < querylength - hardclip_high) {
	/* Print part after hardclip_low */
	stringlength = (querypos + stringlength) - hardclip_low;
	debug1(printf("Case 2: stringlength %d\n",stringlength));
      } else {
	/* Print part between hardclip_low and hardclip_high */
	stringlength = (querylength - hardclip_high) - hardclip_low;
	debug1(printf("Case 3: stringlength %d\n",stringlength));
      }
    } else if (querypos < querylength - hardclip_high) {
      if (querypos + stringlength >= querylength - hardclip_high) {
	/* Print up to hardclip_high */
	stringlength = (querylength - hardclip_high) - querypos;
	debug1(printf("Case 4: stringlength %d\n",stringlength));
      } else {
	/* Print full stringlength */
	debug1(printf("Case 5: stringlength %d\n",stringlength));
      }
    } else {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 6: stringlength 0\n"));
    }

  } else {
    debug1(printf("Comparing querypos %d..%d against %d..%d\n",
		  querypos,querypos - stringlength,hardclip_low,querylength - hardclip_high));
    if (/* querypos >= querylength - hardclip_high && */ querypos - stringlength >= querylength - hardclip_high) {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 1: stringlength 0\n"));
    } else if (querypos >= querylength - hardclip_high) {
      if (querypos - stringlength >= hardclip_low) {
	/* Print part after hardclip_high */
	stringlength = (querylength - hardclip_high) - (querypos - stringlength);
	debug1(printf("Case 2: stringlength %d\n",stringlength));
      } else {
	/* Print part between hardclip_low and hardclip_high */
	stringlength = (querylength - hardclip_high) - hardclip_low;
	debug1(printf("Case 3: stringlength %d\n",stringlength));
      }
    } else if (querypos >= hardclip_low) {
      if (querypos - stringlength < hardclip_low) {
	/* Print up to hardclip_low */
	stringlength = querypos - hardclip_low;
	debug1(printf("Case 4: stringlength %d\n",stringlength));
      } else {
	/* Print full stringlength */
	debug1(printf("Case 5: stringlength %d\n",stringlength));
      }
    } else {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 5: stringlength 0\n"));
    }
  }

  if (stringlength > 0) {
    sprintf(token,"%d%c",stringlength,type);
    debug1(printf("Pushing token %s\n",token));
    tokens = push_token(tokens,token);
  }

  if (lastp == true) {
    debug1(printf("lastp is true\n"));
    if (plusp == true) {
      querypos += stringlength;
      if (querypos < querylength - 1 - hardclip_high) {
	sprintf(token,"%dS",querylength - 1 - hardclip_high - querypos);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (hardclip_high > 0) {
	sprintf(token,"%dH",hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    } else {
      querypos -= stringlength;
      if (querypos > hardclip_low) {
	sprintf(token,"%dS",hardclip_low - querypos);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (hardclip_low > 0) {
	sprintf(token,"%dH",hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  return tokens;
}
#endif


/* Currently used for insertions and deletions */
static List_T
compute_cigar (List_T tokens, char type, int stringlength, int querypos, int querylength,
	       int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  char token[10];
  
  if (plusp == true) {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	sprintf(token,"%d%c",matchlength,type);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }

  } else {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	sprintf(token,"%d%c",matchlength,type);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  return tokens;
}


/* Modified from compute_cigar */
static Intlist_T
compute_cigar_types_only (Intlist_T types, char type, int stringlength, int querypos, int querylength,
			  int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  
  if (plusp == true) {
    debug1(printf("\nEntering compute_cigar_types_only with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	types = Intlist_push(types,type);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
    }

  } else {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	types = Intlist_push(types,type);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	types = Intlist_push(types,'H');
      }
    }
  }

  return types;
}


static void
print_cigar (FILE *fp, char type, int stringlength, int querypos, int querylength,
	     int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  
  if (plusp == true) {
    debug1(printf("\nEntering print_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	fprintf(fp,"%dH",cliplength);
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	fprintf(fp,"%d%c",matchlength,type);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	fprintf(fp,"%dH",cliplength);
      }
    }

  } else {
    debug1(printf("\nEntering print_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	fprintf(fp,"%dH",cliplength);
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	fprintf(fp,"%d%c",matchlength,type);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	fprintf(fp,"%dH",cliplength);
      }
    }
  }

  return;
}


static int
print_md_string (bool *printp, int *nmismatches_refdiff, int *nmismatches_bothdiff,
		 FILE *fp, int matchlength, char *genomicfwd_refdiff, char *genomicfwd_bothdiff,
		 int stringlength, int querypos, int querylength,
		 int hardclip_low, int hardclip_high, bool plusp, bool lastp) {
  int starti, endi, i;
  int local_nmismatches = 0;
  bool hardclip_end_p = false;

  if (plusp == true) {
    debug2(printf("\nEntering md_string with matchlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus: %s ref, %s both\n",
		  matchlength,querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff,genomicfwd_bothdiff));
    if (hardclip_low == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_low > querypos) {
      /* startpos = hardclip_low; */
      starti = hardclip_low - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_low %d - querypos %d\n",
		    starti,hardclip_low,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      /* endpos = querylength - hardclip_high; */
      endi = (querylength - hardclip_high) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_high %d) - querypos %d\n",
		    endi,querylength,hardclip_high,querypos));
    } else {
      /* endpos = querypos + stringlength; */
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd_refdiff == NULL) {
      if (endi > starti) {
	matchlength += (endi - starti);
      }

    } else if (md_lowercase_variant_p == false) {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",genomicfwd_refdiff[i]); /* Leave as lower case */
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;
    }

  } else {
    debug2(printf("\nEntering md_string with matchlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus: %s ref, %s both\n",
		  matchlength,querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff,genomicfwd_bothdiff));
    querypos = querylength - querypos - stringlength;
    debug2(printf("  Revising querypos to be %d\n",querypos));

    if (hardclip_low == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_low > querypos) {
      /* startpos = hardclip_low; */
      starti = hardclip_low - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_low %d - querypos %d\n",
		    starti,hardclip_low,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      /* endpos = querylength - hardclip_high; */
      endi = (querylength - hardclip_high) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_high %d) - querypos %d\n",
		    endi,querylength,hardclip_high,querypos));
    } else {
      /* endpos = querypos + stringlength; */
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd_refdiff == NULL) {
      if (endi > starti) {
	matchlength += (endi - starti);
      }

    } else if (md_lowercase_variant_p == false) {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else {
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",genomicfwd_refdiff[i]); /* Leave as lower case */
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;
    }
  }

  /* Update nmismatches_bothdiff */
  if (genomicfwd_bothdiff == NULL) {
    /* No change to nmismatches_bothdiff */
  } else if (genomicfwd_bothdiff == genomicfwd_refdiff) {
    *nmismatches_bothdiff += local_nmismatches;
  } else {
    for (i = starti; i < endi; i++) {
      if (!isupper(genomicfwd_bothdiff[i])) {
	*nmismatches_bothdiff += 1;
      }
    }
  }

  debug2(printf("  Ending with matchlength %d\n",matchlength));

  if (lastp == false) {
    return matchlength;
  } else if (matchlength > 0) {
    fprintf(fp,"%d",matchlength);
    *printp = true;
    return 0;
  } else {
    return 0;
  }
}


/* Copy also in pair.c for GMAP */
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



static void
print_single (FILE *fp, char *abbrev, Hittype_T hittype, Stage3end_T this, Stage3end_T mate,
	      char *acc1, char *acc2, int pathnum, int npaths,
	      int absmq_score, int first_absmq, int second_absmq, int mapq_score,
	      Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
	      Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip_low, int hardclip_high,
	      Resulttype_T resulttype, bool first_read_p,
	      int npaths_mate, int quality_shift,
	      char *sam_read_group_id, bool invertp, bool invert_mate_p, bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength, substring_start, substring_length;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  bool plusp, printp;

  debug(printf("print_single\n"));

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);
  substring = Stage3end_substring1(this);

  debug(printf("clipdir is %d, hardclip_low %d, hardclip_high %d\n",clipdir,hardclip_low,hardclip_high));

  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(Stage3end_plusp(this),mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,Stage3end_chrlength(this),chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  if (plusp == true) {
    if (hide_soft_clips_p == true && hittype != TERMINAL) {
      print_cigar(fp,/*type*/'M',
		  Substring_querystart(substring) + Substring_match_length(substring) +
		  (querylength - Substring_queryend(substring)),/*querypos*/0,querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		  /*querypos*/Substring_querystart(substring),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		  /*querypos*/Substring_queryend(substring),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (hide_soft_clips_p == true && hittype != TERMINAL) {
      print_cigar(fp,/*type*/'M',
		  (querylength - Substring_queryend(substring)) + 
		  Substring_match_length(substring) + Substring_querystart(substring),
		  /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		  /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		  /*querypos*/Substring_queryend(substring),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		  /*querypos*/Substring_querystart(substring),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
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


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH */
  if (hardclip_low > 0 || hardclip_high > 0) {
    fprintf(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\tMD:Z:");  
  printp = false;

  if (hide_soft_clips_p == true) {
    substring_start = Substring_querystart_orig(substring);
    substring_length = Substring_match_length_orig(substring);
  } else {
    substring_start = Substring_querystart(substring);
    substring_length = Substring_match_length(substring);
  }

  if ((genomicdir_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
    if (plusp == true) {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }

  } else if (plusp == true) {
    genomicdir_refdiff = Substring_genomic_refdiff(substring);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		    fp,/*matchlength*/0,&(genomicdir_refdiff[substring_start]),&(genomicdir_bothdiff[substring_start]),
		    substring_length,/*querypos*/substring_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == genomicdir_bothdiff) {
    genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		    fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		    substring_length,/*querypos*/substring_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREEA(genomicfwd_refdiff);

  } else {
    genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
    genomicfwd_bothdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
    make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		    fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		    substring_length,/*querypos*/substring_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREEA(genomicfwd_bothdiff);
    FREEA(genomicfwd_refdiff);
  }

  if (printp == false) {
    fprintf(fp,"0");
  }
  

  /* 12. TAGS: NH */
  fprintf(fp,"\tNH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  /* fprintf(fp,"\tNM:i:%d",Stage3end_nmismatches_refdiff(this)); */
  fprintf(fp,"\tNM:i:%d",nmismatches_refdiff);

  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\tXW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\tSM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\tXQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XC */
  if (circularp == true) {
    fprintf(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\tXG:Z:A");
  } else if (Stage3end_hittype(this) == TERMINAL) {
    fprintf(fp,"\tXG:Z:T");
  }

  fprintf(fp,"\n");
  return;
}


static bool
check_cigar_single (Hittype_T hittype, Stage3end_T this,
		    int querylength, int clipdir, int hardclip_low, int hardclip_high,
		    bool first_read_p, bool circularp) {
  bool result;
  Intlist_T cigar_types = NULL;
  Substring_T substring;
  bool plusp;

  plusp = Stage3end_plusp(this);
  substring = Stage3end_substring1(this);

  debug1(printf("clipdir is %d, hardclip_low %d, hardclip_high %d\n",clipdir,hardclip_low,hardclip_high));

  if (plusp == true) {
    if (hide_soft_clips_p == true && hittype != TERMINAL) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_querystart(substring) + Substring_match_length(substring) +
					     (querylength - Substring_queryend(substring)),/*querypos*/0,querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring),
					     /*querypos*/Substring_querystart(substring),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring),
					     /*querypos*/Substring_queryend(substring),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (hide_soft_clips_p == true && hittype != TERMINAL) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     (querylength - Substring_queryend(substring)) + 
					     Substring_match_length(substring) + Substring_querystart(substring),
					     /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
					     /*plusp*/false,/*lastp*/true);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring),
					     /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
					     /*plusp*/false,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring),
					     /*querypos*/Substring_queryend(substring),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring),
					     /*querypos*/Substring_querystart(substring),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }
  }

  result = check_cigar_types(cigar_types);

  Intlist_free(&cigar_types);
  return result;
}



static void
print_insertion (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths,
		 int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip_low, int hardclip_high,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring1_start, substring2_start, substring1_length, substring2_length, matchlength, nindels;
  bool plusp, printp;
  List_T cigar_tokens = NULL;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  substring1 = Stage3end_substring1(this);
  substring2 = Stage3end_substring2(this);

  nindels = Stage3end_nindels(this);

  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,Stage3end_chrlength(this),chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  if (plusp == true) {
    if (hide_soft_clips_p == true) {
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',
				   Substring_querystart(substring1) + Substring_match_length(substring1),
				   /*querypos*/0,querylength,hardclip_low,hardclip_high,
				   /*plusp*/true,/*lastp*/false);
    } else {
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'S',Substring_querystart(substring1),
				   /*querypos*/0,querylength,hardclip_low,hardclip_high,
				   /*plusp*/true,/*lastp*/false);
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',Substring_match_length(substring1),
				   /*querypos*/Substring_querystart(substring1),querylength,
				   hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    }

    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'I',nindels,
				 /*querypos*/Substring_queryend(substring1),querylength,
				 hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);

    if (hide_soft_clips_p == true) {
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',
				   Substring_match_length(substring2) + (querylength - Substring_queryend(substring2)),
				   /*querypos*/Substring_querystart(substring2),querylength,
				   hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    } else {
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',Substring_match_length(substring2),
				   /*querypos*/Substring_querystart(substring2),querylength,
				   hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'S',querylength - Substring_queryend(substring2),
				   /*querypos*/Substring_queryend(substring2),querylength,
				   hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (hide_soft_clips_p == true) {
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',
				   (querylength - Substring_queryend(substring2)) +
				   Substring_match_length(substring2),
				   /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
				   /*plusp*/false,/*lastp*/false);
    } else {
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'S',querylength - Substring_queryend(substring2),
				   /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
				   /*plusp*/false,/*lastp*/false);
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',Substring_match_length(substring2),
				   /*querypos*/Substring_queryend(substring2),querylength,
				   hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }

    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'I',nindels,
				 /*querypos*/Substring_querystart(substring2),querylength,
				 hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);

    if (hide_soft_clips_p == true) {
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',
				   Substring_match_length(substring1) +
				   Substring_querystart(substring1),
				   /*querypos*/Substring_queryend(substring1),querylength,
				   hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    } else {
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',Substring_match_length(substring1),
				   /*querypos*/Substring_queryend(substring1),querylength,
				   hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      cigar_tokens = compute_cigar(cigar_tokens,/*type*/'S',Substring_querystart(substring1),
				   /*querypos*/Substring_querystart(substring1),querylength,
				   hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }
  }
  cigar_tokens = Pair_clean_cigar(cigar_tokens,/*watsonp*/true);
  print_tokens_sam(fp,cigar_tokens);
  List_free(&cigar_tokens);


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
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

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH */
  if (hardclip_low > 0 || hardclip_high > 0) {
    fprintf(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\tMD:Z:");
  printp = false;

  if (hide_soft_clips_p == true) {
    substring1_start = Substring_querystart_orig(substring1);
    substring1_length = Substring_match_length_orig(substring1);
    substring2_start = Substring_querystart_orig(substring2);
    substring2_length = Substring_match_length_orig(substring2);
  } else {
    substring1_start = Substring_querystart(substring1);
    substring1_length = Substring_match_length(substring1);
    substring2_start = Substring_querystart(substring2);
    substring2_length = Substring_match_length(substring2);
  }

  if (plusp == true) {
    genomicfwd_refdiff = Substring_genomic_refdiff(substring1);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring1);
    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
				  &(genomicfwd_refdiff[substring1_start]),&(genomicfwd_bothdiff[substring1_start]),
				  substring1_length,/*querypos*/substring1_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);

#if 0
    /* If MD string is supposed to include insertion, then uncomment this */
    matchlength += nindels;
#endif

    genomicfwd_refdiff = Substring_genomic_refdiff(substring2);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring2);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
		    &(genomicfwd_refdiff[substring2_start]),&(genomicfwd_bothdiff[substring2_start]),
		    substring2_length,/*querypos*/substring2_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
  } else {
    genomicdir_refdiff = Substring_genomic_refdiff(substring2);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring2);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				    substring2_length,/*querypos*/substring2_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring2_start]),substring2_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
				    substring2_length,/*querypos*/substring2_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }

#if 0
    /* If MD string is supposed to include insertion, then uncomment this */
    matchlength += nindels;
#endif

    genomicdir_refdiff = Substring_genomic_refdiff(substring1);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring1);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring1_length,/*querypos*/substring1_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring1_start]),substring1_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring1_length,/*querypos*/substring1_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }
  }

  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\tNH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  /* fprintf(fp,"\tNM:i:%d",Stage3end_nmismatches_refdiff(this)); */
  fprintf(fp,"\tNM:i:%d",nmismatches_refdiff + nindels);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\tXW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\tSM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\tXQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XC */
  if (circularp == true) {
    fprintf(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\tXG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}

static bool
check_cigar_insertion (Stage3end_T this, int querylength, int clipdir, int hardclip_low, int hardclip_high,
		       bool first_read_p, bool circularp) {
  bool result;
  Intlist_T cigar_types = NULL;
  Substring_T substring1, substring2;
  bool plusp;
  int nindels;

  plusp = Stage3end_plusp(this);

  substring1 = Stage3end_substring1(this);
  substring2 = Stage3end_substring2(this);

  nindels = Stage3end_nindels(this);

  if (plusp == true) {
    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_querystart(substring1) + Substring_match_length(substring1),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring1),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1),
					     /*querypos*/Substring_querystart(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    }

    cigar_types = compute_cigar_types_only(cigar_types,/*type*/'I',nindels,
					   /*querypos*/Substring_queryend(substring1),querylength,
					   hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);

    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_match_length(substring2) + (querylength - Substring_queryend(substring2)),
					     /*querypos*/Substring_querystart(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2),
					     /*querypos*/Substring_querystart(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring2),
					     /*querypos*/Substring_queryend(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     (querylength - Substring_queryend(substring2)) +
					     Substring_match_length(substring2),
					     /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
					     /*plusp*/false,/*lastp*/false);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring2),
					     /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
					     /*plusp*/false,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2),
					     /*querypos*/Substring_queryend(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }

    cigar_types = compute_cigar_types_only(cigar_types,/*type*/'I',nindels,
					   /*querypos*/Substring_querystart(substring2),querylength,
					   hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);

    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_match_length(substring1) +
					     Substring_querystart(substring1),
					     /*querypos*/Substring_queryend(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1),
					     /*querypos*/Substring_queryend(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring1),
					     /*querypos*/Substring_querystart(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }
  }

  result = check_cigar_types(cigar_types);

  Intlist_free(&cigar_types);
  return result;
}


static void
print_deletion (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		char *acc1, char *acc2, int pathnum, int npaths,
		int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip_low, int hardclip_high,
		Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicfwd_deletion,
    *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring1_start, substring2_start, substring1_length, substring2_length, nindels;
  bool plusp, printp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  substring1 = Stage3end_substring1(this);
  substring2 = Stage3end_substring2(this);


#if 0
  /* These cases are checked below */
  if (hardclip_low >= Substring_querystart(substring2)) {
    nindels = 0;
  } else if (querylength - hardclip_high <= Substring_queryend(substring1)) {
    nindels = 0;
  } else {
    nindels = Stage3end_nindels(this); /* nindels is positive */
  }
#else
  nindels = Stage3end_nindels(this); /* nindels is positive */
#endif


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,Stage3end_chrlength(this),chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  if (hide_soft_clips_p == true) {
    substring1_start = Substring_querystart_orig(substring1);
    substring1_length = Substring_match_length_orig(substring1);
    substring2_start = Substring_querystart_orig(substring2);
    substring2_length = Substring_match_length_orig(substring2);
  } else {
    substring1_start = Substring_querystart(substring1);
    substring1_length = Substring_match_length(substring1);
    substring2_start = Substring_querystart(substring2);
    substring2_length = Substring_match_length(substring2);
  }

  if (plusp == true) {
    if (hide_soft_clips_p == true) {
      if (/*nindels > 0 &&*/ hardclip_low < substring1_start + substring1_length && hardclip_high < querylength - substring2_start) {
	print_cigar(fp,/*type*/'M',
		    Substring_querystart(substring1) +
		    Substring_match_length(substring1),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	fprintf(fp,"%dD",nindels);
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(substring2) +
		    (querylength - Substring_queryend(substring2)),
		    /*querypos*/Substring_querystart(substring2),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      } else {
	print_cigar(fp,/*type*/'M',
		    Substring_querystart(substring1) + 
		    (Substring_match_length(substring1) +
		     Substring_match_length(substring2)) +
		    (querylength - Substring_queryend(substring2)),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/true);
      }


    } else {
      print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      if (/*nindels > 0 &&*/ hardclip_low < substring1_start + substring1_length && hardclip_high < querylength - substring2_start) {
	print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		    /*querypos*/Substring_querystart(substring1),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	fprintf(fp,"%dD",nindels);
	print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		    /*querypos*/Substring_querystart(substring2),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      } else {
	print_cigar(fp,/*type*/'M',Substring_match_length(substring1) + Substring_match_length(substring2),
		    /*querypos*/Substring_querystart(substring1),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      }
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		  /*querypos*/Substring_queryend(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (hide_soft_clips_p == true) {
      if (/*nindels > 0 &&*/ hardclip_low < querylength - substring2_start && hardclip_high < substring1_start + substring1_length) {
	print_cigar(fp,/*type*/'M',
		    (querylength - Substring_queryend(substring2)) + 
		    Substring_match_length(substring2),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false);
	fprintf(fp,"%dD",nindels);
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(substring1) +
		    Substring_querystart(substring1),
		    /*querypos*/Substring_querystart(substring2),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      } else {
	print_cigar(fp,/*type*/'M',
		    (querylength - Substring_queryend(substring2)) + 
		    (Substring_match_length(substring2) + Substring_match_length(substring1)) +
		    Substring_querystart(substring1),		    
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      }

    } else {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      if (/*nindels > 0 &&*/ hardclip_low < querylength - substring2_start && hardclip_high < substring1_start + substring1_length) {
	print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		    /*querypos*/Substring_queryend(substring2),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	fprintf(fp,"%dD",nindels);
	print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		    /*querypos*/Substring_querystart(substring2),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      } else {
	print_cigar(fp,/*type*/'M',Substring_match_length(substring2) + Substring_match_length(substring1),
		    /*querypos*/Substring_queryend(substring2),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      }
      print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		  /*querypos*/Substring_querystart(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
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


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH */
  if (hardclip_low > 0 || hardclip_high > 0) {
    fprintf(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\tMD:Z:");
  printp = false;

  if (plusp == true) {
    genomicfwd_refdiff = Substring_genomic_refdiff(substring1);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring1);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		    &(genomicfwd_refdiff[substring1_start]),&(genomicfwd_bothdiff[substring1_start]),
		    substring1_length,/*querypos*/substring1_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

    debug1(printf("\nhardclip_low %d, hardclip_high %d\n",hardclip_low,hardclip_high));
    debug1(printf("substring1_length %d, substring2_length %d\n",substring1_length,substring2_length));
    debug1(printf("substring1 %d..%d, substring2 %d..%d\n",
		  Substring_querystart(substring1),Substring_queryend(substring1),
		  Substring_querystart(substring2),Substring_queryend(substring2)));
    debug1(printf("trim1: %d..%d, trim2 %d..%d\n",
		  Substring_trim_left(substring1),Substring_trim_right(substring1),
		  Substring_trim_left(substring2),Substring_trim_right(substring2)));
    if (hardclip_low < substring1_start + substring1_length && hardclip_high < querylength - substring2_start) {
      fprintf(fp,"^%s",Stage3end_deletion_string(this));
    }

    genomicfwd_refdiff = Substring_genomic_refdiff(substring2);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring2);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		    &(genomicfwd_refdiff[substring2_start]),&(genomicfwd_bothdiff[substring2_start]),
		    substring2_length,/*querypos*/substring2_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    genomicdir_refdiff = Substring_genomic_refdiff(substring2);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring2);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring2_start]),substring2_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }


    debug1(printf("\nhardclip_low %d, hardclip_high %d\n",hardclip_low,hardclip_high));
    debug1(printf("substring2_length %d, substring1_length %d\n",substring2_length,substring1_length));
    debug1(printf("substring1 %d..%d, substring2 %d..%d\n",
		  Substring_querystart(substring1),Substring_queryend(substring1),
		  Substring_querystart(substring2),Substring_queryend(substring2)));
    debug1(printf("trim1: %d..%d, trim2 %d..%d\n",
		  Substring_trim_left(substring1),Substring_trim_right(substring1),
		  Substring_trim_left(substring2),Substring_trim_right(substring2)));
    if (hardclip_low < querylength - substring2_start && hardclip_high < substring1_start + substring1_length) {
      /* Deletion string: Potential problem if followed by a mismatch, but can be resolved by looking at CIGAR string */
      genomicfwd_deletion = (char *) MALLOCA((nindels+1) * sizeof(char));
      make_complement_buffered(genomicfwd_deletion,Stage3end_deletion_string(this),nindels);
      fprintf(fp,"^%s",genomicfwd_deletion);
      FREEA(genomicfwd_deletion);
    }

    genomicdir_refdiff = Substring_genomic_refdiff(substring1);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring1);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring1_length,/*querypos*/substring1_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring1_start]),substring1_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring1_length,/*querypos*/substring1_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }

  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\tNH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  /* fprintf(fp,"\tNM:i:%d",Stage3end_nmismatches_refdiff(this)); */
  fprintf(fp,"\tNM:i:%d",nmismatches_refdiff + nindels);

  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\tXW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\tSM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\tXQ:i:%d",absmq_score);
  
  /* 12. TAGS: X2 */
  fprintf(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XC */
  if (circularp == true) {
    fprintf(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\tXG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}

static bool
check_cigar_deletion (Stage3end_T this, int querylength, int clipdir, int hardclip_low, int hardclip_high,
		      bool first_read_p, bool circularp) {
  bool result;
  Intlist_T cigar_types = NULL;
  Substring_T substring1, substring2;
  int substring1_start, substring2_start, substring1_length;
  bool plusp;

  plusp = Stage3end_plusp(this);

  substring1 = Stage3end_substring1(this);
  substring2 = Stage3end_substring2(this);

  if (hide_soft_clips_p == true) {
    substring1_start = Substring_querystart_orig(substring1);
    substring1_length = Substring_match_length_orig(substring1);
    substring2_start = Substring_querystart_orig(substring2);
    /* substring2_length = Substring_match_length_orig(substring2); */
  } else {
    substring1_start = Substring_querystart(substring1);
    substring1_length = Substring_match_length(substring1);
    substring2_start = Substring_querystart(substring2);
    /* substring2_length = Substring_match_length(substring2); */
  }

  if (plusp == true) {
    if (hide_soft_clips_p == true) {
      if (/*nindels > 0 &&*/ hardclip_low < substring1_start + substring1_length && hardclip_high < querylength - substring2_start) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_querystart(substring1) +
					       Substring_match_length(substring1),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/false);
	cigar_types = Intlist_push(cigar_types,'D');
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_match_length(substring2) +
					       (querylength - Substring_queryend(substring2)),
					       /*querypos*/Substring_querystart(substring2),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_querystart(substring1) + 
					       (Substring_match_length(substring1) +
						Substring_match_length(substring2)) +
					       (querylength - Substring_queryend(substring2)),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/true);
      }


    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring1),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
      if (/*nindels > 0 &&*/ hardclip_low < substring1_start + substring1_length && hardclip_high < querylength - substring2_start) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1),
					       /*querypos*/Substring_querystart(substring1),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	cigar_types = Intlist_push(cigar_types,'D');
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2),
					       /*querypos*/Substring_querystart(substring2),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1) + Substring_match_length(substring2),
					       /*querypos*/Substring_querystart(substring1),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      }
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring2),
					     /*querypos*/Substring_queryend(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (hide_soft_clips_p == true) {
      if (/*nindels > 0 &&*/ hardclip_low < querylength - substring2_start && hardclip_high < substring1_start + substring1_length) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       (querylength - Substring_queryend(substring2)) + 
					       Substring_match_length(substring2),
					       /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
					       /*plusp*/false,/*lastp*/false);
	cigar_types = Intlist_push(cigar_types,'D');
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_match_length(substring1) +
					       Substring_querystart(substring1),
					       /*querypos*/Substring_querystart(substring2),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       (querylength - Substring_queryend(substring2)) + 
					       (Substring_match_length(substring2) + Substring_match_length(substring1)) +
					       Substring_querystart(substring1),		    
					       /*querypos*/querylength,querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      }

    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring2),
					     /*querypos*/querylength,querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      if (/*nindels > 0 &&*/ hardclip_low < querylength - substring2_start && hardclip_high < substring1_start + substring1_length) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2),
					       /*querypos*/Substring_queryend(substring2),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = Intlist_push(cigar_types,'D');
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1),
					       /*querypos*/Substring_querystart(substring2),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2) + Substring_match_length(substring1),
					       /*querypos*/Substring_queryend(substring2),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      }
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring1),
					     /*querypos*/Substring_querystart(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }
  }

  result = check_cigar_types(cigar_types);

  Intlist_free(&cigar_types);
  return result;
}


static void
halfdonor_dinucleotide (char *donor1, char *donor2, Substring_T donor) {
  bool sensep;
  char *genomic;
  int substring_start, substring_length;

  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  sensep = Substring_chimera_sensep(donor);

  substring_start = Substring_querystart(donor);
  genomic = Substring_genomic_refdiff(donor);

  if (sensep == true) {
    substring_length = Substring_match_length(donor);
    *donor1 = toupper(genomic[substring_start+substring_length]);
    *donor2 = toupper(genomic[substring_start+substring_length+1]);

  } else {			/* sensep == false */
    *donor2 = toupper(complCode[(int) genomic[substring_start-2]]);
    *donor1 = toupper(complCode[(int) genomic[substring_start-1]]);
  }

  return;
}

static void
halfacceptor_dinucleotide (char *acceptor2, char *acceptor1, Substring_T acceptor) {
  bool sensep;
  char *genomic;
  int substring_start, substring_length;

  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  sensep = Substring_chimera_sensep(acceptor);

  substring_start = Substring_querystart(acceptor);
  genomic = Substring_genomic_refdiff(acceptor);

  if (sensep == true) {
    *acceptor2 = toupper(genomic[substring_start-2]);
    *acceptor1 = toupper(genomic[substring_start-1]);

  } else {			/* sensep == false */
    substring_length = Substring_match_length(acceptor);
    *acceptor1 = toupper(complCode[(int) genomic[substring_start+substring_length]]);
    *acceptor2 = toupper(complCode[(int) genomic[substring_start+substring_length+1]]);
  }

  return;
}



static void
print_halfdonor (FILE *fp, char *abbrev, Substring_T donor, Stage3end_T this, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths, int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T concordant_chrpos, Chrpos_T donor_chrpos, Chrpos_T acceptor_chrpos, Chrpos_T mate_chrpos,
		 int clipdir, int hardclip_low, int hardclip_high, Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool use_hardclip_p, bool print_xt_p, char donor_strand, char acceptor_strand,
		 char *donor_chr, char *acceptor_chr, char donor1, char donor2, char acceptor2, char acceptor1,
		 double donor_prob, double acceptor_prob, bool circularp) {
  unsigned int flag = 0U;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring_start, substring_length;
  int transloc_hardclip_low, transloc_hardclip_high;
  bool plusp, printp;
  bool start_ambig, end_ambig;
  int amb_length_start, amb_length_end;
  int n, i;
  Univcoord_T *start_ambcoords, *end_ambcoords, splicecoord;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif


  querylength = Shortread_fulllength(queryseq);
  plusp = Substring_plusp(donor);

  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Substring_chrnum(donor),donor_chrpos,Substring_chrlength(donor),chromosome_iit);
  

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  /* sensedir = Substring_chimera_sensedir(donor); */
  sensep = Substring_chimera_sensep(donor);

  if (use_hardclip_p == true) {
    if (sensep == plusp) {
      transloc_hardclip_low = 0;
      if (plusp == true) {
	/* sensep true */
	transloc_hardclip_high = querylength - Substring_queryend(donor);

      } else {
	/* sensep false */
	transloc_hardclip_high = Substring_querystart(donor);
      }

    } else { /* sensep != Substring_plusp(donor) */
      transloc_hardclip_high = 0;
      if (plusp == true) {
	transloc_hardclip_low = Substring_querystart(donor);

      } else {
	transloc_hardclip_low = querylength - Substring_queryend(donor);
      }
    }

    if (transloc_hardclip_low > hardclip_low) {
      hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > hardclip_high) {
      hardclip_high = transloc_hardclip_high;
    }
  }


  if (sensep == plusp) {
    if (plusp == true) {
      /* sensep true */
      assert(Substring_chimera_pos(donor) == Substring_queryend(donor));
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    Substring_querystart(donor) + 
		    Substring_match_length(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }

    } else {
      /* sensep false */
      assert(Substring_chimera_pos(donor) == Substring_querystart(donor));
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    (querylength - Substring_queryend(donor)) +
		    Substring_match_length(donor),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
		      /*querypos*/Substring_querystart(donor),querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);

      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      }
    }

  } else { /* sensep != Substring_plusp(donor) */
    if (plusp == true) {
      assert(Substring_chimera_pos(donor) == Substring_querystart(donor));
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(donor) + 
		    (querylength - Substring_queryend(donor)),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      } else {
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_querystart(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }

    } else {
      assert(Substring_chimera_pos(donor) == Substring_queryend(donor));
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(donor) +
		    Substring_querystart(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      } else {
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		    /*querypos*/Substring_queryend(donor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		    /*querypos*/Substring_querystart(donor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true);
      }
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  /* For anchor_chrnum, previously used Stage3end_chrnum(this), but this is 0 */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     /*anchor_chrnum*/Substring_chrnum(donor),donor_chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH */
  if (hardclip_low > 0 || hardclip_high > 0) {
    fprintf(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\tMD:Z:");
  printp = false;

  if (hide_soft_clips_p == true) {
    substring_start = Substring_querystart_orig(donor);
    substring_length = Substring_match_length_orig(donor);
  } else {
    substring_start = Substring_querystart(donor);
    substring_length = Substring_match_length(donor);
  }

  if (use_hardclip_p == false) {
    genomicdir_refdiff = Substring_genomic_refdiff(donor);
    genomicdir_bothdiff = Substring_genomic_bothdiff(donor);
    if (plusp == true) {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicdir_refdiff[substring_start]),&(genomicdir_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }

  } else if (sensep == true) {
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(donor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(donor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(donor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(donor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
    }

  } else {			/* sensep == false */
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(donor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(donor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(donor);
      genomicdir_bothdiff = Substring_genomic_refdiff(donor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
    }
  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\tNH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  /* fprintf(fp,"\tNM:i:%d",Substring_nmismatches_refdiff(donor)); */
  fprintf(fp,"\tNM:i:%d",nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\tXW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\tSM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\tXQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
#if 0
  /* Not necessary to compute, because already computed by print_exon_exon */
  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      fprintf(fp,"\tXS:A:+");
    } else {
      fprintf(fp,"\tXS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      fprintf(fp,"\tXS:A:-");
    } else {
      fprintf(fp,"\tXS:A:+");
    }
  } else if (force_xs_direction_p == true) {
    fprintf(fp,"\tXS:A:+");
  } else {
    fprintf(fp,"\tXS:A:?");
  }
#else
  fprintf(fp,"\tXS:A:%c",donor_strand);
#endif



  /* 12. TAGS: XA */
  if ((start_ambig = Stage3end_start_ambiguous_p(this)) == true ||
      (end_ambig = Stage3end_end_ambiguous_p(this)) == true) {
    fprintf(fp,"\tXA:Z:");

    if (plusp == true) {
      if ((n = Stage3end_start_nambcoords(this)) > 0) {
	assert(sensep == false);
	start_ambcoords = Stage3end_start_ambcoords(this);
	splicecoord = Substring_alignstart(donor);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(donor);
	fprintf(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(donor);
	fprintf(fp,"%u",splicecoord - start_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",splicecoord - start_ambcoords[i]);
	}
#endif
      }
      fprintf(fp,"|");
      if ((n = Stage3end_end_nambcoords(this)) > 0) {
	assert(sensep == true);
	end_ambcoords = Stage3end_end_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(donor);
	fprintf(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(donor);
	fprintf(fp,"%u",end_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - splicecoord);
	}
#endif
      }

    } else {
      if ((n = Stage3end_end_nambcoords(this)) > 0) {
	assert(sensep == true);
	end_ambcoords = Stage3end_end_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(donor);
	fprintf(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(donor);
	fprintf(fp,"%u",splicecoord - end_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",splicecoord - end_ambcoords[i]);
	}
#endif
      }
      fprintf(fp,"|");
      if ((n = Stage3end_start_nambcoords(this)) > 0) {
	assert(sensep == false);
	start_ambcoords = Stage3end_start_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(donor);
	fprintf(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(donor);
	fprintf(fp,"%u",start_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",start_ambcoords[i] - splicecoord);
	}
#endif
      }
    }
  }

  /* 12. TAGS: XT */
  if (print_xt_p == true) {
    fprintf(fp,"\tXT:Z:%c%c-%c%c,%.2f,%.2f",donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
    fprintf(fp,",%c%s@%u..%c%s@%u",donor_strand,donor_chr,donor_chrpos,acceptor_strand,acceptor_chr,acceptor_chrpos);
  }

  /* 12. TAGS: XC */
  if (circularp == true) {
    fprintf(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\tXG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}


static bool
check_cigar_halfdonor (Substring_T donor, int querylength, int clipdir, int hardclip_low, int hardclip_high,
		       bool first_read_p, bool circularp) {
  bool result;
  Intlist_T cigar_types = NULL;
  bool plusp, sensep;
  bool use_hardclip_p = false;
  int transloc_hardclip_low, transloc_hardclip_high;

  plusp = Substring_plusp(donor);

  sensep = Substring_chimera_sensep(donor);

  if (use_hardclip_p == true) {
    if (sensep == plusp) {
      transloc_hardclip_low = 0;
      if (plusp == true) {
	/* sensep true */
	transloc_hardclip_high = querylength - Substring_queryend(donor);

      } else {
	/* sensep false */
	transloc_hardclip_high = Substring_querystart(donor);
      }

    } else { /* sensep != Substring_plusp(donor) */
      transloc_hardclip_high = 0;
      if (plusp == true) {
	transloc_hardclip_low = Substring_querystart(donor);

      } else {
	transloc_hardclip_low = querylength - Substring_queryend(donor);
      }
    }

    if (transloc_hardclip_low > hardclip_low) {
      hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > hardclip_high) {
      hardclip_high = transloc_hardclip_high;
    }
  }


  if (sensep == plusp) {
    if (plusp == true) {
      /* sensep true */
      assert(Substring_chimera_pos(donor) == Substring_queryend(donor));
      if (hide_soft_clips_p == true) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_querystart(donor) + 
					       Substring_match_length(donor),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
					       /*querypos*/Substring_queryend(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(donor),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(donor),
					       /*querypos*/Substring_querystart(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
					       /*querypos*/Substring_queryend(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }

    } else {
      /* sensep false */
      assert(Substring_chimera_pos(donor) == Substring_querystart(donor));
      if (hide_soft_clips_p == true) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       (querylength - Substring_queryend(donor)) +
					       Substring_match_length(donor),
					       /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
					       /*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
					       /*querypos*/Substring_querystart(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);

      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(donor),
					       /*querypos*/querylength,querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(donor),
					       /*querypos*/Substring_queryend(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
					       /*querypos*/Substring_querystart(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      }
    }

  } else { /* sensep != Substring_plusp(donor) */
    if (plusp == true) {
      assert(Substring_chimera_pos(donor) == Substring_querystart(donor));
      if (hide_soft_clips_p == true) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_match_length(donor) + 
					       (querylength - Substring_queryend(donor)),
					       /*querypos*/Substring_querystart(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(donor),
					       /*querypos*/Substring_querystart(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(donor),
					       /*querypos*/Substring_queryend(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }

    } else {
      assert(Substring_chimera_pos(donor) == Substring_queryend(donor));
      if (hide_soft_clips_p == true) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
					       /*querypos*/querylength,querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_match_length(donor) +
					       Substring_querystart(donor),
					       /*querypos*/Substring_queryend(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
					       /*querypos*/querylength,querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(donor),
					       /*querypos*/Substring_queryend(donor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(donor),
					       /*querypos*/Substring_querystart(donor),querylength,hardclip_low,hardclip_high,
					       /*plusp*/false,/*lastp*/true);
      }
    }
  }

  result = check_cigar_types(cigar_types);

  Intlist_free(&cigar_types);
  return result;
}



static void
print_halfacceptor (FILE *fp, char *abbrev, Substring_T acceptor, Stage3end_T this, Stage3end_T mate,
		    char *acc1, char *acc2, int pathnum, int npaths, int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		    Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		    Chrpos_T concordant_chrpos, Chrpos_T donor_chrpos, Chrpos_T acceptor_chrpos, Chrpos_T mate_chrpos,
		    int clipdir, int hardclip_low, int hardclip_high, Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		    int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		    bool use_hardclip_p, bool print_xt_p, char donor_strand, char acceptor_strand,
		    char *donor_chr, char *acceptor_chr, char donor1, char donor2, char acceptor2, char acceptor1,
		    double donor_prob, double acceptor_prob, bool circularp) {
  unsigned int flag = 0U;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring_start, substring_length;
  int transloc_hardclip_low, transloc_hardclip_high;
  bool plusp, printp;
  bool start_ambig, end_ambig;
  int amb_length_start, amb_length_end;
  int n, i;
  Univcoord_T *start_ambcoords, *end_ambcoords, splicecoord;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif


  querylength = Shortread_fulllength(queryseq);
  plusp = Substring_plusp(acceptor);

  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Substring_chrnum(acceptor),acceptor_chrpos,Substring_chrlength(acceptor),chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  /* sensedir = Substring_chimera_sensedir(acceptor); */
  sensep = Substring_chimera_sensep(acceptor);

  if (use_hardclip_p == true) {
    if (sensep != plusp) {
      transloc_hardclip_low = 0;
      if (plusp == true) {
	/* sensep false */
	transloc_hardclip_high = querylength - Substring_queryend(acceptor);

      } else {
	/* sensep true */
	transloc_hardclip_high = Substring_querystart(acceptor);
      }

    } else { /*  sensep == Substring_plusp(acceptor) */
      transloc_hardclip_high = 0;
      if (plusp == true) {
	transloc_hardclip_low = Substring_querystart(acceptor);

      } else {
	transloc_hardclip_low = querylength - Substring_queryend(acceptor);
      }
    }

    if (transloc_hardclip_low > hardclip_low) {
      hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > hardclip_high) {
      hardclip_high = transloc_hardclip_high;
    }
  }

  if (sensep != plusp) {
    if (plusp == true) {
      /* sensep false */
      assert(Substring_chimera_pos(acceptor) == Substring_queryend(acceptor));
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    Substring_querystart(acceptor) +
		    Substring_match_length(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      } else {
	print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }

    } else {
      /* sensep true */
      assert(Substring_chimera_pos(acceptor) == Substring_querystart(acceptor));
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/'M',
		    (querylength - Substring_queryend(acceptor)) +
		    Substring_match_length(acceptor),
		    /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true);
      } else {
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true);
      }
    }

  } else { /*  sensep == Substring_plusp(acceptor) */
    if (plusp == true) {
      assert(Substring_chimera_pos(acceptor) == Substring_querystart(acceptor));
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(acceptor) +
		    (querylength - Substring_queryend(acceptor)),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      } else {
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
		    /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }

    } else {
      assert(Substring_chimera_pos(acceptor) == Substring_queryend(acceptor));
      if (hide_soft_clips_p == true) {
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/'M',
		    Substring_match_length(acceptor) +
		    Substring_querystart(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      } else {
	print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
		    /*querypos*/querylength,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		    /*querypos*/Substring_queryend(acceptor),querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		    /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/true);
      }
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  /* For anchor_chrnum, previously used Stage3end_chrnum(this), but this is 0 */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     /*anchor_chrnum*/Substring_chrnum(acceptor),acceptor_chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH */
  if (hardclip_low > 0 || hardclip_high > 0) {
    fprintf(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\tMD:Z:");
  printp = false;

  if (hide_soft_clips_p == true) {
    substring_start = Substring_querystart_orig(acceptor);
    substring_length = Substring_match_length_orig(acceptor);
  } else {
    substring_start = Substring_querystart(acceptor);
    substring_length = Substring_match_length(acceptor);
  }

  if (use_hardclip_p == false) {
    genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
    genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
    if (plusp == true) {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicdir_refdiff[substring_start]),&(genomicdir_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((querylength+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }

  } else if (sensep == false) {
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(acceptor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(acceptor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }

    }

  } else {			/* sensep true */
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(acceptor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(acceptor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
    }
  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\tNH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  /* fprintf(fp,"\tNM:i:%d",Substring_nmismatches_refdiff(acceptor)); */
  fprintf(fp,"\tNM:i:%d",nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\tXW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\tSM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\tXQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
#if 0
  /* Not necessary to compute, because already computed by print_exon_exon */
  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      fprintf(fp,"\tXS:A:+");
    } else {
      fprintf(fp,"\tXS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      fprintf(fp,"\tXS:A:-");
    } else {
      fprintf(fp,"\tXS:A:+");
    }
  } else if (force_xs_direction_p == true) {
    fprintf(fp,"\tXS:A:+");
  } else {
    fprintf(fp,"\tXS:A:?");
  }
#else
  fprintf(fp,"\tXS:A:%c",acceptor_strand);
#endif

  /* 12. TAGS: XA */
  if ((start_ambig = Stage3end_start_ambiguous_p(this)) == true ||
      (end_ambig = Stage3end_end_ambiguous_p(this)) == true) {
    fprintf(fp,"\tXA:Z:");

    if (plusp == true) {
      if ((n = Stage3end_start_nambcoords(this)) > 0) {
	assert(sensep == true);
	start_ambcoords = Stage3end_start_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(acceptor);
	fprintf(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(acceptor);
	fprintf(fp,"%u",splicecoord - start_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",splicecoord - start_ambcoords[i]);
	}
#endif
      }
      fprintf(fp,"|");
      if ((n = Stage3end_end_nambcoords(this)) > 0) {
	assert(sensep == false);
	end_ambcoords = Stage3end_end_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(acceptor);
	fprintf(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(acceptor);
	fprintf(fp,"%u",end_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - splicecoord);
	}
#endif
      }

    } else {
      if ((n = Stage3end_end_nambcoords(this)) > 0) {
	assert(sensep == false);
	end_ambcoords = Stage3end_end_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(acceptor);
	fprintf(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(acceptor);
	fprintf(fp,"%u",splicecoord - end_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",splicecoord - end_ambcoords[i]);
	}
#endif
      }
      fprintf(fp,"|");
      if ((n = Stage3end_start_nambcoords(this)) > 0) {
	assert(sensep == true);
	start_ambcoords = Stage3end_start_ambcoords(this);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(acceptor);
	fprintf(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(acceptor);
	fprintf(fp,"%u",start_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",start_ambcoords[i] - splicecoord);
	}
#endif
      }
    }
  }

  /* 12. TAGS: XT */
  if (print_xt_p == true) {
    fprintf(fp,"\tXT:Z:%c%c-%c%c,%.2f,%.2f",donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
    fprintf(fp,",%c%s@%u..%c%s@%u",donor_strand,donor_chr,donor_chrpos,acceptor_strand,acceptor_chr,acceptor_chrpos);
  }

  /* 12. TAGS: XC */
  if (circularp == true) {
    fprintf(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\tXG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}


static bool
check_cigar_halfacceptor (Substring_T acceptor, int querylength, int clipdir, int hardclip_low, int hardclip_high,
			  bool first_read_p, bool circularp) {
  bool result;
  Intlist_T cigar_types = NULL;
  bool plusp, sensep;
  bool use_hardclip_p = false;
  int transloc_hardclip_low, transloc_hardclip_high;

  plusp = Substring_plusp(acceptor);

  sensep = Substring_chimera_sensep(acceptor);

  if (use_hardclip_p == true) {
    if (sensep != plusp) {
      transloc_hardclip_low = 0;
      if (plusp == true) {
	/* sensep false */
	transloc_hardclip_high = querylength - Substring_queryend(acceptor);

      } else {
	/* sensep true */
	transloc_hardclip_high = Substring_querystart(acceptor);
      }

    } else { /*  sensep == Substring_plusp(acceptor) */
      transloc_hardclip_high = 0;
      if (plusp == true) {
	transloc_hardclip_low = Substring_querystart(acceptor);

      } else {
	transloc_hardclip_low = querylength - Substring_queryend(acceptor);
      }
    }

    if (transloc_hardclip_low > hardclip_low) {
      hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > hardclip_high) {
      hardclip_high = transloc_hardclip_high;
    }
  }

  if (sensep != plusp) {
    if (plusp == true) {
      /* sensep false */
      assert(Substring_chimera_pos(acceptor) == Substring_queryend(acceptor));
      if (hide_soft_clips_p == true) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_querystart(acceptor) +
					       Substring_match_length(acceptor),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
					       /*querypos*/Substring_queryend(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(acceptor),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(acceptor),
					       /*querypos*/Substring_querystart(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
					       /*querypos*/Substring_queryend(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }

    } else {
      /* sensep true */
      assert(Substring_chimera_pos(acceptor) == Substring_querystart(acceptor));
      if (hide_soft_clips_p == true) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       (querylength - Substring_queryend(acceptor)) +
					       Substring_match_length(acceptor),
					       /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		    /*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
					       /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
					       /*plusp*/false,/*lastp*/true);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(acceptor),
					       /*querypos*/querylength,querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(acceptor),
					       /*querypos*/Substring_queryend(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
					       /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
					       /*plusp*/false,/*lastp*/true);
      }
    }

  } else { /*  sensep == Substring_plusp(acceptor) */
    if (plusp == true) {
      assert(Substring_chimera_pos(acceptor) == Substring_querystart(acceptor));
      if (hide_soft_clips_p == true) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_match_length(acceptor) +
					       (querylength - Substring_queryend(acceptor)),
					       /*querypos*/Substring_querystart(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
					       /*querypos*/0,querylength,hardclip_low,hardclip_high,
					       /*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(acceptor),
					       /*querypos*/Substring_querystart(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(acceptor),
					       /*querypos*/Substring_queryend(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
      }

    } else {
      assert(Substring_chimera_pos(acceptor) == Substring_queryend(acceptor));
      if (hide_soft_clips_p == true) {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
					       /*querypos*/querylength,querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					       Substring_match_length(acceptor) +
					       Substring_querystart(acceptor),
					       /*querypos*/Substring_queryend(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      } else {
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
					       /*querypos*/querylength,querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(acceptor),
					       /*querypos*/Substring_queryend(acceptor),querylength,
					       hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(acceptor),
					       /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
					       /*plusp*/false,/*lastp*/true);
      }
    }
  }

  result = check_cigar_types(cigar_types);

  Intlist_free(&cigar_types);
  return result;
}


static void
print_localsplice (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		   char *acc1, char *acc2, int pathnum, int npaths,
		   int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		   Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		   Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip_low, int hardclip_high,
		   Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		   int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		   bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  int sensedir;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring1_start, substring2_start, substring1_length, substring2_length, matchlength;
  bool plusp, printp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  if ((sensedir = Stage3end_sensedir(this)) == SENSE_NULL) {
    sensedir = Stage3end_sensedir(mate);
  }
  sensep = (sensedir == SENSE_FORWARD);

  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,Stage3end_chrlength(this),chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  if (sensep == plusp) {
    substring1 = /* donor */ Stage3end_substring_donor(this);
    substring2 = /* acceptor */ Stage3end_substring_acceptor(this);
  } else {
    substring1 = /* acceptor */ Stage3end_substring_acceptor(this);
    substring2 = /* donor */ Stage3end_substring_donor(this);
  }

  if (plusp == true) {
    if (hide_soft_clips_p == true) {
      print_cigar(fp,/*type*/'M',
		  Substring_querystart(substring1) +
		  Substring_match_length(substring1),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      if (hardclip_low < Substring_queryend(substring1) &&
	  querylength - hardclip_high > Substring_querystart(substring2)) {
	debug1(printf("\ncase 1: hardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
		      hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substring2)));
	fprintf(fp,"%uN",Stage3end_distance(this));
      }
      print_cigar(fp,/*type*/'M',
		  Substring_match_length(substring2) +
		  (querylength - Substring_queryend(substring2)),
		  /*querypos*/Substring_querystart(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		  /*querypos*/Substring_querystart(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      if (hardclip_low < Substring_queryend(substring1) &&
	  querylength - hardclip_high > Substring_querystart(substring2)) {
	debug1(printf("\ncase 1: hardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
		      hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substring2)));
	fprintf(fp,"%uN",Stage3end_distance(this));
      }
      print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		  /*querypos*/Substring_querystart(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		  /*querypos*/Substring_queryend(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (hide_soft_clips_p == true) {
      print_cigar(fp,/*type*/'M',
		  (querylength - Substring_queryend(substring1)) +
		  Substring_match_length(substring1),
		  /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/false);
      if (querylength - hardclip_low > Substring_queryend(substring2) &&
	  hardclip_high < Substring_querystart(substring1)) {
	debug1(printf("\ncase 2: querylength %d - hardclip_low %d > queryend(substring2) %d && hardclip_high %d < querystart(substring1) %d\n",
		      querylength,hardclip_low,Substring_queryend(substring2),hardclip_high,Substring_querystart(substring1)));
	fprintf(fp,"%uN",Stage3end_distance(this));
      }
      print_cigar(fp,/*type*/'M',
		  Substring_match_length(substring2) +
		  Substring_querystart(substring2),
		  /*querypos*/Substring_querystart(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring1),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		  /*querypos*/Substring_queryend(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      if (querylength - hardclip_low > Substring_queryend(substring2) &&
	  hardclip_high < Substring_querystart(substring1)) {
	debug1(printf("\ncase 2: querylength %d - hardclip_low %d > queryend(substring2) %d && hardclip_high %d < querystart(substring1) %d\n",
		      querylength,hardclip_low,Substring_queryend(substring2),hardclip_high,Substring_querystart(substring1)));
	fprintf(fp,"%uN",Stage3end_distance(this));
      }
      print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		  /*querypos*/Substring_querystart(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'S',Substring_querystart(substring2),
		  /*querypos*/Substring_querystart(substring2),querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/true);
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
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


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			   quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				   quality_shift,/*show_chopped_p*/false);
  } 

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH */
  if (hardclip_low > 0 || hardclip_high > 0) {
    fprintf(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\tMD:Z:");
  printp = false;

  if (hide_soft_clips_p == true) {
    substring1_start = Substring_querystart_orig(substring1);
    substring1_length = Substring_match_length_orig(substring1);
    substring2_start = Substring_querystart_orig(substring2);
    substring2_length = Substring_match_length_orig(substring2);
  } else {
    substring1_start = Substring_querystart(substring1);
    substring1_length = Substring_match_length(substring1);
    substring2_start = Substring_querystart(substring2);
    substring2_length = Substring_match_length(substring2);
  }

  if (plusp == true) {
    genomicfwd_refdiff = Substring_genomic_refdiff(substring1);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring1);
    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
				  &(genomicfwd_refdiff[substring1_start]),&(genomicfwd_bothdiff[substring1_start]),
				  substring1_length,/*querypos*/substring1_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    
#if 0
    /* Intron: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd_refdiff = Substring_genomic_refdiff(substring2);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring2);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
		    &(genomicfwd_refdiff[substring2_start]),&(genomicfwd_bothdiff[substring2_start]),
		    substring2_length,/*querypos*/substring2_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    genomicdir_refdiff = Substring_genomic_refdiff(substring1);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring1);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				    substring1_length,/*querypos*/substring1_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring1_start]),substring1_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
				    substring1_length,/*querypos*/substring1_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }

#if 0
    /* Intron: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicdir_refdiff = Substring_genomic_refdiff(substring2);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring2);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring2_start]),substring2_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }
  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\tNH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  /* fprintf(fp,"\tNM:i:%d",Stage3end_nmismatches_refdiff(this)); */
  fprintf(fp,"\tNM:i:%d",nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\tXW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\tSM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\tXQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      fprintf(fp,"\tXS:A:+");
    } else {
      fprintf(fp,"\tXS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      fprintf(fp,"\tXS:A:-");
    } else {
      fprintf(fp,"\tXS:A:+");
    }
  } else if (force_xs_direction_p == true) {
    fprintf(fp,"\tXS:A:+");
  } else {
    fprintf(fp,"\tXS:A:?");
  }

  /* 12. TAGS: XA */
  assert(Stage3end_start_ambiguous_p(this) == false);
  assert(Stage3end_end_ambiguous_p(this) == false);

  /* 12. TAGS: XC */
  if (circularp == true) {
    fprintf(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\tXG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}


static bool
check_cigar_localsplice (Stage3end_T this, Stage3end_T mate, int querylength, int clipdir, int hardclip_low, int hardclip_high,
			 bool first_read_p, bool circularp) {
  bool result;
  Intlist_T cigar_types = NULL;
  Substring_T substring1, substring2;
  bool plusp, sensep;
  int sensedir;

  plusp = Stage3end_plusp(this);

  if ((sensedir = Stage3end_sensedir(this)) == SENSE_NULL) {
    sensedir = Stage3end_sensedir(mate);
  }
  sensep = (sensedir == SENSE_FORWARD);

  if (sensep == plusp) {
    substring1 = /* donor */ Stage3end_substring_donor(this);
    substring2 = /* acceptor */ Stage3end_substring_acceptor(this);
  } else {
    substring1 = /* acceptor */ Stage3end_substring_acceptor(this);
    substring2 = /* donor */ Stage3end_substring_donor(this);
  }

  if (plusp == true) {
    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_querystart(substring1) +
					     Substring_match_length(substring1),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
      if (hardclip_low < Substring_queryend(substring1) &&
	  querylength - hardclip_high > Substring_querystart(substring2)) {
	debug1(printf("\ncase 1: hardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
		      hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substring2)));
	cigar_types = Intlist_push(cigar_types,'N');
      }
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_match_length(substring2) +
					     (querylength - Substring_queryend(substring2)),
					     /*querypos*/Substring_querystart(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring1),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1),
					     /*querypos*/Substring_querystart(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      if (hardclip_low < Substring_queryend(substring1) &&
	  querylength - hardclip_high > Substring_querystart(substring2)) {
	debug1(printf("\ncase 1: hardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
		      hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substring2)));
	cigar_types = Intlist_push(cigar_types,'N');
      }
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2),
					     /*querypos*/Substring_querystart(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring2),
					     /*querypos*/Substring_queryend(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     (querylength - Substring_queryend(substring1)) +
					     Substring_match_length(substring1),
					     /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
					     /*plusp*/false,/*lastp*/false);
      if (querylength - hardclip_low > Substring_queryend(substring2) &&
	  hardclip_high < Substring_querystart(substring1)) {
	debug1(printf("\ncase 2: querylength %d - hardclip_low %d > queryend(substring2) %d && hardclip_high %d < querystart(substring1) %d\n",
		      querylength,hardclip_low,Substring_queryend(substring2),hardclip_high,Substring_querystart(substring1)));
	cigar_types = Intlist_push(cigar_types,'N');
      }
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
		  Substring_match_length(substring2) +
		  Substring_querystart(substring2),
		  /*querypos*/Substring_querystart(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring1),
					     /*querypos*/querylength,querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1),
					     /*querypos*/Substring_queryend(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      if (querylength - hardclip_low > Substring_queryend(substring2) &&
	  hardclip_high < Substring_querystart(substring1)) {
	debug1(printf("\ncase 2: querylength %d - hardclip_low %d > queryend(substring2) %d && hardclip_high %d < querystart(substring1) %d\n",
		      querylength,hardclip_low,Substring_queryend(substring2),hardclip_high,Substring_querystart(substring1)));
	cigar_types = Intlist_push(cigar_types,'N');
      }
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2),
					     /*querypos*/Substring_querystart(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring2),
					     /*querypos*/Substring_querystart(substring2),querylength,hardclip_low,hardclip_high,
					     /*plusp*/false,/*lastp*/true);
    }
  }

  result = check_cigar_types(cigar_types);

  Intlist_free(&cigar_types);
  return result;
}


static void
print_shortexon (FILE *fp, char *abbrev, Stage3end_T shortexon, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths,
		 int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip_low, int hardclip_high,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool circularp) {
  unsigned int flag = 0U;
  /* substring1 is low coordinate on genome, substring2 is high */
  Substring_T substring1, substring2, substringM;
  Chrpos_T distance1, distance2;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  int sensedir;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring1_start, substring2_start, substringM_start,
    substring1_length, substring2_length, substringM_length, matchlength;
  bool plusp, printp;
  bool start_ambig, end_ambig;
  int amb_length_start, amb_length_end;
  int n, i;
  Univcoord_T *start_ambcoords, *end_ambcoords, splicecoord;
#ifdef PRINT_AMBIG_COORDS
  Univcoord_T chroffset;
#endif

  
  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(shortexon);

  if ((sensedir = Stage3end_sensedir(shortexon)) == SENSE_NULL) {
    sensedir = Stage3end_sensedir(mate);
  }
  sensep = (sensedir == SENSE_FORWARD);

  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(shortexon),chrpos,Stage3end_chrlength(shortexon),chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  substringM = Stage3end_substring1(shortexon);

  if (sensep == plusp) {
    substring1 = /* donor */ Stage3end_substringD(shortexon);
    distance1 = Stage3end_shortexonA_distance(shortexon);
    distance2 = Stage3end_shortexonD_distance(shortexon);
    substring2 = /* acceptor */ Stage3end_substringA(shortexon);
  } else {
    substring1 = /* acceptor */ Stage3end_substringA(shortexon);
    distance1 = Stage3end_shortexonD_distance(shortexon);
    distance2 = Stage3end_shortexonA_distance(shortexon);
    substring2 = /* donor */ Stage3end_substringD(shortexon);
  }

  if (substring1 == NULL) {
    if (plusp == true) {
      print_cigar(fp,/*type*/'S',Substring_querystart(substringM),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
    } else {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substringM),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }

  } else if (plusp == true) {
    if (hide_soft_clips_p == true) {
      print_cigar(fp,/*type*/'M',
		  Substring_querystart(substring1) +
		  Substring_match_length(substring1),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
    } else {
      print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		  /*querypos*/Substring_querystart(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    }
    if (hardclip_low < Substring_queryend(substring1) &&
	querylength - hardclip_high > Substring_querystart(substringM)) {
      debug1(printf("\ncase 3: hardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substringM) %d\n",
		    hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substringM)));
      fprintf(fp,"%uN",distance1);
    }

  } else {
    if (hide_soft_clips_p == true) {
      print_cigar(fp,/*type*/'M',
		  (querylength - Substring_queryend(substring1)) +
		  Substring_match_length(substring1),
		  /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/false);
    } else {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring1),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		  /*querypos*/Substring_queryend(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }
    if (querylength - hardclip_low > Substring_queryend(substringM) &&
	hardclip_high < Substring_querystart(substring1)) {
      debug1(printf("\ncase 4: querylength %d - hardclip_low %d > queryend(substringM) %d && hardclip_high %d < querystart(substring1) %d\n",
		    querylength,hardclip_low,Substring_queryend(substringM),hardclip_high,Substring_querystart(substring1)));
      fprintf(fp,"%uN",distance1);
    }
  }

  if (plusp == true) {
    print_cigar(fp,/*type*/'M',Substring_match_length(substringM),
		/*querypos*/Substring_querystart(substringM),querylength,
		hardclip_low,hardclip_high,plusp,/*lastp*/false);
  } else {
    print_cigar(fp,/*type*/'M',Substring_match_length(substringM),
		/*querypos*/Substring_queryend(substringM),querylength,
		hardclip_low,hardclip_high,plusp,/*lastp*/false);
  }

  if (substring2 == NULL) {
    if (plusp == true) {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substringM),
		  /*querypos*/Substring_queryend(substringM),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'S',Substring_querystart(substringM),
		  /*querypos*/Substring_querystart(substringM),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }

  } else if (plusp == true) {
    if (hardclip_low < Substring_queryend(substringM) &&
	querylength - hardclip_high > Substring_querystart(substring2)) {
      debug1(printf("\ncase 5: hardclip_low %d < queryend(substringM) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
		    hardclip_low,Substring_queryend(substringM),querylength,hardclip_high,Substring_querystart(substring2)));
      fprintf(fp,"%uN",distance2);
    }
    if (hide_soft_clips_p == true) {
      print_cigar(fp,/*type*/'M',
		  Substring_match_length(substring2) +
		  (querylength - Substring_queryend(substring2)),
		  /*querypos*/Substring_querystart(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		  /*querypos*/Substring_querystart(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		  /*querypos*/Substring_queryend(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (querylength - hardclip_low > Substring_queryend(substring2) &&
	hardclip_high < Substring_querystart(substringM)) {
      debug1(printf("\ncase 6: querylength %d - hardclip_low %d > queryend(substring2) %d && hardclip_high %d < querystart(substringM) %d\n",
		    querylength,hardclip_low,Substring_queryend(substring2),querylength,Substring_querystart(substringM)));
      fprintf(fp,"%uN",distance2);
    }
    if (hide_soft_clips_p == true) {
      print_cigar(fp,/*type*/'M',
		  Substring_match_length(substring2) +
		  Substring_querystart(substring2),
		  /*querypos*/Substring_queryend(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		  /*querypos*/Substring_queryend(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'S',Substring_querystart(substring2),
		  /*querypos*/Substring_querystart(substring2),querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/true);
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
			     mate_chrpos,Stage3end_chrlength(mate),
			     Stage3end_chrnum(shortexon),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
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


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			   quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				   quality_shift,/*show_chopped_p*/false);
  } 

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH */
  if (hardclip_low > 0 || hardclip_high > 0) {
    fprintf(fp,"\tXH:Z:");
    if (plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\tMD:Z:");
  printp = false;

  if (hide_soft_clips_p == true) {
    substringM_start = Substring_querystart_orig(substringM);
    substringM_length = Substring_match_length_orig(substringM);
  } else {
    substringM_start = Substring_querystart(substringM);
    substringM_length = Substring_match_length(substringM);
  }

  if (substring1 == NULL) {
    substring1_start = 0;
    substring1_length = 0;
  } else if (hide_soft_clips_p == true) {
    substring1_start = Substring_querystart_orig(substring1);
    substring1_length = Substring_match_length_orig(substring1);
  } else {
    substring1_start = Substring_querystart(substring1);
    substring1_length = Substring_match_length(substring1);
  }
  if (substring2 == NULL) {
    substring2_start = 0;
    substring2_length = 0;
  } else if (hide_soft_clips_p == true) {
    substring2_start = Substring_querystart_orig(substring2);
    substring2_length = Substring_match_length_orig(substring2);
  } else {
    substring2_start = Substring_querystart(substring2);
    substring2_length = Substring_match_length(substring2);
  }

  if (plusp == true) {

    if (substring1 == NULL) {
      matchlength = 0;
    } else {
      genomicfwd_refdiff = Substring_genomic_refdiff(substring1);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(substring1);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
				    &(genomicfwd_refdiff[substring1_start]),&(genomicfwd_bothdiff[substring1_start]),
				    substring1_length,/*querypos*/substring1_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    }

#if 0
    /* Intron 1: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd_refdiff = Substring_genomic_refdiff(substringM);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substringM);
    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
				  &(genomicfwd_refdiff[substringM_start]),&(genomicfwd_bothdiff[substringM_start]),
				  substringM_length,/*querypos*/substringM_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);

#if 0
    /* Intron 2: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    if (substring2 == NULL) {
      /* Equivalent: if (matchlength > 0) fprintf(fp,"%d",matchlength); */
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
		      /*substring2_length*/0,/*querypos*/0,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd_refdiff = Substring_genomic_refdiff(substring2);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(substring2);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
		      &(genomicfwd_refdiff[substring2_start]),&(genomicfwd_bothdiff[substring2_start]),
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {

    if (substring1 == NULL) {
      matchlength = 0;
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(substring1);
      genomicdir_bothdiff = Substring_genomic_bothdiff(substring1);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				      substring1_length,/*querypos*/substring1_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring1_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring1_start]),substring1_length);
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
				      substring1_length,/*querypos*/substring1_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
    }

#if 0
    /* Intron 1: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicdir_refdiff = Substring_genomic_refdiff(substringM);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substringM);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) MALLOCA((substringM_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substringM_start]),substringM_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				    substringM_length,/*querypos*/substringM_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREEA(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) MALLOCA((substringM_length+1) * sizeof(char));
      genomicfwd_bothdiff = (char *) MALLOCA((substringM_length+1) * sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substringM_start]),substringM_length);
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substringM_start]),substringM_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
				    substringM_length,/*querypos*/substringM_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREEA(genomicfwd_bothdiff);
      FREEA(genomicfwd_refdiff);
    }

#if 0
    /* Intron 2: Not sure how to handle in MD string */
    fprintf(fp,"^");
#endif

    if (substring2 == NULL) {
      /* Equivalent: if (matchlength > 0) fprintf(fp,"%d",matchlength); */
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
		      /*substring2_length*/0,/*querypos*/0,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(substring2);
      genomicdir_bothdiff = Substring_genomic_bothdiff(substring2);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring2_length,/*querypos*/substring2_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
	genomicfwd_bothdiff = (char *) MALLOCA((substring2_length+1) * sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring2_start]),substring2_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring2_length,/*querypos*/substring2_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREEA(genomicfwd_bothdiff);
	FREEA(genomicfwd_refdiff);
      }
    }
  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\tNH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\tHI:i:%d",pathnum);

  /* 12. TAGS: NM */
  /* fprintf(fp,"\tNM:i:%d",Stage3end_nmismatches_refdiff(shortexon)); */
  fprintf(fp,"\tNM:i:%d",nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\tXW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\tSM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\tXQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\tX2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      fprintf(fp,"\tXS:A:+");
    } else {
      fprintf(fp,"\tXS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      fprintf(fp,"\tXS:A:-");
    } else {
      fprintf(fp,"\tXS:A:+");
    }
  } else if (force_xs_direction_p == true) {
    fprintf(fp,"\tXS:A:+");
  } else {
    fprintf(fp,"\tXS:A:?");
  }

  /* 12. TAGS: XA */
  if ((start_ambig = Stage3end_start_ambiguous_p(shortexon)) == true ||
      (end_ambig = Stage3end_end_ambiguous_p(shortexon)) == true) {
    fprintf(fp,"\tXA:Z:");

    if (plusp == true) {
      if ((n = Stage3end_start_nambcoords(shortexon)) > 0) {
	start_ambcoords = Stage3end_start_ambcoords(shortexon);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(substringM);
	fprintf(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(substringM);
	fprintf(fp,"%u",splicecoord - start_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",splicecoord - start_ambcoords[i]);
	}
#endif
      }
      fprintf(fp,"|");
      if ((n = Stage3end_end_nambcoords(shortexon)) > 0) {
	end_ambcoords = Stage3end_end_ambcoords(shortexon);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(substringM);
	fprintf(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(substringM);
	fprintf(fp,"%u",end_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - splicecoord);
	}
#endif
      }

    } else {
      if ((n = Stage3end_end_nambcoords(shortexon)) > 0) {
	end_ambcoords = Stage3end_end_ambcoords(shortexon);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(substringM);
	fprintf(fp,"%u",end_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",end_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignend(substringM);
	fprintf(fp,"%u",splicecoord - end_ambcoords[0]);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",splicecoord - end_ambcoords[i]);
	}
#endif
      }
      fprintf(fp,"|");
      if ((n = Stage3end_start_nambcoords(shortexon)) > 0) {
	start_ambcoords = Stage3end_start_ambcoords(shortexon);
#ifdef PRINT_AMBIG_COORDS
	chroffset = Substring_chroffset(substringM);
	fprintf(fp,"%u",start_ambcoords[0] - chroffset + 1U);
	for (i = 1; i < n; i++) {
  	  fprintf(fp,",%u",start_ambcoords[i] - chroffset + 1U);
	}
#else
	splicecoord = Substring_alignstart(substringM);
	fprintf(fp,"%u",start_ambcoords[0] - splicecoord);
	for (i = 1; i < n; i++) {
	  fprintf(fp,",%u",start_ambcoords[i] - splicecoord);
	}
#endif
      }
    }
  }

  /* 12. TAGS: XC */
  if (circularp == true) {
    fprintf(fp,"\tXC:A:+");
  }

  /* 12. TAGS: XG */
  if (Stage3end_sarrayp(shortexon) == true) {
    fprintf(fp,"\tXG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}


static bool
check_cigar_shortexon (Stage3end_T shortexon, Stage3end_T mate, int querylength, int clipdir, int hardclip_low, int hardclip_high,
		       bool first_read_p, bool circularp) {
  bool result;
  Intlist_T cigar_types = NULL;
  Substring_T substring1, substring2, substringM;
  bool plusp, sensep;
  int sensedir;

  plusp = Stage3end_plusp(shortexon);

  if ((sensedir = Stage3end_sensedir(shortexon)) == SENSE_NULL) {
    sensedir = Stage3end_sensedir(mate);
  }
  sensep = (sensedir == SENSE_FORWARD);

  substringM = Stage3end_substring1(shortexon);

  if (sensep == plusp) {
    substring1 = /* donor */ Stage3end_substringD(shortexon);
    substring2 = /* acceptor */ Stage3end_substringA(shortexon);
  } else {
    substring1 = /* acceptor */ Stage3end_substringA(shortexon);
    substring2 = /* donor */ Stage3end_substringD(shortexon);
  }

  if (substring1 == NULL) {
    if (plusp == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substringM),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substringM),
					     /*querypos*/querylength,querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }

  } else if (plusp == true) {
    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_querystart(substring1) +
					     Substring_match_length(substring1),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring1),
					     /*querypos*/0,querylength,hardclip_low,hardclip_high,
					     /*plusp*/true,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1),
					     /*querypos*/Substring_querystart(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    }
    if (hardclip_low < Substring_queryend(substring1) &&
	querylength - hardclip_high > Substring_querystart(substringM)) {
      debug1(printf("\ncase 3: hardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substringM) %d\n",
		    hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substringM)));
      cigar_types = Intlist_push(cigar_types,'N');
    }

  } else {
    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     (querylength - Substring_queryend(substring1)) +
					     Substring_match_length(substring1),
					     /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
					     /*plusp*/false,/*lastp*/false);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring1),
					     /*querypos*/querylength,querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring1),
					     /*querypos*/Substring_queryend(substring1),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }
    if (querylength - hardclip_low > Substring_queryend(substringM) &&
	hardclip_high < Substring_querystart(substring1)) {
      debug1(printf("\ncase 4: querylength %d - hardclip_low %d > queryend(substringM) %d && hardclip_high %d < querystart(substring1) %d\n",
		    querylength,hardclip_low,Substring_queryend(substringM),hardclip_high,Substring_querystart(substring1)));
      cigar_types = Intlist_push(cigar_types,'N');
    }
  }

  if (plusp == true) {
    cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substringM),
					   /*querypos*/Substring_querystart(substringM),querylength,
					   hardclip_low,hardclip_high,plusp,/*lastp*/false);
  } else {
    cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substringM),
					   /*querypos*/Substring_queryend(substringM),querylength,
					   hardclip_low,hardclip_high,plusp,/*lastp*/false);
  }

  if (substring2 == NULL) {
    if (plusp == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substringM),
					     /*querypos*/Substring_queryend(substringM),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substringM),
					     /*querypos*/Substring_querystart(substringM),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }

  } else if (plusp == true) {
    if (hardclip_low < Substring_queryend(substringM) &&
	querylength - hardclip_high > Substring_querystart(substring2)) {
      debug1(printf("\ncase 5: hardclip_low %d < queryend(substringM) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
		    hardclip_low,Substring_queryend(substringM),querylength,hardclip_high,Substring_querystart(substring2)));
      cigar_types = Intlist_push(cigar_types,'N');
    }
    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_match_length(substring2) +
					     (querylength - Substring_queryend(substring2)),
					     /*querypos*/Substring_querystart(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2),
					     /*querypos*/Substring_querystart(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',querylength - Substring_queryend(substring2),
					     /*querypos*/Substring_queryend(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {
    if (querylength - hardclip_low > Substring_queryend(substring2) &&
	hardclip_high < Substring_querystart(substringM)) {
      debug1(printf("\ncase 6: querylength %d - hardclip_low %d > queryend(substring2) %d && hardclip_high %d < querystart(substringM) %d\n",
		    querylength,hardclip_low,Substring_queryend(substring2),querylength,Substring_querystart(substringM)));
      cigar_types = Intlist_push(cigar_types,'N');
    }
    if (hide_soft_clips_p == true) {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',
					     Substring_match_length(substring2) +
					     Substring_querystart(substring2),
					     /*querypos*/Substring_queryend(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    } else {
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'M',Substring_match_length(substring2),
					     /*querypos*/Substring_queryend(substring2),querylength,
					     hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      cigar_types = compute_cigar_types_only(cigar_types,/*type*/'S',Substring_querystart(substring2),
					     /*querypos*/Substring_querystart(substring2),querylength,hardclip_low,hardclip_high,
					     /*plusp*/false,/*lastp*/true);
    }
  }

  result = check_cigar_types(cigar_types);

  Intlist_free(&cigar_types);
  return result;
}



/* Distant splicing, including scramble, inversion, translocation */
static void
print_exon_exon (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths,
		 int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T mate_chrpos, int clipdir, int hardclip_low, int hardclip_high,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  Chrpos_T donor_chrpos, acceptor_chrpos, concordant_chrpos;
  Substring_T donor, acceptor;
  char *donor_chr, *acceptor_chr;
  char donor1, donor2, acceptor2, acceptor1;
  double donor_prob, acceptor_prob;
  int circularpos, querylength;
  char donor_strand, acceptor_strand;
  int sensedir;
  bool allocp;

  donor = Stage3end_substring_donor(this);
  acceptor = Stage3end_substring_acceptor(this);

#if 0
  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
  }
#else
  /* Shouldn't have any overlap on a distant splice */
  hardclip_low = hardclip_high = 0;
#endif

  querylength = Shortread_fulllength(queryseq);
  donor_chrpos = SAM_compute_chrpos(hardclip_low,hardclip_high,this,querylength);
  acceptor_chrpos = SAM_compute_chrpos(hardclip_low,hardclip_high,this,querylength);
  if (Stage3end_substring_low(this) == donor) {
    concordant_chrpos = donor_chrpos;
  } else if (Stage3end_substring_low(this) == acceptor) {
    concordant_chrpos = acceptor_chrpos;
  } else {
    fprintf(stderr,"Stage3end_substring_low %p is neither donor %p or acceptor %p\n",
	    Stage3end_substring_low(this),donor,acceptor);
    concordant_chrpos = 0U;
  }

  halfdonor_dinucleotide(&donor1,&donor2,donor);
  halfacceptor_dinucleotide(&acceptor2,&acceptor1,acceptor);
  donor_chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(donor),&allocp);
  acceptor_chr = Univ_IIT_label(chromosome_iit,Substring_chrnum(acceptor),&allocp);
  donor_prob = Substring_chimera_prob(donor);
  acceptor_prob = Substring_chimera_prob(acceptor);

  /* Code taken from that for XS tag for print_halfdonor and print_halfacceptor */
  if ((sensedir = Substring_chimera_sensedir(donor)) == SENSE_FORWARD) {
    if (Substring_plusp(donor) == true) {
      donor_strand = '+';
    } else {
      donor_strand = '-';
    }
  } else if (sensedir == SENSE_ANTI) {
    if (Substring_plusp(donor) == true) {
      donor_strand = '-';
    } else {
      donor_strand = '+';
    }
  } else if (force_xs_direction_p == true) {
    donor_strand = '+';
  } else {
    donor_strand = '?';
  }

  /* Code taken from that for XS tag for print_halfdonor and print_halfacceptor */
  if ((sensedir = Substring_chimera_sensedir(acceptor)) == SENSE_FORWARD) {
    if (Substring_plusp(acceptor) == true) {
      acceptor_strand = '+';
    } else {
      acceptor_strand = '-';
    }
  } else if (sensedir == SENSE_ANTI) {
    if (Substring_plusp(acceptor) == true) {
      acceptor_strand = '-';
    } else {
      acceptor_strand = '+';
    }
  } else if (force_xs_direction_p == true) {
    acceptor_strand = '+';
  } else {
    acceptor_strand = '?';
  }


  if (Stage3end_sensedir(this) == SENSE_FORWARD) {

    /* NEEDS WORK: Need to decide whether to split halfdonor or halfacceptor */
    /* Not sure if circular chromosomes should participate in distant splicing anyway */
    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*donor_chrpos*/1,acceptor_chrpos,mate_chrpos,
		      /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
    } else {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      clipdir,hardclip_low,hardclip_high,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/false);
    }

    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,donor_chrpos,/*acceptor_chrpos*/1,mate_chrpos,
			 /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
    } else {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 clipdir,hardclip_low,hardclip_high,resulttype,first_read_p,
			 npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/false);
    }

  } else if (Stage3end_sensedir(this) == SENSE_ANTI) {
    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,donor_chrpos,/*acceptor_chrpos*/1,mate_chrpos,
			 /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
    } else {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
			 clipdir,hardclip_low,hardclip_high,resulttype,first_read_p,
			 npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor_strand,acceptor_strand,donor_chr,acceptor_chr,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/false);
    }

    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*donor_chrpos*/1,acceptor_chrpos,mate_chrpos,
		      /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
    } else {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,acceptor_chrpos,mate_chrpos,
		      clipdir,hardclip_low,hardclip_high,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor_strand,acceptor_strand,donor_chr,acceptor_chr,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/false);
    }

  } else {
    abort();
  }

  if (allocp == true) {
    FREE(acceptor_chr);
    FREE(donor_chr);
  }

  return;
}



void
SAM_print (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
	   char *acc1, char *acc2, int pathnum, int npaths,
	   int absmq_score, int first_absmq, int second_absmq, int mapq_score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq_mate, int pairedlength, Chrpos_T chrpos, Chrpos_T mate_chrpos,
	   int clipdir, int hardclip5_low, int hardclip5_high, int hardclip3_low, int hardclip3_high,
	   Resulttype_T resulttype, bool first_read_p,
	   int npaths_mate, int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
	   bool merge_samechr_p) {
  Hittype_T hittype;
  Substring_T donor, acceptor;
  bool sensep, normalp;
  unsigned int flag;
  int circularpos, querylength;
  char donor_strand, acceptor_strand;
  int sensedir;


  hittype = Stage3end_hittype(this);
  /* printf("hittype %s, chrpos %u\n",Stage3end_hittype_string(this),chrpos); */

  /* Test for nomapping was chrpos == 0, but we can actually align to chrpos 0 */
  /* Also, can use this test here because --quiet-if-excessive cases go directly to SAM_print_nomapping */
  if (npaths == 0) {
    SAM_print_nomapping(fp,abbrev,queryseq,mate,acc1,acc2,chromosome_iit,resulttype,first_read_p,
			/*npaths*/0,npaths_mate,mate_chrpos,quality_shift,
			sam_read_group_id,invertp,invert_mate_p);

    if (failedinput_root != NULL) {
      if (fastq_format_p == true) {
	if (first_read_p == true) {
	  Shortread_print_query_singleend_fastq(fp_failedinput_1,queryseq,/*headerseq*/queryseq);
	} else {
	  Shortread_print_query_singleend_fastq(fp_failedinput_1,queryseq,/*headerseq*/queryseq_mate);
	}
      } else {
	if (first_read_p == true) {
	  Shortread_print_query_singleend_fasta(fp_failedinput_2,queryseq,/*headerseq*/queryseq);
	} else {
	  Shortread_print_query_singleend_fasta(fp_failedinput_2,queryseq,/*headerseq*/queryseq_mate);
	}
      }
    }

  } else if (hittype == EXACT || hittype == SUB || hittype == TERMINAL) {
    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0 &&
	check_cigar_single(hittype,this,querylength,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			   first_read_p,/*circularp*/true) == true &&
	check_cigar_single(hittype,this,querylength,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			   first_read_p,/*circularp*/true) == true) {
#ifdef CHECK_ASSERTIONS
      if (Stage3end_plusp(this) == true) {
	assert(chrpos-Stage3end_trim_left(this)+circularpos-Stage3end_chrlength(this) == 1);
      } else {
	assert(chrpos-Stage3end_trim_right(this)+circularpos-Stage3end_chrlength(this) == 1);
      }
#endif
      print_single(fp,abbrev,hittype,this,mate,acc1,acc2,pathnum,npaths,
		   absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		   chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		   resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p,/*circularp*/true);
      print_single(fp,abbrev,hittype,this,mate,acc1,acc2,pathnum,npaths,
		   absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		   /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
		   resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p,/*circularp*/true);
    } else if (first_read_p == true) {
      print_single(fp,abbrev,hittype,this,mate,acc1,acc2,pathnum,npaths,
		   absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		   chrpos,mate_chrpos,clipdir,hardclip5_low,hardclip5_high,resulttype,/*first_read_p*/true,
		   npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p,/*circularp*/false);
    } else {
      print_single(fp,abbrev,hittype,this,mate,acc1,acc2,pathnum,npaths,
		   absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		   chrpos,mate_chrpos,clipdir,hardclip3_low,hardclip3_high,resulttype,/*first_read_p*/false,
		   npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p,/*circularp*/false);
    }

  } else if (hittype == INSERTION) {
    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0 &&
	check_cigar_insertion(this,querylength,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			      first_read_p,/*circularp*/true) == true &&
	check_cigar_insertion(this,querylength,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			      first_read_p,/*circularp*/true) == true) {
      print_insertion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/true);
      print_insertion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/true);
    } else if (first_read_p == true) {
      print_insertion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,clipdir,hardclip5_low,hardclip5_high,resulttype,/*first_read_p*/true,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/false);
    } else {
      print_insertion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,clipdir,hardclip3_low,hardclip3_high,resulttype,/*first_read_p*/false,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/false);
    }

  } else if (hittype == DELETION) {
    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0 &&
	check_cigar_deletion(this,querylength,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			     first_read_p,/*circularp*/true) == true &&
	check_cigar_deletion(this,querylength,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			     first_read_p,/*circularp*/true) == true) {
      print_deletion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		     absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		     chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		     resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		     invertp,invert_mate_p,/*circularp*/true);
      print_deletion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		     absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		     /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
		     resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		     invertp,invert_mate_p,/*circularp*/true);
    } else if (first_read_p == true) {
      print_deletion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		     absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		     chrpos,mate_chrpos,clipdir,hardclip5_low,hardclip5_high,resulttype,/*first_read_p*/true,
		     npaths_mate,quality_shift,sam_read_group_id,
		     invertp,invert_mate_p,/*circularp*/false);
    } else {
      print_deletion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		     absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		     chrpos,mate_chrpos,clipdir,hardclip3_low,hardclip3_high,resulttype,/*first_read_p*/false,
		     npaths_mate,quality_shift,sam_read_group_id,
		     invertp,invert_mate_p,/*circularp*/false);
    }

  } else if (hittype == HALFSPLICE_DONOR) {
    donor = Stage3end_substring_donor(this);

    /* Code taken from that for XS tag for print_halfdonor and print_halfacceptor */
    if ((sensedir = Substring_chimera_sensedir(donor)) == SENSE_FORWARD) {
      if (Substring_plusp(donor) == true) {
	donor_strand = '+';
      } else {
	donor_strand = '-';
      }
    } else if (sensedir == SENSE_ANTI) {
      if (Substring_plusp(donor) == true) {
	donor_strand = '-';
      } else {
	donor_strand = '+';
      }
    } else if (force_xs_direction_p == true) {
      donor_strand = '+';
    } else {
      donor_strand = '?';
    }

    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0 &&
	check_cigar_halfdonor(donor,querylength,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			     first_read_p,/*circularp*/true) == true &&
	check_cigar_halfdonor(donor,querylength,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			      first_read_p,/*circularp*/true) == true) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/chrpos,chrpos,/*acceptor_chrpos*/-1U,mate_chrpos,
		      /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		      donor_strand,/*acceptor_strand*/'\0',/*donor_chr*/NULL,/*acceptor_chr*/NULL,
		      /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		      /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/true);
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*chrpos*/1,/*acceptor_chrpos*/-1U,mate_chrpos,
		      /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		      donor_strand,/*acceptor_strand*/'\0',/*donor_chr*/NULL,/*acceptor_chr*/NULL,
		      /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		      /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/true);
    } else if (first_read_p == true) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/chrpos,chrpos,/*acceptor_chrpos*/-1U,mate_chrpos,
		      clipdir,hardclip5_low,hardclip5_high,resulttype,/*first_read_p*/true,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		      donor_strand,/*acceptor_strand*/'\0',/*donor_chr*/NULL,/*acceptor_chr*/NULL,
		      /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		      /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/false);
    } else {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/chrpos,chrpos,/*acceptor_chrpos*/-1U,mate_chrpos,
		      clipdir,hardclip3_low,hardclip3_high,resulttype,/*first_read_p*/false,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		      donor_strand,/*acceptor_strand*/'\0',/*donor_chr*/NULL,/*acceptor_chr*/NULL,
		      /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		      /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/false);
    }

  } else if (hittype == HALFSPLICE_ACCEPTOR) {
    acceptor = Stage3end_substring_acceptor(this);

    /* Code taken from that for XS tag for print_halfdonor and print_halfacceptor */
    if ((sensedir = Substring_chimera_sensedir(acceptor)) == SENSE_FORWARD) {
      if (Substring_plusp(acceptor) == true) {
	acceptor_strand = '+';
      } else {
	acceptor_strand = '-';
      }
    } else if (sensedir == SENSE_ANTI) {
      if (Substring_plusp(acceptor) == true) {
	acceptor_strand = '-';
      } else {
	acceptor_strand = '+';
      }
    } else if (force_xs_direction_p == true) {
      acceptor_strand = '+';
    } else {
      acceptor_strand = '?';
    }

    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0 &&
	check_cigar_halfacceptor(acceptor,querylength,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
				 first_read_p,/*circularp*/true) == true &&
	check_cigar_halfacceptor(acceptor,querylength,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
				 first_read_p,/*circularp*/true) == true) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/chrpos,/*donor_chrpos*/-1U,chrpos,mate_chrpos,
			 /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
			 /*donor_strand*/'\0',acceptor_strand,/*donor_chr*/NULL,/*acceptor_chr*/NULL,
			 /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
			 /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/true);
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,/*donor_chrpos*/-1U,/*chrpos*/1,mate_chrpos,
			 /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
			 /*donor_strand*/'\0',acceptor_strand,/*donor_chr*/NULL,/*acceptor_chr*/NULL,
			 /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
			 /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/true);
    } else if (first_read_p == true) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/chrpos,/*donor_chrpos*/-1U,chrpos,mate_chrpos,
			 clipdir,hardclip5_low,hardclip5_high,resulttype,/*first_read_p*/true,
			 npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
			 /*donor_strand*/'\0',acceptor_strand,/*donor_chr*/NULL,/*acceptor_chr*/NULL,
			 /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
			 /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/false);
    } else {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/chrpos,/*donor_chrpos*/-1U,chrpos,mate_chrpos,
			 clipdir,hardclip3_low,hardclip3_high,resulttype,/*first_read_p*/false,
			 npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
			 /*donor_strand*/'\0',acceptor_strand,/*donor_chr*/NULL,/*acceptor_chr*/NULL,
			 /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
			 /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/false);
    }

  } else if (hittype == SPLICE || hittype == SAMECHR_SPLICE || hittype == TRANSLOC_SPLICE) {
    /* Follows print_splice_distance() in substring.c */
    donor = Stage3end_substring_donor(this);
    acceptor = Stage3end_substring_acceptor(this);

    if (donor == NULL || acceptor == NULL) {
      abort();
    } else if (hittype == TRANSLOC_SPLICE || (hittype == SAMECHR_SPLICE && merge_samechr_p == false)) {
      /* Stage3end_chrnum(this) == 0 || Stage3end_distance(this) == 0U */
      /* distant splice */
      if (first_read_p == true) {
	print_exon_exon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			mate_chrpos,clipdir,hardclip5_low,hardclip5_high,resulttype,/*first_read_p*/true,
			npaths_mate,quality_shift,sam_read_group_id,
			invertp,invert_mate_p);
      } else {
	print_exon_exon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			mate_chrpos,clipdir,hardclip3_low,hardclip3_high,resulttype,/*first_read_p*/false,
			npaths_mate,quality_shift,sam_read_group_id,
			invertp,invert_mate_p);
      }
    } else {
      normalp = true;
      sensep = (Stage3end_sensedir(this) == SENSE_FORWARD);

      if (Substring_plusp(donor) != Substring_plusp(acceptor)) {
	/* inversion */
	normalp = false;
      } else if (Substring_plusp(donor) == true) {
	if (sensep == true) {
	  if (Substring_genomicstart(acceptor) < Substring_genomicstart(donor)) {
	    /* scramble */
	    normalp = false;
	  }
	} else {
	  if (Substring_genomicstart(donor) < Substring_genomicstart(acceptor)) {
	    /* scramble */
	    normalp = false;
	  }
	}
      } else {
	if (sensep == true) {
	  if (Substring_genomicstart(donor) < Substring_genomicstart(acceptor)) {
	    /* scramble */
	    normalp = false;
	  }
	} else {
	  if (Substring_genomicstart(acceptor) < Substring_genomicstart(donor)) {
	    /* scramble */
	    normalp = false;
	  }
	}
      }

      if (normalp == true) {
	querylength = Shortread_fulllength(queryseq);
	if ((circularpos = Stage3end_circularpos(this)) > 0 &&
	    check_cigar_localsplice(this,mate,querylength,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
				    first_read_p,/*circularp*/true) == true &&
	    check_cigar_localsplice(this,mate,querylength,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
				    first_read_p,/*circularp*/true) == true) {
	  print_localsplice(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			    absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			    chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			    resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			    invertp,invert_mate_p,/*circularp*/true);
	  print_localsplice(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			    absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			    /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			    resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			    invertp,invert_mate_p,/*circularp*/true);
	} else if (first_read_p == true) {
	  print_localsplice(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			    absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			    chrpos,mate_chrpos,clipdir,hardclip5_low,hardclip5_high,
			    resulttype,/*first_read_p*/true,npaths_mate,quality_shift,sam_read_group_id,
			    invertp,invert_mate_p,/*circularp*/false);
	} else {
	  print_localsplice(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			    absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			    chrpos,mate_chrpos,clipdir,hardclip3_low,hardclip3_high,
			    resulttype,/*first_read_p*/false,npaths_mate,quality_shift,sam_read_group_id,
			    invertp,invert_mate_p,/*circularp*/false);
	}

      } else {
	if (first_read_p == true) {
	  print_exon_exon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			  absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			  mate_chrpos,clipdir,hardclip5_low,hardclip5_high,
			  resulttype,/*first_read_p*/true,npaths_mate,quality_shift,sam_read_group_id,
			  invertp,invert_mate_p);
	} else {
	  print_exon_exon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			  absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			  mate_chrpos,clipdir,hardclip3_low,hardclip3_high,
			  resulttype,/*first_read_p*/false,npaths_mate,quality_shift,sam_read_group_id,
			  invertp,invert_mate_p);
	}
      }
    }
      
  } else if (hittype == ONE_THIRD_SHORTEXON || hittype == TWO_THIRDS_SHORTEXON || hittype == SHORTEXON) {
    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0 &&
	check_cigar_shortexon(this,mate,querylength,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			      first_read_p,/*circularp*/true) == true &&
	check_cigar_shortexon(this,mate,querylength,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			      first_read_p,/*circularp*/true) == true) {
      print_shortexon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/true);
      print_shortexon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/true);
    } else if (first_read_p == true) {
      print_shortexon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,clipdir,hardclip5_low,hardclip5_high,
		      resulttype,/*first_read_p*/true,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/false);
    } else {
      print_shortexon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,clipdir,hardclip3_low,hardclip3_high,
		      resulttype,/*first_read_p*/false,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/false);
    }

  } else if (hittype == GMAP) {
    /* Note: sam_paired_p must be true because we are calling GMAP only on halfmapping uniq */

    if (mate == NULL) {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/hardclip5_high,
				  this,Shortread_fulllength(queryseq));
      mate_chrpos = 0U;
      hardclip3_low = hardclip3_high = 0;

    } else if (first_read_p == true) {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/hardclip5_low,/*hardclip_high*/hardclip5_high,
				  this,Shortread_fulllength(queryseq));
      mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/hardclip3_low,/*hardclip_high*/hardclip3_high,
				       mate,Shortread_fulllength(queryseq_mate));
    } else {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/hardclip3_low,/*hardclip_high*/hardclip3_high,
				  this,Shortread_fulllength(queryseq));
      mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/hardclip5_low,/*hardclip_high*/hardclip5_high,
				       mate,Shortread_fulllength(queryseq_mate));
    }

    flag = SAM_compute_flag(Stage3end_plusp(this),mate,resulttype,first_read_p,
			    pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			    invertp,invert_mate_p);

    querylength = Shortread_fulllength(queryseq);
    if ((circularpos = Stage3end_circularpos(this)) > 0 &&
	Pair_check_cigar(Stage3end_pairarray(this),Stage3end_npairs(this),querylength,
			 /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,
			 /*watsonp*/Stage3end_plusp(this),Stage3end_cdna_direction(this),
			 first_read_p,/*circularp*/true) == true &&
	Pair_check_cigar(Stage3end_pairarray(this),Stage3end_npairs(this),querylength,
			 /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,
			 /*watsonp*/Stage3end_plusp(this),Stage3end_cdna_direction(this),
			 first_read_p,/*circularp*/true) == true) {
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     /*clipdir*/+1,/*hardclip_low*/0,/*hardclip_high*/querylength-circularpos,Shortread_fulllength(queryseq),
		     /*watsonp*/Stage3end_plusp(this),Stage3end_cdna_direction(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,first_read_p,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,chrpos,Stage3end_chrlength(this),
		     queryseq,resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
		     mate_chrpos,Stage3end_chrlength(mate),/*mate_cdna_direction*/Stage3end_cdna_direction(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/true,/*merged_overlap_p*/false);
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     /*clipdir*/+1,/*hardclip_low*/circularpos,/*hardclip_high*/0,Shortread_fulllength(queryseq),
		     /*watsonp*/Stage3end_plusp(this),Stage3end_cdna_direction(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,first_read_p,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,/*chrpos*/1,Stage3end_chrlength(this),
		     queryseq,resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
		     mate_chrpos,Stage3end_chrlength(mate),/*mate_cdna_direction*/Stage3end_cdna_direction(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/true,/*merged_overlap_p*/false);
    } else if (first_read_p == true) {
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     clipdir,hardclip5_low,hardclip5_high,Shortread_fulllength(queryseq),
		     Stage3end_plusp(this),Stage3end_cdna_direction(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,/*first_read_p*/true,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,chrpos,Stage3end_chrlength(this),
		     queryseq,resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
		     mate_chrpos,Stage3end_chrlength(mate),/*mate_cdna_direction*/Stage3end_cdna_direction(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/false,/*merged_overlap_p*/false);
    } else {
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     clipdir,hardclip3_low,hardclip3_high,Shortread_fulllength(queryseq),
		     Stage3end_plusp(this),Stage3end_cdna_direction(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,/*first_read_p*/false,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,chrpos,Stage3end_chrlength(this),
		     queryseq,resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),
		     mate_chrpos,Stage3end_chrlength(mate),/*mate_cdna_direction*/Stage3end_cdna_direction(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/false,/*merged_overlap_p*/false);
    }
  } else {
    abort();
  }

  return;
}



void
SAM_print_paired (Result_T result, Resulttype_T resulttype,
		  Univ_IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool clip_overlap_p, bool merge_overlap_p,
		  bool merge_samechr_p, int quality_shift, char *sam_read_group_id) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Stage3end_T *stage3array1, *stage3array2, stage3, mate, hit5, hit3;
  Chrpos_T chrpos, chrpos5, chrpos3;
  int npaths, npaths1, npaths2, pathnum;
  int first_absmq, second_absmq, first_absmq1, second_absmq1, first_absmq2, second_absmq2;
  int hardclip5_low = 0, hardclip5_high = 0, hardclip3_low = 0, hardclip3_high = 0, clipdir;
  char *acc1, *acc2;
  Pairtype_T pairtype;
  FILE *fp, *fp_xs;
  char *abbrev, *abbrev_xs;

  struct Pair_T *pairarray;
  int npairs;
  char *queryseq_merged, *quality_merged;
  int querylength_merged;
  int flag;

  acc1 = Shortread_accession(queryseq1);
  acc2 = Shortread_accession(queryseq2); /* NULL, unless --allow-pe-name-mismatch is specified */

  if (resulttype == PAIREDEND_NOMAPPING) {
    if (nofailsp == true) {
      /* No output */
      return;
      
    } else 
      SAM_print_nomapping(fp_nomapping,ABBREV_NOMAPPING_1,queryseq1,/*mate*/(Stage3end_T) NULL,
			  acc1,acc2,chromosome_iit,resulttype,
			  /*first_read_p*/true,/*npaths*/0,/*npaths_mate*/0,
			  /*mate_chrpos*/0U,quality_shift,
			  sam_read_group_id,invert_first_p,invert_second_p);
      SAM_print_nomapping(fp_nomapping,ABBREV_NOMAPPING_2,queryseq2,/*mate*/(Stage3end_T) NULL,
			  acc1,acc2,chromosome_iit,resulttype,
			  /*first_read_p*/false,/*npaths*/0,/*npaths_mate*/0,
			  /*mate_chrpos*/0U,quality_shift,
			  sam_read_group_id,invert_second_p,invert_first_p);
      if (failedinput_root != NULL) {
	if (fastq_format_p == true) {
	  Shortread_print_query_pairedend_fastq(fp_failedinput_1,fp_failedinput_2,queryseq1,queryseq2,
						invert_first_p,invert_second_p);
	} else {
	  Shortread_print_query_pairedend_fasta(fp_failedinput_1,queryseq1,queryseq2,
						invert_first_p,invert_second_p);
	}
      }

  } else {
    if (failsonlyp == true) {
      /* Unwanted success: skip */

    } else if (resulttype == CONCORDANT_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];
      hit5 = Stage3pair_hit5(stage3pair);
      hit3 = Stage3pair_hit3(stage3pair);

      if (Stage3pair_circularp(stage3pair) == true) {
	/* Don't resolve overlaps on a circular alignment */
	clipdir = 0;
	hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;
	fp = fp_concordant_circular;
	abbrev = ABBREV_CONCORDANT_CIRCULAR;

      } else if (clip_overlap_p == false && merge_overlap_p == false) {
	clipdir = 0;
	hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;
	fp = fp_concordant_uniq;
	abbrev = ABBREV_CONCORDANT_UNIQ;

      } else {
	clipdir = Stage3pair_overlap(&hardclip5_low,&hardclip5_high,&hardclip3_low,&hardclip3_high,stage3pair);
	debug3(printf("clipdir %d with hardclip5 = %d..%d, hardclip3 = %d..%d\n",
		      clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high));
	fp = fp_concordant_uniq;
	abbrev = ABBREV_CONCORDANT_UNIQ;
      }

      chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1));
      chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2));

      if (merge_overlap_p == false || clipdir == 0) {
	/* print first end */
	SAM_print(fp,abbrev,hit5,/*mate*/hit3,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		  Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		  Stage3pair_mapq_score(stage3pair),chromosome_iit,
		  /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		  clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		  resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		  quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		  merge_samechr_p);

	/* print second end */
	SAM_print(fp,abbrev,hit3,/*mate*/hit5,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		  Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		  Stage3pair_mapq_score(stage3pair),chromosome_iit,
		  /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		  clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		  resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		  quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		  merge_samechr_p);

      } else {
	/* merge_overlap_p == true and overlap was found */
	pairarray = Stage3pair_merge(&npairs,&querylength_merged,&queryseq_merged,&quality_merged,
				     stage3pair,queryseq1,queryseq2,
				     /*querylength5*/Shortread_fulllength(queryseq1),
				     /*querylength3*/Shortread_fulllength(queryseq2),
				     clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high);
	/* printf("queryseq_merged: %s\n",queryseq_merged); */
	if (clipdir >= 0) {
	  chrpos = chrpos5;
	} else {
	  chrpos = chrpos3;
	}
	/* merging changes resulttype from UNPAIRED_UNIQ to SINGLEEND_UNIQ */
	flag = SAM_compute_flag(Stage3end_plusp(hit5),/*mate*/NULL,/*resulttype*/SINGLEEND_UNIQ,/*first_read_p*/true,
				/*pathnum*/1,/*npaths*/1,/*npaths_mate*/0,
				Stage3pair_absmq_score(stage3pair),first_absmq,/*invertp*/false,
				/*invert_mate_p*/false);
	Pair_print_sam(fp_unpaired_uniq,/*abbrev*/ABBREV_UNPAIRED_UNIQ,pairarray,npairs,
		       acc1,/*acc2*/NULL,Stage3end_chrnum(hit5),chromosome_iit,/*usersegment*/(Sequence_T) NULL,
		       /*queryseq_ptr*/queryseq_merged,/*quality_string*/quality_merged,
		       /*clipdir*/0,/*hardclip_low*/0,/*hardclip_high*/0,/*querylength*/querylength_merged,
		       Stage3end_plusp(hit5),Stage3end_cdna_direction(hit5),
		       /*chimera_part*/0,/*chimera*/NULL,quality_shift,/*first_read_p*/true,
		       /*pathnum*/1,/*npaths*/1,
#if 0
		       Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
#else
		       /*absmq_score*/MAX_QUALITY_SCORE,/*first_absmq*/MAX_QUALITY_SCORE,/*second_absmq*/0,
#endif
		       chrpos,Stage3end_chrlength(hit5),/*queryseq*/NULL,resulttype,flag,
		       /*pair_mapq_score*/MAX_QUALITY_SCORE,/*end_mapq_score*/MAX_QUALITY_SCORE,
		       /*mate_chrnum*/0,/*mate_effective_chrnum*/0,/*mate_chrpos*/0,/*mate_chrlength*/0,
		       /*mate_cdna_direction*/0,/*pairedlength*/0,
		       sam_read_group_id,/*invertp*/false,/*circularp*/false,/*merged_overlap_p*/true);
	if (quality_merged != NULL) {
	  FREE_OUT(quality_merged);
	}
	FREE_OUT(queryseq_merged);
	FREE_OUT(pairarray);
      }

    } else if (resulttype == CONCORDANT_TRANSLOC) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	if (1 || failedinput_root != NULL) {
	  /* Not able to print as input */
	  /* Print as nomapping, but send to fp_concordant_transloc */
	  SAM_print_nomapping(fp_concordant_transloc,ABBREV_CONCORDANT_TRANSLOC,
			      queryseq1,/*mate*/(Stage3end_T) NULL,
			      acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/true,npaths,/*npaths_mate*/npaths,
			      /*mate_chrpos*/0U,quality_shift,
			      sam_read_group_id,invert_first_p,invert_second_p);
	  SAM_print_nomapping(fp_concordant_transloc,ABBREV_CONCORDANT_TRANSLOC,
			      queryseq2,/*mate*/(Stage3end_T) NULL,
			      acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/false,npaths,/*npaths_mate*/npaths,
			      /*mate_chrpos*/0U,quality_shift,
			      sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	/* Note: We are ignoring merge_overlap for concordant_transloc */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);

	  if (Stage3pair_circularp(stage3pair) == true) {
	    /* Don't resolve overlaps on a circular alignment */
	    clipdir = 0;
	    hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  } else if (clip_overlap_p == false) {
	    clipdir = 0;
	    hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  } else {
	    clipdir = Stage3pair_overlap(&hardclip5_low,&hardclip5_high,&hardclip3_low,&hardclip3_high,stage3pair);
	    debug3(printf("clipdir %d with hardclip5 = %d..%d, hardclip3 = %d..%d\n",
			  clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high));
	  }

	  chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1));
	  chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2));

	  /* print first end */
	  SAM_print(fp_concordant_transloc,ABBREV_CONCORDANT_TRANSLOC,
		    hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		    resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp_concordant_transloc,ABBREV_CONCORDANT_TRANSLOC,
		    hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		    resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }
    
    } else if (resulttype == CONCORDANT_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	/* Print as nomapping, but send to fp_concordant_mult_xs */
	SAM_print_nomapping(fp_concordant_mult_xs_1,ABBREV_CONCORDANT_MULT_XS,
			    queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,npaths,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_concordant_mult_xs_1,ABBREV_CONCORDANT_MULT_XS,
			    queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,npaths,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);
	if (failedinput_root != NULL) {
	  if (fastq_format_p == true) {
	    Shortread_print_query_pairedend_fastq(fp_failedinput_1,fp_failedinput_2,queryseq1,queryseq2,
						  invert_first_p,invert_second_p);
	  } else {
	    Shortread_print_query_pairedend_fasta(fp_failedinput_1,queryseq1,queryseq2,
						  invert_first_p,invert_second_p);
	  }
	}

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);

	  if (Stage3pair_circularp(stage3pair) == true) {
	    /* Don't resolve overlaps on a circular alignment */
	    clipdir = 0;
	    hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  } else if (clip_overlap_p == false && merge_overlap_p == false) {
	    clipdir = 0;
	    hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  } else {
	    clipdir = Stage3pair_overlap(&hardclip5_low,&hardclip5_high,&hardclip3_low,&hardclip3_high,stage3pair);
	    debug3(printf("clipdir %d with hardclip5 = %d..%d, hardclip3 = %d..%d\n",
			  clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high));
	  }

	  chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1));
	  chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2));

	  if (merge_overlap_p == false || clipdir == 0) {
	    /* print first end */
	    SAM_print(fp_concordant_mult,ABBREV_CONCORDANT_MULT,
		      hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		      Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		      Stage3pair_mapq_score(stage3pair),chromosome_iit,
		      /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		      clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		      resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		      quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		      merge_samechr_p);

	    /* print second end */
	    SAM_print(fp_concordant_mult,ABBREV_CONCORDANT_MULT,
		      hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		      Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		      Stage3pair_mapq_score(stage3pair),chromosome_iit,
		      /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		      clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		      resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		      quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		      merge_samechr_p);
	    
	  } else {
	    /* merge_overlap_p == true and overlap was found */
	    pairarray = Stage3pair_merge(&npairs,&querylength_merged,&queryseq_merged,&quality_merged,
					 stage3pair,queryseq1,queryseq2,
					 /*querylength5*/Shortread_fulllength(queryseq1),
					 /*querylength3*/Shortread_fulllength(queryseq2),
					 clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high);
	    /* printf("queryseq_merged: %s\n",queryseq_merged); */
	    if (clipdir >= 0) {
	      chrpos = chrpos5;
	    } else {
	      chrpos = chrpos3;
	    }
	    /* merging changes resulttype from UNPAIRED_UNIQ to SINGLEEND_UNIQ */
	    flag = SAM_compute_flag(Stage3end_plusp(hit5),/*mate*/NULL,/*resulttype*/SINGLEEND_UNIQ,/*first_read_p*/true,
				    /*pathnum*/1,/*npaths*/1,/*npaths_mate*/0,
				    Stage3pair_absmq_score(stage3pair),first_absmq,/*invertp*/false,
				    /*invert_mate_p*/false);
	    Pair_print_sam(fp_concordant_mult,ABBREV_CONCORDANT_MULT,pairarray,npairs,
			   acc1,/*acc2*/NULL,Stage3end_chrnum(hit5),chromosome_iit,
			   /*usersegment*/(Sequence_T) NULL,
			   /*queryseq_ptr*/queryseq_merged,/*quality_string*/quality_merged,
			   /*clipdir*/0,/*hardclip_low*/0,/*hardclip_high*/0,/*querylength*/querylength_merged,
			   Stage3end_plusp(hit5),Stage3end_cdna_direction(hit5),
			   /*chimera_part*/0,/*chimera*/NULL,quality_shift,/*first_read_p*/true,pathnum,npaths,
#if 0
			   Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
#else
			   /*absmq_score*/MAX_QUALITY_SCORE,/*first_absmq*/MAX_QUALITY_SCORE,/*second_absmq*/0,
#endif
			   chrpos,Stage3end_chrlength(hit5),/*queryseq*/NULL,resulttype,flag,
			   /*pair_mapq_score*/MAX_QUALITY_SCORE,/*end_mapq_score*/MAX_QUALITY_SCORE,
			   /*mate_chrnum*/0,/*mate_effective_chrnum*/0,/*mate_chrpos*/0,/*mate_chrlength*/0,
			   /*mate_cdna_direction*/0,/*pairedlength*/0,
			   sam_read_group_id,/*invertp*/false,/*circularp*/false,/*merged_overlap_p*/true);
	    if (quality_merged != NULL) {
	      FREE_OUT(quality_merged);
	    }
	    FREE_OUT(queryseq_merged);
	    FREE_OUT(pairarray);
	  }
	}
      }

    } else if (resulttype == PAIRED_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];
      if (Stage3pair_circularp(stage3pair) == true) {
	fp = fp_paired_uniq_circular;
	abbrev = ABBREV_PAIRED_UNIQ_CIRCULAR;
      } else if ((pairtype = Stage3pair_pairtype(stage3pair)) == PAIRED_INVERSION) {
	fp = fp_paired_uniq_inv;
	abbrev = ABBREV_PAIRED_UNIQ_INV;
      } else if (pairtype == PAIRED_SCRAMBLE) {
	fp = fp_paired_uniq_scr;
	abbrev = ABBREV_PAIRED_UNIQ_SCR;
      } else if (pairtype == PAIRED_TOOLONG) {
	fp = fp_paired_uniq_long;
	abbrev = ABBREV_PAIRED_UNIQ_LONG;
      } else {
	fprintf(stderr,"Unexpected pairtype %d\n",pairtype);
	abort();
      }

      hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

      hit5 = Stage3pair_hit5(stage3pair);
      hit3 = Stage3pair_hit3(stage3pair);
      chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1));
      chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2));

      /* print first end */
      SAM_print(fp,abbrev,hit5,/*mate*/hit3,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		/*clipdir*/0,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		merge_samechr_p);

      /* print second end */
      SAM_print(fp,abbrev,hit3,/*mate*/hit5,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		/*clipdir*/0,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		merge_samechr_p);

    } else if (resulttype == PAIRED_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	/* Print as nomapping, but send to fp_concordant_mult */
	SAM_print_nomapping(fp_paired_mult_xs_1,ABBREV_PAIRED_MULT_XS,
			    queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,npaths,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_paired_mult_xs_1,ABBREV_PAIRED_MULT_XS,
			    queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,npaths,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);

	if (failedinput_root != NULL) {
	  if (fastq_format_p == true) {
	    Shortread_print_query_pairedend_fastq(fp_failedinput_1,fp_failedinput_2,queryseq1,queryseq2,
						  invert_first_p,invert_second_p);
	  } else {
	    Shortread_print_query_pairedend_fasta(fp_failedinput_1,queryseq1,queryseq2,
						  invert_first_p,invert_second_p);
	  }
	}

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);
	  chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1));
	  chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2));

	  /* print first end */
	  SAM_print(fp_paired_mult,ABBREV_PAIRED_MULT,
		    hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		    resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp_paired_mult,ABBREV_PAIRED_MULT,
		    hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high,
		    resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      /* Even though they are not related, we should print mate information in this situation */
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

      hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0;

      hit5 = stage3array1[0];
      hit3 = stage3array2[0];
      chrpos5 = SAM_compute_chrpos(hardclip5_low,hardclip5_high,hit5,Shortread_fulllength(queryseq1));
      chrpos3 = SAM_compute_chrpos(hardclip3_low,hardclip3_high,hit3,Shortread_fulllength(queryseq2));

      if (Stage3end_circularpos(hit5) > 0 || Stage3end_circularpos(hit3) > 0) {
	fp = fp_unpaired_circular;
	abbrev = ABBREV_UNPAIRED_CIRCULAR;
      } else {
	fp = fp_unpaired_uniq;
	abbrev = ABBREV_UNPAIRED_UNIQ;
      }

      /* print first end */
      /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      SAM_print(fp,abbrev,hit5,/*mate*/hit3,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3end_absmq_score(stage3array1[0]),first_absmq1,/*second_absmq*/0,
		Stage3end_mapq_score(stage3array1[0]),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		/*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		/*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		resulttype,/*first_read_p*/true,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		invert_first_p,invert_second_p,merge_samechr_p);

      /* print second end */
      /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      SAM_print(fp,abbrev,hit3,/*mate*/hit5,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3end_absmq_score(stage3array2[0]),first_absmq2,/*second_absmq*/0,
		Stage3end_mapq_score(stage3array2[0]),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		/*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		/*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		resulttype,/*first_read_p*/false,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		invert_second_p,invert_first_p,merge_samechr_p);

    } else if (resulttype == UNPAIRED_MULT || resulttype == UNPAIRED_TRANSLOC) {
      if (resulttype == UNPAIRED_MULT) {
	fp = fp_unpaired_mult;
	fp_xs = fp_unpaired_mult_xs_1;
	abbrev = ABBREV_UNPAIRED_MULT;
	abbrev_xs = ABBREV_UNPAIRED_MULT_XS;
      } else {
	fp = fp_xs = fp_unpaired_transloc;
	abbrev = abbrev_xs = ABBREV_UNPAIRED_TRANSLOC;
      }

      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 1) {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 1) {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      }
#endif

      /* print first end results */
      if (npaths2 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else {
	mate = stage3array2[0];
	hardclip3_low = hardclip3_high = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,mate,Shortread_fulllength(queryseq2));
      }

      if (npaths1 == 1) {
	stage3 = stage3array1[0];
	hardclip5_low = hardclip5_high = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1));

	SAM_print(fp,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths1,
		  Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		  /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		  resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		  invert_first_p,invert_second_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	if (1 || failedinput_root != NULL) {
	  /* Just printing one end as nomapping */
	  SAM_print_nomapping(fp_xs,abbrev_xs,queryseq1,mate,acc1,acc2,chromosome_iit,
			      resulttype,/*first_read_p*/true,npaths1,/*npaths_mate*/npaths2,
			      /*mate_chrpos*/chrpos3,
			      quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	}

      } else {
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  hardclip5_low = hardclip5_high = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/hardclip5_high,stage3,Shortread_fulllength(queryseq1));

	  SAM_print(fp,abbrev,stage3,mate,acc1,acc2,pathnum,npaths1,
		    Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p,merge_samechr_p);
	}
      }
			  
      /* print second end results */
      if (npaths1 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else {
	mate = stage3array1[0];
	hardclip5_low = hardclip5_high = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,mate,Shortread_fulllength(queryseq1));
      }

      if (npaths2 == 1) {
	stage3 = stage3array2[0];
	hardclip3_low = hardclip3_high = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq2));

	SAM_print(fp,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths2,
		  Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		  /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		  resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		  invert_second_p,invert_first_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	if (1 || failedinput_root != NULL) {
	  /* Just printing one end as nomapping */
	  SAM_print_nomapping(fp_xs,abbrev_xs,queryseq2,mate,acc1,acc2,chromosome_iit,
			      resulttype,/*first_read_p*/false,npaths2,/*npaths_mate*/npaths1,
			      /*mate_chrpos*/chrpos5,
			      quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else {
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  hardclip3_low = hardclip3_high = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq2));

	  SAM_print(fp,abbrev,stage3,mate,acc1,acc2,pathnum,npaths2,
		    Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p,merge_samechr_p);
	}
      }

    } else {
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

      if (resulttype == HALFMAPPING_UNIQ) {
	if (npaths1 == 1 && Stage3end_circularpos(stage3array1[0]) > 0) {
	  fp = fp_xs = fp_halfmapping_circular;
	  abbrev = abbrev_xs = ABBREV_HALFMAPPING_CIRCULAR;
	} else if (npaths2 == 1 && Stage3end_circularpos(stage3array2[0]) > 0) {
	  fp = fp_xs = fp_halfmapping_circular;
	  abbrev = abbrev_xs = ABBREV_HALFMAPPING_CIRCULAR;
	} else {
	  fp = fp_xs = fp_halfmapping_uniq;
	  abbrev = abbrev_xs = ABBREV_HALFMAPPING_UNIQ;
	}
      } else if (resulttype == HALFMAPPING_TRANSLOC) {
	fp = fp_xs = fp_halfmapping_transloc;
	abbrev = abbrev_xs = ABBREV_HALFMAPPING_TRANSLOC;
      } else if (resulttype == HALFMAPPING_MULT) {
	fp = fp_halfmapping_mult;
	fp_xs = fp_halfmapping_mult_xs_1;
	abbrev = ABBREV_HALFMAPPING_MULT;
	abbrev_xs = ABBREV_HALFMAPPING_MULT_XS;
      } else {
	abort();
      }

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 0) {
	/* Nothing to sort */
      } else if (npaths1 == 1) {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 0) {
	/* Nothing to sort */
      } else if (npaths2 == 1) {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      }
#endif


      /* print first end results */
      if (npaths2 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else {
	mate = stage3array2[0];
	hardclip3_low = hardclip3_high = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,mate,Shortread_fulllength(queryseq2));
      }

      if (npaths1 == 0) {
	if (1 || failedinput_root != NULL) {
	  /* just printing one end as nomapping */
	  /* mate should be non-NULL here */
	  SAM_print_nomapping(fp,abbrev,queryseq1,mate,acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/true,npaths1,/*npaths_mate*/npaths2,
			      /*mate_chrpos*/chrpos3,
			      quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	}

      } else if (npaths1 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array1[0];
	hardclip5_low = hardclip5_high = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1));

	SAM_print(fp,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths1,
		  Stage3end_absmq_score(stage3),first_absmq1,/*second_absmq1*/0,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		  /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		  resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		  invert_first_p,invert_second_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	if (1 || failedinput_root != NULL) {
	  /* Just printing one end as nomapping */
	  /* mate should be NULL here */
	  SAM_print_nomapping(fp_xs,abbrev_xs,queryseq1,mate,acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/true,npaths1,/*npaths_mate*/npaths2,
			      /*mate_chrpos*/chrpos3,
			      quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	}

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  hardclip5_low = hardclip5_high = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq1));

	  SAM_print(fp,abbrev,stage3,mate,acc1,acc2,pathnum,npaths1,
		    Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p,merge_samechr_p);
	}
      }
			  
      /* print second end results */
      if (npaths1 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else {
	mate = stage3array1[0];
	hardclip5_low = hardclip5_high = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,mate,Shortread_fulllength(queryseq1));
      }

      if (npaths2 == 0) {
	if (1 || failedinput_root != NULL) {
	  /* Just printing one end as nomapping */
	  /* mate should be non-NULL here */
	  SAM_print_nomapping(fp,abbrev,queryseq2,mate,acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/false,npaths2,/*npaths_mate*/npaths1,
			      /*mate_chrpos*/chrpos5,
			      quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else if (npaths2 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array2[0];
	hardclip3_low = hardclip3_high = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq2));

	SAM_print(fp,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths2,
		  Stage3end_absmq_score(stage3),first_absmq2,/*second_absmq2*/0,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		  /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		  resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		  invert_second_p,invert_first_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	if (1 || failedinput_root != NULL) {
	  /* Just printing one end as nomapping */
	  /* mate should be NULL here */
	  SAM_print_nomapping(fp_xs,abbrev_xs,queryseq2,mate,acc1,acc2,chromosome_iit,resulttype,
			      /*first_read_p*/false,npaths2,/*npaths_mate*/npaths1,
			      /*mate_chrpos*/chrpos5,
			      quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  hardclip3_low = hardclip3_high = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/0,/*hardclip_high*/0,stage3,Shortread_fulllength(queryseq2));

	  SAM_print(fp,abbrev,stage3,mate,acc1,acc2,pathnum,npaths2,
		    Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,/*hardclip5_low*/0,/*hardclip5_high*/0,/*hardclip3_low*/0,/*hardclip3_high*/0,
		    resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p,merge_samechr_p);
	}
      }

    }
  }

  return;
}



