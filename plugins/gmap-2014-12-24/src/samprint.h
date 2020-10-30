/* $Id: samprint.h 160877 2015-03-13 00:31:23Z twu $ */
#ifndef SAMPRINT_INCLUDED
#define SAMPRINT_INCLUDED

#include <stdio.h>
#include "stage3hr.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "shortread.h"
#include "resulthr.h"
#include "genomicpos.h"
#include "types.h"
#include "substring.h"
#include "bool.h"
#include "intlist.h"

extern void
SAM_setup (bool quiet_if_excessive_p_in, int maxpaths_report_in,
	   char *failedinput_root_in, bool fastq_format_p_in, bool hide_soft_clips_p_in,
	   bool sam_multiple_primaries_p_in,
	   bool force_xs_direction_p_in, bool md_lowercase_variant_p_in, IIT_T snps_iit_in);

extern void
SAM_file_setup_single (FILE *fp_failedinput_in, FILE *fp_nomapping_in,
		       FILE *fp_unpaired_uniq_in, FILE *fp_unpaired_circular_in, FILE *fp_unpaired_transloc_in,
		       FILE *fp_unpaired_mult_in, FILE *fp_unpaired_mult_xs_1_in);

extern void
SAM_file_setup_paired (FILE *failedinput_1_in, FILE *failedinput_2_in, FILE *fp_nomapping_in,
		       FILE *fp_halfmapping_uniq_in, FILE *fp_halfmapping_circular_in, FILE *fp_halfmapping_transloc_in,
		       FILE *fp_halfmapping_mult_in, FILE *fp_halfmapping_mult_xs_1_in, FILE *fp_halfmapping_mult_xs_2_in,
		       FILE *fp_paired_uniq_circular_in, FILE *fp_paired_uniq_inv_in, FILE *fp_paired_uniq_scr_in,
		       FILE *fp_paired_uniq_long_in, FILE *fp_paired_mult_in, FILE *fp_paired_mult_xs_1_in, FILE *fp_paired_mult_xs_2_in,
		       FILE *fp_concordant_uniq_in, FILE *fp_concordant_circular_in, FILE *fp_concordant_transloc_in, 
		       FILE *fp_concordant_mult_in, FILE *fp_concordant_mult_xs_1_in, FILE *fp_concordant_mult_xs_2_in);

extern void
SAM_file_setup_all (FILE *failedinput_1_in, FILE *failedinput_2_in, FILE *fp_nomapping_in,
		    FILE *fp_unpaired_uniq_in, FILE *fp_unpaired_circular_in, FILE *fp_unpaired_transloc_in,
		    FILE *fp_unpaired_mult_in, FILE *fp_unpaired_mult_xs_1_in, FILE *fp_unpaired_mult_xs_2_in,
		    FILE *fp_halfmapping_uniq_in, FILE *fp_halfmapping_circular_in, FILE *fp_halfmapping_transloc_in,
		    FILE *fp_halfmapping_mult_in, FILE *fp_halfmapping_mult_xs_1_in, FILE *fp_halfmapping_mult_xs_2_in,
		    FILE *fp_paired_uniq_circular_in, FILE *fp_paired_uniq_inv_in, FILE *fp_paired_uniq_scr_in,
		    FILE *fp_paired_uniq_long_in, FILE *fp_paired_mult_in, FILE *fp_paired_mult_xs_1_in, FILE *fp_paired_mult_xs_2_in,
		    FILE *fp_concordant_uniq_in, FILE *fp_concordant_circular_in, FILE *fp_concordant_transloc_in, 
		    FILE *fp_concordant_mult_in, FILE *fp_concordant_mult_xs_1_in, FILE *fp_concordant_mult_xs_2_in);

extern Chrpos_T
SAM_compute_chrpos (int hardclip_low, int hardclip_high, Stage3end_T this, int querylength);

extern unsigned int
SAM_compute_flag (bool plusp, Stage3end_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, int npaths_mate,
		  int absmq_score, int first_absmq, bool invertp, bool invert_mate_p);

extern void
SAM_print_nomapping (FILE *fp, char *abbrev, Shortread_T queryseq, Stage3end_T mate, char *acc1, char *acc2,
		     Univ_IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int npaths, int npaths_mate, Chrpos_T mate_chrpos,
		     int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p);

extern void
SAM_print (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
	   char *acc1, char *acc2, int pathnum, int npaths,
	   int absmq_score, int first_absmq, int second_absmq, int mapq_score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq2, int pairedlength, Chrpos_T chrpos, Chrpos_T mate_chrpos,
	   int clipdir, int hardclip5_low, int hardclip5_high, int hardclip3_low, int hardclip3_high,
	   Resulttype_T resulttype, bool first_read_p, int npaths_mate, int quality_shift,
	   char *sam_read_group_id, bool invertp, bool invert_mate_p, bool merge_samechr_p);

extern void
SAM_print_paired (Result_T result, Resulttype_T resulttype,
		  Univ_IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool clip_overlap_p, bool merge_overlap_p,
		  bool merge_samechr_p, int quality_shift, char *sam_read_group_id);

#endif

