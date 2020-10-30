/* $Id: substring.h 154591 2014-12-04 02:00:32Z twu $ */
#ifndef SUBSTRING_INCLUDED
#define SUBSTRING_INCLUDED

#include <stdio.h>
#include "mode.h"
#include "genomicpos.h"
#include "types.h"
#include "chrnum.h"
#include "shortread.h"
#include "genome.h"
#include "compress.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "bool.h"
#include "pairdef.h"

typedef enum {END, INS, DEL, DON, ACC, AMB_DON, AMB_ACC, TERM} Endtype_T;

extern char *
Endtype_string (Endtype_T endtype);

extern void
Substring_setup (bool print_nsnpdiffs_p_in, bool print_snplabels_p_in,
		 bool show_refdiff_p_in, IIT_T snps_iit_in, int *snps_divint_crosstable_in,
		 IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		 IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
		 int donor_typeint_in, int acceptor_typeint_in, int trim_mismatch_score_in,
		 bool novelsplicingp_in, bool knownsplicingp_in,
		 bool output_sam_p_in, Mode_T mode_in, Univcoord_T genomelength_in,
		 int reject_trimlength_in);

#define T Substring_T
typedef struct T *T;

extern void
Substring_alias_circular (T this);
extern void
Substring_unalias_circular (T this);

extern T
Substring_new (int nmismatches_whole, Chrnum_T chrnum, Univcoord_T chroffset,
	       Univcoord_T chrhigh, Chrpos_T chrlength,
	       Univcoord_T left, Univcoord_T genomicstart, Univcoord_T genomicend,
	       Compress_T query_compress, Endtype_T start_endtype, Endtype_T end_endtype,
	       int querystart, int queryend, int querylength,
	       Univcoord_T alignstart, Univcoord_T alignend, int genomiclength,
	       int extraleft, int extraright, bool exactp,
	       bool plusp, int genestrand, bool first_read_p,
	       bool trim_left_p, bool trim_right_p, int minlength);

extern float
Substring_compute_mapq (T this, Compress_T query_compress, char *quality_string, bool trim_terminals_p);

extern int
Substring_display_prep (char **deletion, T this, char *query, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Genome_T genome, int deletion_pos, int deletion_length);

extern bool
Substring_bad_stretch_p (T this, Compress_T query_compress_fwd, Compress_T query_compress_rev);

extern void
Substring_free (T *old);

extern bool
Substring_contains_p (T this, int querypos);
extern void
Substring_print_ends (T this, int chrnum);
extern int
Substring_compare (T substring1, T substring2, int alias1, int alias2, Chrpos_T chrlength1, Chrpos_T chrlength2);
extern bool
Substring_overlap_p (T substring1, T substring2);
extern Chrpos_T
Substring_insert_length (T substring5, T substring3);
extern bool
Substring_overlap_point_trimmed_p (T substring, Univcoord_T endpos);
extern Univcoord_T
Substring_overlap_segment_trimmed (T substring1, T substring2);

extern Univcoord_T
Substring_splicecoord (T this);
extern int
Substring_splicesites_knowni (T this);
extern Univcoord_T
Substring_splicecoord_A (T this);
extern Univcoord_T
Substring_splicecoord_D (T this);

extern bool
Substring_plusp (T this);
extern int
Substring_genestrand (T this);
extern bool
Substring_first_read_p (T this);
extern char *
Substring_genomic_bothdiff (T this);
extern char *
Substring_genomic_refdiff (T this);
extern char *
Substring_genomic_querydir (T this);
extern int
Substring_nmismatches_whole (T this);
extern int
Substring_nmismatches_bothdiff (T this);
extern int
Substring_nmismatches_refdiff (T this);
extern int
Substring_nmatches (T this);
extern int
Substring_nmatches_posttrim (T this);
extern void
Substring_set_nmismatches_terminal (T this, int nmismatches_whole, int nmismatches_bothdiff);
extern Endtype_T
Substring_start_endtype (T this);
extern Endtype_T
Substring_end_endtype (T this);
extern void
Substring_set_endtypes (T this, Endtype_T start_endtype, Endtype_T end_endtype);
extern float
Substring_mapq_loglik (T this);
extern int
Substring_trim_left (T this);
extern int
Substring_trim_right (T this);
extern bool
Substring_trim_left_splicep (T this);
extern bool
Substring_trim_right_splicep (T this);
extern int
Substring_querystart (T this);
extern int
Substring_querystart_orig (T this);
extern int
Substring_queryend (T this);
extern int
Substring_queryend_orig (T this);
extern int
Substring_querylength (T this);
extern int
Substring_match_length (T this);
extern int
Substring_match_length_orig (T this);
extern Chrpos_T
Substring_genomic_alignment_length (T this);

extern Chrnum_T
Substring_chrnum (T this);
extern Univcoord_T
Substring_chroffset (T this);
extern Univcoord_T
Substring_chrhigh (T this);
extern Chrpos_T
Substring_chrlength (T this);
extern Univcoord_T
Substring_alignstart (T this);
extern Univcoord_T
Substring_alignend (T this);
extern Univcoord_T
Substring_alignstart_trim (T this);
extern Univcoord_T
Substring_alignend_trim (T this);
extern Univcoord_T
Substring_left_genomicseg (T this);
extern Univcoord_T
Substring_genomicstart (T this);
extern Univcoord_T
Substring_genomicend (T this);
extern Chrpos_T
Substring_genomiclength (T this);

extern Chrpos_T
Substring_chrstart (T this);
extern Chrpos_T
Substring_chrend (T this);

extern double
Substring_chimera_prob (T this);
extern double
Substring_chimera_prob_2 (T this);
extern int
Substring_chimera_pos (T this);
extern int
Substring_chimera_pos_A (T this);
extern int
Substring_chimera_pos_D (T this);
extern bool
Substring_chimera_knownp (T this);
extern int
Substring_nchimera_known (T this);
extern int
Substring_nchimera_novel (T this);
extern int
Substring_chimera_sensedir (T this);
extern bool
Substring_chimera_sensep (T this);
extern int
Substring_circularpos (T this);


extern T
Substring_copy (T old);

extern T
Substring_new_donor (Univcoord_T donor_coord, int donor_knowni, int donor_pos, int donor_nmismatches,
		     double donor_prob, Univcoord_T left, Compress_T query_compress,
		     int querylength, bool plusp, int genestrand, bool first_read_p, bool sensep,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength);
extern T
Substring_new_acceptor (Univcoord_T acceptor_coord, int acceptor_knowni, int acceptor_pos, int acceptor_nmismatches,
			double acceptor_prob, Univcoord_T left, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, bool first_read_p, bool sensep,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength);
extern T
Substring_new_shortexon (Univcoord_T acceptor_coord, int acceptor_knowni, Univcoord_T donor_coord, int donor_knowni,
			 int acceptor_pos, int donor_pos, int nmismatches,
			 double acceptor_prob, double donor_prob, Univcoord_T left,
			 Compress_T query_compress, int querylength,
			 bool plusp, int genestrand, bool first_read_p, bool sensep,
			 bool acceptor_ambp, bool donor_ambp,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength);

extern List_T
Substring_sort_chimera_halves (List_T hitlist, bool ascendingp);


extern void
Substring_print_m8 (FILE *fp, T substring, Shortread_T headerseq, char *acc_suffix,
		    char *chr, bool invertp);
extern void
Substring_print_single (FILE *fp, T substring, Shortread_T queryseq,
			char *chr, bool invertp);
extern void
Substring_print_insertion_1 (FILE *fp, T substring1, T substring2, int nindels, 
			     Shortread_T queryseq, char *chr, bool invertp);
extern void
Substring_print_insertion_2 (FILE *fp, T substring1, T substring2, int nindels,
			     Shortread_T queryseq, char *chr, bool invertp);
extern void
Substring_print_deletion_1 (FILE *fp, T substring1, T substring2, int nindels, 
			    char *deletion, Shortread_T queryseq, char *chr,
			    bool invertp);
extern void
Substring_print_deletion_2 (FILE *fp, T substring1, T substring2, int nindels, 
			    Shortread_T queryseq, char *chr, bool invertp);
extern void
Substring_print_donor (FILE *fp, T donor, bool sensep, bool invertp, Shortread_T queryseq,
		       Univ_IIT_T chromosome_iit, T acceptor, Chrpos_T chimera_distance);
extern void 
Substring_print_acceptor (FILE *fp, T acceptor, bool sensep, bool invertp, Shortread_T queryseq,
			  Univ_IIT_T chromosome_iit, T donor, Chrpos_T chimera_distance);
extern void
Substring_print_shortexon (FILE *fp, T shortexon, bool sensep, bool invertp, Shortread_T queryseq,
			   Univ_IIT_T chromosome_iit, Chrpos_T distance1, Chrpos_T distance2);

extern void
Substring_print_gmap (FILE *fp, struct Pair_T *pairs, int npairs, int nsegments, bool invertedp,
		      Endtype_T start_endtype, Endtype_T end_endtype,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		      int querylength, bool watsonp, int cdna_direction, int score,
		      int insertlength, int pairscore, int mapq_score, Univ_IIT_T chromosome_iit);

extern bool
Substring_contains_known_splicesite (T this);

extern Overlap_T
Substring_gene_overlap (T this, bool favor_multiexon_p);

extern long int
Substring_tally (T this, IIT_T tally_iit, int *tally_divint_crosstable);

extern bool
Substring_runlength_p (T this, IIT_T runlength_iit, int *runlength_divint_crosstable);


#ifdef USE_OLD_MAXENT
extern void
Substring_assign_donor_prob (T donor, Genome_T genome, Univ_IIT_T chromosome_iit);
extern void
Substring_assign_acceptor_prob (T acceptor, Genome_T genome, Univ_IIT_T chromosome_iit);
extern void
Substring_assign_shortexon_prob (T shortexon, Genome_T genome, Univ_IIT_T chromosome_iit);
#else
extern void
Substring_assign_donor_prob (T donor);
extern void
Substring_assign_acceptor_prob (T acceptor);
extern void
Substring_assign_shortexon_prob (T shortexon);
#endif

extern int
Substring_count_mismatches_region (T this, int trim_left, int trim_right,
				   Compress_T query_compress_fwd, Compress_T query_compress_rev);

extern List_T
Substring_convert_to_pairs (List_T pairs, T substring, Shortread_T queryseq,
			    int clipdir, int hardclip_low, int hardclip_high, bool first_read_p, int queryseq_offset);
extern List_T
Substring_add_insertion (List_T pairs, T substringA, T substringB, int insertionlength, Shortread_T queryseq,
			 int clipdir, int hardclip_low, int hardclip_high, bool first_read_p, int queryseq_offset);
extern List_T
Substring_add_deletion (List_T pairs, T substringA, T substringB, char *deletion, int deletionlength,
			int clipdir, int hardclip_low, int hardclip_high, bool first_read_p, int queryseq_offset);
extern List_T
Substring_add_intron (List_T pairs, T substringA, T substringB,
		      int clipdir, int hardclip_low, int hardclip_high, bool first_read_p, int queryseq_offset);

#undef T
#endif


