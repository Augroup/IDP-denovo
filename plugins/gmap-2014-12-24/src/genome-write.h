/* $Id: genome-write.h 132144 2014-04-02 16:02:28Z twu $ */
#ifndef GENOME_WRITE_INCLUDED
#define GENOME_WRITE_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "types.h"

extern void
Genome_write_comp32 (char *genomesubdir, char *fileroot, FILE *input, 
		     Univ_IIT_T contig_iit, IIT_T altstrain_iit, Univ_IIT_T chromosome_iit,
		     bool uncompressedp, bool rawp, bool writefilep,
		     Univcoord_T genomelength, int index1part, int nmessages);

#endif
