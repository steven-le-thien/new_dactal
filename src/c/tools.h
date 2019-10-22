// File in inc_ml, created by Thien Le in July 2018

#ifndef TOOLS_H
#define TOOLS_H

#include "dactal_19.h"

extern int make_raxml_constraint(
    char * msa_name, 
    char * out_name, 
    DACTAL_GRP * dactal_options,
    int clear_lab
);

extern int fasttree_job(
    DIST_MOD dist_model, 
    char * in_aln, 
    char * out_path
); 

extern int recursive_subroutine(
    int current_counter,
    char * tree_suffix,
    DACTAL_GRP * dactal_options,
    int * node_counter,
    char * out_tree_suffix
);


extern int raxml_job(
    DIST_MOD dist_model, 
    char * in_seq,    
    char * out_pfx);

extern int rm_label_job(char * in_tree, char * out_path);

// Settings
static const char RAXML_GTRGAMMA[] 
  = "GTRGAMMA";
static const char RAXML_JC[]      
  = "GTRCAT -V --JC69";
static const char FT_GTRGAMMA[]     
  = "-gtr -gamma";
static const char FT_JC[]           
  = " ";   
static const char NJ_LOGDET[]     
  = "logDet";
static const char NJ_JC[]           
  = "JC";
static const char IN_LOGDET[]     
  = "logDet";
static const char IN_JC[]           
  = "JC"; 

// Bins
#if 1
static const char centroid_decomposition_bin[]
  = "centroid_decomposition.py";
static const char check_tree_bin[]  
  = "check_tree.py";
static const char FastTree_bin[]        
  = "FastTree"; 
static const char get_label_bin[]
  = "get_label.py";
static const char newick_prune_bin[]
  = "nw_prune";
static const char remove_internal_labels_bin[]
  = "remove_internal_labels.py";
static const char greedy_rfs_bin[]
  = "GreedyRFS.py";
static const char fastrfs_bin[]
  = "FastRFS";
static const char RAxML_bin[]
  = "raxmlHPC-PTHREADS-AVX";
static const char superfine_bin[]
  = "runSuperFine.py";
#else
static const char RAxML_bin[]    = "raxmlHPC-PTHREADS-SSE3";
static const char FastTree_bin[] = "FastTree";
static const char PAUP_bin[]     = "paup4a163_centos64";
static const char constraint_inc_bin[]    = "constraint_inc";
static const char tree_to_dist_bin[]    = "tree_to_dist.py";
static const char fasta_to_phylip_bin[]   = "fasta_to_phylip.py";
static const char fastme_bin[]            = "fastme";
static const char run_upp_bin[]         = "run_upp.py";
static const char build_subsets_bin[]     = "build_subsets_from_tree.py";
static const char nw_prune_bin[]        = "nw_prune";
static const char rm_lbl_bin[]          = "remove_internal_labels.py";
#endif

#endif
