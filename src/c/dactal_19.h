#ifndef DACTAL_H
#define DACTAL_H


typedef enum DIST_MOD{
  D_JC,
  D_LOGDET,
} DIST_MOD;

typedef enum SUPERTREE_M{
  M_FASTRFS,
  M_GREEDYRFS,
  M_SUPERFINE
} SUPERTREE_M;

typedef enum SUBTREE_M{
  M_FASTTREE,
  M_RAXML
} SUBTREE_M;

// Struct to store users' options to the algorithm 
typedef struct dactal_options{
  // Necessary options
  char * in_aln;
  char * out_tree;
  int stopping_t;
  int padding;

  // Other options
  char * tmp_folder; 

  SUPERTREE_M supertree_method;
  SUBTREE_M subtree_method;
  DIST_MOD distance_model;

  int help_mode;

} DACTAL_GRP;

#endif // DACTAL_H