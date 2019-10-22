#include "stdlib.h"
#include "stdio.h"
#include "string.h"


#include "utilities.h"
#include "options.h"
#include "tools.h"

// New DACTAL implementation as of August 2019
// Input:   a set of aligned sequences
// Output:  a tree
int main(int argc, char ** argv){
  DACTAL_GRP dactal_options; 

  int node_counter = 1;

  char starting_tree_path[GENERAL_BUFFER_SIZE]; 
  starting_tree_path[0] = 0;

  printf("Starting DACTAL... Input is an alignment and output is a tree\n");
  printf("Reading in options...\n");
  FCAL_STD(
      read_dactal_cmd_arg(argc, argv, &dactal_options) 
  );

  if(dactal_options.help_mode){
    printf("Usage: dactal -i [input_alignment_file] -o [output_tree_file] "
        "-m [maximum_subset_size] -p [overlapping_size_per_subtree] "
        "-s [supertree_method] -b [subtree_method] -d [distance_model]\n");
    printf("OR dactal -h\n");
    printf("Default: SuperFine + FastTree + JC\n");

    return 0;
  }

  //DEBUG
  printf("Input alignment is %s\n", dactal_options.in_aln);
  printf("Output tree is %s\n", dactal_options.out_tree);
  //End DEBUG

  printf("Running FastTree on the alignment input...\n"); 

  strcpy(starting_tree_path, dactal_options.out_tree);
  starting_tree_path[strlen(dactal_options.out_tree)] = 0; 
  strcat(starting_tree_path, "_starting_tree");
  FCAL_STD(
      fasttree_job(
          D_JC,
          dactal_options.in_aln,
          starting_tree_path
      )
  );

  printf("Doing main recursion loop\n");
  FCAL_STD(
      recursive_subroutine(
          0,
          starting_tree_path, 
          &dactal_options,
          &node_counter,
          dactal_options.out_tree
      )
  );

  // preorder(dc_root);


  return 0;
}