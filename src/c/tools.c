// File in inc_ml, created by Thien Le in July 2018

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>

#include "tools.h"
#include "utilities.h"
#include "dactal_19.h"


// Private functions
int find_prefix_and_dir_from_path(char *, char *, char *, char *);
int set_up_dc_node(
    int * node_counter, 
    char ** naming_subproblem, 
    char * out_tree_suffix,
    int * subproblem_counter
);
int centroid_decomposition_job(
    char * in_tree,
    char * out_tree_prefix,
    int node_counter,
    int padding);
int check_tree_size_job(char * in_tree, int stopping_t, int * stop);
int make_constraint_trees(
    char * decomposed_leaves_path, 
    DACTAL_GRP * dactal_options,
    // char * in_aln,
    char * out_tree);
    // SUBTREE_M subtree_method);
int combine_trees(
    char * out_tree_suffix, 
    char ** sub_problem_out_tree_suffix,
    SUPERTREE_M supertree_method);

// Public implementation
int fasttree_job(DIST_MOD dist_model, char * in_aln, char * out_path)
{
  SYSCAL_STD(
      "%s -nt %s -quiet < %s > %s",
      FastTree_bin,
      dist_model == D_JC ? FT_JC : FT_GTRGAMMA,
      in_aln,
      out_path
  );
  return 0;
}

int recursive_subroutine(
    int current_counter,
    char * tree_suffix,
    DACTAL_GRP * dactal_options,
    int * node_counter,
    char * out_tree_suffix
){  
  int i;
  int stop = -1; 
  int subproblem_counter[4];
  char ** sub_problem_out_tree_suffix;

  // All string local variables
  char in_tree[GENERAL_BUFFER_SIZE];
  char out_prefix[GENERAL_BUFFER_SIZE];

  // Clear all strings
  in_tree[0] = 0;
  out_prefix[0] = 0;
  sub_problem_out_tree_suffix = SAFE_MALLOC(4 * sizeof(char*));
  for(i = 0; i < 4; ++i){
    sub_problem_out_tree_suffix[i] = SAFE_MALLOC(GENERAL_BUFFER_SIZE);
    sub_problem_out_tree_suffix[i][0] = 0;
  }

  // Build all strings local variables
  if(*node_counter == 1)
    sprintf(in_tree, "%s", tree_suffix);
  else{
    sprintf(in_tree, "%s%d", tree_suffix, current_counter);
  }

  strcat(out_prefix, tree_suffix);

  // Check if the current tree is small enough
  FCAL_STD(check_tree_size_job(in_tree, dactal_options->stopping_t, &stop));
  if(stop < 0){
    printf("checking tree size job does not write to 'stop', something is wrong\n");
  } else if(stop) {
    // Do things at the leaves
    FCAL_STD(make_constraint_trees(
        &(in_tree[0]),
        dactal_options,
        // dactal_options->in_aln,
        out_tree_suffix
        // dactal_options->subtree_method
    ));
    return 0;
  }

  // Centroid decomposition
  FCAL_STD(centroid_decomposition_job(
      in_tree, 
      out_prefix, 
      *node_counter,
      dactal_options->padding
  ));

  // Set up the DC tree properly
  FCAL_STD(set_up_dc_node(
      node_counter, 
      sub_problem_out_tree_suffix, 
      tree_suffix,
      &(subproblem_counter[0])
  ));

  // Recurse in a preorder fashion
  for(i = 0; i < 4; i++){
    FCAL_STD(
        recursive_subroutine(
            subproblem_counter[i],
            tree_suffix,
            dactal_options,
            node_counter,
            sub_problem_out_tree_suffix[i]
        )
    );
  }

  // Combine results of subproblem
  FCAL_STD(combine_trees(
      out_tree_suffix, 
      sub_problem_out_tree_suffix,
      dactal_options->supertree_method));

  return 0;
}

// int rm_label_job(char * in_tree, char * out_path)
// {
//   SYSCAL_STD(
//       "%s -t %s -o %s", 
//       rm_lbl_bin,
//       in_tree, 
//       out_path
//   );
//   return 0;
// }

int raxml_job(
    DIST_MOD dist_model, 
    char * in_seq,
    char * out_pfx)
{
  // if(wd_sfx)
  //   SYSCAL_STD(
  //       "%s -T 12 --silent -m %s -j -n %s -s %s -w %s -p 1 > %s", 
  //       RAxML_bin, 
  //       dist_model == D_JC ? RAXML_JC : RAXML_GTRGAMMA,
  //       out_pfx,
  //       in_seq,
  //       wd_sfx,
  //       TMP_FILE1
  //   );
  // else
    SYSCAL_STD(
        "%s -T 12 --silent -m %s -j -n %s -s %s -p 1 > %s",    
        RAxML_bin, 
        dist_model == D_JC ? RAXML_JC : RAXML_GTRGAMMA,
        out_pfx, 
        in_seq,
        TMP_FILE1
    );

  SYSCAL_STD("rm %s", TMP_FILE1); 
  return 0; 
}

// Internal functions

int find_prefix_and_dir_from_path(
    char * path, 
    char * prefix, 
    char * dir, 
    char * name)
{                              
  int i, j;
  // Malloc sequence, assuming path is created
  ASSERT_INPUT(path);

  // Init
  STR_CLR(prefix);
  STR_CLR(dir); 

  // Find the last backslash, this separates the dir from the prefix
  for(i = strlen(path) - 1; i >= 0; i--)
    if(path[i] == '/') break;
  
  // Do a string copy
  for(j = 0; j < strlen(path); j++){
    if(j <= i) dir[j] = path[j];
    else prefix[j - i - 1] = path[j];
  }
  dir[i + 1] = 0;
  prefix[strlen(path) - i - 1] = 0;

  sprintf(name, "%sRAxML_bestTree.%s", dir, prefix);

  return 0; 
}


int set_up_dc_node(
    int * node_counter, 
    char ** naming_subproblem, 
    char * out_tree_suffix,
    int * subproblem_counter 
){
  int i;
  for(i = 0; i < 4; ++i){
    naming_subproblem[i][0] = 0;
    sprintf(
        naming_subproblem[i], 
        "%ssubproblem%d", 
        out_tree_suffix, 
        *node_counter
    ); 
    subproblem_counter[i] = *node_counter;
    (*node_counter)++;

  }

  // printf("current node is %d, children are %d %d %d %d\n", cur_node->node, cur_node->children[0]->node, cur_node->children[1]->node, cur_node->children[2]->node, cur_node->children[3]->node);

  return 0; 
}

int centroid_decomposition_job(
    char * in_tree,
    char * out_tree_prefix,
    int node_counter,
    int padding
){
  int i;
  SYSCAL_STD(
      "%s -t %s -o %s",
      remove_internal_labels_bin,
      in_tree,
      in_tree  
  );

  SYSCAL_STD(
      "%s -t %s -o %s.lab -i %d -p %d",
      centroid_decomposition_bin, 
      in_tree,
      out_tree_prefix,
      node_counter,
      padding
  );
  for(i = 0; i < 4; i++){
    SYSCAL_STD(
        "%s -v %s $(cat %s.lab%d) > %s%d",
        newick_prune_bin,
        in_tree,
        out_tree_prefix, 
        node_counter + i,
        out_tree_prefix, 
        node_counter + i
    );
  }
  return 0;
}

int check_tree_size_job(char * in_tree, int stopping_t, int * stop){
  char n_leaves_file[GENERAL_BUFFER_SIZE];
  FILE * f;
  int n_leaves;

  SYSCAL_STD(
      "%s -t %s",
      check_tree_bin,
      in_tree
  );

  n_leaves_file[0] = 0;
  strcat(n_leaves_file, in_tree);
  strcat(n_leaves_file, "_n_leaves");

  f = SAFE_FOPEN_RD(n_leaves_file);

  fscanf(f, "%d", &n_leaves);

  if(n_leaves <= stopping_t) *stop = 1;
  else *stop = 0;

  return 0;
}

int make_constraint_trees(
    char * decomposed_leaves_path, 
    DACTAL_GRP * dactal_options,
    // char * in_aln ,
    char * out_tree
    // SUBTREE_M subtree_method
){
  FILE *p,*f,*l;
  char taxon_name[GENERAL_BUFFER_SIZE];
  char buffer[GENERAL_BUFFER_SIZE];
  char buffer2[GENERAL_BUFFER_SIZE];
  char tree_lab_path[GENERAL_BUFFER_SIZE];
  char constraints_lab_path[GENERAL_BUFFER_SIZE];
  char c;

  sprintf(tree_lab_path, "%stmp_tree_lab", out_tree);
  sprintf(constraints_lab_path, "%stmp_constraints_lab", out_tree);

  SYSCAL_STD(
      "%s -t %s -o %s", 
      get_label_bin, 
      decomposed_leaves_path, 
      tree_lab_path);

  f = fopen(constraints_lab_path, "w");

  l = SAFE_FOPEN_RD(tree_lab_path);
  while(fscanf(l,"%s", taxon_name) >= 0){
    sprintf(buffer2, ">%s", taxon_name);
    fprintf(f, "%s", buffer2);
    p = SAFE_FOPEN_RD(dactal_options->in_aln);
    while(fscanf(p, "%s", buffer) >= 0){
      if(strcmp(buffer,buffer2) == 0){
        while(fscanf(p, "%c", &c) >= 0 && c != '>'){
          fprintf(f, "%c",c);
        }
        break;
      }
    }
    fclose(p);
  } 
  fclose(l);
  fclose(f);

  switch(dactal_options->subtree_method){
    case M_FASTTREE:
      FCAL_STD(fasttree_job(
          dactal_options->distance_model, 
          constraints_lab_path, 
          out_tree
      ));
      break;
    case M_RAXML:
      FCAL_STD(raxml_job(
          dactal_options->distance_model, 
          constraints_lab_path, 
          out_tree
      ));

      SYSCAL_STD(
          "mv RAxML_result.%s %s",
          out_tree,
          out_tree
      );
      break;
  }
  
  return 0;
}

int combine_trees(
    char * out_tree_suffix, 
    char ** sub_problem_out_tree_suffix,
    SUPERTREE_M supertree_method
){
  char tmp_supertree_in_path[GENERAL_BUFFER_SIZE];
  sprintf(tmp_supertree_in_path, "%stmp_supertree_in", out_tree_suffix);

  SYSCAL_STD(
      "cat %s %s %s %s > %s",
      sub_problem_out_tree_suffix[0],
      sub_problem_out_tree_suffix[1],
      sub_problem_out_tree_suffix[2],
      sub_problem_out_tree_suffix[3],
      tmp_supertree_in_path
  );

  switch(supertree_method){
    case M_GREEDYRFS:
      SYSCAL_STD(
          "%s -t %s -o %s",
          greedy_rfs_bin,
          tmp_supertree_in_path,
          out_tree_suffix
      );
      break;
    case M_FASTRFS:
      SYSCAL_STD(
          "%s -i %s -o %s",
          fastrfs_bin,
          tmp_supertree_in_path,
          out_tree_suffix
      );
      SYSCAL_STD(
          "mv %s.single %s",
          out_tree_suffix,
          out_tree_suffix
      );
      break;
  }

  return 0;
}
