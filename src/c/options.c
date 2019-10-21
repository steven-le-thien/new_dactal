// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "options.h"
#include "utilities.h"

// Internal implementations
int parse_dactal_arg(
    char ** argv, 
    int flag_index, 
    int content_s, 
    int content_e, 
    DACTAL_GRP * dactal_options);

/* This function takes a point to some dactal_options, init its fields and 
 *    set the corresponding values depending on the user input
 * It will also catch errors in incompatible calls
 * Input:       user inputs and an empty dactal_options struct
 * Output:      0 on success, ERROR otherwise
 * Effect:      popularizes fields in dactal_options, may call malloc on stuff
 */ 
int read_dactal_cmd_arg(
      int argc, 
      char ** argv, 
      DACTAL_GRP * dactal_options)
{
  int i;
  int is_flag;
  int flag_index;
  int content_s, content_e; // these are supposed to be left close, right open 

  // Initializing the options to default fields (which assumes that nothing has 
  // been computed before and this is the very first run)
  FCAL_STD(init_dactal_options(dactal_options)); 

  ASSERT_INPUT(
      argv[1][0] == '-'
  );

  // State machine to read in user input
  is_flag = 1; // always start with a flag
  flag_index = 1;
  for(i = 2; i < argc; i++){
    if(is_flag){  // was reading a flag from the previous round, actual flag is 
                  // stored in flag_index
      ASSERT_INPUT(
          argv[1][0] == '-'
      );
      

      // reading the first content, note that there may be a lot more contents 
      // to read
      content_s = i;
      is_flag = 0;

    } else { // have been reading content
      if(argv[i][0] == '-'){ // read a flag here, so no more content to read
        content_e = i;
        is_flag = 1;

        // Read the previous flag into option
        FCAL_STD(
            parse_dactal_arg(
                argv, 
                flag_index, 
                content_s, 
                content_e, 
                dactal_options)
        );
      
        // Reset the machine
        flag_index = i;
        content_s = -1;
        content_e = -1;
      } 
      // Else, reading another content, so we don't need to do anything
    }
  }

  // Read last flag
  FCAL_STD(
      parse_dactal_arg(argv, flag_index, content_s, argc, dactal_options)
  );

  return 0; 
}

int parse_dactal_arg(char ** argv, 
                  int flag_index, 
                  int content_s, 
                  int content_e, 
                  DACTAL_GRP * dactal_options)
{
  char * flag = argv[flag_index];
  // Switch case for the flag index, which should be a single character, 
  // following a dot. 
  // First check that there is exactly one character that follows the dash 
  ASSERT_INPUT(
      strlen(flag) == 2 && 
          flag[0] == '-' && 
          flag[1] >= 'A' && 
          flag[1] <= 'z' &&
          content_e - content_s == 1
  );


  // Actually look up what exact field show we be setting and set said field
  // -i: input alignment              (expect the full path - does not check)
  // -o: output tree                  (expect the full path - does not check)
  // -t: tmp_folder...................(expect the full path - does not check)
  // -m: stopping criteria for CD
  switch(flag[1]){   
    case 'i':
        dactal_options->in_aln = argv[content_s];
        break;
    case 'o':
        dactal_options->out_tree = argv[content_s];
        break;
    case 't':
        dactal_options->tmp_folder = argv[content_s];
        break; 
    case 'm':
        dactal_options->stopping_t = atoi(argv[content_s]);
        break;
    case 'p':
        dactal_options->padding = atoi(argv[content_s]);
        break;
    case 's':
        if(strcmp(argv[content_s], "fastrfs") == 0)
          dactal_options->supertree_method = M_FASTRFS;
        else if(strcmp(argv[content_s], "greedyrfs") == 0)
          dactal_options->supertree_method = M_GREEDYRFS;
    case 'b':
        if(strcmp(argv[content_s], "fasttree") == 0)
          dactal_options->subtree_method = M_FASTTREE;
        else if(strcmp(argv[content_s], "raxml") == 0)
          dactal_options->subtree_method = M_RAXML;
    case 'd':
        if(strcmp(argv[content_s], "jc") == 0)
          dactal_options->distance_model = D_JC;
        else if(strcmp(argv[content_s], "logdet") == 0)
          dactal_options->distance_model = D_LOGDET;
  }
  return 0;
}


int init_dactal_options(DACTAL_GRP * dactal_options){
  // The first two should be initialized properly
  dactal_options->in_aln         = NULL;
  dactal_options->out_tree           = NULL;
  dactal_options->tmp_folder          = NULL;
  return 0;
}
