// File in inc_ml, created by Thien Le in July 2018

#ifndef OPTION_H
#define OPTION_H

#include "dactal_19.h"

// Functions
extern int read_dactal_cmd_arg(int argc,char ** argv, DACTAL_GRP * options);
extern int init_dactal_options(DACTAL_GRP * options);

#endif