### DACTAL (2019 implementation)

DACTAL is described in `Nelesen, S., Liu, K., Wang, L. S., Linder, C. R., & Warnow, T. (2012). DACTAL: divide-and-conquer trees (almost) without alignments. Bioinformatics (Oxford, England), 28(12), i274–i282. doi:10.1093/bioinformatics/bts218`. It is a divide-and-conquer based method for the problem of estimating the model gene tree given a set of alignment.

The implementation design can be summarized as. 1/ Build a fast tree on the full set of alignment; 2/ Decompose the set of alignment into multiple overlapping subalignments with small size using padded centroid decomposition of the initial tree; 3/ Build trees on the small overlapping subalignments; 4/ Combine the subtrees using a supertree method in the same order they were decomposed

## Requirement
DACTAL was built from standard C packages with python3 scripts
1. Make
2. GCC or other C compile (code was developed and tested with GCC) 
3. Python 3.

If you want to use different methods for subset tree computations and supertree then you must also have them on your PATH variable
1. FastTree2 (http://www.microbesonline.org/fasttree/FastTree.c)
2. RAxML (only if needed to run RAxML)
3. Newick Utils (https://github.com/tjunier/newick_utils)
4. Extra scripts in the tools folder
5. FastRFS (only if needed to run FastRFS)
6. GreedyRFS (only if needed to run GreedyRFS)
7. SuperFine (only if needed to run SuperFine) 

Furthermore, the binaries must be named exactly as
1. FastTree2: `FastTree`
2. RAxML: `raxmlHPC-PTHREADS-AVX`
3. FastRFS: `FastRFS`
4. GreedyRFS: `GreedyRFS.py`
5. SuperFine: `runSuperFine.py`
6. All python scripts in the tools folder must have the same name as that in the distribution

If for some reason you would prefer to use another name for any of the binaries, then modify the corresponding string in `tools.h` and recompile with `make clean; make dactal`.  

You will need at least one subset tree computation method and one supertree method installed. 

Make sure all dependencies are on your PATH variable. One quick and dirty way to do this is to type
```
export PATH=$PATH:/path/to/folder/or/file
```
where `/path/to/folder/or/file` is the absolute path to the dependency or the folder containing the dependency file. 

## Installation
Some compilation flag current does not compile with `clang`. Some machines install `clang` under the alias of `gcc` and thus will mess up the current Makefile. To know whether this happens to your machine, try to run the compilation steps below. If you get an error message that comes from `clang` (eg: `clang: error: unsupported option '-fopenmp'`),  do the following steps to use the correct compiler:
1. Install the lastest version of `gcc`. For example, at the time of writing this README, the newest version is `gcc-8`.  
2. Type `export CC=$(which gcc-8)` where `gcc-8` should be replaced with the name of the newest version of `gcc` that you have just installed. 

These steps will affect the `CC` variable in the environment that is used to link the correct compiler in the Makefile

Run `make dactal` to generate the binary `dactal`, used for DACTAL. Also put this binary on your PATH variable. The command for DACTAL is 
```
dactal -i [input_alignment_file] -o [output_tree_file] \
-m [maximum_subset_size] -p [overlapping_size_per_subtree] \
-s [supertree_method] -b [subtree_method] -d [distance_model]
```
All flags are currently mandatory.

## DACTAL options
1. `-i` specifies input alignment. Please specify the full path.
2. `-o` specifies output tree. Only specify the filename (and not the full path) if you are using RAxML for the subtree computation. Any file with the same name will be overwritten.
3. `-m` accepts a positive number that approximates (upperbounds) the subtree size
4. `-p` accepts a positive number that approximates (upperbounds) the padding size. This is per subtree in each step of centroid-edge decomposition. 
5. `-s` accepts either `superfine`,`fastrfs` or `greedyrfs`, which specifies the supertree method used to combine subtrees. 
6. `-b` accepts either `fasttree` or `raxml`, which specifies the subtree method used to build the decomposed subtrees.
7. `-d` accepts either `logdet` or `jc`, which are the different distance models.

## DACTAL outputs
1. The final tree, as specified in `-o`. 
2. The initial tree, at `[output_tree_file]_starting_tree`.
3. All immediate files.

## Format
2. All trees (input / output) are in Newick.
3. The input alignment for `dactal` is in FASTA.

## License and Copyright
See the attached LICENSE and COPYRIGHT file for details. Note that some part of the program (scripts in the tools folder) was distributed under a different license

## Contact
The package is under constant development. Please contact `thienle@mit.edu` for implementation questions/suggestions for the C code. 
