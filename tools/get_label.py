#!/usr/bin/env python3

"""
Modified by Thien Le in 2019 May for HMMDecomposition
"""
"""
Build subsets using centroid decomposition from PASTA

Written by EKM (molloy.erin.k@gmail.com) in Spring 2018.
"""
import dendropy
import argparse
import sys

def leaf_node_names(tree):
    leaves = tree.leaf_nodes()
    return [i.taxon.label for i in leaves]

def main(args):
    # Step 1: Decompose tree
    tree = dendropy.Tree.get(path=args.input_tree_file,
                             schema="newick")
    tree.resolve_polytomies(limit=2,
                            update_bipartitions=True)

    with open(args.output_file, "w") as f:
        f.write("\n".join(leaf_node_names(tree)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-t", "--input_tree_file", type=str,
                        required=True,
                        help="Input tree file")
    parser.add_argument("-o", "--output_file", type=str,
                        required=True,
                        help="Output file")

    main(parser.parse_args())
