#!/usr/bin/env python3
"""
Modified by Thien Le in 2019 May for HMMDecomposition
"""
"""
Build subsets using centroid decomposition from PASTA

Written by EKM (molloy.erin.k@gmail.com) in Spring 2018.
"""

import copy
from dendropy import Tree
from dendropy import Edge
from dendropy import Node
from dendropy import DataSet as Dataset
from dendropy.datamodel.treemodel import _convert_node_to_root_polytomy as convert_node_to_root_polytomy
import dendropy
import argparse
# from pasta.pastaalignerjob import bisect_tree
# from pasta.tree import PhylogeneticTree
import sys

def get_newick_string(tree):
    """
    Creates a newick string with  branch lengths in decimal notation.
    This function was taken from Phybase.py, see
    github.com/ngcrawford/CloudForest/blob/master/cloudforest/phybase.py

    Parameters
    ----------
    tree : dendropy tree

    Returns
    -------
    new_tree : string
               newick representation of tree
    """
    from decimal import Decimal

    str_newick_tree = tree.as_string(schema='newick')

    colon_s = 0
    comma_back_paren_s = 0
#    num = ''
    new_tree = ''
    for count, char in enumerate(str_newick_tree):
        if char == ':':
            colon_s = count
            continue
        if char in (')', ','):
            comma_back_paren_s = 1
#            num = '%1.12f' % Decimal(num)
#            new_tree += ":" + num
            colon_s = 0
#            num = ''
#        if colon_s != 0:
#            num = num + char
        if colon_s == 0:
            new_tree += char
    new_tree = new_tree.strip('\'').strip('\"').strip('\'') + '\n'
    return new_tree

"""
Copied from PASTA
"""
"""SATe/PASTA - Phylogenetic Tree Container, effectively a wrapper of dendropy.Tree"""

# This file is part of PASTA, and is forked form SATe

# PASTA, like SATe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Jiaye Yu and Mark Holder, University of Kansas

def is_valid_tree(t):
    assert t and t
    rc = t.seed_node.child_nodes()
    num_children = len(rc)
    if num_children == 0:
        return True
    elif num_children == 1:
        assert len(rc[0].child_nodes()) != 0
        return True
    elif num_children == 2:
        assert((not list(rc[0].child_nodes())) and (not list(rc[1].child_nodes())))
    return True

def get_centroid_edge(tree,spanning=False):
    """Get centroid edge"""
    root = tree.seed_node
    root_children = root.child_nodes()
    if root_children and (spanning or not hasattr(root_children[0].edge, "num_leaves_below")):

        #Calc split
        n = len(tree.leaf_nodes())
        for i in tree.postorder_edge_iter():
            nd = i.head_node
            if nd.is_leaf():
                i.num_leaves_below = 1
            else:
                i.num_leaves_below = sum([j.edge.num_leaves_below for j in nd.child_nodes()])


        n_leaves = len(tree.leaf_nodes())
        if spanning and len(root_children) == 1:
            n_leaves += 1
    else:
        if root.edge:
            n_leaves = root.edge.num_leaves_below
        else:
            n_leaves = sum([j.edge.num_leaves_below for j in root_children])
    centroid_edge = None
    centroid_imbalance = n_leaves
    half_taxa = n_leaves/2
    if half_taxa == 0:
        half_taxa = 1
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        n_descendants = edge.num_leaves_below
        imbalance = abs(half_taxa - n_descendants)
        if (imbalance < centroid_imbalance):
            centroid_edge = edge
            centroid_imbalance = imbalance
    assert centroid_edge is not None
    return centroid_edge

def leaf_node_names(tree):
    leaves = tree.leaf_nodes()
    return [i.taxon.label for i in leaves]

def bipartition_by_edge(tree, e, padding):
    """Prunes the subtree that attached to the head_node of edge e and returns them as a separate tree."""
    
    t = tree
    n = len(tree.leaf_nodes())
    assert e.tail_node is not None
    assert e.head_node is not None
    
    nr = e.head_node
    nd = e.tail_node
    assert nr.parent_node is nd

    t.reroot_at_edge(e)
    is_valid_tree(t)

    root = t.seed_node
    res = []

    max_leaf_n = -1

    child_array = root.child_nodes()
    assert len(child_array) == 3

    large_i = -1
    for i in range(len(child_array)):
        if(len(child_array[i].leaf_nodes()) > max_leaf_n):
            max_leaf_n = len(child_array[i].leaf_nodes())
            large_i = i

    for j in range(len(child_array)):
        c = child_array[j]
        c.edge.length = None # Length of bisected edge
        c.parent_node = None
        if j != large_i:
            convert_node_to_root_polytomy(c)
            res.append(Tree(seed_node=c))
        else:
            large_tree_child_array = c.child_nodes()
            assert len(large_tree_child_array) == 2
            for lc in large_tree_child_array:
                lc.edge.length = None # Length of bisected edge
                lc.parent_node = None
                convert_node_to_root_polytomy(lc)
                res.append(Tree(seed_node=lc))

    # for i in range(len(res)):
    #     print(len(res[i].leaf_nodes()))

    assert len(res) == 4
    assert len(res[0].leaf_nodes()) + len(res[1].leaf_nodes()) + len(res[2].leaf_nodes()) + len(res[3].leaf_nodes()) == n

    X = [];
    a_arr = [];
    for i in range(len(res)):
        for leaf in res[i].leaf_nodes():
            a_arr.append((leaf.taxon.label, leaf.distance_from_root()))
        a_arr = sorted(a_arr, key = lambda x: x[1])

        for j in range(min(padding, len(a_arr))):
            X.append(a_arr[j][0])
        a_arr = []

    all_lab = []
    for i in range(len(res)): 
        lab = []
        lab += X
        lab += leaf_node_names(res[i])
        all_lab.append(lab)

    return all_lab


def main(args):
    # Step 1: Decompose tree
    tree = dendropy.Tree.get(path=args.input_tree_file,
                             schema="newick")
    tree.resolve_polytomies(limit=2,
                            update_bipartitions=True)

    e = get_centroid_edge(tree)
    ts = bipartition_by_edge(tree,e,args.padding)

    # Step 2: Write out leaf subsets
    for i in range(len(ts)):
        with open(args.output + str(args.indexing + i), "w") as f:
            f.write("\n".join(str(x) for x in ts[i]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-t", "--input_tree_file", type=str,
                        required=True,
                        help="Input tree file")

    parser.add_argument("-o", "--output", type=str,
                        required=True,
                        help="Output prefix")

    parser.add_argument("-i", "--indexing", type=int,
                        required=True,
                        help="Indexing")

    parser.add_argument("-p", "--padding", type=int,
                        required=True,
                        help="Padding")

    main(parser.parse_args())
