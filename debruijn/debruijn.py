#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
#from operator import itemgetter
import random
#from random import randint
import statistics
import networkx as nx
#import matplotlib
random.seed(9001)

__author__ = "Imane MESSAK"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Imane MESSAK"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Imane MESSAK"
__email__ = "imane.messak1@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()



#1. Création du graphe de de Bruijn

def read_fastq(fastq_file):
    """Read a fastq file and return by yield sequences
    Parameters : the fastq file"""
    with open(fastq_file, "r") as fastq:
        for _ in fastq:
            yield next(fastq).strip()
            next(fastq)
            next(fastq)

def cut_kmer(seq, kmer_size):
    """Read all sequences and return kmers
    Parameters: seq :sequence to cut, kmer_size :size of kmers
    Returns : kmers"""
    for i in range(len(seq)-kmer_size+1):
        yield seq[i:i+kmer_size]

def build_kmer_dict(fastq_file, kmer_size):
    """ Build a dictionnary of kmers
    Parameters : fastq_file :fastq file, kmer_size:size of kmers
    Returns: dictionnary of kmers with their occurences """
    l_seq = []
    dic_kmer = {}
    for i in read_fastq(fastq_file):
        l_seq.append(i)
    for j in range(len(l_seq)):
        l_kmer = []
        for i in cut_kmer(l_seq[j], kmer_size):
            l_kmer.append(i)
        for i in range(len(l_kmer)):
            if l_kmer[i] not in dic_kmer:
                dic_kmer[l_kmer[i]] = l_kmer.count(l_kmer[i])
    return dic_kmer

def build_graph(dic_kmer):
    """ Build a tree of kmers prefixes and suffixes
    Parameters : dic_kmer : dictionnary of kmers
    Returns : tree of kmers prefixes and suffixes"""
    tree_graph = nx.DiGraph()
    for key, val in dic_kmer.items():
        tree_graph.add_edge(key[:-1], key[1:], weight=val)
    return tree_graph



#2. Parcours du graphe de de Bruijn

def get_starting_nodes(tree_graph):
    """ Construct a list of starting nodes from the graph
    Parameters : tree_graph : tree of prefixes and suffixes kmers
    Returns : list of starting nodes"""
    start = []
    for node in tree_graph.nodes:
        if not list(tree_graph.predecessors(node)):
            start.append(node)
    return start

def get_sink_nodes(tree_graph):
    """ Construct a list of sink nodes from the tree_graph
    Parameters : tree_graph
    Returns : list of sink nodes """
    sink = []
    for node in tree_graph.nodes:
        if not list(tree_graph.successors(node)):
            sink.append(node)
    return sink

def get_contigs(tree_graph, start, sink):
    """Get the different possible contigs
    Paramters: tree_graph, start: list of starting nodes,sink: list of sink nodes
    Returns: A list of tuples with contig and length of contig"""
    contigs = []
    for start_node in start:
        for sink_node in sink:
            for path in nx.all_simple_paths(tree_graph, start_node, sink_node):
                contig = path[0]
                for i in range(1, len(path)):
                    contig += path[i][-1]
                contigs.append((contig, len(contig)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(tuple_list, output_file):
    """Save the contigs in a text file with fasta form
    Parameters: tuple_list :list of tuples with contig and its length, output_file: output file
    Returns : fasta file with contigs inside"""
    with open(os.path.join('../Results/', output_file), "w+") as filout:
        for i, contig in enumerate(tuple_list):
            filout.write(">contig_"+str(i)+" len="+str(contig[1])+"\n")
            filout.write(fill(contig[0])+"\n")
        filout.close()



#3. Simplification du graphe de de Bruijn

def std(list_val):
    """Calculate the std
    Parameters : list_val: liste of integer
    Returns: the standard deviation of the list of integer"""
    return statistics.stdev(list_val)

def path_average_weight(graphe, path):
    """Caclul the average weights for all paths
    Parameters : graphe, paths : list of paths
    Returns : Mean of the weight for each paths"""
    avg_weights = 0
    for i in range(len(path)-1):
        avg_weights += graphe[path[i]][path[i+1]]["weight"]
        j = i
    return avg_weights/(j+1)

def remove_paths(graphe, paths, delete_entry_node=False, delete_sink_node=False):
    """ Remove a path containing starting or ending nodes
    Parameters : graphe, paths : list of paths, two booleen parameters
    Returns : graphe without indesirable paths"""
    for path in paths:
        if delete_entry_node and delete_sink_node is not True:
            graphe.remove_nodes_from(path[:-1])
        elif delete_entry_node and delete_sink_node:
            graphe.remove_nodes_from(path)
        elif delete_entry_node is not True and delete_sink_node:
            graphe.remove_nodes_from(path[1:])
        elif delete_entry_node is not True and delete_sink_node is not True:
            graphe.remove_nodes_from(path[1:-1])
    return graphe

def select_best_path(graphe, paths, path_len_m, path_weight_m,
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between all paths
    Parameters : graphe, paths : list of paths,
            path_len_m: list of paths lenght
            path_weight_m: list of paths weights
            two boleen parameters
    Returns : graph with only best path"""
    cpath = []
    for path_i in range(len(paths)):
        for path_j in range((path_i+1), len(paths)):
            std_len = std([path_len_m[path_i], path_len_m[path_j]])
            std_weight = std([path_weight_m[path_i], path_weight_m[path_j]])
            if std_len == 0 and std_weight == 0:
                rdm = random.choice([path_i, path_j])
                cpath.append(paths[rdm])
            elif std_len != 0 and std_weight == 0:
                if path_len_m[path_i] > path_len_m[path_j]:
                    cpath.append(paths[path_j])
                else:
                    cpath.append(paths[path_i])
            else:
                if path_weight_m[path_i] > path_weight_m[path_j]:
                    cpath.append(paths[path_j])
                else:
                    cpath.append(paths[path_i])
    remove_paths(graphe, cpath, delete_entry_node, delete_sink_node)
    return graphe

def solve_bubble(graphe, ancestor_node, descendant_node):
    """Remove one bubble in graphe
    Parameters :graphe, ancestor node, descendant node
    Returns: graphe without bubble between these nodes"""
    paths_l = list(nx.all_simple_paths(graphe, ancestor_node, descendant_node))
    path_len = []
    path_weight = []
    for path in paths_l:
        path_len.append(len(path))
        path_weight.append(path_average_weight(graphe, path))
    select_best_path(graphe, paths_l, path_len, path_weight)
    return graphe

def simplify_bubbles(graphe):
    """Remove bubbles in graphe
    Parameters : graphe
    Returns: graph without bubbles"""
    bubbles = []
    for node in graphe.nodes():
        predecessors_list = list(graphe.predecessors(node))
        if len(predecessors_list) > 1:
            ancestor = nx.lowest_common_ancestor(graphe, predecessors_list[0],
                                                 predecessors_list[1])
            bubbles.append([ancestor, node])
    for i in range(len(bubbles)):
        graphe = solve_bubble(graphe, bubbles[i][0], bubbles[i][1])
    return graphe

def solve_entry_tips(graphe, list_node_in):
    """Remove all unwanted starting nodes
    Parameters: graph, list of entry nodes
    Returns: graph without indesirable entry path"""
    path_l = []
    path_len = []
    path_weight = []

    for node in list_node_in:
        for des_node in nx.descendants(graphe, node):
            pred_node = list(graphe.predecessors(des_node))
            if len(pred_node) > 1:
                for path in nx.all_simple_paths(graphe, node, des_node):
                    path_l.append(path)
                    path_len.append(len(path))
                    path_weight.append(path_average_weight(graphe, path))
    graphe = select_best_path(graphe, path_l, path_len, path_weight, True, False)

    return graphe

def solve_out_tips(graphe, list_node_out):
    """Remove all unwanted sinking nodes
    Parameters: graph, list of sink nodes
    Returns: graph without indesirable sink nodes"""
    path_l = []
    path_len = []
    path_weight = []

    for node in list_node_out:
        for anc_node in nx.ancestors(graphe, node):
            succ_node = list(graphe.successors(anc_node))
            if len(succ_node) > 1:
                for path in nx.all_simple_paths(graphe, node, anc_node):
                    path_l.append(path)
                    path_len.append(len(path))
                    path_weight.append(path_average_weight(graphe, path))

    graphe = select_best_path(graphe, path_l, path_len, path_weight, True, False)

    return graphe

#==============================================================
# Main program
#==============================================================
def main():
    """Main program function"""

    args = get_arguments()
    #Lecture du fichier et construction du graphe:
    dic_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(dic_kmer)
    start_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    #Résolution des bulles:
    graph = simplify_bubbles(graph)
    #Resolution des points d'entrée et de sortie
    graph = solve_entry_tips(graph, start_nodes)
    graph = solve_out_tips(graph, sink_nodes)
    #Ecriture des contigs:
    contigs = get_contigs(graph, start_nodes, sink_nodes)
    save_contigs(contigs, args.output_file)

if __name__ == '__main__':
    main()
