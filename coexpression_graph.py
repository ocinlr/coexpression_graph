#!/usr/bin/env python
# coding: utf-8

# coexpression_graph - extraction of a coexpression network
# Nicolas Lopez-Rozo, Apr 01 2022

"""
.. _Module-coexpression:

**********************************************
**Module for Coexpression Network Extraction**
**********************************************

This module takes as input several gene profiles and calculates an undirected
graph based on the coexpression between them, calculated using one of the following metrics:

- Pearson correlation coefficient ('PCC')
- Biweight Midcorrelation ('BICOR')
- Distance correlation ('DCORR')


.. _algorithm:

General algorithm
=================

The general idea behind the module can be summarized as follows:

1. Read expression file and apply log_10 if required
2. Calculate coexpression metric (in parallel using multiprocessing)
3. Calculate bins for quantification
4. Generating auxiliary plots (optional)
5. Finding threshold based on density and marginal addition of nodes
6. Apply threshold to find coexpression graph


.. _input-file:

Reading the input file
======================

The input file for *coexpression_graph* is loaded using *pandas*, which infers automatically if 
the file has headers. However, there is a possibility to specify the file delimiter using the parameter
**sep**, which by default is *None*.

It is also optional to apply a logarithm (base 10) to the expression values before calculating the correlation
metric. The corresponding parameter for this effect is **lognormalize_input**. By default, this option is set to *False*.


.. _threshold:

Calculating the Threshold
=========================

The selection of the threshold is calculated as follows:


1. The selected metric is calculated for all posible pairs of gene profiles.

2. A number of bins is created and the network properties are quantified using those bins. The size of each bin is 0.01, ranging from the minimum to the maximum values calculated previously.

3. Two bins are of importance: the one where the minimum of network density occurs, and the one which has a higher increase in the number of nodes included in the network.

4. The minimum of these two bin values is selected as the threshold.

5. Those pairs whose metric value is greater than or equal to the threshold belong to the coexpression graph. 


The idea behind the threshold selection is based on the following papers:

- van Verk, M. C., Bol, J. F., & Linthorst, H. J. (2011). Prospecting for genes involved in transcriptional regulation of plant defenses, a bioinformatics approach. *BMC plant biology* , 11, 88. https://doi.org/10.1186/1471-2229-11-88
- Zhang, L., Yu, S., Zuo, K., Luo, L., & Tang, K. (2012). Identification of gene modules associated with drought response in rice by network-based analysis. *PloS one* , 7(5), e33748. https://doi.org/10.1371/journal.pone.0033748


.. _save-files:

Saving auxiliary files
======================

It is possible to save some files which can help in visualizing properties such as density, 
number of nodes, number of edges, etc. The parameter **save_files** is set by default to *None*,
but can be set to a string indicating the path where the files should be stored.

The stored files are:

- **Plots**: The saved plots as a function of the metric value are the size of the greatest connected component (*{metric}-CCmaxsize.png*), the network density (*{metric}-Density.png*), the number of edges (*{metric}-Edges.png*) and the number of nodes (*{metric}-Nodes.png*). 
- **Metric values**: Files of the form *{metric}-{i}.txt* are generated, including the correlation information.
- **Coexpression graph**: The file *{metric}trimmed.csv* contains the generated network, after the threshold.
- **Bins**: The file *{metric}-bins.csv* can be used to generate the plots again. It includes a header and a line for each bin, specifying the corresponding values.


.. _api_coex:

|

**API Reference**
=================

"""


# General libraries
import os, sys, random, time
from math import *
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import gc

# Libraries for parallel computing
from concurrent.futures import ProcessPoolExecutor
import itertools as itt

# Metrics
import pingouin as pg


##### VARIABLES #####
nodesize = 10       # gets updated in the code
numsets, par = None, None
setsize = None
exp_mat = None

#####################

# data structure: disjoint sets
def _myparent(x):
    global par
    if par[x] == -1:
        return x
    else:
        a = _myparent(par[x])
        par[x] = a
        return a

def _myunion(x, y):
    global par, numsets, setsize
    a, b = _myparent(x), _myparent(y)
    if a != b:
        par[b] = a
        numsets -= 1
        setsize[a] += setsize[b]
    return a != b

def _mysetsize(x):
    return setsize[_myparent(x)]

def _resetDSets():
    global par, numsets, setsize
    numsets, par = nodesize, [-1 for x in range(nodesize)]
    setsize = [1 for x in range(nodesize)]

def _pcc(tup):
    global exp_mat
    a, b = tup
    return a, b, float(pg.corr(exp_mat[a], exp_mat[b], method='pearson').r)

def _bicor(tup):
    global exp_mat
    a, b = tup
    return a, b, float(pg.corr(exp_mat[a], exp_mat[b], method='bicor').r)

def _dcorr(tup):
    global exp_mat
    a, b = tup
    return a, b, float(pg.distance_corr(exp_mat[a], exp_mat[b]).dcor)


def coexpression_graph(file, metric='PCC', save_files=None, lognormalize_input=False, num_threads=4, sep=None):
    """
    Given a file containing n rows of gene expression profiles with m accessions each, it
    applies the given metric ('PCC', 'BICOR', 'DCORR') to calculate a coexpression matrix.
    Finally, a histogram is calculated with bins of 0.01 and a threshold is calculated to
    create a coexpression graph.

    The threshold is obtained based on two references:
    
        1. The bin which causes minimum network density
        2. The bin which has the highest marginal addition of nodes to the network, that is, whose difference in nodes with the previous bin is maximal.

    |
    
    :param file: Path to the file containing the expression profiles
    :type file: str
    :param metric: Metric to be applied to correlate expression profiles: 'PCC' for Pearson 
        correlation coefficient, 'BICOR' for Biweighted Midcorrelation, and 'DCORR' for Distance 
        correlation. Default: 'PCC'
    :type metric: str
    :param save_files: Wether the intermediate files should be stored in disk. If value is 
        _None_, files are not saved in disk. Otherwise, a string indicating the path should be 
        given. The intermediate files which are stored are:

        - Given a expression matrix with n genes, each with m expression values, ``n-1`` files with the name *{metric}-{i}.txt* are stored containing the upper diagonal coexpression matrix. Each file will have the coexpression values using **metric** for gene *i* with every other *j* genes, for 0 <= i < j < n.
        - File *{metric}-{i}.txt* (e.g. PCC-0.txt) will be a comma-separated file including each *j* and the coexpression value between *i* and *j*.
        - Plots of the form *{metric}-{parameter}.png*, where parameter can be Density, CCmaxsize, Nodes and Edges
        - The trimmed network: *{metric}trimmed.csv*.
        
    :type save_files: None or str
    :param lognormalize_input: Whether a logarithm base 10 will be applied to expression values before 
        calculating the coexpression metric
    :type lognormalize_input: bool
    :param num_threads: Number of processes to be used in parallel when calculating coexpression. 
    :type num_threads: int
    :param sep: The string used to separate values. 
    :type sep: None or str

    :return: E: The coexpression matrix as an edge list. Each value in the list will be a tuple with 
        three values: the first two values are the zero-indexed positions of the genes and the third
        value is the corresponding coexpression value between the gene expression profiles.
    :rtype: List[tuple]
    """
    global exp_mat, nodesize
    
    # 1. Read expression file and apply log if required
    tmp = pd.read_csv(file, sep=sep)
    names = list(tmp[tmp.columns[0]])
    dict_names = {v:u for u, v in enumerate(names)}
    size = len(names)
    nodesize = len(names)
    exp_mat = tmp.drop(columns=tmp.columns[0]).to_numpy()
    del tmp
    if lognormalize_input:
        exp_mat = np.log10(exp_mat)
    result = []
    # 2. Calculate coexpression metric in parallel using multiprocessing
    with ProcessPoolExecutor(max_workers=num_threads) as exe:
        # Maps the method 'cube' with a iterable
        values = list(itt.combinations(range(size), 2))
        if metric == 'PCC':
            result = list(exe.map(_pcc, values, chunksize=num_threads*10))
        elif metric == 'BICOR':
            result = list(exe.map(_bicor, values, chunksize=num_threads*10))
        elif metric == 'DCORR':
            result = list(exe.map(_dcorr, values, chunksize=num_threads*10))
        else:
            raise Exception(f"'{metric}' is not a valid metric name. Valid options are 'BICOR', 'DCORR' or 'PCC'")
    
    if save_files is not None:
        last_i = 0
        save_files = Path(save_files).resolve()
        for i in range(size): 
            with open(save_files / f'{metric}-{last_i}.txt', "w"):
                pass
        f = open(save_files / f'{metric}-{last_i}.txt', "a")
        for k, a in enumerate(result):
            i, j, val = result[k]
            if last_i != i:
                last_i = i
                f.close()
                f = open(save_files / f'{metric}-{last_i}.txt', "a")
            f.write(f"{j},{val:0.6f}\n")
        f.close()

    # 3. Calculate bins measurements

    a, b, vals = list(zip(*result))
    mn, mx = min(vals), max(vals)
    result = list(zip(vals, a, b))
    result.sort(reverse=True)
    gc.collect()
    
    mx = ceil(mx / 0.01)            #mx = ceil(Edges[0][0] / 0.01)
    mn = max(-100, min(floor(mn / 0.01), 0)) #mn = min(floor(Edges[-1][0] / 0.01), 0)
    #print(mn, mx, mx - mn)
    vals = [mn/100 + 0.01*i for i in range(mx - mn)]
    #print(vals)
    #print(vals[:5], vals[-10:-1])
    bins = [0 for x in vals]        # actual number of edges at Threshold
    stat_nodes = [0 for x in vals]  # actual number of nodes (maybe disconnected)
    nodeset = set()                 # unique nodes added to the network
    maxCCsize = [1 for x in vals]   # size of biggest connected component
    density = [0 for x in vals]     # network density at Threshold
    pos = len(bins) - 1
    _resetDSets()
    
    for k in result:
        # move to the correct bin, if necessary, and replicate the values
        while pos > 0 and k[0] < vals[pos]:
            pos -= 1
            bins[pos] = bins[pos + 1]
            stat_nodes[pos] = stat_nodes[pos + 1]
            density[pos + 1] = bins[pos + 1] / \
                            (0.5 * stat_nodes[pos + 1] * (stat_nodes[pos + 1] - 1))
            maxCCsize[pos] = maxCCsize[pos + 1]

        #counting number of edges and nodes
        bins[pos] += 1
        nodeset.add(k[1])
        nodeset.add(k[2])
        stat_nodes[pos] = len(nodeset)
        # if two nodes weren't connected before, then ...
        # print(k); print(_myparent(k[1])); print(_myparent(k[2]))
        if _myunion(k[1], k[2]):
            maxCCsize[pos] = max(maxCCsize[pos], _mysetsize(k[1]))
    
    density[0] = 1.0

    # 4. Generating plots
    if save_files is not None:
        plt.scatter(vals, stat_nodes)
        plt.xlabel(f"{metric} cutoff")
        plt.ylabel("Number of nodes")
        plt.grid()
        plt.savefig(save_files / f"{metric}-Nodes.png", dpi=600)
        plt.close()

        plt.scatter(vals, bins)
        plt.xlabel(f"{metric} cutoff")
        plt.ylabel("Number of edges")
        plt.grid()
        plt.savefig(save_files / f"{metric}-Edges.png", dpi=600)
        plt.close()

        plt.scatter(vals, density)
        plt.xlabel(f"{metric} cutoff")
        plt.ylabel("Network Density")
        plt.grid()
        plt.savefig(save_files / f"{metric}-Density.png", dpi=600)
        plt.close()

        plt.scatter(vals, maxCCsize)
        plt.xlabel(f"{metric} cutoff")
        plt.ylabel("Size of biggest connected component")
        plt.grid()
        plt.savefig(save_files / f"{metric}-CCmaxsize.png", dpi=600)
        plt.close()

        with open(save_files / f'{metric}-bins.csv', 'w') as f:
            f.write("Threshold,Edges,Density,Nodes,maxCCsize\n")
            for v in zip(vals, bins, density, stat_nodes, maxCCsize):
                f.write("{:.2f},{},{:.8f},{},{}\n".format(*v))

    # 5. Finding threshold based on density and marginal addition of nodes
    rho_min_pos = 0
    for i, val in enumerate(density):
        if val < density[rho_min_pos]:
            rho_min_pos = i

    deltas = [0] + [a - b for a, b in zip(stat_nodes[:-1], stat_nodes[1:])]
    delta_nodes_max = 0
    for i, val in enumerate(deltas):
        if val > deltas[delta_nodes_max]:
            delta_nodes_max = i

    print("Position for min density:", rho_min_pos, "->", density[rho_min_pos])
    print("Position for max change in nodes:", delta_nodes_max, "->", deltas[delta_nodes_max])
    #print("deltas =", deltas, flush=True)
    pos_min = min(delta_nodes_max, rho_min_pos)
    print("Selected Threshold:", vals[pos_min], "-> id =", pos_min)
    
    # 6. Apply threshold to find coexpression graph

    Edges = []
    for k in result:
        val, a, b = k
        if val >= vals[pos_min]:
            Edges.append( (a, b, val) )

    if save_files is not None:
        with open(save_files / f"{metric}trimmed.csv", "w") as f:
            for k in Edges:
                f.write("{0},{1},{2:0.6f}\n".format(*k))
    return Edges


if __name__ == '__main__':

    for metric in ['PCC', 'BICOR', 'DCORR']:

        print(time.strftime(f"%Y %m %d %H:%M:%S Coexpression graph extraction\nMETRIC: '{metric}'"), flush=True)
        E = coexpression_graph("data/A0.tsv", metric=metric, save_files="example", lognormalize_input=True, num_threads=os.cpu_count()-1, sep="\t")

    print(time.strftime("%Y %m %d %H:%M:%S Bye!"), flush=True)
