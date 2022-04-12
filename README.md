- Written by: Nicolás López, Miguel Romero, Camilo Rocha and Jorge Finke
- Last updated: Apr 12/2022

---

# coexpression_graph
Python modules for the calculation of coexpression graph and affinity matrix, based on expression profiles

## Table of Contents

- [Module 1: coexpression_graph](#module-1-coexpression_graph)
  - [General algorithm](#general-algorithm)
  - [Reading the input file](#reading-the-input-file)
  - [Calculating the Threshold](#calculating-the-threshold)
  - [Saving auxiliary files](#saving-auxiliary-files)
- [Module 2: affinity_matrix](#module-2-affinity_matrix)
- [Example](#example)
- [API reference](#api-reference)

## Module 1: coexpression_graph

This module takes as input several gene profiles and calculates an undirected
graph based on the coexpression between them, calculated using one of the following metrics:

- Pearson correlation coefficient ('PCC')
- Biweight Midcorrelation ('BICOR')
- Distance correlation ('DCORR')


### General algorithm

The general idea behind the module can be summarized as follows:

1. Read expression file and apply log_10 if required
2. Calculate coexpression metric (in parallel using multiprocessing)
3. Calculate bins for quantification
4. Generating auxiliary plots (optional)
5. Finding threshold based on density and marginal addition of nodes
6. Apply threshold to find coexpression graph


### Reading the input file

The input file for *coexpression_graph* is loaded using *pandas*, which infers automatically if 
the file has headers. However, there is a possibility to specify the file delimiter using the parameter
**sep**, which by default is *None*.

It is also optional to apply a logarithm (base 10) to the expression values before calculating the correlation
metric. The corresponding parameter for this effect is **lognormalize_input**. By default, this option is set to *False*.


### Calculating the Threshold

The selection of the threshold is calculated as follows:

1. The selected metric is calculated for all posible pairs of gene profiles.

2. A number of bins is created and the network properties are quantified using those bins. The size of each bin is 0.01, ranging from the minimum to the maximum values calculated previously.

3. Two bins are of importance: the one where the minimum of network density occurs, and the one which has a higher increase in the number of nodes included in the network.

4. The minimum of these two bin values is selected as the threshold.

5. Those pairs whose metric value is greater than or equal to the threshold belong to the coexpression graph. 


The idea behind the threshold selection is based on the following papers:

- van Verk, M. C., Bol, J. F., & Linthorst, H. J. (2011). Prospecting for genes involved in transcriptional regulation of plant defenses, a bioinformatics approach. *BMC plant biology* , 11, 88. https://doi.org/10.1186/1471-2229-11-88
- Zhang, L., Yu, S., Zuo, K., Luo, L., & Tang, K. (2012). Identification of gene modules associated with drought response in rice by network-based analysis. *PloS one* , 7(5), e33748. https://doi.org/10.1371/journal.pone.0033748


### Saving auxiliary files

It is possible to save some files which can help in visualizing properties such as density, 
number of nodes, number of edges, etc. The parameter **save_files** is set by default to *None*,
but can be set to a string indicating the path where the files should be stored.

The stored files are:

- **Plots**: The saved plots as a function of the metric value are the size of the greatest connected component (*{metric}-CCmaxsize.png*), the network density (*{metric}-Density.png*), the number of edges (*{metric}-Edges.png*) and the number of nodes (*{metric}-Nodes.png*). 
- **Metric values**: Files of the form *{metric}-{i}.txt* are generated, including the correlation information.
- **Coexpression graph**: The file *{metric}trimmed.csv* contains the generated network, after the threshold.
- **Bins**: The file *{metric}-bins.csv* can be used to generate the plots again. It includes a header and a line for each bin, specifying the corresponding values.

## Module 2: affinity_matrix

Create new distance matrix using information from the gene co-expression
network and information from the associations between genes and functions.
The new distance between two genes will be the mean of the weight between 
two nodes and the proportion of shared functions.
Please load the coexpression data without labels. Dataframe between genes 
and functions gene_by_func must have gene ID in the first column and an 
array of functions in the second column.

## Example

Note: Include the files `*.py` in the folder you need to execute your code

- For **coexpression_graph** module, you can use the following snippet: 

```python3
import os
import coexpression_graph

for metric in ['PCC', 'BICOR', 'DCORR']:
    E = coexpression_graph("data/A0.tsv", metric=metric, save_files="example", lognormalize_input=True, num_threads=os.cpu_count()-1, sep="\t")
```

- For **affinity_matrix** module, you can use the following snippet: 

```python3
import os
import affinity_matrix

gene_by_func = pd.read_csv("data/example_gene_term_list.csv")
def array2list(str):
    str = str.replace("\'",'\"').replace(" ", ", ")
    return json.loads(str)
gene_by_func['tList'] = gene_by_func.tList.apply(lambda t: array2list(t))
gene_by_func = dict([(x,y) for x,y in gene_by_func.itertuples(index=False)])

'''Charge the matrix with the data of coexpression between two genes, 
this matrix must have the source in the first column, target in the 
second column and the value of coexpression in the third column.'''

data = pd.read_csv('data/example_gcn.csv')
data = data.astype({'source': int, 'target': int, 'score': np.float64})
affinity_matrix(data, gene_by_func, True, False)

```

## API Reference

The reference can be found in the folder `docs/` as a [html webpage](index.html) or as a pdf [here](docs/API.pdf).
