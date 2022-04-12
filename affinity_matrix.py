"""

|

.. _affinity:

******************************************
**Module for Affinity Matrix Computation**
******************************************

Create new distance matrix using information from the gene co-expression
network and information from the associations between genes and functions.
The new distance between two genes will be the mean of the weight between 
two nodes and the proportion of shared functions.
Please load the coexpression data without labels. Dataframe between genes 
and functions gene_by_func must have gene ID in the first column and an 
array of functions in the second column.

|

.. _api:

**API Reference**
=================

"""

import sys
import json
import pandas as pd
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt

def affinity_matrix(edgelist, gene_by_func, normalize, save):
  """The maximum and minimum values of co-expression are saved and two dictionaries are created to optimize searches
  
  :param edgelist: Coexpression matrix as an edge list. Each value in the list will be a tuple with 
      three values: the first two values are the genes associated and the third value is the corresponding 
      coexpression value between them. *source* and *target* are of type *int*, while *score* is *float* 
  :type edgelist: DataFrame
  :param gene_by_func: The matrix with the data of functions associated with a gene, 
      this matrix must have gene ID in the first column and an array of all its
      functional annotations in the second column identified with the GO term. 
  :type gene_by_func: DataFrame
  :param normalize: The coexpression values given by the edge list are normalized.
  :type normalize: bool
  :param save: The affinity matrix is saved as a csv archive, else return the new variable.
  :type save: bool

  :return: The affinity matrix as a new relationship value between genes is returned.
  :rtype: Matrix

  |
  
  """
  smax = edgelist['score'].max()
  smin = edgelist['score'].min()
  total = len(edgelist)
  edgelist = edgelist.to_dict('records')
  # %%

  data = list()
  for x in tqdm(edgelist, total=total):
    u,v,s = x['source'],x['target'],x['score']
    ns = s
    if normalize == True:
      ''' If parameter normalize is True, data is normalized'''
      ns = (s-smin)/(smax-smin)
    pf = 1
    if u in gene_by_func and v in gene_by_func:
      '''If two genes share biological functions, 
         a relationship factor is assigned.'''
      fu = gene_by_func[u]
      fv = gene_by_func[v]
      cfunc = np.intersect1d(fu, fv)
      afunc = np.union1d(fu, fv)
      pf = 1 - len(cfunc)/len(afunc)
    '''The average value between the co-expression value and the 
       relationship factor is calculated.'''
    w = np.mean([ns, pf])
    data.append([u,v,w])

  '''The calculated value of affinity is saved in a matrix'''
  new_edgelist = pd.DataFrame(data, columns=['source','target','weight'])
  if save == True:
    new_edgelist.to_csv('affinity_edgelist.csv', index=False)
  else:
    return new_edgelist

def test():
  """Charge the matrix with the data of functions associated with a gene, 
  this matrix must have gene ID in the first column and an array of all it's 
  functions in the second column

  |
  
  """

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



if __name__ == '__main__':
  test()
