#Weighted quick-union with path compression.
#This code should produce
""" >>> uf = UnionFind(10)
    >>> for (p, q) in [(3, 4), (4, 9), (8, 0), (2, 3), (5, 6), (5, 9),
    ...                (7, 3), (4, 8), (6, 1)]:
    ...     uf.union(p, q)
    >>> uf._id
    [8, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    >>> uf.find(0, 1)
    True
    >>> uf._id
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    """


import pandas as pd
import numpy as np

class UnionFind:

    def __init__(self, n):
        self._id = list(range(n))
        #uf._sz[i] indicates the size of each node i. That is, if a certain node 
        #class is larger is will subsume the smaller unioned node class
        self._sz = [1] * n

    def _root(self, i):
        j = i
        while (j != self._id[j]):
            self._id[j] = self._id[self._id[j]]
            j = self._id[j]
        return j

    def find(self, p, q):
        return self._root(p) == self._root(q)

    def union(self, p, q):
        i = self._root(p)
        j = self._root(q)
        # Check to see which node class is larger
        if (self._sz[i] < self._sz[j]):
            #if node containing i is smaller, join i to j
            self._id[i] = j
            #and increase the size of the node containing j appropriately
            self._sz[j] += self._sz[i]
        else:
            #else, join j to i
            self._id[j] = i
            #and increase the size of the node containing i appropriately
            self._sz[i] += self._sz[j]

#Read in data
rawNeurondata = pd.read_csv("~/Dropbox/DS/wTrush/rawNeurondata.csv", header = 1)

#Isolate the data for chromosome 1, rename columns
chr1 = rawNeurondata[['start: t##', 'end: t##', 'interaction']]
chr1 = chr1.rename(index = str, 
  columns={"start: t##": "startTAD", "end: t##": "endTAD"})

#Add a column for 'distance'. Note, this is adjustable depending on functionality
chr1['distance'] = (1/chr1['interaction']).astype(float)

#Create a data frame which gives the metric for TADS
dist1=pd.DataFrame()

for i in range(len(chr1)):
    dist1.set_value(chr1['startTAD'][i], chr1['endTAD'][i], chr1['distance'][i])
    
dist1 = dist1.sort_index(axis = 1,ascending = True)
dist1 = dist1.sort_index(axis = 0,ascending = True)

#Ask for threshold distance to union two TADs. Should be a number between 0.0005 and 2.0

d = input("Enter d:  ")


print "\tWhen d = %r, union the TAD pairs...\n" % d
for i in dist1.index.values:
  for j in dist1.columns.values:
    if dist1.loc[i,j] <d:
      uf = UnionFind(len(chr1))
      uf.union(i,j)
      print "\t (%r,%r)" % (i,j)
print "\n"
#Want this to print out the new ID variable with unioned TADs
#But not working right now
print uf._id





"""
uf = UnionFind(10)
for (p, q) in [(3, 4), (4, 9), (8, 0), (2, 3), (5, 6), (5, 9), 
(7, 3), (4, 8), (6, 1)]:
  uf.union(p, q)
print uf._id
print uf.find(0,1)
uf._id
"""