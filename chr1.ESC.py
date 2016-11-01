import pandas as pd
import numpy as np

#Read in data
rawESCdata = pd.read_csv("~/Dropbox/DS/wTrush/rawESCdata.csv", header = 1)

#Isolate the data for chromosome 1, rename columns
chr1 = rawESCdata[['start: t##', 'end: t##', 'interaction']]
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

#Print all TADs which are unioned. Output is a printed list of tuples which are joined
###WRITE THIS A FUNCTION WITH INPUT = d AND OUTPUT = list

graph = []

for i in dist1.index.values:
  for j in dist1.columns.values:
    if dist1.loc[i,j] <d:
      graph.append([i,j])
print "\tWhen d = %r, union the TAD pairs...\n\t\t%r" % (d, graph)

# This is based on depth first search
# An outer loop scans all nodes of graph and starts a search from every node
# Node neighbors are added to the cycle path
# Recursion ends if no more non-visited neighbors can be added
# A new cycle is found if the path is longer than 4 nodes and next neighbor is the start of the path
# To avoid duplicate cycles, the cycles are normalized by rotating the smallest node to the start

# This includes the edit which stops searching edges once the node1 is 'too' large

print "\t This will induce the following unions....\n\t\t"

cycles = []

def main():
    global graph
    global cycles

    # first sort graph so that node1 < node2 for all edges in graph
    for edge in graph:
        edge = edge.sort()

    # order edges in graph by ascending order on node1
    graph = sorted(graph)    
    for edge in graph:
      for vertex in edge:
        findNewCycles([vertex])
    print cycles

def findNewCycles(path):
    start_node = path[0]
    next_node = None
    sub = []

    for edge in graph:
        node1, node2 = edge
        # stop search once you're too far along in the graph
        if node1 <= path[0]:
            if start_node in edge:
                    if node1 == start_node:
                        next_node = node2
                    else:
                        next_node = node1
            if not visited(next_node, path):
		# neighbor node not on path yet
                    sub = [next_node]
                    sub.extend(path)
		# explore extended path
                    findNewCycles(sub);
            elif len(path) > 3  and next_node == path[-1]:
		# cycle found
                    p = rotate_to_smallest(path);
                    inv = invert(p)
                    if isNew(p) and isNew(inv):
                        cycles.append(p)


# next two function is for formating purpose only
def invert(path):
    return rotate_to_smallest(path[::-1])



# rotate cycle path such that it begins with the smallest node
def rotate_to_smallest(path):
    n = path.index(min(path))
    return path[n:]+path[:n]



# check if it is a unique path
def isNew(path):
    return (not path in cycles)



# check if a vertex has been visited
def visited(vertex, path):
    return (vertex in path)

main()












