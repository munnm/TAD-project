import pandas as pd
import numpy as np

import time
from math import log

rawNeurondata = pd.read_csv("~/Dropbox/DS/wTrush/rawNeurondata.csv", header = 1)

chr1 = rawNeurondata[['start: t##', 'end: t##', 'interaction']]
chr1 = chr1.rename(index = str, 
  columns={"start: t##": "startTAD", "end: t##": "endTAD"})
chr1['distance'] = (1/chr1['interaction']).astype(float)
#chr1['distance'] = [log(y,10) for y in (1/chr1['interaction']).astype(float)]


chr2 = rawNeurondata[['start: t##', 'end: t##', 'interaction']]
chr2 = chr2.rename(index = str, 
  columns={"start: t##.1": "startTAD", "end: t##.1": "endTAD"})
chr2['distance'] = (1/chr2['interaction']).astype(float)


chr3 = rawNeurondata[['start: t##.2', 'end: t##.2', 'interaction.2']]
chr3 = chr3.rename(index = str, 
  columns={"start: t##.2": "startTAD", "end: t##.2": "endTAD"})
chr3['distance'] = (1/chr3['interaction.2']).astype(float)


chr4 = rawNeurondata[['start: t##.3', 'end: t##.3', 'interaction.3']]
chr4 = chr4.rename(index = str, 
  columns={"start: t##.3": "startTAD", "end: t##.3": "endTAD"})
chr4['distance'] = (1/chr4['interaction.3']).astype(float)


""" First for chr1"""
#to time the process
start_time = time.time()

dist1=pd.DataFrame()

for i in range(len(chr1)):
    dist1.set_value(chr1['startTAD'][i], chr1['endTAD'][i], chr1['distance'][i])
    
dist1 = dist1.sort_index(axis = 1,ascending = True)
dist1 = dist1.sort_index(axis = 0,ascending = True)

print("---%s seconds ---" % (time.time() - start_time))


#To create the full distance metric on the space of TADs

#to time the process
start_time = time.time()

full_dist = pd.DataFrame()
for i in range(len(chr1)):
    full_dist.set_value(i,i,0)
    full_dist.set_value(chr1['startTAD'][i], chr1['endTAD'][i], chr1['distance'][i]) 
full_dist = full_dist.sort_index(axis = 1,ascending = True)
full_dist = full_dist.sort_index(axis = 0,ascending = True)

print("---%s seconds ---" % (time.time() - start_time))


