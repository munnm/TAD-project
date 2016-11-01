import pandas as pd
import numpy as np

import time
from math import log

rawESCdata = pd.read_csv("~/Dropbox/DS/wTrush/rawESCdata.csv", header = 1)

chr1 = rawESCdata[['start: t##', 'end: t##', 'interaction']]
chr1 = chr1.rename(index = str, 
  columns={"start: t##": "startTAD", "end: t##": "endTAD"})
chr1['distance'] = (1/chr1['interaction']).astype(float)
#chr1['distance'] = [log(y,10) for y in (1/chr1['interaction']).astype(float)]


chr2 = rawESCdata[['start: t##', 'end: t##', 'interaction']]
chr2 = chr2.rename(index = str, 
  columns={"start: t##.1": "startTAD", "end: t##.1": "endTAD"})
chr2['distance'] = (1/chr2['interaction']).astype(float)


chr3 = rawESCdata[['start: t##.2', 'end: t##.2', 'interaction.2']]
chr3 = chr3.rename(index = str, 
  columns={"start: t##.2": "startTAD", "end: t##.2": "endTAD"})
chr3['distance'] = (1/chr3['interaction.2']).astype(float)


chr4 = rawESCdata[['start: t##.3', 'end: t##.3', 'interaction.3']]
chr4 = chr4.rename(index = str, 
  columns={"start: t##.3": "startTAD", "end: t##.3": "endTAD"})
chr4['distance'] = (1/chr4['interaction.3']).astype(float)


""" First for chr1"""
#to time the process
#start_time = time.time()

dist1=pd.DataFrame()

for i in range(len(chr1)):
    dist1.set_value(chr1['startTAD'][i], chr1['endTAD'][i], chr1['distance'][i])
    
dist1 = dist1.sort_index(axis = 1,ascending = True)
dist1 = dist1.sort_index(axis = 0,ascending = True)

chr1TADmetric = pd.DataFrame(dist1)

chr1TADmetric.to_csv("chr1TADmetric.csv")



#print("---%s seconds ---" % (time.time() - start_time))


#To create the full distance metric on the space of TADs

#to time the process
#start_time = time.time()

full_dist1 = pd.DataFrame()
for i in range(len(chr1)):
    full_dist1.set_value(i,i,0)
    full_dist1.set_value(chr1['startTAD'][i], chr1['endTAD'][i], chr1['distance'][i]) 

#re-order TADs into increasing order along columns/rows
full_dist1 = full_dist1.sort_index(axis = 1,ascending = True)
full_dist1 = full_dist1.sort_index(axis = 0,ascending = True)

#remove the column 0 and row 0, not relevant
full_dist1 = full_dist1.drop(full_dist1[[0]], axis = 0)
full_dist1 = full_dist1.drop(full_dist1[[0]], axis = 1)

for i in range(1,len(full_dist1[1])):
  for j in range(1,len(full_dist1[1])):
    if i < j and not pd.isnull(full_dist1[j][i]):
      full_dist1[i][j] = full_dist1[j][i]

chr1fullTADmetric = pd.DataFrame(full_dist1)

chr1fullTADmetric.to_csv("chr1fullTADmetric.csv")


#print("---%s seconds ---" % (time.time() - start_time))


