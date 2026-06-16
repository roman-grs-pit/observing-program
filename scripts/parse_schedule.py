from astropy.table import Table
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

sl = Table.read(
    '../data/schedule-RD26075MP1-run_id_RD26075MP1_R26069HJROLLP1.csv')
t9 = Table.read('../data/994-hlwas-Feb26.sim.ecsv')

# select down to grism and XMM data
selg = t9['BANDPASS'] == 'GRISM'
selx = t9['TARGET_NAME'] == 'XMM-LSS'
selc = t9['TARGET_NAME'] == 'COSMOS'

passgx = np.unique(t9['PASS'][selg & selx])
passgc = np.unique(t9['PASS'][selg & selc])

seggx = np.unique(t9['SEGMENT'][selg & selx])
seggc = np.unique(t9['SEGMENT'][selg & selc])

slseg = sl['visit-id'] % 1000

pseg = sl['visit-id'] % 1000000
slpass = (pseg-slseg)//1000

selpassx = np.isin(slpass, passgx)
selpassc = np.isin(slpass, passgc)

selsegx = np.isin(slseg, seggx)
selsegc = np.isin(slseg, seggc)

yd = []

for i in range(len(sl)):
    yd.append(float(sl['start'][i][:8]))
yd = np.array(yd)
yf = yd % np.array(yd).astype(int)
yf /= .365
yd = yd.astype(int)+yf

plt.hist(yd[selpassx & selsegx], bins=25)
plt.xlabel('Year')
plt.ylabel('Number of visits')
plt.title('XMM-LSS')
plt.show()
plt.hist(yd[selpassc & selsegc], bins=25)
plt.xlabel('Year')
plt.ylabel('Number of visits')
plt.title('COSMOS')
plt.show()
