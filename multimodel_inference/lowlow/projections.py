#!/usr/bin/env python
import moments
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
import sys
infile=sys.argv[1]
popid=[sys.argv[2]]
proj=range(int(sys.argv[3]),int(sys.argv[4]))
dd = moments.Misc.make_data_dict(infile)
maxS=0
maxproj=0
for p in range(len(proj)):
    data = moments.Spectrum.from_data_dict(dd, pop_ids=popid,
        projections=[proj[p]],polarized=True)
    print proj[p],data.S()
    if maxS<data.S():
            maxS=data.S()
            maxproj=proj[p]
print "\n------\n",maxproj,maxS,"\n"
