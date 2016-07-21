#!//anaconda/bin/python

import healpy as hp
import numpy as n
import matplotlib.pyplot as plt
import sys

nFile = sys.argv[1]
sFile = sys.argv[2]
#dipole = sys.argv[3]
#trans = sys.argv[4]

N = hp.read_map(nFile)
S = hp.read_map(sFile)

#ratio = 10*n.log10(N/S)
ratio = N-S
use = ratio[n.where(n.isnan(ratio)==False)[0]]
print use.min(),use.max()
rot=[90,90,0]
#title='Ratio of N/'+dipole+' to S/'+dipole+'\n Transmitter Polarization: '+trans
title = 'GBNorth minus bradley model'
hp.mollview(ratio,rot=rot,cmap='gnuplot',unit='dB',format='%.1f',title=title)#,min=-2,max=2)
plt.show()
