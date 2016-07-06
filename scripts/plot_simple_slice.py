### Script to plot the output of the CST_to_healpix_mod.py Script
###---------------------------------------------------------------
import healpy as hp
import numpy as n
import sys,optparse
from pylab import *
from scipy import optimize
#nFile = sys.argv[1]
#sFile = sys.argv[2]
#dipole = sys.argv[3]
#trans = sys.argv[4]

#NB assumes TX pol is the same as the receiver pol
o = optparse.OptionParser()
o.set_description('make E and H place slices of two maps and take the ratio. optionally correct for transmitter beam')
o.add_option('--trans',type=str,help='transmitter polarization[NS or EW]')
opts,args = o.parse_args(sys.argv[1:])
sFile = args[0]
mFile = args[1]
S = hp.read_map(sFile)
M = hp.read_map(mFile)

s_err_File = sFile.replace('power','rms')
if s_err_File == sFile:
    S_err = n.zeros_like(S)
else:
    S_err = hp.read_map(s_err_File)

#14.9416951725 4.28571428571 8.77551020408
#101.073957306 4.48979591837 8.77551020408

#receiver coordinates
#alt = n.linspace(-n.pi/2,n.pi/2)
alt = n.linspace(-n.pi/2,n.pi/2)
az = n.zeros_like(alt)

#tx_beam = 10*n.log10(1 - n.sin(alt)**2*n.sin(az)**2)
S_slice_E = hp.pixelfunc.get_interp_val(S,alt,az) #- tx_beam
S_slice_E_err = hp.pixelfunc.get_interp_val(S_err,alt,az)

#pi/2
S_slice_H = hp.pixelfunc.get_interp_val(S,alt,az+n.pi/2) #- tx_beam
S_slice_H_err = hp.pixelfunc.get_interp_val(S_err,alt,az+n.pi/2)

M_slice_E = hp.pixelfunc.get_interp_val(M,alt,az) #- tx_beam

#pi/2
M_slice_H = hp.pixelfunc.get_interp_val(M,alt,az+n.pi/2)


figure(figsize=(8,6))
#suptitle('N/S comparison in E and H planes for transmitter orientation: %s'%opts.trans)
#suptitle('\n'.join(sys.argv[1:]))
#suptitle('\n'.join(args)+'\n'+' '.join([l for l in sys.argv[1:] if l not in args]))
#ax = subplot(211)
if opts.trans=='NS':
    EH_plane = ['E','H']
    EH_colors = ['b','g','r','c']
elif opts.trans=='EW':
    EH_plane = ['H','E']
    EH_colors = ['g','b','c','r']
errorbar(alt*180/n.pi,S_slice_E,S_slice_E_err,fmt='.'+EH_colors[1],label='ECHO [E]')
errorbar(alt*180/n.pi,S_slice_H,S_slice_H_err,fmt='.'+EH_colors[3],label='ECHO [H]')
#errorbar(alt*180/n.pi,M_slice_E,fmt='.'+EH_colors[1],label='B [{plane}]'.format(plane=EH_plane[1]))
#errorbar(alt*180/n.pi,M_slice_H,fmt='.'+EH_colors[3],label='B [{plane}]'.format(plane=EH_plane[1]))
plot(alt*180/n.pi,M_slice_E,'b-',label='Model E')
plot(alt*180/n.pi,M_slice_H,'r-',label='Model H')
title(str(sFile.split('/')[4:-1])+'\n'+mFile)
#ax.xaxis.set_ticklabels([])
legend(ncol=2,loc='best')
ylabel('[dB]')
xlabel('elevation [d]')
show()
