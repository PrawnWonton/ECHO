#! /usr/bin/env python

import numpy as n,sys,os
import healpy as hpy
import optparse
from pylab import *

#
######################### Opts #############################

o = optparse.OptionParser()
o.set_description('Convert a CST output txt file into a healpix map')
o.add_option('--nside',type=int,default=64,help='nside of output map')
o.add_option('--rot90',action='store_true',help='rotate the output by 90 degrees')
o.add_option('--rotangle',type=float,help='angle to rotate by')
opts,args = o.parse_args(sys.argv[1:])

if opts.rot90:
    outfile = args[0].replace('.txt','_rot90.fits')
else:
    outfile = args[0].replace('.txt','.fits')


#load the input CST txt file
#find the healpix indices corresponding to the theta/phi
#make a healpix map


D = n.loadtxt(args[0],skiprows=3)
theta = D[:,0]*n.pi/180
#print theta.shape,theta.min(),theta.max()
#theta += n.pi/2
#print theta.min(),theta.max()
#CST will output either elev or phi. phi is zero at beam x=y=0, elev is zero at x=0,phi=0
#if open(args[0]).readlines()[0].startswith('Elev'): theta += n.pi/2

phi = D[:,1]*n.pi/180
beam = D[:,2] #beam amplitude in dB
print theta[n.where(beam==beam.max())[0]]

#phi,beam = (phi[indices],beam[indices])
#print theta.shape,phi.shape,beam.shape

if opts.rot90:
    healpix_indexes = hpy.ang2pix(opts.nside,theta,phi+n.pi/2)
else:
    healpix_indexes = hpy.ang2pix(opts.nside,theta,phi)

if opts.rotangle:
    angle = opts.rotangle/180
    healpix_indexes = hpy.ang2pix(opts.nside,theta,phi+n.pi*angle)
else:
    healpix_indexes = hpy.ang2pix(opts.nside,theta,phi)

#print "healpix index for first point",healpix_indexes[0]
#print "min(theta)",min(theta)*180/n.pi
#print "phi[min(theta)]",phi[n.min(theta)==theta]
#print "value at zenith",beam[theta==n.min(theta)]
hp_map = n.zeros(hpy.nside2npix(opts.nside))

hp_map[healpix_indexes] = beam
hp_map -= hp_map.max()

#receiver coordinates
alt = n.linspace(-n.pi/2,n.pi/2)
az = n.zeros_like(alt)
slice_E = hpy.pixelfunc.get_interp_val(hp_map,alt+n.pi/2,az) #- tx_beam
slice_H = hpy.pixelfunc.get_interp_val(hp_map,alt+n.pi/2,az-n.pi/2)
plot(alt*180/n.pi,slice_E,'b-')
plot(alt*180/n.pi,slice_H,'r-')
xlabel('zenith angle')
ylabel('dB')
show()
hpy.write_map(outfile,hp_map)
#print "Write successfull"
