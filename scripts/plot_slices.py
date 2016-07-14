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
o.add_option('--a_tx_correct',action='store_true',help='correct for the tx beam on the first antenna')
o.add_option('--b_tx_correct',action='store_true',help='correct for the tx beam on the second antenna')
o.add_option('--trans',type=str,help='transmitter polarization[NS or EW]')
o.add_option('--fit_tx',action='store_true',help='fit for a transmitter orientation offset on the first input antenna')
o.add_option('--ratio',type=str,help='load a healpix file with a ratio measurement')
opts,args = o.parse_args(sys.argv[1:])
nFile = args[0]
sFile = args[1]
N = hp.read_map(nFile)
S = hp.read_map(sFile)

n_err_File = nFile.replace('power','rms')
s_err_File = sFile.replace('power','rms')
print "loading error file",n_err_File
if n_err_File == nFile:
    N_err = n.zeros_like(N)
else:
    N_err = hp.read_map(n_err_File)
if s_err_File == sFile:
    S_err = n.zeros_like(S)
else:
    S_err = hp.read_map(s_err_File)
if opts.trans=='NS':
    tx_pol = 0
elif opts.trans=='EW':
    tx_pol = 1
else:
    print "error parsing --trans={trans}, please choose EW or NS".format(trans=opts.trans)
    sys.exit()
#14.9416951725 4.28571428571 8.77551020408
#101.073957306 4.48979591837 8.77551020408

#receiver coordinates
alt = n.linspace(-n.pi/2,n.pi/2)
az = n.zeros_like(alt)

def fit_tx_pointing(pointing,alt,az,N,S,N_err,S_err,tx_pol):
    altoff = pointing[0]
    azoff = pointing[1]
    tx_beam = 10*n.log10(1 - n.sin(alt+altoff*n.pi/180)**2*n.sin(az+azoff*n.pi/180 +n.pi/2*tx_pol)**2) #XXX
    N_slice_E = hp.pixelfunc.get_interp_val(N,alt,az) - tx_beam
    tx_beam = 10*n.log10(1 - n.sin(alt)**2*n.sin(az+n.pi/2*tx_pol )**2) #XXX
    S_slice_E = hp.pixelfunc.get_interp_val(S,alt,az) - tx_beam
    N_slice_E_err = hp.pixelfunc.get_interp_val(N_err,alt,az)
    S_slice_E_err = hp.pixelfunc.get_interp_val(S_err,alt,az)
    
    
    tx_beam = 10*n.log10(1 - n.sin(alt+altoff*n.pi/180)**2*n.sin(az+azoff*n.pi/180+n.pi/2+n.pi/2*tx_pol)**2)
    N_slice_H = hp.pixelfunc.get_interp_val(N,alt,az+n.pi/2) - tx_beam
    tx_beam = 10*n.log10(1 - n.sin(alt)**2*n.sin(az+n.pi/2+n.pi/2*tx_pol)**2)
    S_slice_H = hp.pixelfunc.get_interp_val(S,alt,az+n.pi/2) - tx_beam
    N_slice_H_err = hp.pixelfunc.get_interp_val(N_err,alt,az+n.pi/2)
    S_slice_H_err = hp.pixelfunc.get_interp_val(S_err,alt,az+n.pi/2)

    chisquare = n.ma.mean(n.ma.masked_invalid((N_slice_E - S_slice_E)**2/(N_slice_E_err**2+S_slice_E_err**2) + \
        (N_slice_H - S_slice_H)**2/(N_slice_H_err**2+S_slice_H_err**2)))
    return chisquare   
 
#ALT,AZ = n.meshgrid(n.linspace(0,90),n.linspace(0,180))
#CHI = n.zeros_like(ALT)
#for i,(za,phi) in enumerate(zip(ALT.ravel(),AZ.ravel())):
#    CHI.ravel()[i] = fit_tx_pointing([za,phi],alt,az,N,S,N_err,S_err)
#figure()
#imshow(n.log(CHI),interpolation='nearest',extent=(ALT.min(),ALT.max(),AZ.min(),AZ.max()))
#colorbar()
#show()
if opts.fit_tx:
    print "fitting for a tx pointing offset on these E and H plane slices"
    results = optimize.minimize(fit_tx_pointing,[2.,7.],(alt,az,N,S,N_err,S_err,tx_pol),bounds=([0,90],[0,180]),tol=1e-8)
    print "tilt = ",results.x[0]
    print "rotation = ",results.x[1]
    altoff = results.x[0]
    azoff = results.x[1]
else:
    altoff = 0
    azoff=0
#tx_beam = 10*n.log10(1 - n.sin(alt+altoff*n.pi/180)**2*n.sin(az+azoff*n.pi/180)**2)
N_slice_E = hp.pixelfunc.get_interp_val(N,alt,az) #- tx_beam
#tx_beam = 10*n.log10(1 - n.sin(alt)**2*n.sin(az)**2)
S_slice_E = hp.pixelfunc.get_interp_val(S,alt,az) #- tx_beam
N_slice_E_err = hp.pixelfunc.get_interp_val(N_err,alt,az)
S_slice_E_err = hp.pixelfunc.get_interp_val(S_err,alt,az)
N_slice_H = hp.pixelfunc.get_interp_val(N,alt,az+n.pi/2) #- tx_beam
S_slice_H = hp.pixelfunc.get_interp_val(S,alt,az+n.pi/2) #- tx_beam
N_slice_H_err = hp.pixelfunc.get_interp_val(N_err,alt,az+n.pi/2)
S_slice_H_err = hp.pixelfunc.get_interp_val(S_err,alt,az+n.pi/2)
if opts.a_tx_correct:
    N_E_beam = 10*n.log10(1 - n.sin(alt+altoff*n.pi/180)**2*n.sin(az+azoff*n.pi/180+n.pi/2*tx_pol)**2)
    N_H_beam = 10*n.log10(1 - n.sin(alt+altoff*n.pi/180)**2*n.sin(az+azoff*n.pi/180+n.pi/2+n.pi/2*tx_pol)**2)
else:
    N_E_beam = n.zeros_like(N_slice_E)
    N_H_beam = n.zeros_like(S_slice_E)
if opts.b_tx_correct:
    S_E_beam = 10*n.log10(1 - n.sin(alt)**2*n.sin(az+n.pi/2*tx_pol)**2)
    S_H_beam = 10*n.log10(1 - n.sin(alt)**2*n.sin(az+n.pi/2+n.pi/2*tx_pol)**2)
else:
    S_E_beam = n.zeros_like(N_slice_H)
    S_H_beam = n.zeros_like(S_slice_H)

if not opts.ratio is None:
    A_B_ratio = hp.read_map(opts.ratio)
    A_B_E_slice = hp.pixelfunc.get_interp_val(A_B_ratio,alt,az)
    A_B_H_slice = hp.pixelfunc.get_interp_val(A_B_ratio,alt,az+n.pi/2)


figure(figsize=(8,6))
#suptitle('N/S comparison in E and H planes for transmitter orientation: %s'%opts.trans)
#suptitle('\n'.join(sys.argv[1:]))
suptitle('\n'.join(args)+'\n'+' '.join([l for l in sys.argv[1:] if l not in args]))
ax = subplot(211)
if opts.trans=='NS':
    EH_plane = ['E','H']
    EH_colors = ['b','g','r','c']
elif opts.trans=='EW':
    EH_plane = ['H','E']
    EH_colors = ['g','b','c','r']
errorbar(alt*180/n.pi,N_slice_E,N_slice_E_err,fmt='.'+EH_colors[0],label='A {plane}'.format(plane=EH_plane[0]))
errorbar(alt*180/n.pi,S_slice_E,S_slice_E_err,fmt='.'+EH_colors[1],label='B [{plane}]'.format(plane=EH_plane[1]))
errorbar(alt*180/n.pi,N_slice_H,N_slice_H_err,fmt='.'+EH_colors[2],label='A {plane}'.format(plane=EH_plane[0]))
errorbar(alt*180/n.pi,S_slice_H,S_slice_H_err,fmt='.'+EH_colors[3],label='B [{plane}]'.format(plane=EH_plane[1]))
if opts.a_tx_correct:
    plot(alt*180/n.pi,N_slice_E-N_E_beam,EH_colors[0],label='A [{plane}] (-tx beam)'.format(plane=EH_plane[0]))
    plot(alt*180/n.pi,N_slice_H-N_H_beam,EH_colors[2],label='A [{plane}] (-tx beam)'.format(plane=EH_plane[1]))
if opts.b_tx_correct:
    plot(alt*180/n.pi,S_slice_E-S_E_beam,EH_colors[1],label='B [{plane}] (-tx beam)'.format(plane=EH_plane[0]))
    plot(alt*180/n.pi,S_slice_H-S_H_beam,EH_colors[3],label='B [{plane}] (-tx beam)'.format(plane=EH_plane[1]))
ax.xaxis.set_ticklabels([])
legend(ncol=2,loc='best')
ylabel('[dB]')
subplot(212)
N_S_E_err = n.sqrt(N_slice_E_err**2 + S_slice_E_err**2)
N_S_H_err = n.sqrt(N_slice_H_err**2 + S_slice_H_err**2)

errorbar(alt*180/n.pi,N_slice_E-S_slice_E,N_S_E_err,fmt='.'+EH_colors[0],label='[{plane}]'.format(plane=EH_plane[0]))
errorbar(alt*180/n.pi,N_slice_H-S_slice_H,N_S_H_err,fmt='.'+EH_colors[1],label='[{plane}]'.format(plane=EH_plane[1]))
if opts.a_tx_correct or opts.b_tx_correct:
    plot(alt*180/n.pi,N_slice_E-N_E_beam-(S_slice_E-S_E_beam),EH_colors[0],label='[{plane}] (-tx beam)'.format(plane=EH_plane[0]))
    plot(alt*180/n.pi,N_slice_H-N_H_beam-(S_slice_H-S_H_beam),EH_colors[1],label='[{plane}] (-tx beam)'.format(plane=EH_plane[1]))
if not opts.ratio is None:
    plot(alt*180/n.pi,A_B_E_slice,'--'+EH_colors[0],label = 'ORBCOMM [{plane}]'.format(plane=EH_plane[0]))
    plot(alt*180/n.pi,A_B_H_slice,'--'+EH_colors[1],label = 'ORBCOMM [{plane}]'.format(plane=EH_plane[1]))

ylabel('[dB]')
legend(ncol=2,loc='best')
xlabel('elevation [d]')
subplots_adjust(hspace=0)
show()

