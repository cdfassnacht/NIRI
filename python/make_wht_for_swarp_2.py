"""
Modifies the input _wht files before running the final swarp.  The
files are compared to the median-stacked images, and pixels that differ by
more than n sigma from the median are flagged.
"""

import niri_redux as niri
import glob
import sys

if len(sys.argv)<3:
    print ''
    print 'Usage: make_wht_for_swarp_2.py requires TWO input arguments'
    print '  1. Stacked median file name'
    print '  2. nsigma for flagging'
    print ''
    print 'Example:'
    print 'python make_wht_for_swarp_2.py swarp_median.fits 3.'
    print ''

else:
    outname = sys.argv[1]
    nsig = float(sys.argv[2])
    infiles = glob.glob('f*_sci.resamp.fits')

    niri.make_wht_for_final(infiles,outname,nsig)
