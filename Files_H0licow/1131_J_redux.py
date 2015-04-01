"""
Code to calibrate and coadd NIRI J-band data for 1131

Usage: python 1131_J_redux.py [tile] [redux_step]
"""

import os,sys
import numpy as n
import niri_redux as niri
import astromatic as astrom
import astrom_simple as astsimp
from ccdredux import fixpix_wht

if len(sys.argv)<3:
    print ""
    print "ERROR: 1131_J_redux.py requires two input parameters"
    print "  1. Tile name.  For this data set (1131 J) the only option is 'all'"
    print "  2. Reduction step.  The choices are:"
    print "       calib_1"
    print "       calib_2"
    print "       calib_3"
    print "       calib_4"
    print "       make_cats"
    print "       scamp_1"
    print "       swarp"
    print "       scamp_2"
    print ""
    print "Example: python 1131_J_redux.py all calib_2"
    print ""
    exit()

tile = sys.argv[1]
if tile != 'all':
   print ""
   print "ERROR: Tile must have value 'all' for this system"
   print ""
   exit()

redpass = sys.argv[2]

""" Set up variables that are not tile-dependent"""
rawdir   = '../../Raw/1131_J'
caldir   = '../1131_calib'
rawroot  = 'N20130501S0'
flatfile = '%s/Flat_J.fits' % caldir
bpmfile  = '%s/Flat_J_bpm.pl' % caldir
tilebpm  = 'bpm_from_sky.fits'
astcat   = '1131_J_astrom.cat'
scamp_1  = 'scamp_niri_1131.config'
finalout = '1131_niri_2012B_J.fits'

"""
Description of observation files:
 Tiles 1 and 4
   Observation date: 2013-05-01
   Junk frame:       312
   Frames:           313-334 (2 repetitions of 11 dithered pointings)
   Reference frame:  318
 Tiles 2 and 3
   Observation date: 2013-05-01
   Frames:           347-368 (2 repetitions of 11 dithered pointings)
   Reference frames: 352
"""

""" Set up the tile-dependent variables """
skybpm     = None
t1root   = '1131_J_tile1'
t2root   = '1131_J_tile2'
t3root   = '1131_J_tile3'
t4root   = '1131_J_tile4'
t1and4_frames = n.arange(313,335)
t2and3_frames = n.arange(347,369)
all_frames = n.concatenate((t1and4_frames,t2and3_frames))

if tile=="all":
   sci_frames = all_frames
   outroot    = '1131_J'
   ccdbfile   = 'ccmap_final.txt'
   skybpm     = 'bpm_from_sky.fits'

""" 
Do the nonlinearity correction and create the sky file as the first pass 
"""
if redpass=='calib_1':
   niri.calib_1(sci_frames,rawroot,outroot,bpmfile,rawdir)

""" Calibrate the science frames """
if redpass=='calib_2':
   print ""
   skyfile  = '%s_sky.fits' % outroot
   niri.reduce_sci(sci_frames,outroot,skyfile,flatfile)

""" Make the bad pixel mask from the calibrated science frames """
if redpass=='calib_3':
   if skybpm is not None:
      niri.niri_bpm_from_sky(sci_frames,skybpm,outsky='sky_from_ff.fits')
   else:
      print ""
      print "*** No bad pixel mask requested for this set of input files ***"
      print ""

""" 
Prepare the files for running SExtractor
 The input files will be called ff*fits
 The output files will be called fp*fits
"""
if redpass=='calib_4':
   niri.split_and_fix_ff(sci_frames,badpixfile=tilebpm,fixpix=True)

""" Make the SExtractor catalogs """
if redpass=='make_cats':
   fc_files = []
   for i in sci_frames:
      fc_files.append('fc%d_sci.fits' % i)
   niri.niri_sextractor(fc_files,catformat='ldac')

""" Run scamp the first time """
if redpass=='scamp_1':
    os.system('scamp fc*sci.cat -c %s' % scamp_1)
