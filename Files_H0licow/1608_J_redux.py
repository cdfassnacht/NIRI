"""
Code to calibrate and coadd NIRI J-band data for 1608

Usage: python 1608_J_redux.py [tile] [redux_step]
"""

import os,sys
import numpy as n
import niri_redux as niri
import astromatic as astrom
import astrom_simple as astsimp
from ccdredux import fixpix_wht

if len(sys.argv)<3:
    print ""
    print "ERROR: 1608_J_redux.py requires two input parameters"
    print "  1. Tile name.  For this data set (1608 J) the only option is 'all'"
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
    print "Example: python 1608_J_redux.py all calib_2"
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
rawdir   = '../../Raw/1608_J'
caldir   = '../1608_calib'
rawroot  = 'N20130503S0'
flatfile = '%s/Flat_J.fits' % caldir
bpmfile  = '%s/Flat_J_bpm.pl' % caldir
tilebpm  = 'bpm_from_sky.fits'
astcat   = '1608_J_astrom.cat'
finalout = '1608_niri_2013A_J.fits'
didstep  = False

""" Set up the tile-dependent variables """
skybpm     = None
t1root   = '1608_J_tile1'
t2root   = '1608_J_tile2'
t3root   = '1608_J_tile3'
t4root   = '1608_J_tile4'
all_frames = n.arange(200,254)

if tile=="all":
   sci_frames = all_frames
   outroot    = '1608_J'
   ccdbfile   = 'ccmap_final.txt'
   skybpm     = 'bpm_from_sky.fits'

""" 
Do the nonlinearity correction and create the sky file as the first pass 
"""
if redpass=='calib_1':
   didstep = True
   niri.calib_1(sci_frames,rawroot,outroot,bpmfile,rawdir)

""" Calibrate the science frames """
if redpass=='calib_2':
   didstep = True
   print ""
   skyfile  = '%s_sky.fits' % outroot
   niri.reduce_sci(sci_frames,outroot,skyfile,flatfile)

""" Make the bad pixel mask from the calibrated science frames """
if redpass=='calib_3':
   didstep = True
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
   didstep = True
   niri.split_and_fix_ff(sci_frames,badpixfile=tilebpm,fixpix=True)

""" Make the SExtractor catalogs """
if redpass=='make_cats':
   didstep = True
   fc_files = []
   for i in sci_frames:
      fc_files.append('fc%d_sci.fits' % i)
   niri.niri_sextractor(fc_files)

""" Run first pass of scamp on the files """
if redpass=='scamp_1':
   didstep = True
   os.system('scamp fc*sci.cat -c scamp_niri_1608.config')

""" Run swarp on the files """
if redpass=='swarp':
   didstep = True
   os.system('swarp @good_frames_pass1.txt -c swarp_1608_1.config')

""" If necessary, run second pass of scamp on the files """
if redpass=='scamp_2':
   didstep = True
   addparam = '-ASTREFCAT_NAME swarp_median.cat'
   os.system('scamp fc*sci.cat -c scamp_niri_1608.config %s' % addparam)

""" More error checking """
if didstep == False:
   print ''
   print 'ERROR: Unrecognized reduction step was requested.  Valid values are:'
   print '     calib_1'
   print '     calib_2'
   print '     calib_3'
   print '     calib_4'
   print '     make_cats'
   print '     scamp_1'
   print '     swarp'
   print '     scamp_2'
   print ''
