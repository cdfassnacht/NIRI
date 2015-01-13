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
    print "  1. Tile name.  "
#    print "  1. Tile name.  For this data set (1131 J) the only option is 'all'"
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
#if tile != 'all':
#   print ""
#   print "ERROR: Tile must have value 'all' for this system"
#   print ""
#   exit()

redpass = sys.argv[2]

""" Set up variables that are not tile-dependent"""
rawdir   = '../../Raw/2013_05_01/RXJ1131'
rawroot  = 'N20130501S0'
flatfile = '../1131_calib/Flat_J.fits' 
bpmfile  = '../1131_calib/Flat_J_bpm.pl'
tilebpm  = 'bpm_from_sky.fits'
astcat   = '1131_J_astrom.cat'
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
                     Frames can be split into:
                      2x5 frames for Tile 2 (254-258, 265-269)
                      2x6 frames for Tile 3 (259-264, 270-275)
   Reference frames: 352
"""

""" Set up the tile-dependent variables """
skybpm     = None
t1root   = '1131_J_tile1'
t2root   = '1131_J_tile2'
t3root   = '1131_J_tile3'
t4root   = '1131_J_tile4'
t1_frames  = n.arange(281,291)
t2a_frames = n.arange(254,259)
t2b_frames = n.arange(265,270)
t2_frames  = n.concatenate((t2a_frames,t2b_frames))
t3a_frames = n.arange(259,265)
t3b_frames = n.arange(270,276)
t3_frames  = n.concatenate((t3a_frames,t3b_frames))
t4_frames  = n.arange(202,214)
all_frames = n.concatenate((t1_frames,t2_frames,t3_frames,t4_frames))

if tile=="Tile1":
   sci_frames = t1_frames
   ncoadd     = 10
   outroot    = t1root
   refroot    = 'ff281'
   refcat     = '%s.cat' % refroot
   ccdbfile   = 'ccmap_%s_as_ref.txt' % refroot
if tile=="Tile4":
   sci_frames = t4_frames
   ncoadd     = 12
   outroot    = t4root
   refroot    = 'ff202'
   refcat     = '%s.cat' % refroot
   ccdbfile   = 'ccmap_%s_as_ref.txt' % refroot
if tile=="all":
   sci_frames = all_frames
   outroot    = '1131_J_final'
   ccdbfile   = 'ccmap_final.txt'
if tile=="Tiles1and4":
   sci_frames = n.arange(313,335)
   skybpm     = 'bpm_from_sky.fits'
   outroot    = '1131_J_tiles1and4'
   t2refroot  = 'ff318'
   t2refcat   = '%s.cat' % t2refroot
   t2ccdbfile = 'ccmap_%s_as_ref.txt' % t2refroot
   t3refroot  = 'ff260'
   t3refcat   = '%s.cat' % t3refroot
   t3ccdbfile = 'ccmap_%s_as_ref.txt' % t3refroot
elif tile=="Tiles2and3":
   skybpm     = 'bpm_from_sky.fits'
   outroot    = '1131_J_tiles2and2'
   t2refroot  = 'ff254'
   t2refcat   = '%s.cat' % t2refroot
   t2ccdbfile = 'ccmap_%s_as_ref.txt' % t2refroot
   t3refroot  = 'ff260'
   t3refcat   = '%s.cat' % t3refroot
   t3ccdbfile = 'ccmap_%s_as_ref.txt' % t3refroot

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
   niri.niri_sextractor(fc_files,catformat='ascii')

"""
Do the tile-based pass of the astrometry.
"""

if redpass=='astrom_tile':
   if tile == "Tiles2and3":
      swarplist = '%s_swarp.in'%t2root
      #niri.ccmap_tile(t2_frames,t2refcat,t2ccdbfile,tilebpm,swarplist)
      swarpfile = '@%s'%swarplist
      swarpout  = '%s.fits' % t2root
      astrom.run_swarp(swarpfile,swarpout,'swarp_niri.config')
      astrom.make_cat_niri(swarpout,'%s.cat'%t2root,'%s.reg'%t2root,
                           weight_file='%s_wht.fits'%t2root,texp=42.,ncoadd=10,
                           det_area=70)

      swarplist = '%s_swarp.in'%t3root
      #niri.ccmap_tile(t3_frames,t3refcat,t3ccdbfile,tilebpm,swarplist)
      niri.niri_fixpix(t3_frames)
      swarpfile = '@%s'%swarplist
      swarpout  = '%s.fits' % t3root
      astrom.run_swarp(swarpfile,swarpout,'swarp_niri.config')
      astrom.make_cat_niri(swarpout,'%s.cat'%t3root,'%s.reg'%t3root,
                           weight_file='%s_wht.fits'%t3root,texp=42.,ncoadd=12,
                           det_area=70)
   else:
      swarplist = '%s_swarp.in'%outroot
      niri.ccmap_tile(sci_frames,refcat,ccdbfile,tilebpm,swarplist)
      niri.niri_fixpix(sci_frames)
      swarpfile = '@%s'%swarplist
      swarpout  = '%s.fits' % outroot
      astrom.run_swarp(swarpfile,swarpout,'swarp_niri.config')
      astrom.make_cat_niri(swarpout,'%s.cat'%outroot,'%s.reg'%outroot,
                           weight_file='%s_wht.fits'%outroot,texp=42.,
                           ncoadd=ncoadd,det_area=70)

"""
Align the tiles with each other
"""

if redpass == 'astrom_mosaic':
   """ Align Tile3 to Tile2 """
   fitsfile  = '1131_J_tile3.fits'
   fitscat   = '1131_J_tile3.cat'
   astref    = '1131_J_tile2.cat'
   ccmapfile = '1131_J_tile3.ccmap'
   ccdbfile  = 'ccmap_tile3_to_tile2.txt'
   #astsimp.match_fits_to_ast(fitsfile,fitscat,astref,ccmapfile,max_offset=40.,
   #                          doplot=False)
   #astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)

   """ Align Tile1 to Tile2 """
   fitsfile  = '1131_J_tile1.fits'
   fitscat   = '1131_J_tile1.cat'
   astref    = '1131_J_tile2.cat'
   ccmapfile = '1131_J_tile1.ccmap'
   ccdbfile  = 'ccmap_tile1_to_tile2.txt'
   #astsimp.match_fits_to_ast(fitsfile,fitscat,astref,ccmapfile,max_offset=40.,
   #                          doplot=False)
   #astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)

   """ Coadd the first three tiles """
   swarpin  = '%s.fits %s.fits %s.fits' % (t1root,t2root,t3root)
   swarpout = 'swarp_123.fits'
   ncoadd   = 32
   #astrom.run_swarp(swarpin,swarpout,'swarp_niri.config')
   #astrom.make_cat_niri(swarpout,'swarp_123.cat','swarp_123.reg',
   #                     weight_file='swarp_123_wht.fits',texp=42.,
   #                     ncoadd=ncoadd,det_area=70)

   """ Align Tile4 to the first three tiles """
   fitsfile  = '1131_J_tile4.fits'
   fitscat   = '1131_J_tile4.cat'
   astref    = 'swarp_123.cat'
   ccmapfile = '1131_J_tile4.ccmap'
   ccdbfile  = 'ccmap_tile4_to_tiles123.txt'
   #astsimp.match_fits_to_ast(fitsfile,fitscat,astref,ccmapfile,max_offset=40.,
   #                          doplot=False)
   #astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)

   """ Coadd all four tiles """
   ncoadd   = 42
   os.system('rm -v swarp_123*')
   swarpin  = '%s.fits %s.fits %s.fits %s.fits' % (t1root,t2root,t3root,t4root)
   swarpout = 'swarp_1234.fits'
   astrom.run_swarp(swarpin,swarpout,'swarp_niri.config')
   astrom.make_cat_niri(swarpout,'swarp_1234.cat','swarp_1234.reg',
                        weight_file='swarp_1234_wht.fits',texp=42.,
                        ncoadd=ncoadd,det_area=70)

   """ Align the 4-tile coadd with 2MASS """
   fitsfile  = 'swarp_1234.fits'
   catfile   = 'swarp_1234.cat'
   ccmapfile = 'swarp_1234.ccmap'
   astsimp.match_fits_to_ast(fitsfile,catfile,'1131_2mass_stars.dat',ccmapfile,
                             racol=0,deccol=1,max_offset=40.,doplot=False)
   astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)
   astrom.make_cat_niri(swarpout,'1131_J_astrom.cat','swarp_1234.reg',
                        weight_file='swarp_1234_wht.fits',texp=42.,
                        ncoadd=ncoadd,det_area=100,det_thresh=8.)

"""
Redo the astrometry of all the individual input files onto the full grid
now that we have created it, and then do the final coadd.
"""

if redpass == 'final_coadd':
   """ Set up """
   ff_files  = []
   sci_files = []
   for i in sci_frames:
      ff_files.append('ff%d.fits' % i)
      sci_files.append('ff%d_sci.fits' % i)
   ccdbfile = 'ccmap_final.txt'
   photcat  = '1131_2mass_mags_no_lens.dat'
   photcol  = 3  # Column containing 2MASS J-band magnitudes in photcat
   swarpin  = '1131_J_final_swarp.in'
   f = open(swarpin,'w')

   """ Do the astrometry """
   for i in sci_frames:
      fitsfile  = 'ff%d_sci.fits' % i
      catfile   = 'ff%d.cat' % i
      regfile   = 'ff%d.reg' %i
      ccmapfile = 'ff%d.ccmap' %i
      #secat = astsimp.Secat(catfile)
      #secat.match_fits_to_ast(fitsfile,astcat,ccmapfile,max_offset=40.,
      #                        doplot=False)
      #astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)
      f.write('%s\n'%fitsfile)
   f.close()

   """
   Find the photometric zero points for each file and put that information
   in the fits header
   """
   niri.niri_photom(sci_files,photcat,photcol)

   """ 
   Do the final coadd 
   NOTE: The lens is not quite in the center of the field, so that the
         'center' parameter for run_swarp is not exactly on the lens system.
   """
   swarpfile = '@%s'%swarpin
   imcent = '04:38:14.79,-12:17:14.9'  # *** Not exactly on the lens ***
   astrom.run_swarp(swarpfile,finalout,'swarp_niri.config',centertype='MANUAL',
                    center=imcent,imsize=1725)

   """ Check the photometry of the final coadd """
   niri.niri_photom(finalout,photcat,photcol)
