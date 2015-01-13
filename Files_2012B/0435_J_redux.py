"""
Script to run niri code to reduce the J-band NIRI data
"""

import os
import numpy as n
import niri_redux as niri
import astromatic as astrom
import astrom_simple as astsimp
from ccdredux import fixpix_wht

""" 
**** Choose the tile and pass to run ****
Uncomment one from each of the two groups directly below
"""
#tile = "Tile1"
#tile = "Tiles2and3"
#tile = "Tile4"
tile = "all"

#redpass = 'make_sky'
#redpass = 'final_calib'
#redpass = 'make_cats'
#redpass = 'astrom_tile'
#redpass = 'astrom_mosaic'
redpass = 'final_coadd'

"""
Description of observation files:
 Tile 1
   Observation date: 2012-08-22
   Frames:           281-290 (2 repetitions of 5 dithered pointings)
   Reference frame:  281
 Tiles 2 and 3
   Observation date: 2012-08-26
   Frames:           254-275 (2 repetitions of 11 dithered pointings)
                     Frames can be split into:
                      2x5 frames for Tile 2 (254-258, 265-269)
                      2x6 frames for Tile 3 (259-264, 270-275)
   Reference frames: 254, 260
 Tile 4
   Observation date: 2012-09-16 (NB: Some obs taken on 09-15, but they're junk)
   Frames:           202-213 (2 repetitions of 6 dithered pointings)
"""

""" Set up variables that are not tile-dependent"""
rawdir   = '/Users/cdf/Data/Gemini/2012B/Raw/0435_J'
flatfile = '../0435_calib/Flat_J.fits' 
bpmfile  = '../0435_calib/Flat_J_bpm.pl'
tilebpm = 'bpm_from_sky.fits'
astcat   = '0435_J_astrom.cat'
finalout = '0435_niri_2012B_J.fits'

""" Set up the tile-dependent variables """
skybpm     = None
t1root   = '0435_J_tile1'
t2root   = '0435_J_tile2'
t3root   = '0435_J_tile3'
t4root   = '0435_J_tile4'
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
   rawroot    = 'N20120822S0'
   outroot    = t1root
   refroot    = 'ff281'
   refcat     = '%s.cat' % refroot
   ccdbfile   = 'ccmap_%s_as_ref.txt' % refroot
if tile=="Tile4":
   sci_frames = t4_frames
   ncoadd     = 12
   rawroot    = 'N20120916S0'
   outroot    = t4root
   refroot    = 'ff202'
   refcat     = '%s.cat' % refroot
   ccdbfile   = 'ccmap_%s_as_ref.txt' % refroot
if tile=="all":
   sci_frames = all_frames
   outroot    = '0435_J_final'
   ccdbfile   = 'ccmap_final.txt'
elif tile=="Tiles2and3":
   sci_frames = n.arange(254,276)
   skybpm     = 'bpm_from_sky.fits'
   rawroot    = 'N20120826S0'
   outroot    = '0435_J_tiles2and3'
   t2refroot  = 'ff254'
   t2refcat   = '%s.cat' % t2refroot
   t2ccdbfile = 'ccmap_%s_as_ref.txt' % t2refroot
   t3refroot  = 'ff260'
   t3refcat   = '%s.cat' % t3refroot
   t3ccdbfile = 'ccmap_%s_as_ref.txt' % t3refroot

""" 
Do the nonlinearity correction and Create the sky file as the first pass 
"""
if redpass=='make_sky':
   niri.calib_1(sci_frames,rawroot,outroot,bpmfile,rawdir)

""" Do the final calibration and make the SExtractor catalogs """
if redpass=='final_calib':
   print ""
   skyfile  = '%s_sky.fits' % outroot
   niri.reduce_sci(sci_frames,outroot,skyfile,flatfile)
   ##niri.niri_coadd(outroot) # Don't use this any more
   if skybpm is not None:
      niri.niri_bpm_from_sky(sci_frames,skybpm)

""" Make the SExtractor catalogs """
if redpass=='make_cats':
   ff_files = []
   for i in sci_frames:
      ff_files.append('ff%d.fits' % i)
   niri.niri_sextractor(ff_files,tilebpm,fixpix=True,catformat='ascii')

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
   fitsfile  = '0435_J_tile3.fits'
   fitscat   = '0435_J_tile3.cat'
   astref    = '0435_J_tile2.cat'
   ccmapfile = '0435_J_tile3.ccmap'
   ccdbfile  = 'ccmap_tile3_to_tile2.txt'
   #astsimp.match_fits_to_ast(fitsfile,fitscat,astref,ccmapfile,max_offset=40.,
   #                          doplot=False)
   #astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)

   """ Align Tile1 to Tile2 """
   fitsfile  = '0435_J_tile1.fits'
   fitscat   = '0435_J_tile1.cat'
   astref    = '0435_J_tile2.cat'
   ccmapfile = '0435_J_tile1.ccmap'
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
   fitsfile  = '0435_J_tile4.fits'
   fitscat   = '0435_J_tile4.cat'
   astref    = 'swarp_123.cat'
   ccmapfile = '0435_J_tile4.ccmap'
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
   astsimp.match_fits_to_ast(fitsfile,catfile,'0435_2mass_stars.dat',ccmapfile,
                             racol=0,deccol=1,max_offset=40.,doplot=False)
   astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)
   astrom.make_cat_niri(swarpout,'0435_J_astrom.cat','swarp_1234.reg',
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
   photcat  = '0435_2mass_mags_no_lens.dat'
   photcol  = 3  # Column containing 2MASS J-band magnitudes in photcat
   swarpin  = '0435_J_final_swarp.in'
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
   try:
      niri.niri_photom(sci_files,photcat,photcol)
   except:
      return

   """ 
   Do the final coadd 
   NOTE: The lens is not quite in the center of the field, so that the
         'center' parameter for run_swarp is not exactly on the lens system.
   """
   swarpfile = '@%s'%swarpin
   imcent = '04:38:14.79,-12:17:14.9'  # *** Not exactly on the lens ***
   astrom.run_swarp(swarpfile,finalout,'swarp_niri.config',centertype='MANUAL',
                    center=imcent,imsize=1725,pixscale=0.1165)

   """ Check the photometry of the final coadd """
   niri.niri_photom(finalout,photcat,photcol)
