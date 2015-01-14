"""
Script to run niri code to reduce the Ks-band NIRI data
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
   Frames:           291-295 (5 dithered pointings)
   Reference frame:  293 (291 may be bad)
 Tiles 2 and 3
   Observation date: 2012-08-26
   Frames:           276-291 (NB: Frame 287 appears to be corrupted)
   Reference frame:  276
 Tile 4
   Observation date: 2012-09-16 (NB: Some obs taken on 09-15, but they're junk)
   Frames:           214-225
   Reference frame:  214
"""

""" Set up variables that are not tile-dependent"""
rawdir   = '/Users/cdf/Data/Gemini/2012B/Raw/0435_Ks'
flatfile = '../0435_calib/Flat_Ks.fits' 
bpmfile  = '../0435_calib/Flat_Ks_bpm.pl'
tilebpm  = 'bpm_from_sky.fits'
astcat   = '../0435_J/0435_J_astrom.cat'
finalout = '0435_niri_2012B_Ks.fits'

""" Set up the tile-dependent variables """
skybpm     = None
t1_frames   = n.arange(291,296)
#t23_frames = n.arange(276,292)  # Frame 287 may have problems
t23a_frames = n.arange(276,287)
t23b_frames = n.arange(288,292)
t23_frames  = n.concatenate((t23a_frames,t23b_frames))
t4_frames   = n.arange(214,226)
t1obsdate   = '20120822'
t23obsdate  = '20120826'
t4obsdate   = '20120916'
t1root   = '0435_Ks_tile1'
t2root   = '0435_Ks_tile2'
t3root   = '0435_Ks_tile3'
t23root  = '0435_Ks_tiles23'
t4root   = '0435_Ks_tile4'
all_frames = n.concatenate((t1_frames,t23_frames,t4_frames))

if tile=="Tile1":
   obsdate    = t1obsdate
   rawroot    = 'N%sS0' % obsdate
   sci_frames = t1_frames
   ncoadd     = t1_frames.size
   outroot    = t1root
if tile=="Tiles2and3":
   obsdate    = t23obsdate
   rawroot    = 'N%sS0' % obsdate
   sci_frames = t23_frames
   ncoadd     = t23_frames.size
   outroot    = t23root
   skybpm     = 'bpm_from_sky.fits'
if tile=="Tile4":
   obsdate    = t4obsdate
   rawroot    = 'N%sS0' % obsdate
   sci_frames = t4_frames
   ncoadd     = t4_frames.size
   outroot    = t4root
if tile=="all":
   sci_frames = all_frames
   outroot    = '0435_Ks_final'
   ccdbfile   = 'ccmap_final.txt'

fullnames = []
if tile == 'all':
   for i in t1_frames:
      fullnames.append('%s_%d' % (t1obsdate,i))
   for i in t23_frames:
      fullnames.append('%s_%d' % (t23obsdate,i))
   for i in t4_frames:
      fullnames.append('%s_%d' % (t4obsdate,i))
else:
   for i in sci_frames:
      fullnames.append('%s_%d' % (obsdate,i))

""" 
Do the nonlinearity correction and Create the sky file as the first pass 
"""
if redpass=='make_sky':
   niri.calib_1(sci_frames,rawroot,outroot,bpmfile,rawdir,obsdate)

""" Do the final calibration and make the SExtractor catalogs """
if redpass=='final_calib':
   print ""
   skyfile  = '%s_sky.fits' % outroot
   niri.reduce_sci(fullnames,outroot,skyfile,flatfile)
   if skybpm is not None:
      niri.niri_bpm_from_sky(fullnames,skybpm)

""" Make the SExtractor catalogs """
if redpass=='make_cats':
   ff_files = []
   for i in fullnames:
      ff_files.append('ff%s.fits' % i)
   niri.niri_sextractor(ff_files,tilebpm,fixpix=True,catformat='ascii')

"""
Redo the astrometry of all the individual input files onto the full grid
now that we have created it, and then do the final coadd.
"""

if redpass == 'final_coadd':
   """ Set up """
   ff_files  = []
   sci_files = []
   for i in fullnames:
      ff_files.append('ff%s.fits' % i)
      sci_files.append('ff%s_sci.fits' % i)
   ccdbfile = 'ccmap_final.txt'
   photcat  = '0435_2mass_mags_no_lens.dat'
   photcol  = 5  # Column containing 2MASS Ks-band magnitudes in photcat
   swarpin  = '0435_Ks_final_swarp.in'
   f = open(swarpin,'w')

   """ Do the astrometry """
   for i in fullnames:
      fitsfile  = 'ff%s_sci.fits' % i
      catfile   = 'ff%s.cat' % i
      regfile   = 'ff%s.reg' %i
      ccmapfile = 'ff%s.ccmap' %i
      f.write('%s\n'%fitsfile)
      #secat = astsimp.Secat(catfile)
      #secat.match_fits_to_ast(fitsfile,astcat,ccmapfile,max_offset=40.,
      #                        doplot=False)
      #astsimp.rscale_ccmap(ccmapfile,ccdbfile,fitsfile,interactive=True)
   f.close()

   """
   Find the photometric zero points for each file and put that information
   in the fits header
   """
   #niri.niri_photom(sci_files,photcat,photcol)

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
