#DO FIRST on command line
#export PATH=$PATH:$HOME/Downloads/montage/bin

import os
import pyfits
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import montage_wrapper as montage
from astropy.io import fits
import matplotlib.cm as cm
import aplpy
import matplotlib.lines as ln


#========================================================HERE============================================

input_fits_image_band1='./input/shift_i_gama12_ks_2h_20160225_astro_2mass_0p25rms.fits'
input_fits_image_band2='./input/shift_i_gama12_i_band_ACAM_53min_20160226_astro_0p241rms.fits'
input_fits_image_band3='./input/shift_gama12_j_1h12m_20160225_astro_2mass_0p22rms.fits'
input_fits_image_band4='./input/G12_cutout_herschel_250.fits'
#input_fits_image_band5='./input/G12v230_IRAC_Mosaic_36.fits'
#input_fits_image_band6='./input/G12v230_IRAC_Mosaic_45.fits'

band4_present=1  #=================================================================1:yes 0:no (I band)
band5_present=1  #=================================================================1:yes 0:no (3.6 band)
band6_present=1  #=================================================================1:yes 0:no (4.5 band)



#montage.mSubimage(input_fits_image_band5, './G12_3p6_trimmed.fits', ra=, dec=, xsize=0.0023)






aplpy.make_rgb_cube([input_fits_image_band1, input_fits_image_band3, input_fits_image_band2], 'G12_cube.fits')   #R, G, B
aplpy.make_rgb_image('G12_cube.fits','G12_cube.png', vmin_b=500.0, vmax_b=3000.0, vmax_r=73.0, vmin_r=-20.0)
#aplpy.make_rgb_image('G12_cube.fits','G12_cube.png', vmin_g=-20.0, vmin_r=-10.0, vmin_b=-10.0)
#aplpy.make_rgb_image('NGP7_cube.fits','NGP7_cube.png', vmin_b=-30.0, vmax_b=3000.0, vmax_r=73.0)
#aplpy.make_rgb_image('G12_JHK_cube.fits', 'G12_JHK_cube.png', pmin_r=0., pmax_r=80., pmin_g=0., pmax_g=80., pmin_b=0., pmax_b=80.)




rgb_image = aplpy.FITSFigure('./G12_cube_2d.fits')
rgb_image.show_rgb('G12_cube.png')
rgb_image.show_contour(input_fits_image_band4, levels=5, colors='white', linewidths=0.5)

rgb_image.recenter(176.65864, -0.19688142, radius=0.035)
#rgb_image.recenter(176.65903,-0.19764512, width=0.066, height=0.07)

rgb_image.add_scalebar(1.0/60.0)
rgb_image.scalebar.set_label("1'")
rgb_image.scalebar.set_color('white')
rgb_image.add_grid()
rgb_image.grid.set_alpha(1.0)
rgb_image.grid.set_linewidth(0.1)
rgb_image.set_title('G12 Field (Red: Ks, Green: J, Blue: I, Contour: Herschel 250 micron)')     #===========================================================HERE
rgb_image.show_regions('./vla_region_wcs.reg')

rgb_image.save('G12_final.png')
