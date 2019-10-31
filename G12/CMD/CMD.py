import pyfits
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from astropy.io import fits
import matplotlib.cm as cm


names_cutoff = pyfits.open('./G12_IRAC_cat.fits')                            #========================================================HERE
names_cutoff_data = names_cutoff[1].data
#J_Ks= names_cutoff_data.field('J-Ks')
#J_Ks_err= names_cutoff_data.field('J-Ks_err')
#M_calibrated_K = names_cutoff_data.field('M_calibrated_K') 
#M_calibrated_err_K = names_cutoff_data.field('M_calibrated_err_K')
MAG_ISOCOR_3p6_4p5_color = names_cutoff_data.field('MAG_ISOCOR_3p6_4p5_color')
MAG_ISOCOR_3p6_4p5_color_err = names_cutoff_data.field('MAG_ISOCOR_3p6_4p5_color_err')
MAG_ISOCOR_4p5= names_cutoff_data.field('MAG_ISOCOR_4p5')
MAGERR_ISOCOR_4p5= names_cutoff_data.field('MAGERR_ISOCOR_4p5')
#I_3p6= names_cutoff_data.field('I_3p6')
#I_3p6_err= names_cutoff_data.field('I_3p6_err')
#MAG_ISOCOR_3p6= names_cutoff_data.field('MAG_ISOCOR_3p6')
#MAGERR_ISOCOR_3p6= names_cutoff_data.field('MAGERR_ISOCOR_3p6')
#X_WORLD_K= names_cutoff_data.field('X_WORLD_K')
#Y_WORLD_K= names_cutoff_data.field('Y_WORLD_K')
#NUMBER_K= names_cutoff_data.field('NUMBER_K')









plt.errorbar(MAG_ISOCOR_4p5, MAG_ISOCOR_3p6_4p5_color, xerr=MAGERR_ISOCOR_4p5, yerr=MAG_ISOCOR_3p6_4p5_color_err, color='r',fmt='o') #==========================================================HERE
plt.scatter(MAG_ISOCOR_4p5, MAG_ISOCOR_3p6_4p5_color, color='b') #==========================================================HERE


#plt.ylim([-10.0,20.0])
plt.grid()
plt.xlabel('4.5 mircon [mag]')#==========================================================HERE
plt.ylabel('3.6 mircon - 4.5 micron [mag]')#==========================================================HERE
plt.rc('font', size=30)  #==========================================================HERE (font size)
plt.show()
