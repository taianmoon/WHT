import pyfits
import matplotlib.pyplot as plt
import numpy as np
import math




#====================================Plot WHT colors as scatter plot======================================================================


names_wht = pyfits.open('../NGP7/Colour_NGP7_160520.fits')                            #========================================================HERE
names_wht_data = names_wht[1].data
wht_j_ks= names_wht_data.field('J-Ks')
wht_h_ks = names_wht_data.field('H-Ks')
wht_j_h = names_wht_data.field('J-H')
wht_j_ks_err = names_wht_data.field('J-Ks_err')
wht_h_ks_err = names_wht_data.field('H-Ks_err')
wht_j_h_err = names_wht_data.field('J-H_err')
wht_ks_err=names_wht_data.field('M_calibrated_err_K')
wht_j_err=names_wht_data.field('M_calibrated_err_J')
wht_h_err=names_wht_data.field('M_calibrated_err_H')



wht_j_ks_k5sig=[]
wht_j_ks_err_k5sig=[]
wht_h_ks_k5sig=[]
wht_h_ks_err_k5sig=[]

wht_j_ks_k5sig_nonJ=[]
wht_j_ks_err_k5sig_nonJ=[]
wht_h_ks_k5sig_nonJ=[]
wht_h_ks_err_k5sig_nonJ=[]
wht_ks_err_k5sig_nonJ=[]

wht_j_ks_k5sig_nonH=[]
wht_j_ks_err_k5sig_nonH=[]
wht_h_ks_k5sig_nonH=[]
wht_h_ks_err_k5sig_nonH=[]
wht_ks_err_k5sig_nonH=[]

wht_j_ks_k5sig_nonJH=[]
wht_j_ks_err_k5sig_nonJH=[]
wht_h_ks_k5sig_nonJH=[]
wht_h_ks_err_k5sig_nonJH=[]
wht_ks_err_k5sig_nonJH=[]







for counter_wht in range(len(wht_j_ks)):
    
    #wht_ks_err_loop=wht_ks_err[counter_wht]
    if wht_ks_err[counter_wht] <= 0.198 and wht_j_ks[counter_wht]<10 and wht_j_ks[counter_wht]>-10 and wht_h_ks[counter_wht]<10 and wht_h_ks[counter_wht]>-10 and wht_j_err[counter_wht]<=0.198 and wht_h_err[counter_wht]<=0.198: #detected in all JHKs
        wht_j_ks_k5sig.append(wht_j_ks[counter_wht])
        wht_j_ks_err_k5sig.append(wht_j_ks_err[counter_wht])
        wht_h_ks_k5sig.append(wht_h_ks[counter_wht])
        wht_h_ks_err_k5sig.append(wht_h_ks_err[counter_wht])

    elif wht_ks_err[counter_wht] <= 0.198 and wht_j_ks[counter_wht]<10 and wht_j_ks[counter_wht]>-10 and wht_h_ks[counter_wht]<10 and wht_h_ks[counter_wht]>-10 and wht_j_err[counter_wht]>=0.198 and wht_h_err[counter_wht]<=0.198: #non-dection in J
        wht_j_ks_k5sig_nonJ.append(wht_j_ks[counter_wht])
        wht_j_ks_err_k5sig_nonJ.append(wht_j_ks_err[counter_wht])
        wht_h_ks_k5sig_nonJ.append(wht_h_ks[counter_wht])
        wht_h_ks_err_k5sig_nonJ.append(wht_h_ks_err[counter_wht])
        wht_ks_err_k5sig_nonJ.append(wht_ks_err[counter_wht])

    elif wht_ks_err[counter_wht] <= 0.198 and wht_j_ks[counter_wht]<10 and wht_j_ks[counter_wht]>-10 and wht_h_ks[counter_wht]<10 and wht_h_ks[counter_wht]>-10 and wht_j_err[counter_wht]<=0.198 and wht_h_err[counter_wht]>=0.198: #non-dection in H
        wht_j_ks_k5sig_nonH.append(wht_j_ks[counter_wht])
        wht_j_ks_err_k5sig_nonH.append(wht_j_ks_err[counter_wht])
        wht_h_ks_k5sig_nonH.append(wht_h_ks[counter_wht])
        wht_h_ks_err_k5sig_nonH.append(wht_h_ks_err[counter_wht])
        wht_ks_err_k5sig_nonH.append(wht_ks_err[counter_wht])

    elif wht_ks_err[counter_wht] <= 0.198 and wht_j_ks[counter_wht]<10 and wht_j_ks[counter_wht]>-10 and wht_h_ks[counter_wht]<10 and wht_h_ks[counter_wht]>-10 and wht_j_err[counter_wht]>=0.198 and wht_h_err[counter_wht]>=0.198: #non-dection in J and H
        wht_j_ks_k5sig_nonJH.append(wht_j_ks[counter_wht])
        wht_j_ks_err_k5sig_nonJH.append(wht_j_ks_err[counter_wht])
        wht_h_ks_k5sig_nonJH.append(wht_h_ks[counter_wht])
        wht_h_ks_err_k5sig_nonJH.append(wht_h_ks_err[counter_wht])
        wht_ks_err_k5sig_nonJH.append(wht_ks_err[counter_wht])

        
#print min(wht_j_ks_k5sig-0.5)
plt.scatter(wht_j_ks_k5sig, wht_h_ks_k5sig, color='0.5')
plt.scatter(wht_j_ks_k5sig_nonJ, wht_h_ks_k5sig_nonJ, color='0.5')
plt.scatter(wht_j_ks_k5sig_nonH, wht_h_ks_k5sig_nonH, color='0.5')
plt.scatter(wht_j_ks_k5sig_nonJH, wht_h_ks_k5sig_nonJH, color='0.5')
plt.errorbar(wht_j_ks_k5sig, wht_h_ks_k5sig, xerr=wht_j_ks_err_k5sig, yerr=wht_h_ks_err_k5sig, ls='none', color='0.5')
plt.errorbar(wht_j_ks_k5sig_nonJ, wht_h_ks_k5sig_nonJ, xerr=[wht_ks_err_k5sig_nonJ, np.full(len(wht_ks_err_k5sig_nonJ), 0.25)], yerr=wht_h_ks_err_k5sig_nonJ, xuplims=True, ls='none', color='0.5')
plt.errorbar(wht_j_ks_k5sig_nonH, wht_h_ks_k5sig_nonH, xerr=wht_j_ks_err_k5sig_nonH, yerr=[wht_ks_err_k5sig_nonH, np.full(len(wht_ks_err_k5sig_nonH), 0.25)], uplims=True, ls='none', color='0.5')
plt.errorbar(wht_j_ks_k5sig_nonJH, wht_h_ks_k5sig_nonJH, xerr=[wht_ks_err_k5sig_nonJH, np.full(len(wht_ks_err_k5sig_nonJH), 0.25)], yerr=[wht_ks_err_k5sig_nonJH, np.full(len(wht_ks_err_k5sig_nonJH), 0.25)], uplims=True, xuplims=True, ls='none', color='0.5')
plot_xlim=np.hstack((wht_j_ks_k5sig,wht_j_ks_k5sig_nonJ,wht_j_ks_k5sig_nonH,wht_j_ks_k5sig_nonJH))
plot_ylim=np.hstack((wht_h_ks_k5sig,wht_h_ks_k5sig_nonJ,wht_h_ks_k5sig_nonH,wht_h_ks_k5sig_nonJH))
plt.xlim((min(plot_xlim)-0.5, max(plot_xlim)+0.5))
plt.ylim((min(plot_ylim)-0.5, max(plot_ylim)+0.5))
plt.xlabel('J-Ks [mag]')
plt.ylabel('H-Ks [mag]')
plt.grid(True)








#=============================Overlaid templates onto the scatter plot================================================

template_list=['Elliptical', 'Ly_break', 'Red_SF_glx_2'] #====================================================================================HERE

names_template = pyfits.open('./all_templates_160505.fits')                            #========================================================HERE
names_template_data = names_template[1].data

for counter_template in range(len(template_list)):
    template_j_ks= names_template_data.field(template_list[counter_template]+'_j_ks')
    template_h_ks= names_template_data.field(template_list[counter_template]+'_h_ks')
    template_j_h= names_template_data.field(template_list[counter_template]+'_j_h')
    template_redshift= names_template_data.field(template_list[counter_template]+'_redshift')

    plt.plot(template_j_ks, template_h_ks, label=str(template_list[counter_template]), marker='o')
    for i, txt in enumerate(template_redshift):
        plt.annotate(str(txt), (template_j_ks[i],template_h_ks[i]))
    plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0., fontsize="small")

plt.title('NGP7')                                                                       #========================================================HERE
#plt.show()


#plt.savefig('./JKs_HKs_redshift_160505.jpg')                                                        




















#============================Overlaid EzGal templates to the scatter plot========================================================


import ezgal
from pylab import *

model_list=['cb07_burst_0.1_z_0.008_salp.model','c09_exp_1.0_z_0.002_salp.model','basti_ssp_n_0.4_z_0.008_krou.model']              #===============================HERE
#try_name='c09 model, exp SFH 0.1 Gyr, different IMF'





for elements in range(len(model_list)):

#============================Load in Ezgal model================================================


    # load ezgal model file
    model = ezgal.model( model_list[elements] )
    # desired formation redshift
    zf = 4.5
    # fetch an array of redshifts out to given formation redshift
    zs = model.get_zs( zf )
    # plot magnitude evolution versus redshift for three filters
    
    model.set_vega_output()
    j_mag=model.get_apparent_mags( zf, filters='j', zs=zs )
    h_mag=model.get_apparent_mags( zf, filters='h', zs=zs )
    ks_mag=model.get_apparent_mags( zf, filters='ks', zs=zs )
    
    #===============Now we have zs, j_mag, h_mag, ks_mag============================================
    
    
    
    #plot( zs, j_mag, 'k-', label='J' )
    #plot( zs, h_mag, 'r--', label='H' )
    #plot( zs, ks_mag, 'b:', label='Ks' )
    # and set labels
    #xlabel( 'z' )
    #ylabel( 'Apparent Mag' )
    #title('cb07_burst_0.1_z_0.008_salp.model')
    
    
    # how about a legend?
    #legend(loc=4)
    # all done
    #show()
    
    
    j_ks=j_mag-ks_mag
    j_h=j_mag-h_mag
    h_ks=h_mag-ks_mag
    
    #print zs[::2]
    #zs_part=zs[::2]
    
    plot( j_ks, h_ks, 's-', label=str(model_list[elements]) )
    for i, txt in enumerate(zs[::10]):
        annotate(txt, (j_ks[10*i],h_ks[10*i]))
    

    
#xlabel( 'J-Ks' )
#ylabel( 'H-Ks' )
#title(try_name)
legend(loc=1)
show()














