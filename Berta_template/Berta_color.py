import matplotlib.pyplot as plt
import glob
import numpy as np
import os.path
from scipy import interpolate
from scipy import integrate
import math

#==============================HERE=================================================

bandwidth_j=3.0e14/1.17 - 3.0e14/1.33 #UKIRT: 3.0e14/1.17 - 3.0e14/1.33; 2MASS: 2.0e14/0.162; unit of Hz
flux_density_j=1570                   #UKIRT: 1570; 2MASS: 1594; unit of Jy
bandwidth_h=3.0e14/1.49 - 3.0e14/1.78 #UKIRT: 3.0e14/1.49 - 3.0e14/1.78; 2MASS: 2.0e14/0.251; unit of Hz
flux_density_h=1020                   #UKIRT: 1020; 2MASS: 1024; unit of Jy
bandwidth_ks=3.0e14/2.03 - 3.0e14/2.37 #UKIRT: 3.0e14/2.03 - 3.0e14/2.37; 2MASS: 2.0e14/0.262; unit of Hz
flux_density_ks=636                    #UKIRT: 636; 2MASS: 666.8; unit of Jy
                                      #UKIRT in K, not Ks

redshift_loop=np.arange(0,5,0.5)  
file_array = glob.glob("./templates_berta/*") 

    #==============Reproduce Berta's SEDs============================================

color_color_redshift_array=[]
file_array_total=[]

for elements in range(0, len(file_array)):


    filename=str(file_array[elements])
    print filename
    sed_data = open(filename, 'r')

    figure_title = os.path.splitext(filename)[0]

    x=[]
    y=[]
    z=[]

    for line in sed_data:
        columns = line.split()
        x.append(columns[0])
        y.append(columns[1])
    
    x = np.array(x)
    y = np.array(y)

    x=np.array([float(i) for i in x])
    y=np.array([float(i) for i in y])

    z=x*x*y*1.0e12*1.0e-10/3.0e8
    z=z*3.846e52  #change from L_solar/m^2Hz flux to Jy
    x=x*1.0e-10*1.0e6 #change from angstrom to micron

    sed_data.close()


    #=================Add a transmission function of LIRIS filter to Berta's template==========
    
    transmission_j=open("cold_j.txt", 'r')
    transmission_h=open("cold_h.txt", 'r')
    transmission_ks=open("cold_ks.txt", 'r')

    x_liris_j=[]
    x_liris_h=[]
    x_liris_ks=[]

    T_liris_j=[]
    T_liris_h=[]
    T_liris_ks=[]

    for line in transmission_j:
        columns_liris = line.split()
        x_liris_j.append(columns_liris[0])
        T_liris_j.append(columns_liris[1])

    for line in transmission_h:
        columns_liris = line.split()
        x_liris_h.append(columns_liris[0])
        T_liris_h.append(columns_liris[1])

    for line in transmission_ks:
        columns_liris = line.split()
        x_liris_ks.append(columns_liris[0])
        T_liris_ks.append(columns_liris[1])

    x_liris_j = np.array(x_liris_j)
    T_liris_j = np.array(T_liris_j)
    x_liris_h = np.array(x_liris_h)
    T_liris_h = np.array(T_liris_h)
    x_liris_ks = np.array(x_liris_ks)
    T_liris_ks = np.array(T_liris_ks)

    x_liris_j=np.array([float(i) for i in x_liris_j])
    T_liris_j=np.array([float(i) for i in T_liris_j])
    x_liris_h=np.array([float(i) for i in x_liris_h])
    T_liris_h=np.array([float(i) for i in T_liris_h])
    x_liris_ks=np.array([float(i) for i in x_liris_ks])
    T_liris_ks=np.array([float(i) for i in T_liris_ks])

    transmission_j.close()
    transmission_h.close()
    transmission_ks.close()


    interpolate_function_transm_j = interpolate.interp1d(x_liris_j, T_liris_j, bounds_error=False, fill_value=1)  #assume linear interpolation
    interpolate_function_transm_h = interpolate.interp1d(x_liris_h, T_liris_h, bounds_error=False, fill_value=1)  #assume linear interpolation
    interpolate_function_transm_ks = interpolate.interp1d(x_liris_ks, T_liris_ks, bounds_error=False, fill_value=1)  #assume linear interpolation

    z_transm=[]




    for elements_t in range(0, len(x)):

        z_modified=z[elements_t]*interpolate_function_transm_j(x[elements_t])*interpolate_function_transm_h(x[elements_t])*interpolate_function_transm_ks(x[elements_t])
        z_transm.append(z_modified)

    #now, we have x[micron] and z_transm[Jy].





    #=================Interpolate SEDs to extract fluxes at certain wavelengths==========


    #interpolate_function = interpolate.interp1d(x, z)  #assume linear interpolation

    #flux_J=interpolate_function(60)
    #flux_H=interpolate_function(100)
    #flux_K=interpolate_function(100)
    
    #print "flux at 60 microm= "+str(flux_60)+" (Jy)"
    #print "flux at 100 microm= "+str(flux_100)+" (Jy)"

    x_Hz =3.0e14/x #transform from micron (wavelength) to Hz (frequency)
    x_Hz_re=x_Hz[::-1] #reverse the order of x_Hz and z_transm to get increasing x_Hz in frequency
    z_transm_re=z_transm[::-1]
    z_int = integrate.cumtrapz(z_transm_re, x_Hz_re, initial=0) #integrate flux density (z): output cumulative integrated value (array) along x_Hz or z.


    interpolate_function_int=interpolate.interp1d(x_Hz_re, z_int, bounds_error=False, fill_value=0)  #assume linear interpolation



    flux_j=interpolate_function_int(3.0e14/1.17)-interpolate_function_int(3.0e14/1.33)  #in unit of 1.0e-26*W/m^2
    flux_h=interpolate_function_int(3.0e14/1.49)-interpolate_function_int(3.0e14/1.78)
    flux_ks=interpolate_function_int(3.0e14/1.99)-interpolate_function_int(3.0e14/2.31)
    
    #print flux_j, flux_h, flux_ks

    #========================transform to magnitudes, and then colours=============================================


    mag_j=-2.5*math.log10(flux_j)+2.5*math.log10(bandwidth_j*flux_density_j) #taking into account the reference flux
    mag_h=-2.5*math.log10(flux_h)+2.5*math.log10(bandwidth_h*flux_density_h)
    mag_ks=-2.5*math.log10(flux_ks)+2.5*math.log10(bandwidth_ks*flux_density_ks)

    j_ks=mag_j-mag_ks
    h_ks=mag_h-mag_ks
    j_h=mag_j-mag_h

    print mag_j, mag_h, mag_ks
    #print j_ks, h_ks, j_h









    #===============================================================================================
    #===============================================================================================
    #===============================================================================================
    #========================================REDSHIFTED==========================================

    #redshift_loop=np.arange(0,5,0.5)  #=========================================================HERE
    j_ks_z=[]
    h_ks_z=[]
    j_h_z=[]
    redshift=[]

    for i in range(0, len(redshift_loop)):

        redshift.append(redshift_loop[i])
        x_z=x*(1+redshift_loop[i])
        print redshift_loop[i]

        #now, we have x_z[micron] and z[Jy].

        #===========add a transmission function===================================================
        #most parts are inherent from above

        z_transm_z=[]

        for elements_t in range(0, len(x_z)):

            z_modified_z=z[elements_t]*interpolate_function_transm_j(x_z[elements_t])*interpolate_function_transm_h(x_z[elements_t])*interpolate_function_transm_ks(x_z[elements_t])
            z_transm_z.append(z_modified_z)

        #now, we have x_z[micron] and z_transm_z[Jy].


        #==============================do flux, magnitude, colours===================================

        x_Hz_z =3.0e14/x_z #transform from micron (wavelength) to Hz (frequency)
        x_Hz_re_z=x_Hz_z[::-1] #reverse the order of x_Hz_z and z_transm to get increasing x_Hz in frequency
        z_transm_re_z=z_transm_z[::-1] 
        z_int_z = integrate.cumtrapz(z_transm_re_z, x_Hz_re_z, initial=0) #integrate flux density (z): output cumulative integrated value (array) along x_Hz or z.


        interpolate_function_int_z=interpolate.interp1d(x_Hz_re_z, z_int_z, bounds_error=False, fill_value=0)  #assume linear interpolation


 
        flux_j_z=interpolate_function_int_z(3.0e14/1.17)-interpolate_function_int_z(3.0e14/1.33)  #in unit of 1.0e-26*W/m^2
        flux_h_z=interpolate_function_int_z(3.0e14/1.49)-interpolate_function_int_z(3.0e14/1.78)
        flux_ks_z=interpolate_function_int_z(3.0e14/1.99)-interpolate_function_int_z(3.0e14/2.31)
        

        #print flux_j_z, flux_h_z, flux_ks_z

        mag_j_z=-2.5*math.log10(flux_j_z)+2.5*math.log10(bandwidth_j*flux_density_j) #taking into account the reference flux
        mag_h_z=-2.5*math.log10(flux_h_z)+2.5*math.log10(bandwidth_h*flux_density_h)
        mag_ks_z=-2.5*math.log10(flux_ks_z)+2.5*math.log10(bandwidth_ks*flux_density_ks)

        j_ks_z.append(mag_j_z-mag_ks_z)
        h_ks_z.append(mag_h_z-mag_ks_z)
        j_h_z.append(mag_j_z-mag_h_z)

        print mag_j_z-mag_ks_z, mag_h_z-mag_ks_z, mag_j_z-mag_h_z


    #now, we have j_ks_z, h_ks_z, j_h_z arrays (colours at each redshift)


    #=========================Plot colour-colour v.s. redshift========================================================

    np.savetxt(filename+'_j_ks_redshift_0_2_1.txt',j_ks_z)
    np.savetxt(filename+'_h_ks_redshift_0_2_1.txt',h_ks_z)
    np.savetxt(filename+'_j_h_redshift_0_2_1.txt',j_h_z)
    np.savetxt(filename+'_redshift_0_2_1.txt', redshift)

    if elements==0:
        color_color_redshift_array=np.column_stack((j_ks_z, h_ks_z, j_h_z, redshift))
    else:
        color_color_redshift_array=np.column_stack((color_color_redshift_array,j_ks_z, h_ks_z, j_h_z, redshift))

    #file_array_total.append(filename+'_j_ks')
    #file_array_total.append(filename+'_h_ks')
    #file_array_total.append(filename+'_j_h')
    #file_array_total.append(filename+'_redshift')

    if elements==0:
        file_array_total=np.column_stack((filename+'_j_ks',filename+'_h_ks',filename+'_j_h',filename+'_redshift'))
    else:
        file_array_total=np.column_stack((file_array_total,filename+'_j_ks',filename+'_h_ks',filename+'_j_h',filename+'_redshift'))



    fig1 = plt.figure(1)
    plt.plot(j_h_z, j_ks_z, label=str(figure_title))
    #plt.xscale('log')
    #plt.yscale('log')    
    plt.xlabel('J-H [mag]')
    plt.ylabel('J-Ks [mag]')
    plt.title('Colour-colour-redshift relation')
    #plt.legend(loc=4)
    plt.grid(True)
    #for i, txt in enumerate(redshift):
    #    plt.annotate('z='+str(txt), (j_h_z[i],j_ks_z[i]))
    plt.legend(bbox_to_anchor=(0.8, 1), loc=2, borderaxespad=0., fontsize="xx-small")
    plt.xlim((0.5,3.1))
    #plt.savefig(filename+'_JH_JKs_redshift.jpg')
    #plt.close()
    #plt.show()
    #plt.hold(True)


    fig2 = plt.figure(2)
    plt.plot(j_ks_z, h_ks_z, label=str(figure_title))
    #plt.xscale('log')
    #plt.yscale('log')    
    plt.xlabel('J-Ks [mag]')
    plt.ylabel('H-Ks [mag]')
    plt.title('Colour-colour-redshift relation')
    #plt.legend(loc=4)
    plt.grid(True)
    #for i, txt in enumerate(redshift):
    #    plt.annotate('z='+str(txt), (j_ks_z[i],h_ks_z[i]))
    plt.legend(bbox_to_anchor=(0.8, 1), loc=2, borderaxespad=0., fontsize="xx-small")
    plt.xlim((1,5.5))
    #plt.savefig(filename+'_JKs_HKs_redshift.jpg')
    #plt.close()
    #plt.show()    
    #plt.hold(True)


    fig3 = plt.figure(3)
    plt.plot(j_h_z, h_ks_z, label=str(figure_title))
    #plt.xscale('log')
    #plt.yscale('log')    
    plt.xlabel('J-H [mag]')
    plt.ylabel('H-Ks [mag]')
    plt.title('Colour-colour-redshift relation')
    #plt.legend(loc=4)
    plt.grid(True)
    #for i, txt in enumerate(redshift):
    #    plt.annotate('z='+str(txt), (j_h_z[i],h_ks_z[i]))
    plt.legend(bbox_to_anchor=(0.8, 1), loc=2, borderaxespad=0., fontsize="xx-small")
    plt.xlim((0.5,3))
    #plt.savefig(filename+'_JH_HKs_redshift.jpg')
    #plt.close()
    #plt.show()
    #plt.hold(True)


fig1.savefig('./JH_JKs_redshift.jpg')
fig2.savefig('./JKs_HKs_redshift.jpg')
fig3.savefig('./JH_HKs_redshift.jpg')

np.savetxt('./color_color_redshift_array.txt',color_color_redshift_array)
#file_array_total_transpose=np.transpose(file_array_total)
#file_array_total.reshape(-1, 1)
#file_array_total.shape = (1,len(redshift_loop))
np.savetxt('./file_array_total.txt', file_array_total, delimiter="\t", fmt="%s") 



