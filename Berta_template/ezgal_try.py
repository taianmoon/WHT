import ezgal
from pylab import *


#model_list=['cb07_burst_0.1_z_0.008_salp.model', 'cb07_const_1.0_tV_0.2_z_0.008_salp.model', 'cb07_exp_0.1_z_0.008_salp.model', 'cb07_ssp_z_0.008_salp.model']
#try_name='cb07 model, different SFH'

#model_list=['basti_ssp_n_0.4_z_0.008_krou.model','bc03_exp_0.1_z_0.008_salp.model','c09_exp_0.1_z_0.0096_salp.model','cb07_exp_0.1_z_0.008_salp.model','m05_exp_0.1_z_0.02_salp.model']
#try_name='exp or ssp SFH, different models'


#model_list=['c09_exp_0.1_z_0.002_salp.model','c09_exp_0.5_z_0.002_salp.model','c09_exp_1.0_z_0.002_salp.model','c09_exp_10.0_z_0.002_salp.model']
#try_name='c09 model, mettalicity 0.002, different e-folding time for exponential decay SFH'


#model_list=['c09_exp_0.1_z_0.002_salp.model','c09_exp_0.1_z_0.0049_salp.model','c09_exp_0.1_z_0.0096_salp.model','c09_exp_0.1_z_0.019_salp.model','c09_exp_0.1_z_0.024_salp.model','c09_exp_0.1_z_0.03_salp.model']
#try_name='c09 model, exp SFH 0.1 Gyr, different mettalicity'

model_list=['c09_exp_0.1_z_0.002_chab.model','c09_exp_0.1_z_0.002_krou.model','c09_exp_0.1_z_0.002_salp.model']
try_name='c09 model, exp SFH 0.1 Gyr, different IMF'


#'p2_ssp_z_0.0004_salp.model' N/A


#===========================J-H v.s. J-Ks===================================================



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
    
    plot( j_h, j_ks, 's-', label=str(model_list[elements]) )
    #plot(j_h, j_ks, 's')
    for i, txt in enumerate(zs[::4]):
        annotate(txt, (j_h[4*i],j_ks[4*i]))
    

    
xlabel( 'J-H' )
ylabel( 'J-Ks' )
title(try_name)
legend(loc=4)
show()










































#===========================J-Ks v.s. H-Ks===================================================




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
    for i, txt in enumerate(zs[::4]):
        annotate(txt, (j_ks[4*i],h_ks[4*i]))
    

    
xlabel( 'J-Ks' )
ylabel( 'H-Ks' )
title(try_name)
legend(loc=4)
show()









































#===========================J-H v.s. H-Ks===================================================




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
    
    plot( j_h, h_ks, 's-', label=str(model_list[elements]) )
    for i, txt in enumerate(zs[::4]):
        annotate(txt, (j_h[4*i],h_ks[4*i]))
    

    
xlabel( 'J-H' )
ylabel( 'H-Ks' )
title(try_name)
legend(loc=4)
show()















