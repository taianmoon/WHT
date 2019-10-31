#import ezgal
import matplotlib.pyplot as plt
#import glob
import numpy as np
import os.path
#from scipy import interpolate
#from scipy import integrate
import math



photzout_file = open("z_m_column.txt", "r")    #========================================================HERE
lines = photzout_file.readlines()
#print lines[2:]
#print len(lines)
photzout_file.close()
#lines_t=np.transpose(lines)

#lines = np.array(lines)

#z_a=[]

#for i in range(0, len(lines)):
#    z_a_elements=lines[i]
#    z_a_elements.split()
#    z_a.append(z_a_elements[2])
    #print z_a_elements

#z_a=lines[:][2]
#z_m=lines[:][2]

#print z_a_elements

#z_m=[]


lines = np.array(lines)

lines=np.array([float(i) for i in lines])

z_m=[]

for i in range(0, len(lines)):
    if lines[i]>0.04:         #========================================================HERE (exclude z~0 ones)
        z_m.append(lines[i])





#plt.hist(z_m, bins=np.arange(min(z_m), max(z_m) + 0.2, 0.2))                     #========================================================HERE (binwidth is 0.2)
plt.hist(z_m, bins=70)                     #========================================================HERE (binwidth is 0.2)
#plt.title("Marginalized Photometric Redshift Estimation for H12-00 in G12 with WHT (I,J,H,Ks) Filters")   
plt.xlabel("Photometric Redshift")
plt.ylabel("Number of Sources")
#plt.arrow( 3.26, 30.0, 0.0, -5.0, fc="k", ec="k",head_width=0.2, head_length=2.0 ) #========================================================HERE (redshift is 3.26)
plt.annotate('Redshift of H12-00 lensed source', xy=(3.26, 10.0), xytext=(3.26, 20.0), arrowprops=dict(facecolor='black', shrink=0.05, width=1.0, frac=0.1, headwidth=17.0))

#plt.rc('font', size=40)  #==========================================================HERE (font size)

plt.show()





