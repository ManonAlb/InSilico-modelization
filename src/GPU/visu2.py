import matplotlib.pyplot as plt
import numpy as np
import abtem.measure
import os


from ase.io import read
from abtem.visualize import show_atoms
from abtem.potentials import Potential
from abtem.waves import PlaneWave
from abtem.measure import block_zeroth_order_spot



def remove_nmax(array, N):
    alpha_sorted = np.sort(array) #sorted from min to max
    clip = alpha_sorted[-N]
    result = np.where( array >= clip, clip, array) #when > n biggest, replace by

    return(result)


ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..')) + '/DIFFRACTION/'

def diffracto(path, title,nmax):

    files = os.listdir(ROOT+path)
    multi = 0
    c=0
    for file in files :

        if file.endswith(('.hdf5')):
                print(file)
                c+=1
                measurement = abtem.measure.Measurement.read(os.path.join(ROOT +path, file))

                alphax = measurement.calibrations[-2].coordinates(measurement.array.shape[-2])  
                size = alphax.shape[0]
                m_array = measurement.array[0:size,int(size/2)]
                multi+=m_array

    average = multi/c
    print(c)
    corrected = remove_nmax(average, nmax)


    plt.xlim(-50,50)
    plt.xlabel('2alpha x [mrad]')
    plt.ylabel('Intesity')
    plt.title(title, y=1.05)

    return(plt.plot(2*alphax, corrected), plt.savefig('diffractogram_{}'.format(path),dpi=300, bbox_inches='tight'))

diffracto('Fixed_protein_in_water', 'Fixed ubiquitin in water, 1000 frames diffractogramm, 1ns simulation', 2)
plt.show()               
