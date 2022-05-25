import matplotlib.pyplot as plt
import numpy as np
import abtem.measure
import os

from ase.io import read
from abtem.visualize import show_atoms
from abtem.potentials import Potential
from abtem.waves import PlaneWave
from abtem.measure import block_zeroth_order_spot

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..')) + '/DIFFRACTION/'

def remove_nmax(array, N):
    if N <= 1:
        return array

    alpha_sorted = np.sort(array) #sorted from min to max
    clip = alpha_sorted[-N]
    result = np.where( array >= clip, clip, array) #when > n biggest, replace by

    return(result)

def diffracto(path, title, nmax):

    # list files in directory
    files = os.listdir(ROOT+path)
    multi = 0
    c=0
    for file in files :

        if file.endswith(('.hdf5')) :
                print(file)
                c+=1

                measurement = abtem.measure.Measurement.read(os.path.join(ROOT + path, file))

                alphax = measurement.calibrations[-2].coordinates(measurement.array.shape[-2])

                # take the slice
                image = None
                if len(measurement.array.shape) == 3:
                    image = measurement.array[0]
                elif len(measurement.array.shape) == 2:
                    image = measurement.array
                else:
                    raise

                i_size, j_size = image.shape[0], image.shape[1]
                i_mid, j_mid = int(i_size/2), int(j_size/2)

                size = alphax.shape[0] #int(np.sqrt(i_size**2 + j_size**2))
                m_array = np.zeros(size)
                count = np.zeros(size)

                middle = int(size/2)

                # add rest of the circle
                for i in range(i_size):
                    for j in range(j_size):
                        # radius
                        r = np.sqrt((i - i_mid)**2 + (j - j_mid)**2)

                        # for sign
                        sign = 0
                        if (i - i_mid) != 0:
                            sign = (i - i_mid)/abs(i - i_mid)
                        elif (j - j_mid) != 0:
                            sign = (j - j_mid)/abs(j - j_mid)

                        # add to average
                        pixel = middle + sign*r
                        pixel_left = int(np.floor(pixel))
                        pixel_right = int(np.ceil(pixel))

                        weight_left = pixel - pixel_left
                        weight_right = 1 - weight_left

                        if pixel_left >= 0 and pixel_left < size:
                            m_array[pixel_left] += image[i, j]*weight_left
                            count[pixel_left] += 1

                        if pixel_right >= 0 and pixel_right < size:
                            m_array[pixel_right] += image[i, j]*weight_right
                            count[pixel_right] += 1

                # compute average
                for i in range(size):
                    m_array[i] /= count[i]



                multi+=m_array

    average = multi/c
    print(c)
    corrected = remove_nmax(average, nmax)


    plt.xlim(0,50)
    plt.xlabel('q')
    plt.ylabel('Intesity')

    plt.yscale('log')

    return(plt.plot((alphax/2, corrected), plt.savefig('diffractogram_{}'.format(path),dpi=300, bbox_inches='tight'))

diffracto('Protein', 'Ubiquitin in vacuum for 1 ns, I in function of momentum transfer q', 1)
plt.show()
