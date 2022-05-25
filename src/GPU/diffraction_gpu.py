"""
Generate diffraction data using a Multi-Slice algorithm
"""

#MODULES
import matplotlib.pyplot as plt
import os
import abtem

from ase.io import read
from abtem.visualize import show_atoms
from abtem.potentials import Potential
from abtem.waves import PlaneWave
from abtem.measure import block_zeroth_order_spot


#DIRECTORY

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
water_path = ROOT + '/DYNAMIC/Fixed_protein_in_water/FRAMES'
diffraction_path = ROOT + '/DIFFRACTION/Fixed_protein_in_water/'

files = os.listdir(water_path)
for file in files:
    if file.endswith(('.pdb')):

        system = read(water_path + '/' + file )
        system.center() #show_atoms(system, plane='yz')
        #WAVE CHARACTERISTICS

        wave = PlaneWave(energy=3e5, sampling=0.05, device='gpu')
        potential = Potential(system, parametrization='lobato', sampling=0.05, slice_thickness=2, projection='finite', device='gpu')
        wave.match_grid(potential)
#RUN

        exit_wave = wave.multislice(potential)
        pw_diffraction_pattern = exit_wave.diffraction_pattern()

#KEEP IN MEMORY
        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 3))
        pw_diffraction_pattern.write(diffraction_path + '/'+file.strip('.pdb') + '.hdf5', mode='w', format='hdf5')
        block_zeroth_order_spot(pw_diffraction_pattern).show(ax=ax1, power=0.2, cmap='jet') #blockin zero order for output image

        plt.savefig(diffraction_path + file.strip('.pdb') + '.png',bbox_inches='tight', dpi=300)
        plt.close(fig)
            
alphax = pw_diffraction_pattern.calibrations[-2].coordinates(pw_diffraction_pattern.array.shape[-2]) #alpha x in mrad

size = alphax.shape[0]
