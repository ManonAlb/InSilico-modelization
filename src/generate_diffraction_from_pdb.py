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
protein_path = ROOT + '/DYNAMIC/Protein/FRAMES' # Where I keep all of my pdbs files
diffraction_path = ROOT + '/DIFFRACTION/Protein/'

files = os.listdir(protein_path)
for file in files:
    if file.endswith(('.pdb')):

        system = read(protein_path + '/' + file )
        system.center() #show_atoms(system, plane='yz')

        #WAVE CHARACTERISTICS

        wave = PlaneWave(energy=3e5, sampling=0.05, device='gpu')
        potential = Potential(system, parametrization='lobato', sampling=0.05, slice_thickness=1, projection='finite', device='gpu')
        wave.match_grid(potential)
#RUN

        exit_wave = wave.multislice(potential)
        pw_diffraction_pattern = exit_wave.diffraction_pattern()

#KEEP IN MEMORY THE .PNG AND .HDF5 FILES WITHOUT ZERO ORDER BLOCKED

        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 3))
        pw_diffraction_pattern.write(diffraction_path + '/'+file.strip('.pdb') + '.hdf5', mode='w', format='hdf5')
        block_zeroth_order_spot(pw_diffraction_pattern).show(ax=ax1, power=0.2, cmap='jet') #blockin zero order for output image

        plt.savefig(diffraction_path + file.strip('.pdb') + '.png',bbox_inches='tight', dpi=300)
        plt.close(fig)
