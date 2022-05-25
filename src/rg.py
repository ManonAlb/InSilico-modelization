import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np


t = md.load_pdb('/home/manon/Documents/Protein_in_vacuum.pdb')
tw = md.load('/home/manon/Documents/moving_protein_in_water.pdb').remove_solvent()
t.save_dcd('/home/manon/Documents/ubiq.dcd')

traj = md.load('/home/manon/Documents/ubiq.dcd', top = '/home/manon/Documents/Protein_in_vacuum.pdb')



rgw = md.compute_rg(tw)
rg=md.compute_rg(traj)
plt.plot(np.arange(1,1000,1),rg,  '-')
plt.plot( np.arange(1,1000,1),rgw, '-')
plt.savefig('rg.png',dpi=300)
plt.show()
