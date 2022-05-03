#Dynamic

""" 1s simulation ( time_steps*steps ) of a 15 nm water box only at room temperature with an equilibration phase of approximately
10ps ( verify by plotting Temperature as a function of time )"""


#MODULES
import os
import openmm
import numpy as np
import openmm.app as app

from openmm import LangevinIntegrator
from openmm.app import *
from openmm import *
from openmm.unit import *

from sys import stdout

#DIRECTORY


ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__)))

"""SIMULATION"""

#parameters

box_edge=15*nanometers
cutoff=1*nanometers
steps = 10000
intervall = 500
time_step = 2.0*femtosecond
friction = 1/picosecond
temperature =  298*kelvin

pdb_title ='/Water_box'


#Forcefields

force_field=app.ForceField('tip3p.xml')

#empty
top = app.Topology()
pos = Quantity((), angstroms)

# Create new Modeller instance.
m = app.Modeller(top, pos)

# Add solvent to specified box dimensions.
boxSize = Quantity(np.ones([3]) * box_edge/box_edge, box_edge)
m.addSolvent(force_field, boxSize=boxSize)


# New positions
newtop = m.getTopology()
newpos = m.getPositions()

 # Convert positions to numpy.
#positions = Quantity(np.array(newpos / newpos), newpos)

# Create a system
system = force_field.createSystem(newtop, nonbondedCutoff=cutoff, constraints=app.HBonds,nonbondedMethod=app.PME)


# Get new topology and coordinates.
n_atoms = system.getNumParticles()


#water_box = [system,positions]

integrator = LangevinIntegrator( temperature, friction,time_step)
simulation = Simulation(newtop, system, integrator)
simulation.context.setPositions(newpos)
simulation.minimizeEnergy()

# Create a single PDB file that reports the atoms positions every intervall steps

simulation.reporters.append(PDBReporter(ROOT + pdb_title + '.pdb', intervall, enforcePeriodicBox=True))

simulation.reporters.append(StateDataReporter(stdout, intervall, step=True,potentialEnergy=True, temperature=True))

#Run
simulation.step(steps)
