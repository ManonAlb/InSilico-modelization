"""
Molecular Dynamics for the water box only, 15nm and 20ps simu AND water + fixed protein
"""

#MODULES

import os
import simtk.openmm as openmm
import simtk.openmm.app as app


from simtk.openmm import LangevinIntegrator
import numpy as np

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
#from openmmtools import testsystems

from sys import stdout

#DIRECTORY

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))

#SIMULATION PARAMETERS

temperature = 298*kelvin
pdb_title ='/Water_box'
time_step = 1.0*femtosecond
friction = 1/picosecond
cutoff=1*nanometers
box_edge=15*nanometers
step_number = 1000000
intervall = 1000

#PLATFORM
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'} 
properties["DeviceIndex"] = "0"


#Forcefields
force_field=app.ForceField('tip3p.xml')

#empty
top = app.Topology()
pos = Quantity((), angstroms)

# Create new Modeller instance.
m = app.Modeller(top, pos)
# Add solvent to specified box dimensions.
boxSize = Quantity(np.ones([3]) * box_edge/nanometers, nanometers)
m.addSolvent(force_field, boxSize=boxSize)


# New positions
newtop = m.getTopology()
newpos = m.getPositions()


# Create a system
system = force_field.createSystem(newtop, nonbondedCutoff=cutoff, constraints=app.HBonds,nonbondedMethod=app.PME)

#SIMULATION
integrator = LangevinIntegrator( temperature, friction,time_step)
simulation = Simulation(newtop, system, integrator, platform, properties)
simulation.context.setPositions(newpos)
simulation.minimizeEnergy()

#SAVE EVERY SNAPSHOTS

simulation.reporters = [StateDataReporter(stdout, intervall, step=True, potentialEnergy=True, temperature=True)]

n_saved = 0
for i in range(step_number):
    simulation.step(1)

    if i % intervall == 0 or i == step_number - 1:
        reporter = PDBReporter(ROOT + '/DYNAMIC/Water/FRAMES' + pdb_title + '_{}.pdb'.format(n_saved), 1)
        state = simulation.context.getState(getPositions=True)
        reporter.report(simulation, state)

        n_saved += 1
