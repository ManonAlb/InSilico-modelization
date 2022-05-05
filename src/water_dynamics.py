"""
Molecular Dynamics for the water box only, 15nm and 20ps simu AND water + fixed protein
"""

#MODULES

import os
import simtk.openmm as openmm

from simtk.openmm import LangevinIntegrator
import numpy as np

from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app import *
from openmmtools import testsystems

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
steps = 1000000
intervall = 1000

#SYSTEM CREATION

waterbox = testsystems.WaterBox(box_edge, cutoff)
system = waterbox.system
positions = waterbox.positions
topology = waterbox.topology

#INTEGRATOR

integrator = LangevinIntegrator( temperature, friction,time_step)
#SIMULATION

simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy()

#SAVE EVERY SNAPSHOTS

for i in range (1,int(steps/intervall)+1):
    simulation.reporters = [PDBReporter( ROOT + '/DYNAMIC/Fixed_protein_in_water/FRAMES'+ pdb_title + '_{}.pdb'.format(i), intervall)]

    #REPORTER

    simulation.reporters.append(StateDataReporter(stdout, intervall, step=True,
            potentialEnergy=True, temperature=True))

    simulation.step(intervall)
    #print(f"saved to {'{}'.format(DATA) + pdb_title + '{}.pdb'.format(i)}")


#simulation.reporters.append(PDBReporter(ROOT + pdb_title + '{}.pdb'.format(step), intervall)) #./DATA
