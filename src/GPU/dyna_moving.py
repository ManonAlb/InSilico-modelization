
import os
import pdbfixer
from sys import stdout
import simtk.openmm
import simtk.openmm.app as app
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
pdb_path = ROOT + '/DYNAMIC/Moving_protein_in_water/1ubq.pdb'


def prepare_protein(
    pdb_file, ignore_missing_residues=True, ignore_terminal_missing_residues=True, ph=7.0
):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.
    """
    fixer = pdbfixer.PDBFixer(str(pdb_file))
    fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer
    fixer.findMissingResidues()  # identify missing residues, needed for identification of missing atoms

    # if missing terminal residues shall be ignored, remove them from the dictionary
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    # if all missing residues shall be ignored ignored, clear the dictionary
    if ignore_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()  # find non-standard residue
    fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues
    fixer.addMissingHydrogens(ph)  # add missing hydrogens
    return fixer

# prepare protein and build only missing non-terminal residues
prepared_protein = prepare_protein(pdb_path, ignore_missing_residues=False)




pdbid='1ubq'
step_number=1000000 #duration time = time_step*step number =4ps
intervall=1000 #save every intervall
time_step = 1.0*femtosecond
friction = 1/picosecond
temperature =  298*kelvin


pdb_title='/{}_in_water_for_{}fs'.format(pdbid, str(step_number*time_step))  # Name of the output 

#Force field parameters

""" From left to right, protein force field and water force field """


forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")


#Modeller

modeller = app.Modeller(prepared_protein.topology,prepared_protein.positions)

modeller.deleteWater()
modeller.addHydrogens(forcefield)
modeller.addSolvent(forcefield, boxSize=Vec3(15, 15, 15)* nanometers) #box dimension 15nm


#In  Water


#SIMULATION PARAMETERS


#pdb_title='/{}_in_water_for_{}s'.format(pdbid, str(step_number/intervall))  # Name of the output pdb


#Periodic boundary conditions are used as well as Hbonds constrainst during the system creation

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer,constraints=HBonds)

#Langevin integrator

integrator = LangevinIntegrator(temperature, friction ,time_step)

#Energy minimization

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=300)

#keep the minimization in memory
""""position = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, position, open('{}'.format(DATA) + '/minimized.pdb', 'w'))"""


simulation.reporters = [StateDataReporter(stdout, intervall, step=True, potentialEnergy=True, temperature=True)]


n_saved = 0
for i in range(step_number):
    simulation.step(1)

    if i % intervall == 0 or i == step_number - 1:
        reporter = PDBReporter(ROOT + '/DYNAMIC/Moving_protein_in_water/FRAMES' + pdb_title + '_{}.pdb'.format(n_saved), 1)
        state = simulation.context.getState(getPositions=True)
        reporter.report(simulation, state)

        n_saved += 1
