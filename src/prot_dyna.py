import os
import pdbfixer
from sys import stdout
import openmm
from openmm.app import *
from openmm import *
from openmm.unit import *

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
pdb_path = ROOT + '/DYNAMIC/Protein/1ubq.pdb'


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

pdb_title='/{}_in_vacuum_for_{}fs'.format(pdbid, str(step_number/intervall))  # Name of the output pdb

#Force field parameters and modellers

forcefield = app.ForceField("amber14-all.xml")

#modeller

modeller = app.Modeller(prepared_protein.topology,prepared_protein.positions)
modeller.deleteWater()
modeller.addHydrogens(forcefield)

#system

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer,constraints=HBonds)

integrator = LangevinIntegrator( 298*kelvin, 1/picosecond,1.0*femtosecond)

#simulation

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=300)



for i in range (1,int(step_number/intervall)+1):
    simulation.reporters = [PDBReporter( ROOT + '/DYNAMIC/Protein/FRAMES'+ pdb_title + '_{}.pdb'.format(i), intervall)]

    #REPORTER

    simulation.reporters.append(StateDataReporter(stdout, intervall, step=True,
            potentialEnergy=True, temperature=True))

    simulation.step(intervall)

#simulation.saveState('{}'.format(HERE) + '/data/' + 'simulation_vacuum.xml')
