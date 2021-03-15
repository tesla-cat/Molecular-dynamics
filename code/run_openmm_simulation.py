
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

def comparePDB(pdb1, pdb2):
  print(pdb1.topology)
  print(pdb2.topology)
  residues1 = [r.name for r in pdb1.topology.residues()]
  residues2 = [r.name for r in pdb2.topology.residues()]
  for i in range(len(residues1)):
    if residues1[i] != residues2[i]:
      print(i, residues1[i], '=>', residues2[i])

def runSim(
  pdb = r'D:\MD\ASAP-sim\models-charmm-gui\4g7y-VSD.pdb',
  reportSteps = 1000,
  steps = 10000,
):
  pdb = PDBFile(pdb)
  forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
  system = forcefield.createSystem(
    pdb.topology, 
    nonbondedMethod=CutoffNonPeriodic,
    constraints=HBonds
  )
  T, friction, dt = 300*kelvin, 1/picosecond, 0.004*picoseconds
  integrator = LangevinMiddleIntegrator(T, friction, dt)
  simulation = Simulation(pdb.topology, system, integrator)
  simulation.context.setPositions(pdb.positions)
  print('minimize Energy')
  simulation.minimizeEnergy()
  print('run simulation')
  simulation.reporters.extend([
    StateDataReporter(
      stdout, reportSteps, 
      step=True, potentialEnergy=True, temperature=True
    ),
    DCDReporter('%s.dcd' % pdb, reportSteps),
  ])
  simulation.step(steps)