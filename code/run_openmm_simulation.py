
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import mdtraj, pandas
import matplotlib.pyplot as plt

def runSim(pdbFile, steps = 10000, reportSteps = 100):
  pdb = PDBFile(pdbFile)
  forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
  system = forcefield.createSystem(
    pdb.topology, nonbondedMethod=CutoffNonPeriodic, constraints=HBonds
  )
  # add Coulomb force along x axis, initially off
  force = CustomExternalForce("-q*Ex*x")
  force.addGlobalParameter("Ex", 0)
  force.addPerParticleParameter("q")
  nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
  for i in range(system.getNumParticles()):
    charge, sigma, epsilon = nonbonded.getParticleParameters(i)
    force.addParticle(i, [charge])
  system.addForce(force)
  # simulation
  T, friction, dt = 300*kelvin, 1/picosecond, 0.004*picoseconds
  integrator = LangevinMiddleIntegrator(T, friction, dt)
  simulation = Simulation(pdb.topology, system, integrator)
  simulation.context.setPositions(pdb.positions)
  print('minimize energy')
  simulation.minimizeEnergy()
  print('run simulation')
  simulation.reporters.extend([
    StateDataReporter(
      '%s.csv' % pdbFile, reportSteps, 
      time=True, potentialEnergy=True, temperature=True, 
    ),
    StateDataReporter(stdout, reportSteps, step=True),
    DCDReporter('%s.dcd' % pdbFile, reportSteps),
  ])
  simulation.step(int(steps/2))
  simulation.context.setParameter("Ex", 10)
  simulation.step(int(steps/2))

def analyzeSim(pdbFile):
  df = pandas.read_csv('%s.csv' % pdbFile, index_col=0)
  df.plot()
  plt.show()

pdbFile = r'D:\GitHub\MD\models\ASAP3\step3_pbcsetup.pdb'
if 1:
  runSim(pdbFile)
else:
  analyzeSim(pdbFile)