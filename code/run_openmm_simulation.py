
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import mdtraj
import numpy as np
import matplotlib.pyplot as plt
from success.comparePDB import comparePDB

class MySimulation:
  def __init__(self, pdbFile, dt):
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
    T, friction, dt = 30*kelvin, 1/picosecond, dt*picoseconds
    integrator = LangevinMiddleIntegrator(T, friction, dt)
    simulation = Simulation(pdb.topology, system, integrator)
    self.simulation, self.pdb, self.pdbFile = simulation, pdb, pdbFile
    
  def run(self, simId, ExList, stepsPerEx, reportInterval):
    simulation, pdb, pdbFile = self.simulation, self.pdb, self.pdbFile
    fileName = pdbFile.replace('.pdb', '-'+simId)
    print(fileName)
    simulation.reporters = [
      StateDataReporter(
        fileName+'.csv', reportInterval, 
        time=True, potentialEnergy=True, kineticEnergy=True, 
      ),
      StateDataReporter(stdout, reportInterval, step=True),
      DCDReporter(fileName+'.dcd', reportInterval),
    ]
    simulation.context.setPositions(pdb.positions)
    print('minimize energy')
    simulation.minimizeEnergy()
    print('run simulation')
    for Ex in ExList:
      simulation.context.setParameter("Ex", Ex)
      simulation.step(stepsPerEx)

def experiment1():
  ASAP3 = r'D:\GitHub\MD\models\ASAP3\ASAP3.pdb'
  ASAP4 = r'D:\GitHub\MD\models\ASAP4\ASAP4.pdb'
  pdbFiles = [ASAP3, ASAP4]
  amplitudes = np.linspace(150, 300, 4).astype(int)
  
  numSteps = 10000
  numEx = 400
  numFrames = 100
  dt = 0.004 

  stepsPerEx = int(numSteps/numEx)
  reportInterval = int(numSteps/numFrames)
  t = np.arange(numEx) * stepsPerEx * dt
  T = 1/2 * numSteps * dt
  Ex = np.sin(2*np.pi*(1/T)*t)
  print('len(Ex), stepsPerEx, reportInterval', len(Ex), stepsPerEx, reportInterval)
  plt.step(t, Ex)
  plt.show()
  print(amplitudes)
  if 1:
    for pdbFile in pdbFiles:
      mySim = MySimulation(pdbFile, dt)
      for amp in amplitudes:
        mySim.run(str(amp), amp*Ex, stepsPerEx, reportInterval)

experiment1()