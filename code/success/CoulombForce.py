
import mdtraj, json
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

def findCharges(pdbFile, log=False):
  """
  http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
  Aspartic acid 天冬氨酸 - OD1(2) OD2(2)
  Glutamic acid	谷氨酸   - OE1(2) OE2(2)
  Lysine        赖氨酸   + NZ(3)
  Arginine      精氨酸	 + NE     NH1(2) NH2(2)
  Histidine     组氨酸	 + ND1    NE2
  """
  pdb = mdtraj.load_pdb(pdbFile)
  topo = pdb.topology
  charges = {
    'ASP': { 'OD1': -2, 'OD2': -2 },
    'GLU': { 'OE1': -2, 'OE2': -2 },
    'LYS': { 'NZ' :  3 },
    'ARG': { 'NE' :  1, 'NH1':  2, 'NH2': 2 },
    'HIS': { 'ND1':  1, 'NE2':  1 },
  }
  for resname, names in charges.items():
    for name, charge in names.items():
      atoms = topo.select('resname %s and name %s'%(resname, name))
      charges[resname][name] = {
        'charge': charge, 'atoms': atoms.tolist()
      }
  if log: print(json.dumps(charges, indent=4))
  return charges

def findChargesCorrectly(pdbFile, log=False):
  pdb = PDBFile(pdbFile)
  pdb2 = mdtraj.load_pdb(pdbFile)
  protein = pdb2.topology.select('protein')
  forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
  system = forcefield.createSystem(
    pdb.topology, 
    nonbondedMethod=CutoffNonPeriodic,
    constraints=HBonds
  )
  nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
  for i in range(system.getNumParticles()):
    charge, sigma, epsilon = nonbonded.getParticleParameters(i)
    if i in protein: 
      print(pdb2.topology.atom(i), charge)

def addForce(pdbFile):
  potential = "-q*Ex*x"
  force = CustomExternalForce(potential)
  force.addPerParticleParameter("q")
  force.addGlobalParameter("Ex", 0)
  charges = findCharges(pdbFile)
  for resname, names in charges.items():
    for name, obj in names.items():
      for atomIndex in obj['atoms']:
        print(resname, name, atomIndex, obj['charge'])
        force.addParticle(atomIndex, [obj['charge']])
  return force

if __name__=='__main__':
  pdbFile = r'D:\GitHub\MD\models\ASAP3\step3_pbcsetup.pdb'
  #findCharges(pdbFile, log=True)
  #addForce(pdbFile)
  findChargesCorrectly(pdbFile)