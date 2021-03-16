
from simtk.openmm.app import PDBFile

def comparePDB(pdb1, pdb2):
  print(pdb1.topology)
  print(pdb2.topology)
  residues1 = [r.name for r in pdb1.topology.residues()]
  residues2 = [r.name for r in pdb2.topology.residues()]
  for i in range(len(residues1)):
    if residues1[i] != residues2[i]:
      print(i, residues1[i], '=>', residues2[i])

if __name__=='__main__':
  pdb1 = PDBFile(r'D:\GitHub\MD\models-raw\ASAP3.pdb')
  pdb2 = PDBFile(r'D:\GitHub\MD\models-raw\ASAP4.pdb')
  comparePDB(pdb1, pdb2)
