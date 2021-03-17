
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.signal import savgol_filter
import mdtraj
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

def getResAtoms(topo, seq):
  atoms = topo.select("protein and resSeq %d" % seq)
  res = topo.atom(atoms[0]).residue
  return res, atoms

def resDist(traj, topo, seq1, seq2):
  res1, atoms1 = getResAtoms(topo, seq1)
  res2, atoms2 = getResAtoms(topo, seq2)
  center1 = np.mean(traj.xyz[:, atoms1, :], axis=1)
  canter2 = np.mean(traj.xyz[:, atoms2, :], axis=1)
  dist = center1 - canter2
  dist = (dist[:,0]**2 + dist[:,1]**2 + dist[:,2]**2)**0.5
  return res1, res2, dist

def forceFinder(topoFile, trajFiles):
  data = []
  for trajFile in trajFiles:  
    traj = mdtraj.load(trajFile, top=topoFile)
    topo = traj.topology
    # pymol: select chain X and (resi 1 or resi 365) 
    res1, res2, dist = resDist(traj, topo, 1, 365)  
    dist = savgol_filter(dist, 61, 2)
    v = np.gradient(dist)
    a = np.gradient(v)
    data.append([dist, v, a])
  plot(data, res1, res2)

def plot(data, res1, res2):
  fig, axs = plt.subplots(3, len(data), figsize=(16,8))
  for n in range(len(data)):
    dist, v, a = data[n]
    axs[0, n].set_title('traj %d' % (n+1))
    axs[0, n].set_xlabel('frame')
    axs[0, n].set_ylabel('dist [nm]') 
    axs[0, n].plot(dist)
    axs[1, n].set_ylabel('velocity')
    axs[1, n].plot(v)
    axs[2, n].set_ylabel('acceleration')
    axs[2, n].plot(a)
  fig.suptitle('distance between %s and %s' % (res1, res2))
  plt.tight_layout(); plt.show()
    
if __name__=='__main__':
  pdbFile = r'D:\GitHub\MD\models\ASAP3\step3_pbcsetup.pdb'
  forceFinder(pdbFile, ['%s.dcd' % pdbFile]*2)