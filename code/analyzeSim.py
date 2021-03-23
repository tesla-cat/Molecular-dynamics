
import pandas, mdtraj, h5py
import numpy as np
import matplotlib.pyplot as plt
from success.forceFinder import resDist

def analyzeSim(pdbFile, simId):
  pdbFile = pdbFile.replace('.pdb', '-'+simId)
  df = pandas.read_csv(pdbFile+'.csv', index_col=0)
  df.plot(subplots=True, layout=(2,1))
  plt.show()

def analysis1():
  ASAP3 = r'D:\GitHub\MD\models\ASAP3\ASAP3.pdb'
  ASAP4 = r'D:\GitHub\MD\models\ASAP4\ASAP4.pdb'
  pdbNames = ['ASAP3', 'ASAP4']
  pdbFiles = [ASAP3, ASAP4]
  amplitudes = np.linspace(50, 100, 6).astype(int)
  fig, axs = plt.subplots(4, len(amplitudes), figsize=(len(amplitudes)*3, 4*3))
  for i,pdbFile in enumerate(pdbFiles):
    for j,amp in enumerate(amplitudes):
      # handle csv data
      simId = str(amp)
      fileName = pdbFile.replace('.pdb', '-'+simId)
      df = pandas.read_csv(fileName+'.csv', index_col=0)
      t = df.index
      t -= t[0]
      T = 1/2 * t[-1]
      fGHz = 1/(T*1e-12) / 1e9
      Ex = amp*np.sin(2*np.pi*(1/T)*t)
      Ep = 'Potential Energy (kJ/mole)'
      Ek = 'Kinetic Energy (kJ/mole)'
      axs[0,j].set_title('Ex=A sin(2pi f t), f=%d GHz'%fGHz)
      axs[0,j].set_xlabel('Time [ps]')
      axs[0,j].set_ylim([-120, 120])
      axs[0,j].plot(t, Ex)
      axs[1,j].set_title(Ep)
      axs[1,j].set_ylim([-1.5 * 1e6, -1.2 * 1e6])
      axs[1,j].plot(t, df[Ep], label=pdbNames[i])
      axs[1,j].legend()
      axs[2,j].set_title(Ek)
      axs[2,j].set_ylim([2.5 * 1e5, 3 * 1e5])
      axs[2,j].plot(t, df[Ek], label=pdbNames[i])
      axs[2,j].legend()
      # handle dcd data
      resPair = (99,357) # ASP- ARG+
      f = h5py.File(fileName+'.h5', "r")
      dist = f['dist'+str(resPair)]
      axs[3,j].set_title('resDist%s [nm]'% str(resPair))
      #axs[2,j].set_ylim([0.5, 0.9])
      axs[3,j].plot(t, dist, label=pdbNames[i])
      axs[3,j].legend()
  plt.tight_layout()
  plt.show()

def analyzeDist():
  ASAP3 = r'D:\GitHub\MD\models\ASAP3\ASAP3.pdb'
  ASAP4 = r'D:\GitHub\MD\models\ASAP4\ASAP4.pdb'
  pdbFiles = [ASAP3, ASAP4]
  amplitudes = np.linspace(50, 100, 6).astype(int)
  for i,pdbFile in enumerate(pdbFiles):
    for j,amp in enumerate(amplitudes):
      print(i,j)
      simId = str(amp)
      fileName = pdbFile.replace('.pdb', '-'+simId)
      traj = mdtraj.load(fileName+'.dcd', top=pdbFile)
      resPair = (99,357) # ASP- ARG+
      # http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/
      # http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/#MLformula
      atomPair = [
        traj.topology.select("protein and resSeq 99 and name OD1")[0], 
        traj.topology.select("protein and resSeq 357 and name NH1")[0], 
      ]
      print('atomPair', atomPair)
      dist = mdtraj.compute_distances(traj, [atomPair])
      f = h5py.File(fileName+'.h5', "w")
      f.create_dataset('dist'+str(resPair), data=dist)
      f.close()

#analyzeDist()
analysis1()
