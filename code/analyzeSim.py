
import pandas, mdtraj, h5py
import numpy as np
import matplotlib.pyplot as plt
from success.forceFinder import resDist

def analyzeDist(amplitudes):
  ASAP3 = r'D:\GitHub\MD\models\ASAP3\ASAP3.pdb'
  ASAP4 = r'D:\GitHub\MD\models\ASAP4\ASAP4.pdb'
  pdbFiles = [ASAP3, ASAP4]
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

def analysis1(amplitudes):
  ASAP3 = r'D:\GitHub\MD\models\ASAP3\ASAP3.pdb'
  ASAP4 = r'D:\GitHub\MD\models\ASAP4\ASAP4.pdb'
  pdbNames = ['ASAP3', 'ASAP4']
  pdbFiles = [ASAP3, ASAP4]
  fig, axs = plt.subplots(4, len(amplitudes), figsize=(len(amplitudes)*3, 4*3))
  lims = [None] * 4
  data = []
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
      # handle dcd data
      resPair = (99,357) # ASP- ARG+
      f = h5py.File(fileName+'.h5', "r")
      dist = f['dist'+str(resPair)][:].flatten()
      axs[0,j].set_title('Ex= %d sin(2pi %dGHz t)'%(amp, fGHz))
      axs[0,j].set_xlabel('Time [ps]')
      axs[0,j].plot(t, Ex)
      axs[1,j].set_title(Ep)
      axs[1,j].plot(t, df[Ep], label=pdbNames[i])
      axs[1,j].legend()
      axs[2,j].set_title(Ek)
      axs[2,j].plot(t, df[Ek], label=pdbNames[i])
      axs[2,j].legend()
      axs[3,j].set_title('resDist%s [nm]'% str(resPair))
      axs[3,j].plot(t, dist, label=pdbNames[i])
      axs[3,j].legend()
      data.append([Ex, df[Ep], df[Ek], dist])
  data = np.array(data, dtype=object)
  print(data.shape)
  customLim = np.full([4, 2], None)
  customLim[2,0] = 26000
  for k in range(4):
    ymin = np.min(data, axis=(0,2))[k] if customLim[k,0] is None else customLim[k,0]
    ymax = np.max(data, axis=(0,2))[k] if customLim[k,1] is None else customLim[k,1]
    for j,amp in enumerate(amplitudes):
      axs[k,j].set_ylim([ymin, ymax])
  plt.tight_layout()
  plt.show()

amplitudes = np.linspace(60, 100, 5).astype(int)
#analyzeDist(amplitudes)
analysis1(amplitudes)
