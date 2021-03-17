
import os, subprocess

if 0:
  for root, subdirs, files in os.walk('.'):
    for f in files:
      if f.endswith('.dcd'):
        filePath = os.path.join(root, f)
        print(subprocess.run([
          "mdconvert", 
          filePath, 
          "-o",
          filePath.replace('.dcd', '.h5')
        ]))
else:
  print(subprocess.run([
    "mdconvert", 
    r'D:\GitHub\MD\models\ASAP3\step3_pbcsetup.pdb.dcd', 
    "-o",
    r'D:\GitHub\MD\models\ASAP3\step3_pbcsetup.xtc',
    #"-t",
    #r'D:\GitHub\MD\models\ASAP3\step3_pbcsetup.pdb',
  ]))