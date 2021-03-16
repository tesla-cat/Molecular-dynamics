
import os, subprocess

for root, subdirs, files in os.walk('.'):
  for f in files:
    if f.endswith('.nc'):
      filePath = os.path.join(root, f)
      print(subprocess.run([
        "mdconvert", 
        filePath, 
        "-o",
        filePath.replace('.nc', '.dcd')
      ]))