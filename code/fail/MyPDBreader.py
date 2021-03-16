
class MyPDBreader:
  def __init__(self, pdbFile):
    with open(pdbFile) as f:
      lines = f.readlines()
    for line in lines:
      if line.startswith('ATOM'):
        self.addAtom(line)

  def addAtom(self, line):
    recordName = line[0:6]
    atomSerialNumber = line[6:11]
    atomName = line[12:16]
    alternateLocationIndicator = line[16]
    residueName = line[17:20]
    chainIdentifier = line[21]
    residueSequenceNumber = line[22:26]
    codeForInsertionOfResidues = line[26]
    x = line[30:38]
    y = line[38:46]
    z = line[46:54]
    occupancy = line[54:60]
    temperatureFactor = line[60:66]
    segmentIdentifier = line[72:76]
    elementSymbol = line[76:78]
    charge = line[78:80]
    print(atomSerialNumber, residueName, x, y, z)