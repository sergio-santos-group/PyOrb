import numpy as np
from gaussianset import  *

class GaussianFCHK:
  def __init__(self, fname):

    # self.atomicNumbers = []
    # self.xyz = []
    #
    # fin = open(fname, 'r')
    # self.parse_line = self.parse_field
    # for line in fin:
    #   self.parse_line(line)
    # fin.close()
    self.parsedFields = {}

    self.fin = open(fname, 'r')

    line = self.fin.readline()
    while line:
      self.parse_line(line)
      line = self.fin.readline()

    self.fin.close()

  def parse_line(self, line):
    if line.startswith('Number of electrons'):
      self.parsedFields['Number of electrons'] = 'nElectrons'
      self.nElectrons = int(line.split()[4])
    elif line.startswith('Number of basis functions'):
      self.parsedFields['Number of basis functions'] = 'nBasisFcns'
      self.nBasisFcns = int(line.split()[5])
    elif line.startswith('Atomic numbers'):
      self.parsedFields['Atomic numbers'] = 'atomicNumbers'
      self.atomicNumbers = self.read_int_array(int(line.split()[4]))

    elif line.startswith('Current cartesian coordinates'):
      self.parsedFields['Current cartesian coordinates'] = 'coordinates'
      self.coordinates = self.read_float_array(int(line.split()[5])).reshape((-1,3))

    elif line.startswith('Shell types'):
      self.parsedFields['Shell types'] = 'shellTypes'
      self.shellTypes = self.read_int_array(int(line.split()[4]))

    elif line.startswith('Number of primitives per shell'):
      self.parsedFields['Number of primitives per shell'] = 'shellNums'
      self.shellNums = self.read_int_array(int(line.split()[7]))

    elif line.startswith('Shell to atom map'):
      self.parsedFields['Shell to atom map'] = 'shellToAtom'
      self.shellToAtom = self.read_int_array(int(line.split()[6]))

    #--------------------------------------------------
    elif line.startswith('Primitive exponents'):
      self.parsedFields['Primitive exponents'] = 'a'
      self.a = self.read_float_array(int(line.split()[4]))

    elif line.startswith('Contraction coefficients'):
      self.parsedFields['Contraction coefficientss'] = 'c'
      self.c = self.read_float_array(int(line.split()[4]))

    elif line.startswith('P(S=P) Contraction coefficients'):
      self.parsedFields['P(S=P) Contraction coefficients'] = 'csp'
      self.csp = self.read_float_array(int(line.split()[5]))

    elif line.startswith('Alpha Orbital Energies'):
      self.parsedFields['Alpha Orbital Energies'] = 'orbitalEnergy'
      self.orbitalEnergy = self.read_float_array(int(line.split()[5]))

    elif line.startswith('Alpha MO coefficients'):
      self.parsedFields['Alpha MO coefficients'] = 'MOCoefficients'
      self.MOCoefficients = self.read_float_array(int(line.split()[5]))

    elif line.startswith('Total SCF Density'):
      self.parsedFields['Total SCF Density'] = 'density'
      self.density = self.read_float_array(int(line.split()[5]))


  def read_int_array(self, size):
    n = 0
    array = []
    while n < size:
      tmp = self.fin.readline().split()
      array.extend(tmp)
      n += len(tmp)
    return np.array(array, dtype=int)

  def read_float_array(self, size):
    n = 0
    array = []
    while n < size:
      tmp = self.fin.readline().split()
      array.extend(tmp)
      n += len(tmp)
    return np.array(array, dtype=float)


  shellType2OrbType = {
     0: OrbitalType.S,
     1: OrbitalType.P,
    #-1: OrbitalType.P,
     2: OrbitalType.D,
    -2: OrbitalType.D5,
     3: OrbitalType.F,
    -3: OrbitalType.F7,
     4: OrbitalType.G,
    -4: OrbitalType.G9,
     5: OrbitalType.H,
    -5: OrbitalType.H11,
     6: OrbitalType.I,
    -6: OrbitalType.I13,
  }

  def load(self):
    gaussianSet = GaussianSet(self.nElectrons, self.atomicNumbers, self.coordinates)

    nGTO = 0

    for i,shellType in enumerate(self.shellTypes):
      if shellType == -1:
        tmpGTO = nGTO
        s = gaussianSet.add_basis(self.shellToAtom[i]-1, OrbitalType.S)
        for j in xrange(self.shellNums[i]):
          gaussianSet.add_GTO(s, self.c[nGTO], self.a[nGTO])
          nGTO += 1
        p = gaussianSet.add_basis(self.shellToAtom[i]-1, OrbitalType.P)
        for j in xrange(self.shellNums[i]):
          gaussianSet.add_GTO(p, self.csp[tmpGTO], self.a[tmpGTO])
          tmpGTO += 1
      else:

        orbtype = self.shellType2OrbType.get(shellType, None)
        if (orbtype is not None):
          b = gaussianSet.add_basis(self.shellToAtom[i]-1, orbtype)
          for j in xrange(self.shellNums[i]):
            gaussianSet.add_GTO(b, self.c[nGTO], self.a[nGTO])
            nGTO += 1
    if gaussianSet.is_valid():
      if self.MOCoefficients.size > 0:
        gaussianSet.set_coefficients(self.MOCoefficients)
      if self.density.size > 0:
        gaussianSet.set_density(self.density)

    return gaussianSet


  def __str__(self):
    s = ''
    for key,val in self.parsedFields.iteritems():
      attr = getattr(self, val)
      if isinstance(attr, np.ndarray):
        N = attr.size
        t = 'R' if attr.dtype.name.startswith('float') else 'I'
      elif isinstance(attr, int):
        N = attr
        t = 'I'
      s += '%-42s %s   N=%12d\n'%(key,t,N)
    return s

  def __repr__(self):
     return self.__str__()





if __name__ == '__main__':
  #g = GaussianFCHK('/Users/ssantos/PycharmProjects/cube/water/water.fchk')
  g = GaussianFCHK('/Users/ssantos/terpiridinas/2Py/2Py.fchk')
  print g
  gs = g.load()
  for n in ['homo-1','homo','lumo','lumo+1']:
    mo = gs.get_cube(n, spacing = 0.25)
    if mo is not None:
      mo.write('ccube_'+n+'.cube')