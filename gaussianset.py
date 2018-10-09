import numpy as np
from sys import stdout
#import _cube as ccube
import ccube


# https://github.com/cryos/avogadro/blob/master/libavogadro/src/extensions/surfaces/openqube/gaussianset.cpp

class OrbitalType:
  S, SP, P, D, D5, F, F7, G, G9, H, H11, I, I13, UU = range(14)




class GaussianSet:

  type2numMO = {
      OrbitalType.S   : 1,
      OrbitalType.P   : 3,
      OrbitalType.SP  : 4,
      OrbitalType.D   : 6,
      OrbitalType.D5  : 5,
      OrbitalType.F   : 8,
      OrbitalType.F7  : 7 }

  def __init__(self, nElectrons, atmNum, atmCrd):
    self.nElectrons = nElectrons
    self.atmNum = atmNum.copy()
    self.atmCrd = atmCrd.copy()
    self.numMOs = 0

    self.symmetry = []
    self.atomIndices = []
    self.gtoIndices = []
    self.gtoA = []
    self.gtoC = []
    self.initialized = False


  def add_basis(self, atom, orbType):
    self.numMOs += GaussianSet.type2numMO.get(orbType, 0)
    self.atomIndices.append(atom)
    self.symmetry.append(orbType)

  def add_GTO(self, basisID, c, a):
    if len(self.gtoIndices) < len(self.atomIndices):
      self.gtoIndices.append(len(self.gtoA))
    self.gtoA.append(a)
    self.gtoC.append(c)

  def set_coefficients(self, MOcoefficients):
    self.moMatrix = MOcoefficients.reshape((self.numMOs, self.numMOs)).T


  def init(self):
    if self.initialized:
      return


    gtoIndices = self.gtoIndices
    gtoA = self.gtoA
    gtoC = self.gtoC
    indexMO = 0

    cIndices = self.cIndices = []
    gtoCN = self.gtoCN = []

    moIndices = self.moIndices = np.zeros(len(self.symmetry), dtype=int)
    PI3 = np.pi**3

    gtoIndices.append(len(self.gtoA))
    for i,symmetry in enumerate(self.symmetry):
      size_gtoCN = len(gtoCN)
      if symmetry == OrbitalType.S:
        moIndices[i] = indexMO
        indexMO += 1
        #cIndices.append(size_gtoCN)
        #size_gtoCN += indexMO
        cIndices.append(size_gtoCN)
        for j in xrange(gtoIndices[i], gtoIndices[i+1]):
          gtoCN.append(gtoC[j] * np.power(gtoA[j], 0.75) * 0.71270547)

      elif symmetry == OrbitalType.P:
        moIndices[i] = indexMO
        indexMO += 3
        cIndices.append(size_gtoCN)
        for j in xrange(gtoIndices[i], gtoIndices[i+1]):
          p = gtoC[j] * np.power(gtoA[j], 1.25) * 1.425410941
          gtoCN.extend([p,p,p])

      elif symmetry == OrbitalType.D:
        moIndices[i] = indexMO
        indexMO += 6
        cIndices.append(size_gtoCN)
        for j in xrange(gtoIndices[i], gtoIndices[i+1]):
          p = gtoC[j] * np.power(gtoA[j], 1.75) * 1.645922781
          gtoCN.extend([p,p,p])
          p = gtoC[j] * np.power(gtoA[j], 1.75) * 2.850821881
          gtoCN.extend([p,p,p])

      elif symmetry == OrbitalType.D5:
        moIndices[i] = indexMO
        indexMO += 5
        cIndices.append(size_gtoCN)
        for j in xrange(gtoIndices[i], gtoIndices[i+1]):
          f = np.power(gtoA[j], 7.0) / PI3
          gtoCN.append(gtoC[j] * np.power(2048 * f/9, 0.25))
          gtoCN.append(gtoC[j] * np.power(2048 * f,   0.25))
          gtoCN.append(gtoC[j] * np.power(2048 * f,   0.25))
          gtoCN.append(gtoC[j] * np.power( 128 * f,   0.25))
          gtoCN.append(gtoC[j] * np.power(2048 * f,   0.25))

      elif symmetry == OrbitalType.F:
        moIndices[i] = indexMO
        indexMO += 10
        cIndices.append(size_gtoCN)
        for j in xrange(gtoIndices[i], gtoIndices[i+1]):
          f = gtoC[j] * np.power(gtoA[j], 2.25)
          gtoCN.append(f * 1.472158089299094)
          gtoCN.append(f * 3.291845561298979)
          gtoCN.append(f * 3.291845561298979)
          gtoCN.append(f * 3.291845561298979)
          gtoCN.append(f * 5.701643762839922)
          gtoCN.append(f * 3.291845561298979)
          gtoCN.append(f * 1.472158089299094)
          gtoCN.append(f * 3.291845561298979)
          gtoCN.append(f * 3.291845561298979)
          gtoCN.append(f * 1.472158089299094)

      elif symmetry == OrbitalType.F7:
        moIndices[i] = indexMO
        indexMO += 7
        cIndices.append(size_gtoCN)
        for j in xrange(gtoIndices[i], gtoIndices[i+1]):
          f = gtoC[j] * np.power(gtoA[j], 2.25) * 1.4721580892990935
          gtoCN.extend(7*[f])

      # elif symmetry == OrbitalType.F:
      #   moIndices[i] = indexMO
      #   indexMO += 10
      #   cIndices.append(size_gtoCN)
      #   size_gtoCN += indexMO
      # elif symmetry == OrbitalType.F7:
      #   moIndices[i] = indexMO
      #   indexMO += 7
      #   cIndices.append(size_gtoCN)
      #   size_gtoCN += indexMO
      # elif symmetry == OrbitalType.G:
      #   moIndices[i] = indexMO
      #   indexMO += 15
      #   cIndices.append(size_gtoCN)
      #   size_gtoCN += indexMO
      # elif symmetry == OrbitalType.G9:
      #   moIndices[i] = indexMO
      #   indexMO += 9
      #   cIndices.append(size_gtoCN)
      #   size_gtoCN += indexMO
      # elif symmetry == OrbitalType.H:
      #   moIndices[i] = indexMO
      #   indexMO += 21
      #   cIndices.append(size_gtoCN)
      #   size_gtoCN += indexMO
      # elif symmetry == OrbitalType.H11:
      #   moIndices[i] = indexMO
      #   indexMO += 11
      #   cIndices.append(size_gtoCN)
      #   size_gtoCN += indexMO
      # elif symmetry == OrbitalType.I:
      #   moIndices[i] = indexMO
      #   indexMO += 28
      #   cIndices.append(size_gtoCN)
      #   size_gtoCN += indexMO
      # elif symmetry == OrbitalType.I13:
      #   moIndices[i] = indexMO
      #   indexMO += 13
      #   cIndices.append(size_gtoCN)
      #   size_gtoCN += indexMO
      else:
        print 'Unknown basis set'

    self.gtoCN = np.array(self.gtoCN)
    self.gtoA  = np.array(self.gtoA)
    self.gtoC  = np.array(self.gtoC)
    self.initialized = True


  def get_cube(self, state, *args, **kargs):
    if isinstance(state, str):
      state = state.lower()
      if state.startswith('homo'):
        s = self.nElectrons/2
      elif state.startswith('lumo'):
        s = self.nElectrons/2+1
      else:
        print 'Cube.get_cube: unknown orbital %s.'%state
        return None
      if '-' in state:
        s -= int(state.split('-')[1])
      elif '+' in state:
        s += int(state.split('+')[1])
      state = s

    self.init()

    cube = Cube(self.atmNum, self.atmCrd, *args, **kargs)

    status = ccube.get_cube(cube.get_points().flatten(),
                   cube.data,
                   self.atmCrd.flatten(),
                   self.moMatrix[:,state-1].copy(),
                   self.gtoA,
                   self.gtoCN,
                   self.symmetry,
                   self.atomIndices,
                   self.moIndices,
                   self.gtoIndices,
                   self.cIndices)
    # print status

    # cube = Cube(self.atmNum, self.atmCrd, *args, **kargs)
    # points = cube.get_points()
    # nPoints = points.shape[0]/50
    #
    # stdout.write('Orbital %s: '%state)
    # stdout.flush()
    # for n,pt in enumerate(points):
    #   if n%nPoints == 0:
    #     stdout.write('.')
    #     stdout.flush()
    #   r = pt - self.atmCrd
    #   rSq = np.square(r).sum(1)
    #   tmp = 0.0
    #   for i,symmetry in enumerate(self.symmetry):
    #     idx = self.atomIndices[i]
    #     if   symmetry == OrbitalType.S:  tmp += self.pointS( i,         rSq[idx], state-1)
    #     elif symmetry == OrbitalType.P:  tmp += self.pointP( i, r[idx], rSq[idx], state-1)
    #     elif symmetry == OrbitalType.D:  tmp += self.pointD( i, r[idx], rSq[idx], state-1)
    #     elif symmetry == OrbitalType.F:  tmp += self.pointF( i, r[idx], rSq[idx], state-1)
    #     elif symmetry == OrbitalType.D5: tmp += self.pointD5(i, r[idx], rSq[idx], state-1)
    #     elif symmetry == OrbitalType.F7: tmp += self.pointF7(i, r[idx], rSq[idx], state-1)
    #     else:
    #       print 'unknown orbtype'
    #   cube.data[n] = tmp
    #
    # stdout.write(' Done!\n')
    # stdout.flush()

    return cube

  def pointS(self, moIndex, rSq, indexMO):
    coef = self.moMatrix[self.moIndices[moIndex], indexMO]
    if abs(coef) < 1e-15:
      return 0.0

    tmp = 0.0
    cIndex = self.cIndices[moIndex]
    for i in xrange(self.gtoIndices[moIndex],self.gtoIndices[moIndex+1]):
      tmp += self.gtoCN[cIndex] * np.exp(-self.gtoA[i] * rSq)
      cIndex += 1
    return tmp * coef

  def pointP(self, moIndex, r, rSq, indexMO):
    baseIndex = self.moIndices[moIndex]
    tmp = np.zeros(3)
    cIndex = self.cIndices[moIndex]
    for i in xrange(self.gtoIndices[moIndex],self.gtoIndices[moIndex+1]):
      tmp += self.gtoCN[cIndex:cIndex+3] * r * np.exp(-self.gtoA[i] * rSq)
      cIndex += 3

    return tmp.dot(self.moMatrix[baseIndex:baseIndex+3,indexMO])

  def pointD(self, moIndex, r, rSq, indexMO):
    baseIndex = self.moIndices[moIndex]
    tmp = np.zeros(6)
    cIndex = self.cIndices[moIndex]
    for i in xrange(self.gtoIndices[moIndex],self.gtoIndices[moIndex+1]):
      tmp += self.gtoCN[cIndex:cIndex+6] * np.exp(-self.gtoA[i] * rSq)
      cIndex += 6

    return tmp.dot(self.moMatrix[baseIndex:baseIndex+6,indexMO] * r[[0,1,2,0,0,1]] * r[[0,1,2,1,2,2]])


  def pointF(self, moIndex, r, rSq, indexMO):
    baseIndex = self.moIndices[moIndex]
    tmp = np.zeros(10)
    cIndex = self.cIndices[moIndex]
    for i in xrange(self.gtoIndices[moIndex],self.gtoIndices[moIndex+1]):
      tmp += self.gtoCN[cIndex:cIndex+10] * np.exp(-self.gtoA[i] * rSq)
      cIndex += 10

    return tmp.dot(self.moMatrix[baseIndex:baseIndex+10,indexMO] * r[[0,0,0,0,0,0,1,1,1,2]] *
                                                                   r[[0,0,0,1,1,2,1,1,2,2]] *
                                                                   r[[0,1,2,1,2,2,1,2,2,2]] )

  def pointD5(self, moIndex, r, rSq, indexMO):
    baseIndex = self.moIndices[moIndex]
    tmp = np.zeros(5)
    cIndex = self.cIndices[moIndex]
    for i in xrange(self.gtoIndices[moIndex],self.gtoIndices[moIndex+1]):
      tmp += self.gtoCN[cIndex:cIndex+5] * np.exp(-self.gtoA[i] * rSq)
      cIndex += 5

    x,y,z = r
    coefs = self.moMatrix[baseIndex:baseIndex+5, indexMO]
    coefs[0] *= (3*z*z - rSq)
    coefs[1] *= x*z
    coefs[2] *= y*z
    coefs[3] *= (x*x - y*y)
    coefs[4] *= x*y

    return tmp.dot(coefs)


  def pointF7(self, moIndex, r, rSq, indexMO):
    baseIndex = self.moIndices[moIndex]
    tmp = np.zeros(7)
    cIndex = self.cIndices[moIndex]
    for i in xrange(self.gtoIndices[moIndex],self.gtoIndices[moIndex+1]):
      tmp += self.gtoCN[cIndex:cIndex+7] * np.exp(-self.gtoA[i] * rSq)
      cIndex += 7

    x,y,z = r
    #root6 = 2.449489742783178
    #root60 = 7.745966692414834
    #root360 = 18.973665961010276
    coefs = self.moMatrix[baseIndex:baseIndex+7, indexMO]
    coefs[0] *= (z*z*z - 3.0/2.0 * (x*x*z + y*y*z))
    coefs[1] *= ((6.0 * x*z*z - 3.0/2.0 * (x*x*x + x*y*y))/2.449489742783178)
    coefs[2] *= ((6.0 * y*z*z - 3.0/2.0 * (x*x*y + y*y*y))/2.449489742783178)
    coefs[3] *= ((15.0 * (x*x*z - y*y*z))/7.745966692414834)
    coefs[4] *= ((30.0 * x*y*z)/7.745966692414834)
    coefs[4] *= ((15.0 * x*x*x - 45.0 * x*y*y)/18.973665961010276)
    coefs[4] *= ((45.0 * x*x*y - 15.0 * y*y*y)/18.973665961010276)

    return tmp.dot(coefs)

  def set_density(self, density):
    pass

  def is_valid(self):
    return True



class Cube:
  def __init__(self, atomNums, atomCoords, padding=3.0, spacing=0.25):
    self.padding = padding
    self.spacing = spacing

    self.lMin = atomCoords.min(0)-padding
    self.lMax = atomCoords.max(0)+padding
    self.atomCoords = atomCoords.copy()
    self.atomNums = atomNums.copy()


  def get_points(self):
    from itertools import product

    self.size = 1+((self.lMax-self.lMin)/self.spacing).astype(int)
    x = np.linspace(self.lMin[0], self.lMax[0], self.size[0])
    y = np.linspace(self.lMin[1], self.lMax[1], self.size[1])
    z = np.linspace(self.lMin[2], self.lMax[2], self.size[2])

    self.effSpacing = np.array([x[1]-x[0], y[1]-y[0], z[1]-z[0]])
    self.points = np.array([pt for pt in product(x,y,z)])
    self.data = np.zeros(self.size.prod(), dtype=float)
    return self.points# / 0.529177249


  def as_cube(self):
    s = 'Molecular orbital\nMO coefficients\n%5d%12.6f%12.6f%12.6f\n'
    s = s%(-self.atomNums.size, self.lMin[0], self.lMin[1], self.lMin[2])

    for i,v in enumerate(np.diag(self.effSpacing)):
      s += '%5d%12.6f%12.6f%12.6f\n'%(self.size[i],v[0],v[1],v[2])

    for z,crd in zip(self.atomNums, self.atomCoords):
      s += '%5d%12.6f%12.6f%12.6f%12.6f\n'%(z,z,crd[0],crd[1],crd[2])

    s += '%5d%5d\n'%(1,1)

    nz = self.size[2]
    
    k = 0
    for n,pt in enumerate(self.data):
      k += 1
      if k%6 == 0:
        s += '%13.5E\n'%pt
        k = 0
      elif (n+1)%nz == 0:
        s += '%13.5E\n'%pt
        k = 0
      else:
        s += '%13.5E'%pt
    
    # class local:
    #   k = 0
    # def str_builder(n, value):
    #   local.k += 1
    #   if ((local.k%6) == 0) or ((n+1)%nz == 0):
    #     local.k = 0
    #     return '%13.5E\n'%value
    #   return '%13.5E'%value
    # s += ''.join([str_builder(i,val) for i,val in enumerate(self.data)])
    return s
    #with open(fname, 'w') as fout:
    #  fout.write(s)
