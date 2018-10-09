
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "cube.h"


#define ROOT6    2.449489742783178
#define ROOT60   7.745966692414834
#define ROOT360 18.973665961010276


// check https://github.com/cryos/avogadro/blob/master/libavogadro/src/extensions/surfaces/openqube/gaussianset.cpp

void get_cube(double *grid, double *cube, int64_t nCubePts,
       double *coords, int64_t nCoords,
       int64_t *symmetry, int64_t nSymmetry,
       int64_t *atomIndices,
       double *coefs, double *gtoA, double *gtoCN, int64_t *moIndices, int64_t *gtoIndices, int64_t *cIndices)
{
  int64_t i, j, k, idx;
  double **r, *rSq;
  double sum;
  
  int nProcs = omp_get_num_procs();
  omp_set_num_threads(nProcs);
  fprintf(stdout, "FOUND %d OMP PROCS\nUSING UP TO %d OMP THREADS\n", nProcs , omp_get_max_threads());
  
  /*
  rSq = (double  *)malloc((int)nCoords * sizeof(double));
  r   = (double **)malloc((int)nCoords * sizeof(double*));
  for(i = 0; i < nCoords; i++)
    r[i] = (double *)malloc(3 * sizeof(double));
  */
  #pragma omp parallel shared(grid,cube,coords,symmetry,atomIndices,coefs,gtoA,gtoCN,moIndices,gtoIndices,cIndices) private(r,rSq,i,j,k,idx,sum)
  {
    
    rSq = (double  *)malloc((int)nCoords * sizeof(double));
    r   = (double **)malloc((int)nCoords * sizeof(double*));
    for(i = 0; i < nCoords; i++)
      r[i] = (double *)malloc(3 * sizeof(double));
    
    #pragma omp for
    for(i=0; i<nCubePts; i++) {
      if ((omp_get_thread_num() == 0) && (i%(nCubePts/(nProcs*50)) == 0)) {
        fprintf(stdout, ".");
        fflush(stdout);
      }
      // for (j=0; j<nCoords; j+=3) {
      //   k = i*3;
      //   r[j  ] = grid[k  ] - coords[j  ];
      //   r[j+1] = grid[k+1] - coords[j+1];
      //   r[j+2] = grid[k+2] - coords[j+2];
      //   rSq[j] = pow(r[j][0],2) + pow(r[j][1],2) + pow(r[j][2],2);
      // }
      for (j=0; j<nCoords; j++) {
        r[j][0] = grid[3*i  ] - coords[3*j  ];
        r[j][1] = grid[3*i+1] - coords[3*j+1];
        r[j][2] = grid[3*i+2] - coords[3*j+2];
        rSq[j] = pow(r[j][0],2) + pow(r[j][1],2) + pow(r[j][2],2);
      }
      
      
      // for (j=0; j<nCoords; j++) {
      //   r[j][0] = grid[i][0] - coords[j][0];
      //   r[j][1] = grid[i][1] - coords[j][1];
      //   r[j][2] = grid[i][2] - coords[j][2];
      //   rSq[j] = pow(r[j][0],2) + pow(r[j][1],2) + pow(r[j][2],2);
      // }
      
      sum = 0.0;
      for(k=0; k<nSymmetry; k++) {
        idx = atomIndices[k];
        switch (symmetry[k]) {
          case 0: sum += pointS( k, r[idx], rSq[idx], coefs, gtoA, gtoCN, moIndices, gtoIndices, cIndices); break;
          case 2: sum += pointP( k, r[idx], rSq[idx], coefs, gtoA, gtoCN, moIndices, gtoIndices, cIndices); break;
          case 3: sum += pointD( k, r[idx], rSq[idx], coefs, gtoA, gtoCN, moIndices, gtoIndices, cIndices); break;
          case 4: sum += pointD5(k, r[idx], rSq[idx], coefs, gtoA, gtoCN, moIndices, gtoIndices, cIndices); break;
          case 5: sum += pointF( k, r[idx], rSq[idx], coefs, gtoA, gtoCN, moIndices, gtoIndices, cIndices); break;
          case 6: sum += pointF7(k, r[idx], rSq[idx], coefs, gtoA, gtoCN, moIndices, gtoIndices, cIndices); break;
        }
      }
      cube[i] = sum;
    }
    for(i = 0; i < nCoords; i++)
      free(r[i]);
    free(r);
    free(rSq);
    
    if (omp_get_thread_num() == 0) {
      fprintf(stdout, " Done!\n");
      fflush(stdout);
    }
  }


}



inline double pointS(int64_t moIndex, double *r, double rSq,
                     double *coefs, double *gtoA, double *gtoCN, int64_t *moIndices, int64_t *gtoIndices, int64_t *cIndices) {
  double sum=0.0, coef = coefs[moIndices[moIndex]];
  int64_t i, cIndex;
  
  if (fabs(coef) < 1e-15)
    return 0.0;
  
  for(i=gtoIndices[moIndex], cIndex=cIndices[moIndex]; i<gtoIndices[moIndex+1]; i++)
    sum += gtoCN[cIndex++] * exp(-rSq * gtoA[i]);
    
  return (sum*coef);
}


inline double pointP(int64_t moIndex, double *r, double rSq,
                     double *coefs, double *gtoA, double *gtoCN, int64_t *moIndices, int64_t *gtoIndices, int64_t *cIndices) {
  double x=0.0, y=0.0, z=0.0, tmpGTO;
  int64_t i, cIndex, baseIndex=moIndices[moIndex];
  
  for(i=gtoIndices[moIndex], cIndex=cIndices[moIndex]; i<gtoIndices[moIndex+1]; i++) {
    tmpGTO = exp(-rSq * gtoA[i]);
    x += gtoCN[cIndex++] * r[0] * tmpGTO;
    y += gtoCN[cIndex++] * r[1] * tmpGTO;
    z += gtoCN[cIndex++] * r[2] * tmpGTO;
  }
  
  return (x*coefs[baseIndex  ] +
          y*coefs[baseIndex+1] +
          z*coefs[baseIndex+2]);
}


inline double pointD(int64_t moIndex, double *r, double rSq,
                     double *coefs, double *gtoA, double *gtoCN, int64_t *moIndices, int64_t *gtoIndices, int64_t *cIndices) {
  double xx = 0.0, yy = 0.0, zz = 0.0, xy = 0.0, xz = 0.0, yz = 0.0, tmpGTO;
  int64_t i, cIndex, baseIndex=moIndices[moIndex];
  double rx = r[0], ry = r[1], rz = r[2];
  
  for(i=gtoIndices[moIndex], cIndex=cIndices[moIndex]; i<gtoIndices[moIndex+1]; i++) {
    tmpGTO = exp(-rSq * gtoA[i]);
    xx += gtoCN[cIndex++] * tmpGTO; // Dxx
    yy += gtoCN[cIndex++] * tmpGTO; // Dyy
    zz += gtoCN[cIndex++] * tmpGTO; // Dzz
    xy += gtoCN[cIndex++] * tmpGTO; // Dxy
    xz += gtoCN[cIndex++] * tmpGTO; // Dxz
    yz += gtoCN[cIndex++] * tmpGTO; // Dyz
  }
  
  return (xx * rx * rx * coefs[baseIndex  ] +
          yy * ry * ry * coefs[baseIndex+1] +
          zz * rz * rz * coefs[baseIndex+2] +
          xy * rx * ry * coefs[baseIndex+3] +
          xz * rx * rz * coefs[baseIndex+4] +
          yz * ry * rz * coefs[baseIndex+5]);
}

inline double pointD5(int64_t moIndex, double *r, double rSq,
                     double *coefs, double *gtoA, double *gtoCN, int64_t *moIndices, int64_t *gtoIndices, int64_t *cIndices) {
  double d0 = 0.0, d1p = 0.0, d1n = 0.0, d2p = 0.0, d2n = 0.0, tmpGTO;
  int64_t i, cIndex, baseIndex=moIndices[moIndex];
  double rx = r[0], ry = r[1], rz = r[2];
  
  for(i=gtoIndices[moIndex], cIndex=cIndices[moIndex]; i<gtoIndices[moIndex+1]; i++) {
    tmpGTO = exp(-rSq * gtoA[i]);
    d0  += gtoCN[cIndex++] * tmpGTO;
    d1p += gtoCN[cIndex++] * tmpGTO;
    d1n += gtoCN[cIndex++] * tmpGTO;
    d2p += gtoCN[cIndex++] * tmpGTO;
    d2n += gtoCN[cIndex++] * tmpGTO;
  }
  
  double xx = rx * rx;
  double yy = ry * ry;
  double zz = rz * rz;
  double xy = rx * ry;
  double xz = rx * rz;
  double yz = ry * rz;
  
  return (d0  * coefs[baseIndex  ] * (3*zz - rSq) +
          d1p * coefs[baseIndex+1] * xz +
          d1n * coefs[baseIndex+2] * yz +
          d2p * coefs[baseIndex+3] * (xx - yy) +
          d2n * coefs[baseIndex+4] * xy);
}



inline double pointF(int64_t moIndex, double *r, double rSq,
                     double *coefs, double *gtoA, double *gtoCN, int64_t *moIndices, int64_t *gtoIndices, int64_t *cIndices) {
  int64_t i, cIndex, baseIndex=moIndices[moIndex];
  double tmpGTO, rx = r[0], ry = r[1], rz = r[2];
  double xxx = 0.0;
  double xxy = 0.0;
  double xxz = 0.0;
  double xyy = 0.0;
  double xyz = 0.0;
  double xzz = 0.0;
  double yyy = 0.0;
  double yyz = 0.0;
  double yzz = 0.0;
  double zzz = 0.0;
  
  for(i=gtoIndices[moIndex], cIndex=cIndices[moIndex]; i<gtoIndices[moIndex+1]; i++) {
    tmpGTO = exp(-rSq * gtoA[i]);
    xxx += gtoCN[cIndex++] * tmpGTO;
    xxy += gtoCN[cIndex++] * tmpGTO;
    xxz += gtoCN[cIndex++] * tmpGTO;
    xyy += gtoCN[cIndex++] * tmpGTO;
    xyz += gtoCN[cIndex++] * tmpGTO;
    xzz += gtoCN[cIndex++] * tmpGTO;
    yyy += gtoCN[cIndex++] * tmpGTO;
    yyz += gtoCN[cIndex++] * tmpGTO;
    yzz += gtoCN[cIndex++] * tmpGTO;
    zzz += gtoCN[cIndex++] * tmpGTO;
  }
  
  return (xxx * rx * rx * rx * coefs[baseIndex  ] +
          xxy * rx * rx * ry * coefs[baseIndex+1] +
          xxz * rx * rx * rz * coefs[baseIndex+2] +
          xyy * rx * ry * ry * coefs[baseIndex+3] +
          xyz * rx * ry * rz * coefs[baseIndex+4] +
          xzz * rx * rz * rz * coefs[baseIndex+5] +
          yyy * ry * ry * ry * coefs[baseIndex+6] +
          yyz * ry * ry * rz * coefs[baseIndex+7] +
          yzz * ry * rz * rz * coefs[baseIndex+8] +
          zzz * rz * rz * rz * coefs[baseIndex+9]);
}


inline double pointF7(int64_t moIndex, double *r, double rSq,
                     double *coefs, double *gtoA, double *gtoCN, int64_t *moIndices, int64_t *gtoIndices, int64_t *cIndices) {
  int64_t i, cIndex, baseIndex=moIndices[moIndex];
  double tmpGTO, rx = r[0], ry = r[1], rz = r[2];
  double f0 = 0.0, f1p = 0.0, f1n = 0.0, f2p = 0.0, f2n = 0.0, f3p = 0.0, f3n = 0.0;
  
  for(i=gtoIndices[moIndex], cIndex=cIndices[moIndex]; i<gtoIndices[moIndex+1]; i++) {
    tmpGTO = exp(-rSq * gtoA[i]);
    f0  += gtoCN[cIndex++] * tmpGTO;
    f1p += gtoCN[cIndex++] * tmpGTO;
    f1n += gtoCN[cIndex++] * tmpGTO;
    f2p += gtoCN[cIndex++] * tmpGTO;
    f2n += gtoCN[cIndex++] * tmpGTO;
    f3p += gtoCN[cIndex++] * tmpGTO;
    f3n += gtoCN[cIndex++] * tmpGTO;
  }
  
  double xxx = rx * rx * rx;
  double xxy = rx * rx * ry;
  double xxz = rx * rx * rz;
  double xyy = rx * ry * ry;
  double xyz = rx * ry * rz;
  double xzz = rx * rz * rz;
  double yyy = ry * ry * ry;
  double yyz = ry * ry * rz;
  double yzz = ry * rz * rz;
  double zzz = rz * rz * rz;

  // double root6 = 2.449489742783178;
  // double root60 = 7.745966692414834;
  // double root360 = 18.973665961010276;
  
  return (f0  * coefs[baseIndex  ] * (       zzz - 3.0/2.0 * (xxz + yyz)) +
          f1p * coefs[baseIndex+1] * ((6.0 * xzz - 3.0/2.0 * (xxx + xyy))/ROOT6) +
          f1n * coefs[baseIndex+2] * ((6.0 * yzz - 3.0/2.0 * (xxy + yyy))/ROOT6) +
          f2p * coefs[baseIndex+3] * ((15.0 * (xxz - yyz))/ROOT60) +
          f2n * coefs[baseIndex+4] * ((30.0 * xyz)/ROOT60) +
          f3p * coefs[baseIndex+3] * ((15.0 * xxx - 45.0 * xyy)/ROOT360) +
          f3p * coefs[baseIndex+3] * ((45.0 * xxy - 15.0 * yyy)/ROOT360));
}
// int main() {
//   return 0;
// }
