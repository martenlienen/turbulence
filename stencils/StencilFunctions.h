#ifndef _DERIVATIVES_H_
#define _DERIVATIVES_H_

#include <math.h>
#include "../Definitions.h"
#include "../Parameters.h"
#include "../FlowFieldTurbA.h"
#include <iostream>
#include <iostream>
using namespace std;

inline void loadLocalNu2D(FlowFieldTurbA& flowField, FLOAT* const localNu,
                          int i, int j) {
  for (int row = -1; row <= 1; row++) {
    for (int column = -1; column <= 1; column++) {
      FLOAT point = flowField.getNu(i + column, j + row);
      localNu[39 + 9 * row + 3 * column] = point;  // x-component
    }
  }
}

inline void loadLocalNu3D(FlowFieldTurbA& flowField, FLOAT* const localNu,
                          int i, int j, int k) {
  for (int layer = -1; layer <= 1; layer++) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        FLOAT point = flowField.getNu(i + column, j + row, k + layer);
        localNu[39 + 27 * layer + 9 * row + 3 * column] = point;
      }
    }
  }
}

// Load the local velocity cube with relevant velocities of the 2D plane
inline void loadLocalVelocity2D(FlowField& flowField,
                                FLOAT* const localVelocity, int i, int j) {
  for (int row = -1; row <= 1; row++) {
    for (int column = -1; column <= 1; column++) {
      const FLOAT* const point =
          flowField.getVelocity().getVector(i + column, j + row);
      localVelocity[39 + 9 * row + 3 * column] = point[0];      // x-component
      localVelocity[39 + 9 * row + 3 * column + 1] = point[1];  // y-component
    }
  }
}

// Load the local velocity cube with surrounding velocities
inline void loadLocalVelocity3D(FlowField& flowField,
                                FLOAT* const localVelocity, int i, int j,
                                int k) {
  for (int layer = -1; layer <= 1; layer++) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        const FLOAT* const point =
            flowField.getVelocity().getVector(i + column, j + row, k + layer);
        localVelocity[39 + 27 * layer + 9 * row + 3 * column] =
            point[0];  // x-component
        localVelocity[39 + 27 * layer + 9 * row + 3 * column + 1] =
            point[1];  // y-component
        localVelocity[39 + 27 * layer + 9 * row + 3 * column + 2] =
            point[2];  // z-component
      }
    }
  }
}

// load local meshsize for 2D -> same as loadLocalVelocity2D, but invoking call
// to meshsize-ptr
inline void loadLocalMeshsize2D(const Parameters& parameters,
                                FLOAT* const localMeshsize, int i, int j) {
  for (int row = -1; row <= 1; row++) {
    for (int column = -1; column <= 1; column++) {
      localMeshsize[39 + 9 * row + 3 * column] =
          parameters.meshsize->getDx(i + column, j + row);
      localMeshsize[39 + 9 * row + 3 * column + 1] =
          parameters.meshsize->getDy(i + column, j + row);
    }
  }
}

// load local meshsize for 3D
inline void loadLocalMeshsize3D(const Parameters& parameters,
                                FLOAT* const localMeshsize, int i, int j,
                                int k) {
  for (int layer = -1; layer <= 1; layer++) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localMeshsize[39 + 27 * layer + 9 * row + 3 * column] =
            parameters.meshsize->getDx(i + column, j + row, k + layer);
        localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 1] =
            parameters.meshsize->getDy(i + column, j + row, k + layer);
        localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 2] =
            parameters.meshsize->getDz(i + column, j + row, k + layer);
      }
    }
  }
}

// Maps an index and a component to the corresponding value in the cube.
inline int mapd(int i, int j, int k, int component) {
  return 39 + 27 * k + 9 * j + 3 * i + component;
}

// Derivative functions. They are applied to a cube of 3x3x3 cells. lv stands
// for the local velocity, lm represents the local mesh sizes
// dudx <-> first derivative of u-component of velocity field w.r.t. x-direction
inline FLOAT dudx(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv [mapd(0,0,0,0)] - lv [mapd(-1,0,0,0)] ) /
  // GeometricParameters::dx;

  // evaluate dudx in the cell center by a central difference
  const int index0 = mapd(0, 0, 0, 0);
  const int index1 = mapd(-1, 0, 0, 0);
  return (lv[index0] - lv[index1]) / lm[index0];
  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "dudx");}

  return tmp2;*/
}

/*
 * 2015.11.07-pm:
 * calculate dudy and dudz in the centre of the cell
 *
 *
 *
 *
 */
inline FLOAT dudy(const FLOAT* const lv, const FLOAT* const lm) {
  return ((lv[mapd(+0, +1, +0, +0)] + lv[mapd(-1, +1, +0, +0)]) / 2 -
          (lv[mapd(+0, -1, +0, +0)] + lv[mapd(-1, -1, +0, +0)]) / 2) /
         (lm[mapd(+0, -1, +0, +1)] * 0.5 + lm[mapd(+0, +0, +0, +1)] +
          lm[mapd(+0, +1, +0, +1)] * 0.5);
}

inline FLOAT dvdy(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv [mapd(0,0,0,1)] - lv [mapd(0,-1,0,1)] ) /
  // GeometricParameters::dy;
  const int index0 = mapd(0, 0, 0, 1);
  const int index1 = mapd(0, -1, 0, 1);
  return (lv[index0] - lv[index1]) / lm[index0];

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "dvdy");}

  return tmp2;*/
}

/*
 * 2015.11.07-pm:
 * calculate dvdx and dvdz in the centre of the cell
 *
 *
 *
 *
 */
inline FLOAT dvdx(const FLOAT* const lv, const FLOAT* const lm) {
  return ((lv[mapd(+1, +0, +0, +1)] + lv[mapd(+1, -1, +0, +1)]) / 2 -
          (lv[mapd(-1, +0, +0, +1)] + lv[mapd(-1, -1, +0, +1)]) / 2) /
         (lm[mapd(-1, +0, +0, +0)] * 0.5 + lm[mapd(+0, +0, +0, +0)] +
          lm[mapd(+1, +0, +0, +0)] * 0.5);
}

inline FLOAT dwdz(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv [mapd(0,0,0,2)] - lv [mapd(0,0,-1,2)] ) /
  // GeometricParameters::dz;
  const int index0 = mapd(0, 0, 0, 2);
  const int index1 = mapd(0, 0, -1, 2);
  return (lv[index0] - lv[index1]) / lm[index0];

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "dwdz");}

  return tmp2;*/
}

// second derivative of u-component w.r.t. x-direction, evaluated at the
// location of the u-component
inline FLOAT d2udx2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv[mapd(1,0,0,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(-1,0,0,0)]
  // )
  //    / ( GeometricParameters::dx * GeometricParameters::dx );

  // evaluate the second derivative at the location of the u-component of the
  // velocity field;
  // we therefore use the two neighbouring u-components and assume arbitrary
  // mesh sizes in both
  // directions -> the formula arises from a straight-forward taylor expansion.
  // -> for equal meshsizes, we obtain the usual [1 -2 1]-like stencil
  const int index_M1 = mapd(-1, 0, 0, 0);
  const int index_0 = mapd(0, 0, 0, 0);
  const int index_P1 = mapd(1, 0, 0, 0);

  const FLOAT dx0 = lm[index_0];
  const FLOAT dx1 = lm[index_P1];
  const FLOAT dxSum = dx0 + dx1;
  return 2.0 * (lv[index_P1] / (dx1 * dxSum) - lv[index_0] / (dx1 * dx0) +
                lv[index_M1] / (dx0 * dxSum));

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2udx2");}

  return tmp2;*/
}

inline FLOAT d2udy2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1=( lv[mapd(0,1,0,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(0,-1,0,0)]
  // )
  //    / ( GeometricParameters::dy * GeometricParameters::dy );
  // average mesh sizes, since the component u is located in the middle of the
  // cell's face
  const FLOAT dy_M1 = lm[mapd(0, -1, 0, 1)];
  const FLOAT dy_0 = lm[mapd(0, 0, 0, 1)];
  const FLOAT dy_P1 = lm[mapd(0, 1, 0, 1)];
  const FLOAT dy0 = 0.5 * (dy_0 + dy_M1);
  const FLOAT dy1 = 0.5 * (dy_0 + dy_P1);
  const FLOAT dySum = dy0 + dy1;
  return 2.0 * (lv[mapd(0, 1, 0, 0)] / (dy1 * dySum) -
                lv[mapd(0, 0, 0, 0)] / (dy1 * dy0) +
                lv[mapd(0, -1, 0, 0)] / (dy0 * dySum));

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2udy2");}

  return tmp2;*/
}

inline FLOAT d2udz2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv[mapd(0,0,1,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(0,0,-1,0)]
  // )
  //    / ( GeometricParameters::dz * GeometricParameters::dz );
  const FLOAT dz_M1 = lm[mapd(0, 0, -1, 2)];
  const FLOAT dz_0 = lm[mapd(0, 0, 0, 2)];
  const FLOAT dz_P1 = lm[mapd(0, 0, 1, 2)];
  const FLOAT dz0 = 0.5 * (dz_0 + dz_M1);
  const FLOAT dz1 = 0.5 * (dz_0 + dz_P1);
  const FLOAT dzSum = dz0 + dz1;
  return 2.0 * (lv[mapd(0, 0, 1, 0)] / (dz1 * dzSum) -
                lv[mapd(0, 0, 0, 0)] / (dz1 * dz0) +
                lv[mapd(0, 0, -1, 0)] / (dz0 * dzSum));
  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2udz2");}

  return tmp2;*/
}

/** second derivative of the v-component, evaluated at the location of the
 * v-component */
inline FLOAT d2vdx2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv[mapd(1,0,0,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(-1,0,0,1)]
  // )
  //    / ( GeometricParameters::dx * GeometricParameters::dx );
  const FLOAT dx_M1 = lm[mapd(-1, 0, 0, 0)];
  const FLOAT dx_0 = lm[mapd(0, 0, 0, 0)];
  const FLOAT dx_P1 = lm[mapd(1, 0, 0, 0)];
  const FLOAT dx0 = 0.5 * (dx_0 + dx_M1);
  const FLOAT dx1 = 0.5 * (dx_0 + dx_P1);
  const FLOAT dxSum = dx0 + dx1;
  return 2.0 * (lv[mapd(1, 0, 0, 1)] / (dx1 * dxSum) -
                lv[mapd(0, 0, 0, 1)] / (dx1 * dx0) +
                lv[mapd(-1, 0, 0, 1)] / (dx0 * dxSum));

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2vdx2");}

  return tmp2;*/
}

inline FLOAT d2vdy2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv[mapd(0,1,0,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(0,-1,0,1)]
  // )
  //    / ( GeometricParameters::dy * GeometricParameters::dy );
  const int index_M1 = mapd(0, -1, 0, 1);
  const int index_0 = mapd(0, 0, 0, 1);
  const int index_P1 = mapd(0, 1, 0, 1);

  const FLOAT dy0 = lm[index_0];
  const FLOAT dy1 = lm[index_P1];
  const FLOAT dySum = dy0 + dy1;
  return 2.0 * (lv[index_P1] / (dy1 * dySum) - lv[index_0] / (dy1 * dy0) +
                lv[index_M1] / (dy0 * dySum));

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2vdy2");}

  return tmp2;*/
}

inline FLOAT d2vdz2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv[mapd(0,0,1,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(0,0,-1,1)]
  // )
  //    / ( GeometricParameters::dz * GeometricParameters::dz );
  const FLOAT dz_M1 = lm[mapd(0, 0, -1, 2)];
  const FLOAT dz_0 = lm[mapd(0, 0, 0, 2)];
  const FLOAT dz_P1 = lm[mapd(0, 0, 1, 2)];
  const FLOAT dz0 = 0.5 * (dz_0 + dz_M1);
  const FLOAT dz1 = 0.5 * (dz_0 + dz_P1);
  const FLOAT dzSum = dz0 + dz1;
  return 2.0 * (lv[mapd(0, 0, 1, 1)] / (dz1 * dzSum) -
                lv[mapd(0, 0, 0, 1)] / (dz1 * dz0) +
                lv[mapd(0, 0, -1, 1)] / (dz0 * dzSum));

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2vdz2");}

  return tmp2;*/
}

/** second derivative of the w-component, evaluated at the location of the
 * w-component */
inline FLOAT d2wdx2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv[mapd(1,0,0,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(-1,0,0,2)]
  // )
  //    / ( GeometricParameters::dx * GeometricParameters::dx );
  const FLOAT dx_M1 = lm[mapd(-1, 0, 0, 0)];
  const FLOAT dx_0 = lm[mapd(0, 0, 0, 0)];
  const FLOAT dx_P1 = lm[mapd(1, 0, 0, 0)];
  const FLOAT dx0 = 0.5 * (dx_0 + dx_M1);
  const FLOAT dx1 = 0.5 * (dx_0 + dx_P1);
  const FLOAT dxSum = dx0 + dx1;
  return 2.0 * (lv[mapd(1, 0, 0, 2)] / (dx1 * dxSum) -
                lv[mapd(0, 0, 0, 2)] / (dx1 * dx0) +
                lv[mapd(-1, 0, 0, 2)] / (dx0 * dxSum));

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2wdx2");}

  return tmp2;*/
}

inline FLOAT d2wdy2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv[mapd(0,1,0,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(0,-1,0,2)]
  // )
  //    / ( GeometricParameters::dy * GeometricParameters::dy );
  const FLOAT dy_M1 = lm[mapd(0, -1, 0, 1)];
  const FLOAT dy_0 = lm[mapd(0, 0, 0, 1)];
  const FLOAT dy_P1 = lm[mapd(0, 1, 0, 1)];
  const FLOAT dy0 = 0.5 * (dy_0 + dy_M1);
  const FLOAT dy1 = 0.5 * (dy_0 + dy_P1);
  const FLOAT dySum = dy0 + dy1;
  return 2.0 * (lv[mapd(0, 1, 0, 2)] / (dy1 * dySum) -
                lv[mapd(0, 0, 0, 2)] / (dy1 * dy0) +
                lv[mapd(0, -1, 0, 2)] / (dy0 * dySum));

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2wdy2");}

  return tmp2;*/
}

inline FLOAT d2wdz2(const FLOAT* const lv, const FLOAT* const lm) {
  // double tmp1= ( lv[mapd(0,0,1,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(0,0,-1,2)]
  // )
  //    / ( GeometricParameters::dz * GeometricParameters::dz );
  const int index_M1 = mapd(0, 0, -1, 2);
  const int index_0 = mapd(0, 0, 0, 2);
  const int index_P1 = mapd(0, 0, 1, 2);

  const FLOAT dz0 = lm[index_0];
  const FLOAT dz1 = lm[index_P1];
  const FLOAT dzSum = dz0 + dz1;
  return 2.0 * (lv[index_P1] / (dz1 * dzSum) - lv[index_0] / (dz1 * dz0) +
                lv[index_M1] / (dz0 * dzSum));

  /*if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "d2wdz2");}

  return tmp2;*/
}

/** first-derivative of product (u*v), evaluated at the location of the
 * v-component */
inline FLOAT duvdx(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*
      const FLOAT tmp1= 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv
     [mapd(0,1,0,0)] ) *
                           ( lv [mapd(0,0,0,1)] + lv [mapd(1,0,0,1)] ) ) -
                         ( ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,1,0,0)] ) *
                           ( lv [mapd(-1,0,0,1)] + lv [mapd(0,0,0,1)] ) ) )
        + parameters.solver.gamma *( ( fabs ( lv [mapd(0,0,0,0)] + lv
     [mapd(0,1,0,0)] ) *
                                ( lv [mapd(0,0,0,1)] - lv [mapd(1,0,0,1)] ) ) -
                         ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,1,0,0)] ) *
                                ( lv [mapd(-1,0,0,1)] - lv [mapd(0,0,0,1)] ) ) )
                         ) / lm[mapd(0,0,0,0)];
  */

  const FLOAT hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];  // distance of corner
                                                     // points in x-direction
                                                     // from center v-value
  const FLOAT hxLong0 =
      0.5 *
      (lm[mapd(0, 0, 0, 0)] +
       lm[mapd(-1, 0, 0, 0)]);  // distance between center and west v-value
  const FLOAT hxLong1 =
      0.5 * (lm[mapd(0, 0, 0, 0)] +
             lm[mapd(1, 0, 0, 0)]);  // distance between center and east v-value
  const FLOAT hyShort =
      0.5 * lm[mapd(0, 0, 0,
                    1)];  // distance of center u-value from upper edge of cell
  const FLOAT hyLong =
      0.5 * (lm[mapd(0, 0, 0, 1)] +
             lm[mapd(0, 1, 0, 1)]);  // distance of north and center u-value

  const FLOAT u00 = lv[mapd(0, 0, 0, 0)];
  const FLOAT u01 = lv[mapd(0, 1, 0, 0)];
  const FLOAT v00 = lv[mapd(0, 0, 0, 1)];
  const FLOAT v10 = lv[mapd(1, 0, 0, 1)];

  const FLOAT uM10 = lv[mapd(-1, 0, 0, 0)];
  const FLOAT uM11 = lv[mapd(-1, 1, 0, 0)];
  const FLOAT vM10 = lv[mapd(-1, 0, 0, 1)];

  // this a central difference expression for the first-derivative. We therefore
  // linearly interpolate u*v onto the surface of the
  // current cell (in 2D: upper left and upper right corner) and then take the
  // central difference
  const FLOAT secondOrder =
      (((hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01) *
           ((hxLong1 - hxShort) / hxLong1 * v00 + hxShort / hxLong1 * v10) -
       ((hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11) *
           ((hxLong0 - hxShort) / hxLong0 * v00 + hxShort / hxLong0 * vM10)) /
      (2.0 * hxShort);

  // this is a forward-difference in donor-cell style. We apply donor cell and
  // again interpolate the velocity values (u-comp.)
  // onto the surface of the cell. We then apply the standard donor cell scheme.
  // This will, however, result in non-equal
  // mesh spacing evaluations (in case of stretched meshes)
  const FLOAT kr = (hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01;
  const FLOAT kl = (hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11;

  const FLOAT firstOrder = 1.0 / (4.0 * hxShort) *
                           (kr * (v00 + v10) - kl * (vM10 + v00) +
                            fabs(kr) * (v00 - v10) - fabs(kl) * (vM10 - v00));

  // return linear combination of central and donor-cell difference
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;

  //    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1, "Error duv_dx"); }
  return tmp2;
}

/** evaluates first derivative w.r.t. y for u*v at location of u-component. For
 * details on implementation, see duvdx */
inline FLOAT duvdy(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*    const FLOAT tmp1 = 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv
     [mapd(1,0,0,1)] ) *
                           ( lv [mapd(0,0,0,0)] + lv [mapd(0,1,0,0)] ) ) -
                         ( ( lv [mapd(0,-1,0,1)] + lv [mapd(1,-1,0,1)] ) *
                           ( lv [mapd(0,-1,0,0)] + lv [mapd(0,0,0,0)] ) ) ) +
        parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv
     [mapd(1,0,0,1)] ) *
                                ( lv [mapd(0,0,0,0)] - lv [mapd(0,1,0,0)] ) ) -
                         ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(1,-1,0,1)] ) *
                                ( lv [mapd(0,-1,0,0)] - lv [mapd(0,0,0,0)] ) ) )
     ) /
                         lm[mapd(0,0,0,1)];
  */
  const FLOAT hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];  // distance of corner
                                                     // points in x-direction
                                                     // from center v-value
  const FLOAT hyLong0 =
      0.5 *
      (lm[mapd(0, 0, 0, 1)] +
       lm[mapd(0, -1, 0, 1)]);  // distance between center and west v-value
  const FLOAT hyLong1 =
      0.5 * (lm[mapd(0, 0, 0, 1)] +
             lm[mapd(0, 1, 0, 1)]);  // distance between center and east v-value
  const FLOAT hxShort =
      0.5 * lm[mapd(0, 0, 0,
                    0)];  // distance of center u-value from upper edge of cell
  const FLOAT hxLong =
      0.5 * (lm[mapd(0, 0, 0, 0)] +
             lm[mapd(1, 0, 0, 0)]);  // distance of north and center u-value

  const FLOAT v00 = lv[mapd(0, 0, 0, 1)];
  const FLOAT v10 = lv[mapd(1, 0, 0, 1)];
  const FLOAT u00 = lv[mapd(0, 0, 0, 0)];
  const FLOAT u01 = lv[mapd(0, 1, 0, 0)];

  const FLOAT v0M1 = lv[mapd(0, -1, 0, 1)];
  const FLOAT v1M1 = lv[mapd(1, -1, 0, 1)];
  const FLOAT u0M1 = lv[mapd(0, -1, 0, 0)];

  const FLOAT secondOrder =
      (((hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10) *
           ((hyLong1 - hyShort) / hyLong1 * u00 + hyShort / hyLong1 * u01) -
       ((hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1) *
           ((hyLong0 - hyShort) / hyLong0 * u00 + hyShort / hyLong0 * u0M1)) /
      (2.0 * hyShort);

  const FLOAT kr = (hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10;
  const FLOAT kl = (hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1;

  const FLOAT firstOrder = 1.0 / (4.0 * hyShort) *
                           (kr * (u00 + u01) - kl * (u0M1 + u00) +
                            fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;
  //    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"Error duvdy"); }
  return tmp2;
}

/** evaluates first derivative w.r.t. x for u*w at location of w-component. For
 * details on implementation, see duvdx */
inline FLOAT duwdx(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*    const FLOAT tmp1 = 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv
     [mapd(0,0,1,0)] ) *
                           ( lv [mapd(0,0,0,2)] + lv [mapd(1,0,0,2)] ) ) -
                         ( ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,0,1,0)] ) *
                           ( lv [mapd(-1,0,0,2)] + lv [mapd(0,0,0,2)] ) ) ) +
        parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,0)] + lv
     [mapd(0,0,1,0)] ) *
                                ( lv [mapd(0,0,0,2)] - lv [mapd(1,0,0,2)] ) ) -
                         ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,0,1,0)] ) *
                                ( lv [mapd(-1,0,0,2)] - lv [mapd(0,0,0,2)] ) ) )
     ) /
                         lm[mapd(0,0,0,0)];
  */
  const FLOAT hxShort = 0.5 * lm[mapd(0, 0, 0, 0)];  // distance of corner
                                                     // points in x-direction
                                                     // from center v-value
  const FLOAT hxLong0 =
      0.5 *
      (lm[mapd(0, 0, 0, 0)] +
       lm[mapd(-1, 0, 0, 0)]);  // distance between center and west v-value
  const FLOAT hxLong1 =
      0.5 * (lm[mapd(0, 0, 0, 0)] +
             lm[mapd(1, 0, 0, 0)]);  // distance between center and east v-value
  const FLOAT hzShort =
      0.5 * lm[mapd(0, 0, 0,
                    2)];  // distance of center u-value from upper edge of cell
  const FLOAT hzLong =
      0.5 * (lm[mapd(0, 0, 0, 2)] +
             lm[mapd(0, 0, 1, 2)]);  // distance of north and center u-value

  const FLOAT u00 = lv[mapd(0, 0, 0, 0)];
  const FLOAT u01 = lv[mapd(0, 0, 1, 0)];
  const FLOAT w00 = lv[mapd(0, 0, 0, 2)];
  const FLOAT w10 = lv[mapd(1, 0, 0, 2)];

  const FLOAT uM10 = lv[mapd(-1, 0, 0, 0)];
  const FLOAT uM11 = lv[mapd(-1, 0, 1, 0)];
  const FLOAT wM10 = lv[mapd(-1, 0, 0, 2)];

  const FLOAT secondOrder =
      (((hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01) *
           ((hxLong1 - hxShort) / hxLong1 * w00 + hxShort / hxLong1 * w10) -
       ((hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11) *
           ((hxLong0 - hxShort) / hxLong0 * w00 + hxShort / hxLong0 * wM10)) /
      (2.0 * hxShort);

  const FLOAT kr = (hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01;
  const FLOAT kl = (hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11;

  const FLOAT firstOrder = 1.0 / (4.0 * hxShort) *
                           (kr * (w00 + w10) - kl * (wM10 + w00) +
                            fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;
  //    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"Error duwdx");}
  return tmp2;
}

/** evaluates first derivative w.r.t. z for u*w at location of u-component. For
 * details on implementation, see duvdx */
inline FLOAT duwdz(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*    const FLOAT tmp1= 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv
     [mapd(1,0,0,2)] ) *
                           ( lv [mapd(0,0,0,0)] + lv [mapd(0,0,1,0)] ) ) -
                         ( ( lv [mapd(0,0,-1,2)] + lv [mapd(1,0,-1,2)] ) *
                           ( lv [mapd(0,0,-1,0)] + lv [mapd(0,0,0,0)] ) ) ) +
        parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv
     [mapd(1,0,0,2)] ) *
                                ( lv [mapd(0,0,0,0)] - lv [mapd(0,0,1,0)] ) ) -
                         ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(1,0,-1,2)] ) *
                                ( lv [mapd(0,0,-1,0)] - lv [mapd(0,0,0,0)] ) ) )
     ) /
                         lm[mapd(0,0,0,2)];
  */
  const FLOAT hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];  // distance of corner
                                                     // points in x-direction
                                                     // from center v-value
  const FLOAT hzLong0 =
      0.5 *
      (lm[mapd(0, 0, 0, 2)] +
       lm[mapd(0, 0, -1, 2)]);  // distance between center and west v-value
  const FLOAT hzLong1 =
      0.5 * (lm[mapd(0, 0, 0, 2)] +
             lm[mapd(0, 0, 1, 2)]);  // distance between center and east v-value
  const FLOAT hxShort =
      0.5 * lm[mapd(0, 0, 0,
                    0)];  // distance of center u-value from upper edge of cell
  const FLOAT hxLong =
      0.5 * (lm[mapd(0, 0, 0, 0)] +
             lm[mapd(1, 0, 0, 0)]);  // distance of north and center u-value

  const FLOAT w00 = lv[mapd(0, 0, 0, 2)];
  const FLOAT w10 = lv[mapd(1, 0, 0, 2)];
  const FLOAT u00 = lv[mapd(0, 0, 0, 0)];
  const FLOAT u01 = lv[mapd(0, 0, 1, 0)];

  const FLOAT w0M1 = lv[mapd(0, 0, -1, 2)];
  const FLOAT w1M1 = lv[mapd(1, 0, -1, 2)];
  const FLOAT u0M1 = lv[mapd(0, 0, -1, 0)];

  const FLOAT secondOrder =
      (((hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10) *
           ((hzLong1 - hzShort) / hzLong1 * u00 + hzShort / hzLong1 * u01) -
       ((hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1) *
           ((hzLong0 - hzShort) / hzLong0 * u00 + hzShort / hzLong0 * u0M1)) /
      (2.0 * hzShort);

  const FLOAT kr = (hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10;
  const FLOAT kl = (hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1;

  const FLOAT firstOrder = 1.0 / (4.0 * hzShort) *
                           (kr * (u00 + u01) - kl * (u0M1 + u00) +
                            fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;

  //    if (fabs(tmp1-tmp2)> 1.0e-12){handleError(1,"duwdz");}
  return tmp2;
}

/** evaluates first derivative w.r.t. y for v*w at location of w-component. For
 * details on implementation, see duvdx */
inline FLOAT dvwdy(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*    const FLOAT tmp1 =  1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv
     [mapd(0,0,1,1)] ) *
                           ( lv [mapd(0,0,0,2)] + lv [mapd(0,1,0,2)] ) ) -
                         ( ( lv [mapd(0,-1,0,1)] + lv [mapd(0,-1,1,1)] ) *
                           ( lv [mapd(0,-1,0,2)] + lv [mapd(0,0,0,2)] ) ) ) +
        parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv
     [mapd(0,0,1,1)] ) *
                                ( lv [mapd(0,0,0,2)] - lv [mapd(0,1,0,2)] ) ) -
                         ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(0,-1,1,1)] ) *
                                ( lv [mapd(0,-1,0,2)] - lv [mapd(0,0,0,2)] ) ) )
     ) /
                         lm[mapd(0,0,0,1)];
  */
  const FLOAT hyShort = 0.5 * lm[mapd(0, 0, 0, 1)];  // distance of corner
                                                     // points in x-direction
                                                     // from center v-value
  const FLOAT hyLong0 =
      0.5 *
      (lm[mapd(0, 0, 0, 1)] +
       lm[mapd(0, -1, 0, 1)]);  // distance between center and west v-value
  const FLOAT hyLong1 =
      0.5 * (lm[mapd(0, 0, 0, 1)] +
             lm[mapd(0, 1, 0, 1)]);  // distance between center and east v-value
  const FLOAT hzShort =
      0.5 * lm[mapd(0, 0, 0,
                    2)];  // distance of center u-value from upper edge of cell
  const FLOAT hzLong =
      0.5 * (lm[mapd(0, 0, 0, 2)] +
             lm[mapd(0, 0, 1, 2)]);  // distance of north and center u-value

  const FLOAT v00 = lv[mapd(0, 0, 0, 1)];
  const FLOAT v01 = lv[mapd(0, 0, 1, 1)];
  const FLOAT w00 = lv[mapd(0, 0, 0, 2)];
  const FLOAT w10 = lv[mapd(0, 1, 0, 2)];

  const FLOAT vM10 = lv[mapd(0, -1, 0, 1)];
  const FLOAT vM11 = lv[mapd(0, -1, 1, 1)];
  const FLOAT wM10 = lv[mapd(0, -1, 0, 2)];

  const FLOAT secondOrder =
      (((hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01) *
           ((hyLong1 - hyShort) / hyLong1 * w00 + hyShort / hyLong1 * w10) -
       ((hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11) *
           ((hyLong0 - hyShort) / hyLong0 * w00 + hyShort / hyLong0 * wM10)) /
      (2.0 * hyShort);

  const FLOAT kr = (hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01;
  const FLOAT kl = (hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11;

  const FLOAT firstOrder = 1.0 / (4.0 * hyShort) *
                           (kr * (w00 + w10) - kl * (wM10 + w00) +
                            fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;
  //    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dvwdy");}
  return tmp2;
}

/** evaluates first derivative w.r.t. z for v*w at location of v-component. For
 * details on implementation, see duvdx */
inline FLOAT dvwdz(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*    const FLOAT tmp1 = 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv
     [mapd(0,1,0,2)] ) *
                           ( lv [mapd(0,0,0,1)] + lv [mapd(0,0,1,1)] ) ) -
                         ( ( lv [mapd(0,0,-1,2)] + lv [mapd(0,1,-1,2)] ) *
                           ( lv [mapd(0,0,-1,1)] + lv [mapd(0,0,0,1)] ) ) ) +
        parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv
     [mapd(0,1,0,2)] ) *
                                ( lv [mapd(0,0,0,1)] - lv [mapd(0,0,1,1)] ) ) -
                         ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(0,1,-1,2)] ) *
                                ( lv [mapd(0,0,-1,1)] - lv [mapd(0,0,0,1)] ) ) )
     ) /
                         lm[mapd(0,0,0,2)];
  */
  const FLOAT hzShort = 0.5 * lm[mapd(0, 0, 0, 2)];  // distance of corner
                                                     // points in x-direction
                                                     // from center v-value
  const FLOAT hzLong0 =
      0.5 *
      (lm[mapd(0, 0, 0, 2)] +
       lm[mapd(0, 0, -1, 2)]);  // distance between center and west v-value
  const FLOAT hzLong1 =
      0.5 * (lm[mapd(0, 0, 0, 2)] +
             lm[mapd(0, 0, 1, 2)]);  // distance between center and east v-value
  const FLOAT hyShort =
      0.5 * lm[mapd(0, 0, 0,
                    1)];  // distance of center u-value from upper edge of cell
  const FLOAT hyLong =
      0.5 * (lm[mapd(0, 0, 0, 1)] +
             lm[mapd(0, 1, 0, 1)]);  // distance of north and center u-value

  const FLOAT w00 = lv[mapd(0, 0, 0, 2)];
  const FLOAT w10 = lv[mapd(0, 1, 0, 2)];
  const FLOAT v00 = lv[mapd(0, 0, 0, 1)];
  const FLOAT v01 = lv[mapd(0, 0, 1, 1)];

  const FLOAT w0M1 = lv[mapd(0, 0, -1, 2)];
  const FLOAT w1M1 = lv[mapd(0, 1, -1, 2)];
  const FLOAT v0M1 = lv[mapd(0, 0, -1, 1)];

  const FLOAT secondOrder =
      (((hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10) *
           ((hzLong1 - hzShort) / hzLong1 * v00 + hzShort / hzLong1 * v01) -
       ((hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1) *
           ((hzLong0 - hzShort) / hzLong0 * v00 + hzShort / hzLong0 * v0M1)) /
      (2.0 * hzShort);

  const FLOAT kr = (hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10;
  const FLOAT kl = (hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1;

  const FLOAT firstOrder = 1.0 / (4.0 * hzShort) *
                           (kr * (v00 + v01) - kl * (v0M1 + v00) +
                            fabs(kr) * (v00 - v01) - fabs(kl) * (v0M1 - v00));
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;
  //    if (fabs(tmp1-tmp2) > 1.0e-12){std::cout << tmp1 << ", " << tmp2 <<
  //    std::endl;handleError(1,"dvwdz");}
  return tmp2;
}

/** first derivative of u*u w.r.t. x, evaluated at location of u-component. */
inline FLOAT du2dx(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*    const FLOAT tmp1 = 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv
     [mapd(1,0,0,0)] ) *
                           ( lv [mapd(0,0,0,0)] + lv [mapd(1,0,0,0)] ) ) -
                         ( ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) *
                           ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) ) ) +
        parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,0)] + lv
     [mapd(1,0,0,0)] ) *
                                ( lv [mapd(0,0,0,0)] - lv [mapd(1,0,0,0)] ) ) -
                         ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) *
                                ( lv [mapd(-1,0,0,0)] - lv [mapd(0,0,0,0)] ) ) )
     ) /
                         lm[mapd(0,0,0,0)];
  */
  const FLOAT dxShort = 0.5 * lm[mapd(0, 0, 0, 0)];
  const FLOAT dxLong0 = 0.5 * (lm[mapd(-1, 0, 0, 0)] + lm[mapd(0, 0, 0, 0)]);
  const FLOAT dxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);

  const FLOAT u0 = lv[mapd(0, 0, 0, 0)];
  const FLOAT uM1 = lv[mapd(-1, 0, 0, 0)];
  const FLOAT u1 = lv[mapd(1, 0, 0, 0)];

  const FLOAT kr = (dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1;
  const FLOAT kl = (dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1;

  // central difference expression which is second-order accurate for uniform
  // meshes. We interpolate u half-way between
  // neighboured u-component values and afterwards build the central difference
  // for u*u
  const FLOAT secondOrder =
      (((dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1) *
           ((dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1) -
       ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1) *
           ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1)) /
      (2.0 * dxShort);

  // donor-cell like derivative expression. We evaluate u half-way between
  // neighboured u-components and use this as a prediction
  // of the transport direction
  const FLOAT firstOrder =
      1.0 / (4.0 * dxShort) * (kr * (u0 + u1) - kl * (uM1 + u0) +
                               fabs(kr) * (u0 - u1) - fabs(kl) * (uM1 - u0));

  // return linear combination of central- and upwind difference
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;
  //    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"du2dx");}
  return tmp2;
}

/** first derivative of v*v w.r.t. y, evaluated at location of v-component; for
 * details, see du2dx */
inline FLOAT dv2dy(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*    const FLOAT tmp1= 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv
     [mapd(0,1,0,1)] ) *
                           ( lv [mapd(0,0,0,1)] + lv [mapd(0,1,0,1)] ) ) -
                         ( ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) *
                           ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) ) ) +
        parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv
     [mapd(0,1,0,1)] ) *
                                ( lv [mapd(0,0,0,1)] - lv [mapd(0,1,0,1)] ) ) -
                         ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) *
                                ( lv [mapd(0,-1,0,1)] - lv [mapd(0,0,0,1)] ) ) )
     ) /
                         lm[mapd(0,0,0,1)];
  */
  const FLOAT dyShort = 0.5 * lm[mapd(0, 0, 0, 1)];
  const FLOAT dyLong0 = 0.5 * (lm[mapd(0, -1, 0, 1)] + lm[mapd(0, 0, 0, 1)]);
  const FLOAT dyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);

  const FLOAT v0 = lv[mapd(0, 0, 0, 1)];
  const FLOAT vM1 = lv[mapd(0, -1, 0, 1)];
  const FLOAT v1 = lv[mapd(0, 1, 0, 1)];

  const FLOAT kr = (dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1;
  const FLOAT kl = (dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1;

  const FLOAT secondOrder =
      (((dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1) *
           ((dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1) -
       ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1) *
           ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1)) /
      (2.0 * dyShort);
  const FLOAT firstOrder =
      1.0 / (4.0 * dyShort) * (kr * (v0 + v1) - kl * (vM1 + v0) +
                               fabs(kr) * (v0 - v1) - fabs(kl) * (vM1 - v0));
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;
  //    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dv2dy");}
  return tmp2;
}

/** first derivative of w*w w.r.t. z, evaluated at location of w-component; for
 * details, see du2dx */
inline FLOAT dw2dz(const FLOAT* const lv, const Parameters& parameters,
                   const FLOAT* const lm) {
  /*    const FLOAT tmp1= 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv
     [mapd(0,0,1,2)] ) *
                           ( lv [mapd(0,0,0,2)] + lv [mapd(0,0,1,2)] ) ) -
                         ( ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) *
                           ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) ) ) +
        parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv
     [mapd(0,0,1,2)] ) *
                                ( lv [mapd(0,0,0,2)] - lv [mapd(0,0,1,2)] ) ) -
                         ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) *
                                ( lv [mapd(0,0,-1,2)] - lv [mapd(0,0,0,2)] ) ) )
     ) /
                         lm[mapd(0,0,0,2)];
  */
  const FLOAT dzShort = 0.5 * lm[mapd(0, 0, 0, 2)];
  const FLOAT dzLong0 = 0.5 * (lm[mapd(0, 0, -1, 2)] + lm[mapd(0, 0, 0, 2)]);
  const FLOAT dzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);

  const FLOAT w0 = lv[mapd(0, 0, 0, 2)];
  const FLOAT wM1 = lv[mapd(0, 0, -1, 2)];
  const FLOAT w1 = lv[mapd(0, 0, 1, 2)];

  const FLOAT kr = (dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1;
  const FLOAT kl = (dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1;

  const FLOAT secondOrder =
      (((dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1) *
           ((dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1) -
       ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1) *
           ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1)) /
      (2.0 * dzShort);
  const FLOAT firstOrder =
      1.0 / (4.0 * dzShort) * (kr * (w0 + w1) - kl * (wM1 + w0) +
                               fabs(kr) * (w0 - w1) - fabs(kl) * (wM1 - w0));
  const FLOAT tmp2 = (1.0 - parameters.solver.gamma) * secondOrder +
                     parameters.solver.gamma * firstOrder;
  //    if (fabs(tmp1-tmp2) > 1.0e-12){handleError(1,"dw2dz");}
  return tmp2;
}

inline FLOAT computeF2D(const FLOAT* const localVelocity,
                        const FLOAT* const localMeshsize,
                        const Parameters& parameters, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 0)] +
         dt * (1 / parameters.flow.Re * (d2udx2(localVelocity, localMeshsize) +
                                         d2udy2(localVelocity, localMeshsize)) -
               du2dx(localVelocity, parameters, localMeshsize) -
               duvdy(localVelocity, parameters, localMeshsize) +
               parameters.environment.gx);
}

inline FLOAT computeG2D(const FLOAT* const localVelocity,
                        const FLOAT* const localMeshsize,
                        const Parameters& parameters, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 1)] +
         dt * (1 / parameters.flow.Re * (d2vdx2(localVelocity, localMeshsize) +
                                         d2vdy2(localVelocity, localMeshsize)) -
               duvdx(localVelocity, parameters, localMeshsize) -
               dv2dy(localVelocity, parameters, localMeshsize) +
               parameters.environment.gy);
}

inline FLOAT computeF3D(const FLOAT* const localVelocity,
                        const FLOAT* const localMeshsize,
                        const Parameters& parameters, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 0)] +
         dt * (1 / parameters.flow.Re * (d2udx2(localVelocity, localMeshsize) +
                                         d2udy2(localVelocity, localMeshsize) +
                                         d2udz2(localVelocity, localMeshsize)) -
               du2dx(localVelocity, parameters, localMeshsize) -
               duvdy(localVelocity, parameters, localMeshsize) -
               duwdz(localVelocity, parameters, localMeshsize) +
               parameters.environment.gx);
}

inline FLOAT computeG3D(const FLOAT* const localVelocity,
                        const FLOAT* const localMeshsize,
                        const Parameters& parameters, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 1)] +
         dt * (1 / parameters.flow.Re * (d2vdx2(localVelocity, localMeshsize) +
                                         d2vdy2(localVelocity, localMeshsize) +
                                         d2vdz2(localVelocity, localMeshsize)) -
               dv2dy(localVelocity, parameters, localMeshsize) -
               duvdx(localVelocity, parameters, localMeshsize) -
               dvwdz(localVelocity, parameters, localMeshsize) +
               parameters.environment.gy);
}

inline FLOAT computeH3D(const FLOAT* const localVelocity,
                        const FLOAT* const localMeshsize,
                        const Parameters& parameters, FLOAT dt) {
  return localVelocity[mapd(0, 0, 0, 2)] +
         dt * (1 / parameters.flow.Re * (d2wdx2(localVelocity, localMeshsize) +
                                         d2wdy2(localVelocity, localMeshsize) +
                                         d2wdz2(localVelocity, localMeshsize)) -
               dw2dz(localVelocity, parameters, localMeshsize) -
               duwdx(localVelocity, parameters, localMeshsize) -
               dvwdy(localVelocity, parameters, localMeshsize) +
               parameters.environment.gz);
}

/*
 * 2015.11.07-pm:
 * calculate nut
 *
 *
 *
 *
 */
inline FLOAT computeNUT2D(const FLOAT* const localVelocity,
                          const FLOAT* const localMeshsize,
                          const Parameters& parameters, const FLOAT& h) {
  const FLOAT S11 = 2 * dudx(localVelocity, localMeshsize);
  const FLOAT S22 = 2 * dvdy(localVelocity, localMeshsize);
  const FLOAT S12 =
      dudy(localVelocity, localMeshsize) + dvdx(localVelocity, localMeshsize);
  const FLOAT grad = sqrt(2 * (S11 * S11 + S22 * S22 + 2 * S12 * S12));
  const FLOAT lm = 0.41 * h;
  return lm * lm * grad;
  //    cout << S11 << " " << S22 << " " << S12 << " " <<grad << endl;
  //    return 0.0;
}

static int X = 0;
static int U = 0;
static int Y = 1;
static int V = 1;
static int Z = 2;
static int W = 2;
static int shift = 0;

/*
 * 2015.11.07-pm:
 * interpolate nu
 *
 *
 *
 *
 */
inline FLOAT interpolateXY(const FLOAT* const lm, const FLOAT* const ln,
                           const int x, const int y) {
  if (shift == 0) {
    return 0.0;
  } else {
    return ((lm[mapd(+0, +0, +0, X)] * ln[mapd(+x, +0, +0, +0)] +
             lm[mapd(+x, +0, +0, X)] * ln[mapd(+0, +0, +0, +0)]) /
                (lm[mapd(+0, +0, +0, X)] + lm[mapd(+x, +0, +0, X)]) *
                lm[mapd(+0, +y, +0, Y)] +
            (lm[mapd(+0, +y, +0, X)] * ln[mapd(+x, +y, +0, +0)] +
             lm[mapd(+x, +y, +0, X)] * ln[mapd(+0, +y, +0, +0)]) /
                (lm[mapd(+0, +y, +0, X)] + lm[mapd(+x, +y, +0, X)]) *
                lm[mapd(+0, +0, +0, Y)]) /
           (lm[mapd(+0, +0, +0, Y)] + lm[mapd(+0, +y, +0, Y)]);
  }
}

inline FLOAT interpolateXZ(const FLOAT* const lm, const FLOAT* const ln,
                           const int x, const int z) {
  if (shift == 0) {
    return 0.0;
  } else {
    return ((lm[mapd(+0, +0, +0, X)] * ln[mapd(+x, +0, +0, +0)] +
             lm[mapd(+x, +0, +0, X)] * ln[mapd(+0, +0, +0, +0)]) /
                (lm[mapd(+0, +0, +0, X)] + lm[mapd(+x, +0, +0, X)]) *
                lm[mapd(+0, +0, +z, Z)] +
            (lm[mapd(+0, +0, +z, X)] * ln[mapd(+x, +0, +z, +0)] +
             lm[mapd(+x, +0, +z, X)] * ln[mapd(+0, +0, +z, +0)]) /
                (lm[mapd(+0, +0, +z, X)] + lm[mapd(+x, +0, +z, X)]) *
                lm[mapd(+0, +0, +0, Z)]) /
           (lm[mapd(+0, +0, +0, Z)] + lm[mapd(+0, +0, +z, Z)]);
  }
}

inline FLOAT interpolateYZ(const FLOAT* const lm, const FLOAT* const ln,
                           const int y, const int z) {
  if (shift == 0) {
    return 0.0;
  } else {
    return ((lm[mapd(+0, +0, +0, Y)] * ln[mapd(+0, +y, +0, +0)] +
             lm[mapd(+0, +y, +0, Y)] * ln[mapd(+0, +0, +0, +0)]) /
                (lm[mapd(+0, +0, +0, Y)] + lm[mapd(+0, +y, +0, Y)]) *
                lm[mapd(+0, +0, +z, Z)] +
            (lm[mapd(+0, +0, +z, Y)] * ln[mapd(+0, +y, +z, +0)] +
             lm[mapd(+0, +y, +z, Y)] * ln[mapd(+0, +0, +z, +0)]) /
                (lm[mapd(+0, +0, +z, Y)] + lm[mapd(+0, +y, +z, Y)]) *
                lm[mapd(+0, +0, +0, Z)]) /
           (lm[mapd(+0, +0, +0, Z)] + lm[mapd(+0, +0, +z, Z)]);
  }
}

/*
 * 2015.11.07-pm:
 * F-stencils for algebric turbulence model
 *
 *
 *
 *
 */
inline FLOAT Fdnududxdx(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + ln[mapd(+1, +0, +0, 0)];
  const FLOAT nu2 = 1 / parameters.flow.Re + ln[mapd(+0, +0, +0, 0)];
  return 2.0 / (lm[mapd(+0, +0, +0, X)] + lm[mapd(+1, +0, +0, X)]) *
         (nu1 * (lv[mapd(+1, +0, +0, U)] - lv[mapd(+0, +0, +0, U)]) /
              (lm[mapd(1, 0, 0, X)]) -
          nu2 * (lv[mapd(+0, +0, +0, U)] - lv[mapd(-1, +0, +0, U)]) /
              (lm[mapd(0, 0, 0, X)]));
}

inline FLOAT Fdnududydy(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateXY(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateXY(lm, ln, +1, -1);
  return 1 / (lm[mapd(+0, +0, +0, Y)]) *
         (nu1 * 2 * (lv[mapd(+0, +1, +0, U)] - lv[mapd(+0, +0, +0, U)]) /
              (lm[mapd(+0, +0, +0, Y)] + lm[mapd(+0, +1, +0, Y)]) -
          nu2 * 2 * (lv[mapd(+0, +0, +0, U)] - lv[mapd(+0, -1, +0, U)]) /
              (lm[mapd(+0, -1, +0, Y)] + lm[mapd(+0, +0, +0, Y)]));
}

inline FLOAT Fdnududzdz(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateXZ(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateXZ(lm, ln, +1, -1);
  return 1 / (lm[mapd(+0, +0, +0, Z)]) *
         (nu1 * 2 * (lv[mapd(+0, +0, +1, U)] - lv[mapd(+0, +0, +0, U)]) /
              (lm[mapd(+0, +0, +0, Z)] + lm[mapd(+0, +0, +1, Z)]) -
          nu2 * 2 * (lv[mapd(+0, +0, +0, U)] - lv[mapd(+0, +0, -1, U)]) /
              (lm[mapd(+0, +0, -1, Z)] + lm[mapd(+0, +0, +0, Z)]));
}

inline FLOAT Fdnudvdydx(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateXY(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateXY(lm, ln, +1, -1);
  return 1 / (lm[mapd(+0, +0, +0, Y)]) *
         (nu1 * 2 * (lv[mapd(+1, +0, +0, V)] - lv[mapd(+0, +0, +0, V)]) /
              (lm[mapd(+0, +0, +0, X)] + lm[mapd(+1, +0, +0, X)]) -
          nu2 * 2 * (lv[mapd(+1, -1, +0, V)] - lv[mapd(+0, -1, +0, V)]) /
              (lm[mapd(+0, -1, +0, X)] + lm[mapd(+1, -1, +0, X)]));
}

inline FLOAT Fdnudwdzdx(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateXZ(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateXZ(lm, ln, +1, -1);
  return 1 / (lm[mapd(+0, +0, +0, Z)]) *
         (nu1 * 2 * (lv[mapd(+1, +0, +0, W)] - lv[mapd(+0, +0, +0, W)]) /
              (lm[mapd(+0, +0, +0, X)] + lm[mapd(+1, +0, +0, X)]) -
          nu2 * 2 * (lv[mapd(+1, +0, -1, W)] - lv[mapd(+0, +0, -1, W)]) /
              (lm[mapd(+0, +0, -1, X)] + lm[mapd(+1, +0, -1, X)]));
}

/*
 * 2015.11.07-pm:
 * G-stencils for algebric turbulence model
 *
 *
 *
 *
 */
inline FLOAT Gdnudvdxdx(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateXY(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateXY(lm, ln, +1, -1);
  return 1 / (lm[mapd(+0, +0, +0, X)]) *
         (nu1 * 2 * (lv[mapd(+1, +0, +0, V)] - lv[mapd(+0, +0, +0, V)]) /
              (lm[mapd(+0, +0, +0, X)] + lm[mapd(+1, +0, +0, X)]) -
          nu2 * 2 * (lv[mapd(+0, +0, +0, V)] - lv[mapd(-1, +0, +0, V)]) /
              (lm[mapd(-1, +0, +0, X)] + lm[mapd(+0, +0, +0, X)]));
}

inline FLOAT Gdnudvdydy(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + ln[mapd(+0, +1, +0, +0)];
  const FLOAT nu2 = 1 / parameters.flow.Re + ln[mapd(+0, +0, +0, +0)];
  return 2.0 / (lm[mapd(+0, +0, +0, Y)] + lm[mapd(+0, +1, +0, Y)]) *
         (nu1 * (lv[mapd(+0, +1, +0, V)] - lv[mapd(+0, +0, +0, V)]) /
              (lm[mapd(+0, +1, +0, Y)]) -
          nu2 * (lv[mapd(+0, +0, +0, V)] - lv[mapd(+0, -1, +0, V)]) /
              (lm[mapd(+0, +0, +0, Y)]));
}

inline FLOAT Gdnudvdzdz(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateYZ(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateYZ(lm, ln, +1, -1);
  return 1 / (lm[mapd(+0, +0, +0, Z)]) *
         (nu1 * 2 * (lv[mapd(+0, +0, +1, V)] - lv[mapd(+0, +0, +0, V)]) /
              (lm[mapd(+0, +0, +0, Z)] + lm[mapd(+0, +0, +1, Z)]) -
          nu2 * 2 * (lv[mapd(+0, +0, +0, V)] - lv[mapd(+0, +0, -1, V)]) /
              (lm[mapd(+0, +0, -1, Z)] + lm[mapd(+0, +0, +0, Z)]));
}

inline FLOAT Gdnududxdy(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateXY(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateXY(lm, ln, -1, +1);
  return 1 / (lm[mapd(+0, +0, +0, X)]) *
         (nu1 * 2 * (lv[mapd(+0, +1, +0, U)] - lv[mapd(+0, +0, +0, U)]) /
              (lm[mapd(+0, +0, +0, Y)] + lm[mapd(+0, +1, +0, Y)]) -
          nu2 * 2 * (lv[mapd(-1, +1, +0, U)] - lv[mapd(-1, +0, +0, U)]) /
              (lm[mapd(-1, +0, +0, Y)] + lm[mapd(-1, +1, +0, Y)]));
}

inline FLOAT Gdnudwdzdy(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateYZ(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateYZ(lm, ln, +1, -1);
  return 1 / (lm[mapd(+0, +0, +0, Z)]) *
         (nu1 * 2 * (lv[mapd(+0, +1, +0, W)] - lv[mapd(+0, +0, +0, W)]) /
              (lm[mapd(+0, +0, +0, Y)] + lm[mapd(+0, +1, +0, Y)]) -
          nu2 * 2 * (lv[mapd(+0, +1, -1, W)] - lv[mapd(+0, +0, -1, W)]) /
              (lm[mapd(+0, +0, -1, Y)] + lm[mapd(+0, +1, -1, Y)]));
}

/*
 * 2015.11.07-pm:
 * H-stencils for algebric turbulence model
 *
 *
 *
 *
 */
inline FLOAT Hdnudwdxdx(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateXZ(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateXZ(lm, ln, -1, +1);
  return 1 / (lm[mapd(+0, +0, +0, X)]) *
         (nu1 * 2 * (lv[mapd(+1, +0, +0, W)] - lv[mapd(+0, +0, +0, W)]) /
              (lm[mapd(+0, +0, +0, X)] + lm[mapd(+1, +0, +0, X)]) -
          nu2 * 2 * (lv[mapd(+0, +0, +0, W)] - lv[mapd(-1, +0, +0, W)]) /
              (lm[mapd(-1, +0, +0, X)] + lm[mapd(+0, +0, +0, X)]));
}

inline FLOAT Hdnudwdydy(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateYZ(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateYZ(lm, ln, +1, -1);
  return 1 / (lm[mapd(+0, +0, +0, Y)]) *
         (nu1 * 2 * (lv[mapd(+0, +1, +0, W)] - lv[mapd(+0, +0, +0, W)]) /
              (lm[mapd(+0, +0, +0, Y)] + lm[mapd(+0, +1, +0, Y)]) -
          nu2 * 2 * (lv[mapd(+0, +0, +0, W)] - lv[mapd(+0, -1, +0, W)]) /
              (lm[mapd(+0, -1, +0, Y)] + lm[mapd(+0, +0, +0, Y)]));
}

inline FLOAT Hdnudwdzdz(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + ln[mapd(+0, +0, +1, +0)];
  const FLOAT nu2 = 1 / parameters.flow.Re + ln[mapd(+0, +0, +0, +0)];
  return 2.0 / (lm[mapd(+0, +0, +0, Z)] + lm[mapd(+0, +0, +1, Z)]) *
         (nu1 * (lv[mapd(+0, +0, +1, W)] - lv[mapd(+0, +0, +0, W)]) /
              (lm[mapd(+0, +0, +1, Z)]) -
          nu2 * (lv[mapd(+0, +0, +0, W)] - lv[mapd(+0, +0, -1, W)]) /
              (lm[mapd(+0, +0, +0, Z)]));
}

inline FLOAT Hdnududxdz(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateXZ(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateXZ(lm, ln, -1, +1);
  return 1 / (lm[mapd(+0, +0, +0, X)]) *
         (nu1 * 2 * (lv[mapd(+0, +0, +1, W)] - lv[mapd(+0, +0, +0, W)]) /
              (lm[mapd(+0, +0, +0, Z)] + lm[mapd(+0, +0, +1, Z)]) -
          nu2 * 2 * (lv[mapd(-1, +0, +1, W)] - lv[mapd(-1, +0, +0, W)]) /
              (lm[mapd(-1, +0, +0, Z)] + lm[mapd(-1, +0, +1, Z)]));
}

inline FLOAT Hdnudvdydz(const Parameters& parameters, const FLOAT* const lv,
                        const FLOAT* const lm, const FLOAT* const ln) {
  const FLOAT nu1 = 1 / parameters.flow.Re + interpolateYZ(lm, ln, +1, +1);
  const FLOAT nu2 = 1 / parameters.flow.Re + interpolateYZ(lm, ln, -1, +1);
  return 1 / (lm[mapd(+0, +0, +0, Y)]) *
         (nu1 * 2 * (lv[mapd(+0, +0, +1, V)] - lv[mapd(+0, +0, +0, V)]) /
              (lm[mapd(+0, +0, +0, Z)] + lm[mapd(+0, +0, +1, Z)]) -
          nu2 * 2 * (lv[mapd(+0, -1, +1, V)] - lv[mapd(+0, -1, +0, V)]) /
              (lm[mapd(+0, -1, +0, Z)] + lm[mapd(+0, -1, +1, Z)]));
}

/*
 * 2015.11.07-pm:
 * stencils for calculating F, G - 2D
 *
 *
 *
 *
 */
inline FLOAT computeF2DT(const FLOAT* const localVelocity,
                         const FLOAT* const localMeshsize,
                         const FLOAT* const localNu,
                         const Parameters& parameters, FLOAT dt) {
  cout << "d2udx2" << endl;
  cout << 1 / parameters.flow.Re * d2udx2(localVelocity, localMeshsize) << endl;
  cout << Fdnududxdx(parameters, localVelocity, localMeshsize, localNu) << endl;
  cout << "d2udy2" << endl;
  cout << 1 / parameters.flow.Re * d2udy2(localVelocity, localMeshsize) << endl;
  cout << Fdnududydy(parameters, localVelocity, localMeshsize, localNu) << endl;
  cout << "rest" << endl;
  cout << Fdnudvdydx(parameters, localVelocity, localMeshsize, localNu) << endl;
  //            cout <<
  //            "-----------------------------------------------------------" <<
  //            endl;

  return localVelocity[mapd(0, 0, 0, 0)] +
         dt * (
                  //            1 / parameters.flow.Re *
                  (
                      // + d2udx2(localVelocity, localMeshsize)
                      // + d2udy2(localVelocity, localMeshsize)
                      +Fdnududxdx(parameters, localVelocity, localMeshsize,
                                  localNu) +
                      Fdnududydy(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Fdnudvdydx(parameters, localVelocity, localMeshsize,
                                 localNu)) -
                  du2dx(localVelocity, parameters, localMeshsize) -
                  duvdy(localVelocity, parameters, localMeshsize) +
                  parameters.environment.gx);
}

inline FLOAT computeG2DT(const FLOAT* const localVelocity,
                         const FLOAT* const localMeshsize,
                         const FLOAT* const localNu,
                         const Parameters& parameters, FLOAT dt) {
  cout << "d2vdx2" << endl;
  cout << 1 / parameters.flow.Re * d2vdx2(localVelocity, localMeshsize) << endl;
  cout << Gdnudvdxdx(parameters, localVelocity, localMeshsize, localNu) << endl;
  cout << "d2vdy2" << endl;
  cout << 1 / parameters.flow.Re * d2vdy2(localVelocity, localMeshsize) << endl;
  cout << Gdnudvdydy(parameters, localVelocity, localMeshsize, localNu) << endl;
  cout << "rest" << endl;
  cout << Gdnududxdy(parameters, localVelocity, localMeshsize, localNu) << endl;
  cout << "-----------------------------------------------------------" << endl;

  return localVelocity[mapd(0, 0, 0, 1)] +
         dt * (
                  //            1 / parameters.flow.Re *
                  (
                      // + d2vdx2(localVelocity, localMeshsize)
                      // + d2vdy2(localVelocity, localMeshsize)
                      +Gdnudvdxdx(parameters, localVelocity, localMeshsize,
                                  localNu) +
                      Gdnudvdydy(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Gdnududxdy(parameters, localVelocity, localMeshsize,
                                 localNu)) -
                  duvdx(localVelocity, parameters, localMeshsize) -
                  dv2dy(localVelocity, parameters, localMeshsize) +
                  parameters.environment.gy);
}

/*
 * 2015.11.07-pm:
 * stencils for calculating F, G, H - 3D
 *
 *
 *
 *
 */
inline FLOAT computeF3DT(const FLOAT* const localVelocity,
                         const FLOAT* const localMeshsize,
                         const FLOAT* const localNu,
                         const Parameters& parameters, FLOAT dt) {
  //        cout << 1 / parameters.flow.Re *d2udx2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Fdnududxdx(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << 1 / parameters.flow.Re *d2udy2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Fdnududydy(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << 1 / parameters.flow.Re *d2udz2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Fdnududzdz(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << Fdnudvdydx(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << Fdnudwdzdx(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout <<
  //        "-----------------------------------------------------------" <<
  //        endl;

  return localVelocity[mapd(0, 0, 0, 0)] +
         dt * (
                  // 1 / parameters.flow.Re *
                  (
                      // d2udx2(localVelocity, localMeshsize) +
                      // d2udy2(localVelocity, localMeshsize) +
                      // d2udz2(localVelocity, localMeshsize)
                      +Fdnududxdx(parameters, localVelocity, localMeshsize,
                                  localNu) +
                      Fdnududydy(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Fdnududzdz(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Fdnudvdydx(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Fdnudwdzdx(parameters, localVelocity, localMeshsize,
                                 localNu)) -
                  du2dx(localVelocity, parameters, localMeshsize) -
                  duvdy(localVelocity, parameters, localMeshsize) -
                  duwdz(localVelocity, parameters, localMeshsize) +
                  parameters.environment.gx);
}

inline FLOAT computeG3DT(const FLOAT* const localVelocity,
                         const FLOAT* const localMeshsize,
                         const FLOAT* const localNu,
                         const Parameters& parameters, FLOAT dt) {
  //        cout << 1 / parameters.flow.Re *d2vdx2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Gdnudvdxdx(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << 1 / parameters.flow.Re *d2vdy2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Gdnudvdydy(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << 1 / parameters.flow.Re *d2vdz2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Gdnudvdzdz(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << Gdnududxdy(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << Gdnudwdzdy(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout <<
  //        "-----------------------------------------------------------" <<
  //        endl;

  return localVelocity[mapd(0, 0, 0, 1)] +
         dt * (
                  // 1 / parameters.flow.Re *
                  (
                      // d2vdx2(localVelocity, localMeshsize) +
                      // d2vdy2(localVelocity, localMeshsize) +
                      // d2vdz2(localVelocity, localMeshsize)
                      +Gdnudvdxdx(parameters, localVelocity, localMeshsize,
                                  localNu) +
                      Gdnudvdydy(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Gdnudvdzdz(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Gdnududxdy(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Gdnudwdzdy(parameters, localVelocity, localMeshsize,
                                 localNu)) -
                  dv2dy(localVelocity, parameters, localMeshsize) -
                  duvdx(localVelocity, parameters, localMeshsize) -
                  dvwdz(localVelocity, parameters, localMeshsize) +
                  parameters.environment.gy);
}

inline FLOAT computeH3DT(const FLOAT* const localVelocity,
                         const FLOAT* const localMeshsize,
                         const FLOAT* const localNu,
                         const Parameters& parameters, FLOAT dt) {
  //        cout << 1 / parameters.flow.Re *d2wdx2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Hdnudwdxdx(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << 1 / parameters.flow.Re *d2wdy2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Hdnudwdydy(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << 1 / parameters.flow.Re *d2wdz2(localVelocity, localMeshsize)
  //        << endl;
  //        cout << Hdnudwdzdz(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << Hdnududxdz(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout << Hdnudvdydz(parameters, localVelocity, localMeshsize,
  //        localNu) << endl;
  //        cout <<
  //        "-----------------------------------------------------------" <<
  //        endl;

  return localVelocity[mapd(0, 0, 0, 2)] +
         dt * (
                  // 1 / parameters.flow.Re *
                  (
                      // d2wdx2(localVelocity, localMeshsize) +
                      // d2wdy2(localVelocity, localMeshsize) +
                      // d2wdz2(localVelocity, localMeshsize)
                      +Hdnudwdxdx(parameters, localVelocity, localMeshsize,
                                  localNu) +
                      Hdnudwdydy(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Hdnudwdzdz(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Hdnududxdz(parameters, localVelocity, localMeshsize,
                                 localNu) +
                      Hdnudvdydz(parameters, localVelocity, localMeshsize,
                                 localNu)) -
                  dw2dz(localVelocity, parameters, localMeshsize) -
                  duwdx(localVelocity, parameters, localMeshsize) -
                  dvwdy(localVelocity, parameters, localMeshsize) +
                  parameters.environment.gz);
}

#endif
