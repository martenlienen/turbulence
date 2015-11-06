#ifndef _MESHSIZE_H_
#define _MESHSIZE_H_

#include "Definitions.h"
#include <cmath>

// forward declaration of Parameters
class Parameters;

enum MeshsizeType{Uniform=0,TanhStretching=1};

/** defines the local mesh size.
 *  @author Philipp Neumann
 */
class Meshsize {
  public:
    virtual ~Meshsize(){}
    // returns the meshsize of cell i,j or i,j,k, respectively
    virtual FLOAT getDx(int i, int j) const = 0;
    virtual FLOAT getDy(int i, int j) const = 0;

    virtual FLOAT getDx(int i, int j, int k) const = 0;
    virtual FLOAT getDy(int i, int j, int k) const = 0;
    virtual FLOAT getDz(int i, int j, int k) const = 0;

    // returns the global geometric position in x-,y-,z-direction
    // of the lower/left/front corner of the local cell at (i,j,k)
    virtual FLOAT getPosX(int i, int j, int k) const = 0;
    virtual FLOAT getPosY(int i, int j, int k) const = 0;
    virtual FLOAT getPosZ(int i, int j, int k) const = 0;

    virtual FLOAT getPosX(int i, int j) const = 0;
    virtual FLOAT getPosY(int i, int j) const = 0;

    // returns the min. meshsize used in this simulation
    // -> required for adaptive time stepping
    virtual FLOAT getDxMin() const = 0;
    virtual FLOAT getDyMin() const = 0;
    virtual FLOAT getDzMin() const = 0;
};


/** implements a uniform, equidistant grid spacing */
class UniformMeshsize: public Meshsize {
  public:
    UniformMeshsize(const Parameters &parameters);
    virtual ~UniformMeshsize();

    virtual FLOAT getDx(int i, int j) const { return _dx; }
    virtual FLOAT getDy(int i, int j) const { return _dy; }

    virtual FLOAT getDx(int i, int j, int k) const { return _dx; }
    virtual FLOAT getDy(int i, int j, int k) const { return _dy; }
    virtual FLOAT getDz(int i, int j, int k) const { return _dz; }

    virtual FLOAT getPosX(int i,int j, int k) const { return _dx*(_firstCornerX-2+i); }
    virtual FLOAT getPosY(int i,int j, int k) const { return _dy*(_firstCornerY-2+j); }
    virtual FLOAT getPosZ(int i,int j, int k) const { return _dz*(_firstCornerZ-2+k); }
    virtual FLOAT getPosX(int i,int j) const { return getPosX(i,j,0); }
    virtual FLOAT getPosY(int i,int j) const { return getPosY(i,j,0); }


    virtual FLOAT getDxMin() const {return _dx;}
    virtual FLOAT getDyMin() const {return _dy;}
    virtual FLOAT getDzMin() const {return _dz;}

  private:
    const FLOAT _dx;
    const FLOAT _dy;
    const FLOAT _dz;
    const int _firstCornerX;
    const int _firstCornerY;
    const int _firstCornerZ;
};


/** implements a stretched mesh for e.g. channel flow. For each dimension, a stretching of the mesh can be introduced towards
 *  the outer boundaries, i.e. if stretchX is true (in constructor), then the mesh will be finer close to the left and right
 *  boundary. The stretching is based on a formular involving tanh-functions, as e.g. used in the dissertation by Tobias Neckel,
 *  chair of scientific computing (TUM).
 *  For non-stretched meshes, the UniformMeshsize implementation is used to create a uniform mesh.
 *  @author Philipp Neumann
 */
class TanhMeshStretching: public Meshsize {
  public:
    TanhMeshStretching(const Parameters & parameters,bool stretchX, bool stretchY, bool stretchZ);
    virtual ~TanhMeshStretching();

    virtual FLOAT getDx(int i,int j) const {
      if (_stretchX){
        return getMeshsize(i,_firstCornerX,_sizeX,_lengthX,_dxMin);
      } else { return _uniformMeshsize.getDx(i,j); }
    }
    virtual FLOAT getDy(int i,int j) const {
      if (_stretchY){
        return getMeshsize(j,_firstCornerY,_sizeY,_lengthY,_dyMin);
      } else { return _uniformMeshsize.getDy(i,j); }
    }
    virtual FLOAT getDx(int i,int j,int k) const {
      return getDx(i,j);
    }
    virtual FLOAT getDy(int i,int j,int k) const {
      return getDy(i,j);
    }
    virtual FLOAT getDz(int i,int j,int k) const {
      if (_stretchZ){
        return getMeshsize(k,_firstCornerZ,_sizeZ,_lengthZ,_dzMin);
      } else { return _uniformMeshsize.getDz(i,j,k); }
    }

    virtual FLOAT getPosX(int i,int j, int k) const {
      if (_stretchX){
        return computeCoordinate(i,_firstCornerX,_sizeX,_lengthX,_dxMin);
      } else { return _uniformMeshsize.getPosX(i,j,k); }
    }
    virtual FLOAT getPosY(int i,int j, int k) const {
      if (_stretchY){ return computeCoordinate(j,_firstCornerY,_sizeY,_lengthY,_dyMin);
      } else { return _uniformMeshsize.getPosY(i,j,k); }
    }
    virtual FLOAT getPosZ(int i,int j, int k) const {
      if (_stretchZ){ return computeCoordinate(k,_firstCornerZ,_sizeZ,_lengthZ,_dzMin);
      } else { return _uniformMeshsize.getPosZ(i,j,k); }
    }
    virtual FLOAT getPosX(int i,int j) const { return getPosX(i,j,0); }
    virtual FLOAT getPosY(int i,int j) const { return getPosY(i,j,0); }


    virtual FLOAT getDxMin() const {return _dxMin; }
    virtual FLOAT getDyMin() const {return _dyMin; }
    virtual FLOAT getDzMin() const {return _dzMin; }


  private:
    // computes the coordinate of the lower/left/front corner of the 1D-cell at index i w.r.t. having "size" cells along
    // an interval of length "length". We refer to local indexing, so "firstCorner" denotes the first non-ghost cell index
    // of this process. We use a stretched mesh for all nodes inside the comput. bounding box, and a regular mesh outside this box,
    // using the meshsize of the next inner cell.
    FLOAT computeCoordinate(int i, int firstCorner,int size, FLOAT length, FLOAT dxMin) const {
      const int index = i-2+firstCorner;
      // equidistant mesh on lower/left part
      if (index < 0){
        return dxMin*index;
      // equidistant mesh on upper/right part
      } else if (index > size-1){
        return length+dxMin*(index-size);
      } else {
        // stretched mesh on lower half of channel -> we check if we are in lower 50% and then use stretching for 2.0*p
        FLOAT p = ((FLOAT) index)/size;
        if (p<0.5){
          return 0.5*length*(1.0 + tanh(_deltaS*(2.0*p-1.0))/_tanhDeltaS);
        // stretched mesh on upper half of channel -> we mirror the stretching
        } else {
          p = ((FLOAT) size-index)/size;
          return length-0.5*length*(1.0 + tanh(_deltaS*(2.0*p-1.0))/_tanhDeltaS);
        }
      }
    }

    // returns the meshsize based on vertex coordinates that span the respective 1D-cell
    FLOAT getMeshsize(int i,int firstCorner,int size, FLOAT length,FLOAT dxMin) const {
      const FLOAT pos0 = computeCoordinate(i,firstCorner,size,length,dxMin);
      const FLOAT pos1 = computeCoordinate(i+1,firstCorner,size,length,dxMin);
      #ifdef DEBUG
      if (pos1-pos0<1.0e-12){handleError(1,"Error TanhMeshStretching::getMeshsize(): dx < 1.0e-12!");}
      #endif
      return pos1-pos0;
    }

    const UniformMeshsize _uniformMeshsize;
    const FLOAT _lengthX;
    const FLOAT _lengthY;
    const FLOAT _lengthZ;
    const int _sizeX;
    const int _sizeY;
    const int _sizeZ;
    const int _firstCornerX;
    const int _firstCornerY;
    const int _firstCornerZ;
    const bool _stretchX;
    const bool _stretchY;
    const bool _stretchZ;
    const FLOAT _deltaS;
    const FLOAT _tanhDeltaS;
    const FLOAT _dxMin;
    const FLOAT _dyMin;
    const FLOAT _dzMin;
};
#endif // _MESHSIZE_H_
