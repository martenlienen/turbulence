#include "Meshsize.h"

#include "Parameters.h"

UniformMeshsize::UniformMeshsize(
  const Parameters &parameters
): Meshsize(),
_dx(parameters.geometry.lengthX/parameters.geometry.sizeX),
_dy(parameters.geometry.lengthY/parameters.geometry.sizeY),
_dz(parameters.geometry.dim==3 ? parameters.geometry.lengthZ/parameters.geometry.sizeZ : 0.0),
_firstCornerX(parameters.parallel.firstCorner[0]),
_firstCornerY(parameters.parallel.firstCorner[1]),
_firstCornerZ(parameters.geometry.dim==3 ? parameters.parallel.firstCorner[2] : 0)
{
  if (_dx<= 0.0){ handleError(1,"_dx<=0.0!"); }
  if (_dy<= 0.0){ handleError(1,"_dy<=0.0!"); }
  if (parameters.geometry.dim==3){
    if (_dz<=0.0){ handleError(1,"_dz<=0.0!"); }
  }
}

UniformMeshsize::~UniformMeshsize(){}




TanhMeshStretching::TanhMeshStretching(
  const Parameters & parameters,bool stretchX, bool stretchY, bool stretchZ
): Meshsize(), _uniformMeshsize(parameters),
   _lengthX(parameters.geometry.lengthX), _lengthY(parameters.geometry.lengthY),
   _lengthZ(parameters.geometry.dim==3 ? parameters.geometry.lengthZ : 0.0),
   _sizeX(parameters.geometry.sizeX), _sizeY(parameters.geometry.sizeY),
   _sizeZ(parameters.geometry.dim==3 ? parameters.geometry.sizeZ : 1),
   _firstCornerX(parameters.parallel.firstCorner[0]), _firstCornerY(parameters.parallel.firstCorner[1]),
   _firstCornerZ(parameters.geometry.dim==3 ? parameters.parallel.firstCorner[2] : 0),
   _stretchX(stretchX), _stretchY(stretchY), _stretchZ(stretchZ),
   _deltaS(2.7), _tanhDeltaS(tanh(2.7)), // this parameters is chosen as 2.7 as used also in the dissertation by Tobias Neckel
   _dxMin(stretchX ? 0.5*parameters.geometry.lengthX*(1.0 + tanh(_deltaS*(2.0/_sizeX-1.0))/_tanhDeltaS) : _uniformMeshsize.getDx(0,0) ),
   _dyMin(stretchY ? 0.5*parameters.geometry.lengthY*(1.0 + tanh(_deltaS*(2.0/_sizeY-1.0))/_tanhDeltaS) : _uniformMeshsize.getDy(0,0) ),
   _dzMin(stretchZ ? 0.5*parameters.geometry.lengthZ*(1.0 + tanh(_deltaS*(2.0/_sizeZ-1.0))/_tanhDeltaS) : _uniformMeshsize.getDz(0,0,0) )
{
}

TanhMeshStretching::~TanhMeshStretching(){}



