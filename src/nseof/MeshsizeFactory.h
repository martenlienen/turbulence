#ifndef _MESHSIZEFACTORY_H_
#define _MESHSIZEFACTORY_H_

#include <memory>

#include "Definitions.h"
#include "Parameters.h"
#include "MemoizedMesh.h"

namespace nseof {

/** initialises the meshsize in the Parameters. Must be called after configuring
 * (Configuration and PetscParallelConfiguration).
 *  We therefore make use of the singleton/factory pattern.
 *  @author Philipp Neumann
 */
class MeshsizeFactory {
 public:
  static MeshsizeFactory& getInstance() {
    static MeshsizeFactory singleton;
    return singleton;
  }

  void initMeshsize(Parameters& parameters) {
    Meshsize* mesh;

    // initialise meshsize
    switch (parameters.geometry.meshsizeType) {
      // uniform meshsize
      case Uniform:
        mesh = new UniformMeshsize(parameters);
        break;
      // tanh-stretched mesh
      case TanhStretching:
        mesh = new TanhMeshStretching(
            parameters, (bool)parameters.geometry.stretchX,
            (bool)parameters.geometry.stretchY,
            (bool)parameters.geometry.stretchZ);
        break;
      default:
        handleError(1, "Unknown meshsize type!");
        break;
    }

    // check that meshsize is initialised
    if (mesh == NULL) {
      handleError(1, "parameters.meshsize==NULL!");
    } else {
      //parameters.meshsize = mesh;
      parameters.meshsize = new MemoizedMesh(parameters,
                                             std::unique_ptr<Meshsize>(mesh));
    }
  }

 private:
  MeshsizeFactory() {}
  ~MeshsizeFactory() {}
};
}

#endif  // _MESHSIZEFACTORY_H_
