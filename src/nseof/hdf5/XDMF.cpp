#include <sstream>

#include "XDMF.h"

namespace nseof {

namespace hdf5 {

std::unique_ptr<XDMF> XDMF::fromParameters(Parameters& params, std::string hdf,
                                           int nranks) {
  GeometricParameters& gp = params.geometry;
  ParallelParameters& pp = params.parallel;

  int cellsX = pp.localSize[0];
  int cellsY = pp.localSize[1];
  int cellsZ = gp.dim == 3 ? pp.localSize[2] : 1;

  int pointsX = cellsX + 1;
  int pointsY = cellsY + 1;
  int pointsZ = gp.dim == 3 ? cellsZ + 1 : 1;

  return std::unique_ptr<XDMF>(
      new XDMF(hdf, cellsX, cellsY, cellsZ, pointsX, pointsY, pointsZ, nranks));
}

XDMF::XDMF(std::string hdf, int cellsX, int cellsY, int cellsZ, int pointsX,
           int pointsY, int pointsZ, int nranks)
    : hdf(hdf), nranks(nranks), cellsX(cellsX), cellsY(cellsY), cellsZ(cellsZ) {
  int npoints = pointsX * pointsY * pointsZ;
  std::ostringstream topologyShape;

  // If you give a Z-dimension of 1 in the 2D case, paraview will get confused
  // because there are no real 3D cells to attach the attribute data to
  if (pointsZ > 1) {
    topologyShape << pointsZ << " " << pointsY << " " << pointsX;
  } else {
    topologyShape << pointsY << " " << pointsX;
  }

  tinyxml2::XMLElement* xdmf = this->doc.NewElement("Xdmf");
  xdmf->SetAttribute("Version", "2.0");

  tinyxml2::XMLElement* domain = this->doc.NewElement("Domain");
  xdmf->InsertFirstChild(domain);

  tinyxml2::XMLElement* topology = this->doc.NewElement("Topology");
  topology->SetAttribute("TopologyType", "3DSMesh");
  topology->SetAttribute("Dimensions", topologyShape.str().c_str());
  domain->InsertEndChild(topology);

  for (int i = 0; i < nranks; i++) {
    std::ostringstream name;
    name << "Rank-" << i;

    std::ostringstream hdfPath;
    hdfPath << hdf << ":/Geometries/Rank-" << i;

    // 3 because the points are 3-dimensional
    std::ostringstream geomDim;
    geomDim << npoints << " 3";

    tinyxml2::XMLElement* geometry = this->doc.NewElement("Geometry");
    geometry->SetAttribute("GeometryType", "XYZ");
    geometry->SetAttribute("Name", name.str().c_str());

    tinyxml2::XMLElement* dataItem = this->doc.NewElement("DataItem");
    dataItem->SetAttribute("Format", "HDF");
    dataItem->SetAttribute("Dimensions", geomDim.str().c_str());
    tinyxml2::XMLText* text = this->doc.NewText(hdfPath.str().c_str());
    dataItem->InsertFirstChild(text);

    geometry->InsertFirstChild(dataItem);
    domain->InsertEndChild(geometry);
  }

  this->timesteps = this->doc.NewElement("Grid");
  this->timesteps->SetAttribute("GridType", "Collection");
  this->timesteps->SetAttribute("CollectionType", "Temporal");
  domain->InsertEndChild(this->timesteps);

  this->doc.InsertFirstChild(xdmf);
}

void XDMF::addTimestep(int timestep,
                       const std::vector<std::unique_ptr<nseof::plotting::Reader>>& readers) {
  tinyxml2::XMLElement* grids = this->doc.NewElement("Grid");
  grids->SetAttribute("GridType", "Collection");
  grids->SetAttribute("CollectionType", "Spacial");
  this->timesteps->InsertEndChild(grids);

  tinyxml2::XMLElement* time = this->doc.NewElement("Time");
  time->SetAttribute("Value", timestep);
  grids->InsertEndChild(time);

  for (int i = 0; i < this->nranks; i++) {
    tinyxml2::XMLElement* grid = this->doc.NewElement("Grid");
    grids->InsertEndChild(grid);

    tinyxml2::XMLElement* topology = this->doc.NewElement("Topology");
    topology->SetAttribute("Reference", "/Xdmf/Domain/Topology");
    grid->InsertEndChild(topology);

    std::ostringstream geomRef;
    geomRef << "/Xdmf/Domain/Geometry[@Name='Rank-" << i << "']";
    tinyxml2::XMLElement* geometry = this->doc.NewElement("Geometry");
    geometry->SetAttribute("Reference", geomRef.str().c_str());
    grid->InsertEndChild(geometry);

    for (auto& reader : readers) {
      tinyxml2::XMLElement* attribute = this->doc.NewElement("Attribute");
      attribute->SetAttribute("Name", reader->getName().c_str());
      attribute->SetAttribute("Center", "Cell");

      if (reader->getDim() == 1) {
        attribute->SetAttribute("AttributeType", "Scalar");
      } else {
        attribute->SetAttribute("AttributeType", "Vector");
      }

      grid->InsertEndChild(attribute);

      std::ostringstream dimensions;

      if (this->cellsZ > 1) {
        dimensions << this->cellsZ << " ";
      }

      dimensions << this->cellsY << " " << this->cellsX;

      if (reader->getDim() > 1) {
        dimensions << " " << reader->getDim();
      }

      tinyxml2::XMLElement* dataItem = this->doc.NewElement("DataItem");
      dataItem->SetAttribute("Format", "HDF");
      dataItem->SetAttribute("Dimensions", dimensions.str().c_str());

      std::ostringstream dataPath;
      dataPath << this->hdf << ":/Data/Time-" << timestep << "/Rank-" << i
               << "/" << reader->getName();
      tinyxml2::XMLText* text = this->doc.NewText(dataPath.str().c_str());
      dataItem->InsertEndChild(text);

      attribute->InsertEndChild(dataItem);
    }
  }
}

void XDMF::write(std::string filename) { this->doc.SaveFile(filename.c_str()); }
}
}
