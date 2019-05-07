#include <DGtal/io/readers/MeshReader.h>
#include "NormalAccumulator.h"
#include "GeodesicGraphComputer.h"

typedef typename DGtal::Z3i::RealPoint P3d;
typedef typename DGtal::Z3i::Point P3;
typedef typename DGtal::Z3i::Domain Dom3d;
typedef DGtal::ImageContainerBySTLVector<Dom3d, double> Image3Dd;
typedef DGtal::DigitalSetBySTLSet<Dom3d> AcceptedPointSet;

int main(int argc, char *const *argv)
{

  DGtal::trace.info() << "------------------------------------ "<< std::endl;
  DGtal::trace.info() << "Step 1: Reading input mesh ... ";

  // 1) Reading input mesh:
  DGtal::Mesh<P3d> aMesh;
  aMesh << "example.off";
  DGtal::trace.info() << " [done] " << std::endl;  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;  

  // 2) Init accumulator:
  DGtal::trace.info() << "Step 2: Init accumulator ... ";
  NormalAccumulator acc(7);
  acc.initFromMesh(aMesh);
    
  DGtal::trace.info() << " [done] " << std::endl;  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;  

  
  // 3) Compute accumulation and confidence
  DGtal::trace.info() << "Step 3: Compute accumulation/confidence ... ";
  acc.computeAccumulation();
  acc.computeConfidence();
  
  Image3Dd imageConfidence = acc.getConfidenceImage();
  
  // 4) Apply graph reconstruction from confidence image.
  // 4a) get initial sef of points:
  GeodesicGraphComputer::TSet aConfidenceSet(imageConfidence.domain());
  for (auto const &p: imageConfidence.domain())
      if(imageConfidence(p)>= 0.5)
          aConfidenceSet.insert(p);
  
  P3 p0 = acc.getMaxAccumulationPoint();   
  
  DGtal::trace.info() << "Starting point: " << p0 << std::endl;
  GeodesicGraphComputer gg(4, aConfidenceSet, 3,  acc.getDomain(), p0);
  gg.computeGraphFromGeodesic();
  
  // export vertex
  std::ofstream outV, outE;
  outV.open("vertex.dat");
  outE.open("edges.dat");
  std::vector<P3> vertices = gg.getGraphVerticesRepresentant();
  for (auto const &p: vertices)
    outV << p[0] << " " << p[1] << " " << p[2] << std::endl;

  std::vector<std::pair<unsigned int, unsigned int> > edges;
  edges = gg.myGraph.getEdges();
  for (auto const &e: edges)
    outE << e.first << " " << e.second << std::endl;  
  outV.close(); outE.close();
  DGtal::trace.info() << "] [Done]."<< std::endl ;  
  
  return 0;
}


